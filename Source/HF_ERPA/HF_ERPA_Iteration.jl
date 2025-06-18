function HF_ERPA_Iteration(Int_Params::Vector{Any},NN_File::String,NNN_File::String,Params::Vector{Any},N_nu::Matrix{Int64},Orb_Phonon::Matrix{Vector{Int64}},Phonon::Vector{PhState},N_Particle::pnInteger,Particle::pnSVector,N_Hole::pnInteger,Hole::pnSVector,Orb::Vector{NOrb},Orb_NN::NNOrb,VNN::NNInt,TrOp::TranOper)
    # Read parameters ...
    A = Params[1]
    Z = Params[2]
    HbarOmega = Int_Params[1]
    N_max = Params[4]
    N_2max = 2*N_max
    N_3max = Int_Params[4]
    Output_File = Params[8]
    a_max = div((N_max+1)*(N_max+2),2)

    # Iteration parameters ...
    Eta = 1.0
    Eps = 1e-6
    Iter = 0
    Iter_max = 100

    # Bare orbitals ...
    Orb_bare = Make_Orbitals(Params[1],Params[2],Int_Params[2])
    
    # Import HF basis transformation matrices ...
    pU_Import = Output_File * "/Bin/pU.bin"
    pU = Matrix{Float64}(undef,a_max,a_max)
    open(pU_Import, "r") do Read_File
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                ME = read(Read_File, Float64)
                pU[a,b] = ME
            end
        end
    end

    nU_Import = Output_File * "/Bin/nU.bin"
    nU = Matrix{Float64}(undef,a_max,a_max)
    open(nU_Import, "r") do Read_File
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                ME = read(Read_File, Float64)
                nU[a,b] = ME
            end
        end
    end

    # Initialize density matrices ...
    pRho, nRho = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)
    pRho_0, nRho_0 = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)

    # Initialize mean-field 1-body Hamiltonians ...
    pH, nH = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)

    # Initialize parameters for bare NN+NNN interactions ...
    Params_bare = Any[]
    Calc_Params = Any[A, Z, N_max, N_2max, N_3max, "CMS1+2B", true, "Bin"]
    for k in 1:length(Int_Params)
        Params_bare = push!(Params_bare, Int_Params[k])
    end
    for k in eachindex(Calc_Params)
        Params_bare = push!(Params_bare,Calc_Params[k])
    end
    Params_bare = push!(Params_bare, "A" * string(A) * "_Z" * string(Z) *
                    "_hw" * string(HbarOmega) * "_Nmax" * string(N_max) *
                    "_N2max" * string(N_2max) * "_N3max" * string(N_3max) *
                    "_CMS1+2B")

    # Make 1-body kinetic operator ...
    @time T = T1B(N_max,Orb,HbarOmega)

    # 2-body NN interaction & Orbitals ...
    @time VNN_bare, Orb_NN_bare = V2B_Read(NN_File,Params_bare,Orb_bare)

    # 3-body NNN interaction & Orbitals ...
    @time VNNN_bare, Orb_NNN_bare = V3B_NO2B_Read(NNN_File,Params_bare,Orb_bare)
    
    # Initialize iteration parameters ...
    Params_Iter = Any[Params[1],Params[2],Params[3],Params[4],false,Params[6],Params[7],Params[8]]
    
    # Standard RPA calculation - ERPA(0) ...

    # Allocate matrices A & B ...
    @time A, B = HF_ERPA_0_Allocate(Params_Iter,N_nu,Orb_Phonon,Phonon,Particle,Hole,Orb,Orb_NN,VNN)

    # Solve RPA generalized-eigenvalue problem ...
    @time E_RPA, X_RPA, Y_RPA = HF_ERPA_0_Diagonalize(Params_Iter,A,B,N_nu)

    # ERPA iteration ...
    @time while (Eta > Eps) && (Iter < Iter_max)

        # Construct OBDM ...
        @time pRho, nRho = HF_ERPA_OBDM(Params_Iter,N_nu,N_Particle,Particle,N_Hole,Hole,Orb,X_RPA,Y_RPA)

        pH, nH = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)

        # Make 1-body Mean-field Hamiltonian ... in the LHO basis ...
        pH, nH = HF_Iteration(Params_bare,pH,nH,pU * pRho * pU', nU * nRho * nU',Orb_bare,Orb_NN_bare,Orb_NNN_bare,T,VNN_bare,VNNN_bare)

        # Transformation to the HF basis ...
        pH, nH =  pU' * pH * pU, nU' * nH * nU

        VNN = deepcopy(VNN_bare)

        # Update residual interaction
        @time V2B_Res_Density(Params_bare,Orb_bare,Orb_NN_bare,Orb_NNN_bare,VNN,VNNN_bare,pU * pRho * pU', nU * nRho * nU')

        # Transform residual interaction to the HF basis ...
        @time VNN, Orb_NN = V2B_Res(Params_bare,Orb_bare,Orb_NN_bare,VNN)

        # Allocate A & B matrices ...
        @time A, B = HF_ERPA_Allocate(Params_Iter,N_nu,Orb_Phonon,Phonon,Particle,Hole,Orb,Orb_NN,VNN,pRho,nRho,pH,nH)

        # Solve RPA eqs. ...
        @time E_RPA, X_RPA, Y_RPA = HF_ERPA_I_Diagonalize(Params_Iter,A,B,N_nu)

        Sum = 0.0
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                ME = (abs(pRho[a,b] - pRho_0[a,b]) + abs(nRho[a,b] - nRho_0[a,b])) / Float64(2*a_max^2)
                Sum += ME
            end
        end

        VNN_bare, Orb_NN_bare = V2B_Read(NN_File,Params_bare,Orb_bare)

        Eta = Sum
        Iter += 1

        println("\nERPA iteration:   Iteration Number = " * string(Iter) * ",\tEta = " * string(Eta))

        pRho_0, nRho_0 = deepcopy(pRho), deepcopy(nRho)

    end

    println("\nERPA iteration has finished ...")

    # Final calculation with/without Orthogonalization ...

    # Allocate A & B matrices ...
    @time A, B = HF_ERPA_Allocate(Params,N_nu,Orb_Phonon,Phonon,Particle,Hole,Orb,Orb_NN,VNN,pRho,nRho,pH,nH)

    # Solve RPA eqs. ...
    @time E_RPA, X_RPA, Y_RPA = HF_ERPA_Diagonalize(Params,A,B,N_nu,Orb_Phonon,Phonon,Particle,Hole,TrOp,pRho,nRho,X_RPA,Y_RPA)
    
    # Compute the new Mean-Field energy ... ?
    #E_HF = HF_Energy(Params_bare,pU * pRho * pU', nU * nRho * nU',Orb,Orb_NN_bare,Orb_NNN_bare,T,VNN_bare,VNNN_bare)

    return E_RPA, X_RPA, Y_RPA, pRho, nRho
end
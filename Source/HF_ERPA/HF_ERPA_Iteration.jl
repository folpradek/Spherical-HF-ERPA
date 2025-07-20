function HF_ERPA_Iteration(Params::Parameters,N_nu::Matrix{Int64},Orb_Phonon::Matrix{Vector{Int64}},Phonon::Vector{PhState},N_Particle::pnInteger,Particle::pnSVector,N_Hole::pnInteger,Hole::pnSVector,Orb::Vector{NOrb},Orb_NN_bare::NNOrb,Orb_NNN_bare::NNNOrb,T::Matrix{Float64},VNN_bare::NNInt,VNNN_bare::Array{Vector{Vector{Float32}},4})
    # Read parameters ...
    N_max = Params.Calc.Nmax
    a_max = div((N_max+1)*(N_max+2),2)

    # Iteration parameters ...
    Eta = 1.0
    Eps = 1e-6
    Iter = 0
    Iter_max = 100

    # Import residual 2-body interaction ... HF basis ...
    println("\nImporting Residual 2-body interaction ...")
    @time VNN_res, Orb_NN_res = V2B_Res_Import(Params,Orb)

    # Transformation matrices LHO -> HF bases ...
    U_LHO_HF = pnMatrix(Read_Transformation_Matrix(Params,"IO/" * Params.Calc.Path * "/Bin/pU_HF.bin"),
                        Read_Transformation_Matrix(Params,"IO/" * Params.Calc.Path * "/Bin/nU_HF.bin"))

    U_HF_ERPA = pnMatrix(diagm(ones(Float64,a_max)),diagm(ones(Float64,a_max)))

    # Initialize density matrices ...
    Rho = pnMatrix(zeros(Float64,a_max,a_max),zeros(Float64,a_max,a_max))

    Rho_LHO = pnMatrix(zeros(Float64,a_max,a_max),zeros(Float64,a_max,a_max))

    pRho_HF, nRho_HF = zeros(Float64,a_max,a_max),zeros(Float64,a_max,a_max)
    @inbounds for a in 1:a_max
        pRho_HF[a,a] = Orb[a].pO
        nRho_HF[a,a] = Orb[a].nO
    end
    Rho_HF = pnMatrix(pRho_HF,nRho_HF)

    Rho_0 = pnMatrix(deepcopy(pRho_HF), deepcopy(nRho_HF))

    # Initialize 1-body Hamiltonian matrix ...
    H = pnMatrix(zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max))
    
    # Standard RPA calculation - ERPA(0) ...

    # Allocate matrices A & B ...
    println("\nAllocating RPA matrices A & B ...")
    @time A, B = HF_ERPA_0_Allocate(Params,N_nu,Orb_Phonon,Phonon,Particle,Hole,Orb,Orb_NN_res,VNN_res)

    # Solve RPA generalized-eigenvalue problem ...
    println("\nSolving the RPA Generalized-Eigenvalue Problem ...")
    @time E_RPA, X_RPA, Y_RPA = HF_ERPA_0_Diagonalize(Params,A,B,N_nu)

    # ERPA iteration ...
    @time while (Eta > Eps) && (Iter < Iter_max)

        @time Rho, U = HF_ERPA_OBDM_Iteration(Params,N_nu,N_Particle,Particle,N_Hole,Hole,X_RPA,Y_RPA,Rho_HF,Rho_0)

        U_HF_ERPA = pnMatrix(U_HF_ERPA.p * U.p, U_HF_ERPA.n * U.n)

        Rho_0 = pnMatrix(U.p' * Rho_0.p * U.p, U.n' * Rho_0.n * U.n)

        pRho_LHO = U_LHO_HF.p * U_HF_ERPA.p * Rho.p * U_HF_ERPA.p' * U_LHO_HF.p'
        nRho_LHO = U_LHO_HF.n * U_HF_ERPA.n * Rho.n * U_HF_ERPA.n' * U_LHO_HF.n'

        pRho_LHO .= 0.5 .* (pRho_LHO .+ pRho_LHO')
        nRho_LHO .= 0.5 .* (nRho_LHO .+ nRho_LHO')

        Rho_LHO = pnMatrix(pRho_LHO,nRho_LHO)

        # Make 1-body Mean-field Hamiltonian ... in the LHO basis ...
        if Params.Calc.ERPA.ScOBH == true
            H = HF_ERPA_H_MF(Params,Rho_LHO,Orb,Orb_NN_bare,Orb_NNN_bare,T,VNN_bare,VNNN_bare)
        elseif Params.Calc.ERPA.ScOBH == false
            H = HF_ERPA_H_MF(Params,Rho_HF,Orb,Orb_NN_bare,Orb_NNN_bare,T,VNN_bare,VNNN_bare)
        end

        # Transformation to the HF-ERPA basis ...
        H = pnMatrix(deepcopy(U_HF_ERPA.p' * U_LHO_HF.p' * H.p * U_LHO_HF.p * U_HF_ERPA.p), deepcopy(U_HF_ERPA.n' * U_LHO_HF.n' * H.n * U_LHO_HF.n * U_HF_ERPA.n))
        H = pnMatrix(0.5 .* (H.p + H.p'), 0.5 .* (H.n + H.n'))

        println("\nUpdating the residual interaction ...")
        if Params.Calc.ERPA.Sc3N== true
            @time VNN_res = HF_ERPA_VNN(Params,Orb,Orb_NN_bare,Orb_NNN_bare,VNN_bare,VNNN_bare,Rho_LHO,pnMatrix(U_LHO_HF.p * U_HF_ERPA.p, U_LHO_HF.n * U_HF_ERPA.n))
        elseif Params.Calc.ERPA.Sc3N== false && Params.Calc.ERPA.OBDM == "Full"
            @time VNN_res = HF_ERPA_Transform_VNN(Params,Orb,Orb_NN_res,VNN_res,pnMatrix(U.p, U.n))
        end

        # Allocate A & B matrices ...
        println("\nAllocating RPA matrices A & B ...")
        @time A, B = HF_ERPA_Allocate(Params,N_nu,Orb_Phonon,Phonon,Particle,Hole,Orb,Orb_NN_res,VNN_res,Rho,H)

        # Solve RPA eqs. ...
        println("\nSolving the RPA Generalized-Eigenvalue Problem ...")
        @time E_RPA, X_RPA, Y_RPA = HF_ERPA_I_Diagonalize(Params,A,B,N_nu)

        Sum = 0.0
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                ME = (abs(Rho.p[a,b] - Rho_0.p[a,b]) + abs(Rho.n[a,b] - Rho_0.n[a,b])) / Float64(2*a_max^2)
                Sum += ME
            end
        end
        
        if abs(Sum - Eta) < 1e-6
            println("\nHF-ERPA iteration finished - cycled degenerate solution found ...")
            break
        elseif (abs(Sum) < 1e-5) && (Iter > 50)
            println("\nERPA iteration has finished ... Solution converged with lesser precisions 1e-5 ...")
            break
        elseif (abs(Sum) < 1e-4) && (Iter > 75)
            println("\nERPA iteration has finished ... Solution converged with lesser precisions 1e-4 ...")
            break
        end

        Eta = Sum
        Iter += 1

        println("\nERPA iteration:   Iteration Number = " * string(Iter) * ",\tEta = " * string(Eta))

        Rho_0 = pnMatrix(deepcopy(Rho.p),deepcopy(Rho.n))
    end

    println("\nERPA iteration has terminated ...")

    # Final calculation with/without Orthogonalization ...

    # Allocate A & B matrices ...
    println("\nAllocating RPA matrices A & B ...")
    @time A, B = HF_ERPA_Allocate(Params,N_nu,Orb_Phonon,Phonon,Particle,Hole,Orb,Orb_NN_res,VNN_res,Rho,H)

    # Initialize transition operators ...
    @time TrOp = Transition_Operators_Ini(Params,Orb)

    @time TrOp = Transition_Operators_Transform(Params,TrOp,Orb,pnMatrix(U_LHO_HF.p * U_HF_ERPA.p, U_LHO_HF.n * U_HF_ERPA.n))

    # Solve RPA eqs. ...
    println("\nSolving the RPA Generalized-Eigenvalue Problem ...")
    @time E_RPA, X_RPA, Y_RPA = HF_ERPA_Diagonalize(Params,A,B,N_nu,Orb_Phonon,Phonon,Particle,Hole,TrOp,Rho,X_RPA,Y_RPA)

    return E_RPA, X_RPA, Y_RPA, Rho, pnMatrix(U_LHO_HF.p * U_HF_ERPA.p, U_LHO_HF.n * U_HF_ERPA.n)
end
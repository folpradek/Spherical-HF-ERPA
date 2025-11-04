function HF_Solver(Params::Parameters)
    # Make single-particle orbitals - NuHamil ordering
    Orb = Make_Orbitals(Params.Calc.A,Params.Calc.Z,Params.Int.Nmax)

    # Load 1-body kinetic operator ...
    T = T1B(Params.Calc.Nmax,Orb,Params.Int.hw)

    # 2-body bare NN interaction & Orbitals ...
    @time VNN, Orb_NN = V2B_Read(Params,Orb)

    # 3-body bare NNN interaction & Orbitals ...
    @time VNNN, Orb_NNN = V3B_NO2B_Read(Params,Orb)

    # Solve HF equations ...
    @time SPE, C, Rho, h, Iteration = HF_Solve(Params,Orb,Orb_NN,Orb_NNN,T,VNN,VNNN,1e-6)

    # HF energy calculation ...
    E_HF = HF_Energy(Params,Rho,Orb,Orb_NN,Orb_NNN,T,VNN,VNNN)

    # Calculate the total HF kinetic energy ...
    T_HF = Kinetic_Energy(Params,Rho,Orb,T)
        # Older, non-general calculation of kinetic energy & kinetic operator t ...
    #T_HF, t = HF_Kinetic_Energy(Params,Rho,Orb,T)

    # Calculation summary ...
    HF_Summary(Params,E_HF,T_HF,1e-6,Iteration)

    # Calculate & export radial HF densities & radii ...
    Summary_File = "IO/" * Params.Calc.Path * "/HF/HF_Summary.dat"
    Densities_File = "IO/" * Params.Calc.Path * "/HF/Densities/HF_Radial_Densities.dat"
    OBDM_Export(Params,Summary_File,Densities_File,Rho,C,Orb)
        # Older export function ... not general ...
    #HF_Radial_Density(Params,Orb,C)

    # Calculate & export radial HF potential ...
        # To be refined ... export non-local V^HF_lj 
    #HF_Radial_Potential(Params,Orb,C,pnMatrix((h.p .- t.p),(h.n .- t.n)))


    # Single-Particle States export ...
    HF_SPS_Summary(Params,SPE,Orb)

    # Export of single-particle energies & density matrices Rho ...
    HF_Export(Params,C,SPE)

    # Perform Beyond-HF MBPT calculations & export residual interation ...
    if Params.Calc.BMF == true

        # Include NO2B NNN interaction to NN component ...
        @time VNN = V2B_Res_Density(Params,Orb,Orb_NN,Orb_NNN,VNN,VNNN,Rho)

        # Create residual NN (density-dependent) interaction ... in HF basis ...
        @time VNN_res, Orb_NN_res = V2B_Res(Params,Orb,Orb_NN,VNN,C)

        # Perform HF-MBPT calculations ...
        @time HF_MBPT(Params,Orb,Orb_NN_res,VNN_res)

        # Export residual 2-body interaction ...
        @time V2B_Res_Export(Params,Orb_NN_res,VNN_res)

        # Additional export of single-particle orbitals ...
        if Params.Calc.Format == "HRBin" || Params.Calc.Format == "HR"
            @time Orbitals_Export(Params,Orb)
        end
    end

    return
end

function HF_Solve(Params::Parameters,Orb::Vector{NOrb},Orb_NN::NNOrb,Orb_NNN::NNNOrb,T::Matrix{Float64},VNN::NNInt,VNNN::Array{Vector{Vector{Float32}},4},epsilon::Float64)
    # Read calculation parameters ...
    hw = Params.Int.hw
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1) * (N_max + 2), 2)

    # Set some parameters for iteration ...
    Iteration, Iteraction_max = 0, 100
    delta = 1.0

    # Preallocate arrays ...
        # Vectors for single-particle energies ...
    SPE = pnVector(zeros(Float64,a_max),zeros(Float64,a_max))
    SPE_old = pnVector(zeros(Float64,a_max),zeros(Float64,a_max))

        # HF single-particle Hamiltonian matrix ...
    h = pnMatrix(zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max))

        # Initial guess on single-particle states & 1-body density matrix ... LHO orbitals ...
    C = pnMatrix(diagm(ones(Float64,a_max)), diagm(ones(Float64,a_max)))
    Rho = HF_Density_Operator(Orb,a_max,C)

    # Iteration of spherical HF equations ...
    @time while (Iteration < Iteraction_max) && (delta > epsilon)

        # Perform iteration of HF equations ... fills the HF Hamiltonian ...
        h = HF_Allocate(Params,Rho,Orb,Orb_NN,Orb_NNN,T,VNN,VNNN)

        # Diagonalize the HF Hamiltonians ...
        pSPE, pC = eigen(Symmetric(h.p), sortby=+)
        nSPE, nC = eigen(Symmetric(h.n), sortby=+)

        # Allocate single-particle energies & HF orbitals ...
        SPE = pnVector(pSPE,nSPE)
        C = pnMatrix(pC,nC)

        # Reorder HF orbitals & SPEs ...
        C, SPE = HF_Orbital_Ordering(Orb,a_max,C,SPE)

        # Generate new HF density matrix ...
        Rho = HF_Density_Operator(Orb,a_max,C)
    
        # Check on convergence of HF SPEs ...
        delta = (sum(abs.(SPE.p .- SPE_old.p )) + sum(abs.(SPE.n .- SPE_old.n))) / Float64(2 * a_max)

        SPE_old  = pnVector(deepcopy(SPE.p), deepcopy(SPE.n))
        Iteration += 1

        println("Iteration number:   " * string(Iteration) * "   Energy difference:   " * string(round(delta, sigdigits=8)) * " MeV")
    end

    println("\nIteration of HF eqs. with NO2B NN+NNN interaction finished ...")

    return SPE, C, Rho, h, Iteration
end
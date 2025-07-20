function HF_RPA(Params::Parameters,Orb::Vector{NOrb},Orb_NN::NNOrb,VNN::NNInt)
    # Import LHO -> HF transformation matrices ...
    U_HF = pnMatrix(Read_Transformation_Matrix(Params,"IO/" * Params.Calc.Path * "/Bin/pU_HF.bin"),
                Read_Transformation_Matrix(Params,"IO/" * Params.Calc.Path * "/Bin/nU_HF.bin"))

    # Prepare Particle & Hole orbitals ...
    @time N_Particle, Particle, N_Hole, Hole = Make_ParticleHole_Orbitals(Params.Calc.Nmax,Params.Calc.Path,Orb)

    # Prepare 1p-1h phonon states ...
    @time N_Phonon, Phonon = Make_Phonon_Orbitals(N_Particle,Particle,N_Hole,Hole)

    # Initialize transition operators ...
    @time TrOp = Transition_Operators_Ini(Params,Orb)

    @time TrOp = Transition_Operators_Transform(Params,TrOp,Orb,U_HF)

    # Count & pre-index all phonon states in JP subspaces ...
    @time N_nu, Orb_Phonon = HF_RPA_Phonon_Count(Params,N_Phonon,Phonon,Particle,Hole)

    # Allocate matrices A & B ...
    @time A, B = HF_RPA_Allocate(Params,N_nu,Orb_Phonon,Phonon,Particle,Hole,Orb,Orb_NN,VNN)

    # Solve TDA eigenvalue problem and RPA generalized-eigenvalue problem ...
    @time E_TDA, X_TDA, E_RPA, X_RPA, Y_RPA = HF_RPA_Diagonalize(Params,A,B,N_nu,Orb_Phonon,Phonon,Particle,Hole,TrOp)

    # Calculate RPA One-Body Density Matrix (OBDM) ...
    @time Rho_RPA = HF_RPA_OBDM(Params,Orb,Orb_Phonon,Phonon,Particle,Hole,N_nu,Y_RPA)

    @time CI_TDA, CI_RPA = HF_RPA_Collectivity(Params,N_nu,X_TDA,X_RPA,Y_RPA)

    # Calculate RPA correlation energy ...
    @time E_Corr_RPA = HF_RPA_Corr_Energy(Params,N_nu,E_RPA,Y_RPA)

    # Elmag. reduced multipole operators ...
    @time rM_TDA, rM_RPA = HF_RPA_rM(Params,N_nu,Orb_Phonon,Phonon,Particle,Hole,X_TDA,X_RPA,Y_RPA,TrOp)

    # Elmag. reduced transition intensities ...
    @time rB_TDA = HF_RPA_rB(Params,N_nu,rM_TDA)
    @time rB_RPA = HF_RPA_rB(Params,N_nu,rM_RPA)

    # RPA & TDA complete solutions export ...
    @time HF_RPA_Export(Params,Orb,N_nu,E_Corr_RPA,E_TDA,E_RPA,X_RPA,Y_RPA,CI_TDA,CI_RPA,rB_TDA,rB_RPA,Rho_RPA)

    # RPA transition radial densities export ...
    # Requires manual control in corresponding function in HF_RPA_Transitions.jl file
    #   !!! One has to choose the phonons to export !!!
    #@time HF_RPA_Transition_Densities_Export(Params,Orb,N_nu,Orb_Phonon,Phonon,Particle,Hole,X_RPA,Y_RPA,TrOp)

    return
end
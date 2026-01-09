function HF_ScRPA_solver(Params::Parameters)
    # Make s.p. orbitals - NuHamil ordering convention ...
    @time Orb = orbitals_make(Params)

    # Import Residual 2-body NN interaction ...
    println("\nImporting Residual 2-body interaction ...")
    @time V_NN, Orb_NN = O2b_import(Params,Orb,"IO/" * Params.Calc.Path * "/Bin/V2B_HF.bin",Make_Orb_NN=true)

    # Import LHO -> HF transformation matrices ...
    @time C_HF = O1b_import(Params,Orb,"IO/" * Params.Calc.Path * "/Bin/C_HF.bin")

    # Prepare Particle & Hole orbitals ...
    @time N_Particle, Particle, N_Hole, Hole = orbitals_ph_make(Params,Orb)

    # Prepare 1p-1h phonon states ...
    @time N_Phonon, Phonon = orbitals_one_phonon_make(N_Particle,Particle,N_Hole,Hole)

    # Initialize transition operators ...
    @time TrOp = Tr1b_initialize(Params,Orb)

    @time TrOp = Tr1b_transformation(Params,TrOp,Orb,C_HF)

    # Count & pre-index all phonon states in JP subspaces ...
    @time N_nu, Orb_Phonon = HF_RPA_phonon_count(Params,N_Phonon,Phonon)

        # Allocate h_N ... 1-body Hamiltonian ...
        h_N = O1b_import(Params,Orb,"IO/" * Params.Calc.Path * "/Bin/h_HF.bin")

        # Allocate HF-MBPT2 TBDM ...
        Rho_NN = HF_MBPT2_TBDM(Params,h_N,V_NN,Orb,Orb_NN)

        # 3-body NNN interaction & Orbitals ...
        @time V_NNN, Orb_NNN = V3b_no2b_read(Params,Orb)

        # Include 2-body MBPT correlations to h_N ...
        @time h_N = h1b_correlations(Params,h_N,Rho_NN,V_NNN,Orb,Orb_NN,Orb_NNN)

        # Free memory ...
        V_NNN = nothing
        Orb_NNN = nothing
        GC.gc()

    # Allocate matrices A, B, N & M ...
    @time A, B, N, M = HF_ScRPA_allocate(Params,N_nu,Orb_Phonon,Phonon,Particle,Hole,h_N,V_NN,Rho_NN,Orb,Orb_NN)

    # Solve TDA eigenvalue problem and RPA generalized-eigenvalue problem ...
    @time E_TDA, X_TDA, E_RPA, X_RPA, Y_RPA = HF_ScRPA_diagonalize(Params,A,B,N,M,N_nu,Orb_Phonon,Phonon,Particle,Hole,TrOp)

        # Try to allocate the 2-body correlation matrix Sigma ...
        @time Sigma_NN = HF_ScRPA_TBDM_allocate(Params,N_nu,Orb_Phonon,Phonon,N_Particle,Particle,N_Hole,Hole,Orb,Orb_NN,X_RPA,Y_RPA)

    # Calculate RPA One-Body Density Matrix (OBDM) ...
    @time Rho_RPA = HF_RPA_OBDM(Params,Orb,Orb_Phonon,Phonon,Particle,Hole,N_nu,Y_RPA)

    # Transform RPA OBDM to the LHO basis ...
    Rho_RPA = O1B(C_HF.p * Rho_RPA.p * C_HF.p', C_HF.n * Rho_RPA.n * C_HF.n')

    @time CI_TDA, CI_RPA = HF_RPA_collectivity(Params,N_nu,X_TDA,X_RPA,Y_RPA)

    # Evaluation of the RPA correlation energy ... for testing purposes!!!
    HF_RPA_energy_density(Params,N_nu,N_Particle,Particle,N_Hole,Hole,Orb,Orb_NN,V_NN,E_RPA,Y_RPA,Rho_RPA)

    HF_RPA_energy_bosonic(Params,N_nu,E_RPA,Y_RPA)
    #_____________________________________________________________________________________________________

    # Calculate RPA correlation energy ...
    @time E_RPA_corr = HF_RPA_energy(Params,N_nu,E_RPA,Y_RPA)

    # Electromagnetic reduced multipole operators ...
    @time rM_TDA, rM_RPA = HF_RPA_rM(Params,N_nu,Orb_Phonon,Phonon,Particle,Hole,X_TDA,X_RPA,Y_RPA,TrOp)

    # Electromagnetic reduced transition intensities ...
    @time rB_TDA = HF_RPA_rB(Params,N_nu,rM_TDA)
    @time rB_RPA = HF_RPA_rB(Params,N_nu,rM_RPA)

    # Export of RPA & TDA solutions ...
    @time HF_RPA_export(Params,Orb,N_nu,E_RPA_corr,E_TDA,E_RPA,X_RPA,Y_RPA,CI_TDA,CI_RPA,rB_TDA,rB_RPA,Rho_RPA,C_HF)

    return
end
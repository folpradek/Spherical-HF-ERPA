function HF_RPA(Params::Vector{Any},Orb::Vector{NOrb},Orb_NN::NNOrb,VNN::NNInt)
    # Calculation parameters review ...
    # A, Z, HbarOmega, N_max = Params[1], Params[2], Params[3], Params[4]
    # N_2max = 2*N_max
    # J_max = N_2max + 1
    # Orthogonalization = Params[5]
    # Input_File = Params[8]

    # Prepare Particle & Hole orbitals ...
    @time N_Particle, Particle, N_Hole, Hole = Make_ParticleHole_Orbitals(Params[4],Params[8],Orb)

    # Prepare 1p-1h phonon states ...
    @time N_Phonon, Phonon = Make_Phonon_Orbitals(N_Particle,Particle,N_Hole,Hole)

    # Initialize transition operators ...
    @time TrOp = Transition_Operators_Ini(Params[3],Params[4],Orb)

    @time TrOp = Transition_Operators_Transform(Params[4],Params[8],TrOp,Orb)

    # Count & pre-index all phonon states in JP subspaces ...
    @time N_nu, Orb_Phonon = HF_RPA_Phonon_Count(Params,N_Phonon,Phonon,Particle,Hole)

    # Allocate matrices A & B ...
    @time A, B = HF_RPA_Allocate(Params,N_nu,Orb_Phonon,Phonon,Particle,Hole,Orb,Orb_NN,VNN)

    # Solve TDA eigenvalue problem and RPA generalized-eigenvalue problem ...
    @time E_TDA, X_TDA, E_RPA, X_RPA, Y_RPA = HF_RPA_Diagonalize(Params,A,B,N_nu,Orb_Phonon,Phonon,Particle,Hole,TrOp)

    # Calculate RPA One-Body Density Matrix (OBDM) ...
    @time pRho_RPA, nRho_RPA = HF_RPA_OBDM(Params,Orb,Orb_Phonon,Phonon,Particle,Hole,N_nu,Y_RPA)

    @time s_TDA, s_RPA = HF_RPA_Collectivity(Params,N_nu,X_TDA,X_RPA,Y_RPA)

    # Calculate RPA correlation energy ...
    @time E_Corr_RPA = HF_RPA_Corr_Energy(Params,N_nu,E_RPA,Y_RPA)

    # Elmag. reduced multipole operators ...
    @time rM_TDA, rM_RPA = HF_RPA_rM(Params,N_nu,Orb_Phonon,Phonon,Particle,Hole,X_TDA,X_RPA,Y_RPA,TrOp)

    # Elmag. reduced transition intensities ...
    @time rB_TDA = HF_RPA_rB(Params,N_nu,rM_TDA)
    @time rB_RPA = HF_RPA_rB(Params,N_nu,rM_RPA)

    # RPA & TDA complete solutions export ...
    @time HF_RPA_Export(Params,Orb,N_nu,E_Corr_RPA,E_TDA,E_RPA,X_RPA,Y_RPA,s_TDA,s_RPA,rB_TDA,rB_RPA,pRho_RPA,nRho_RPA)

    # RPA transition radial densities export ...
    # Requires manual control in corresponding function in HF_RPA_Transitions.jl file
    #   !!! One has to choose the phonons to export !!!
    #@time HF_RPA_Transition_Densities_Export(Params,Orb,N_nu,Orb_Phonon,Phonon,Particle,Hole,X_RPA,Y_RPA,TrOp)

    return
end
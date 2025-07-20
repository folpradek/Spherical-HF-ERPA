function HF_ERPA(Params::Parameters,Orb::Vector{NOrb},Orb_NN_bare::NNOrb,Orb_NNN_bare::NNNOrb,T::Matrix{Float64},VNN_bare::NNInt,VNNN_bare::Array{Vector{Vector{Float32}},4})
    # Prepare Particle & Hole orbitals ...
    N_Particle, Particle, N_Hole, Hole = Make_ParticleHole_Orbitals(Params.Calc.Nmax,Params.Calc.Path,Orb)

    # Prepare 1p-1h phonon states ...
    N_Phonon, Phonon = Make_Phonon_Orbitals(N_Particle,Particle,N_Hole,Hole)

    # Count & pre-index all phonon states in JP subspaces ...
    N_nu, Orb_Phonon = HF_RPA_Phonon_Count(Params,N_Phonon,Phonon,Particle,Hole)

    # Start ERPA iteration ...
    @time E_RPA, X_RPA, Y_RPA, Rho, U = HF_ERPA_Iteration(Params,N_nu,Orb_Phonon,Phonon,N_Particle,Particle,N_Hole,Hole,Orb,Orb_NN_bare,Orb_NNN_bare,T,VNN_bare,VNNN_bare)

    # Initialize transition operators ...
    TrOp = Transition_Operators_Ini(Params,Orb)

    TrOp = Transition_Operators_Transform(Params,TrOp,Orb,U)

    # Elmag. reduced multipole operators ...
    rM = HF_ERPA_rM(Params,N_nu,Orb_Phonon,Phonon,Particle,Hole,X_RPA,Y_RPA,TrOp,Rho)

    # Elmag. reduced transition intensities ...
    rB = HF_ERPA_rB(Params,N_nu,rM)

    # Collectiviy of 1-phonon ERPA states ...
    CI = HF_ERPA_Collectivity(Params,N_nu,X_RPA,Y_RPA)

    # ERPA ground state energy ...
    E_0_RPA = HF_ERPA_Energy(Params,N_nu,N_Particle,Particle,N_Hole,Hole,E_RPA,Y_RPA)

    # ERPA solutions export ...
    @time HF_ERPA_Export(Params,Orb,N_nu,E_0_RPA,E_RPA,X_RPA,Y_RPA,CI,rB,Rho,U)

    # ERPA transition radial densities export ...
    # Requires manual control in corresponding function in HF_ERPA_Transitions.jl file
    #   !!! One has to manually choose the phonon transition to export !!!
    #@time HF_ERPA_Transition_Densities_Export(Params,Orb,N_nu,Orb_Phonon,Phonon,Particle,Hole,X_RPA,Y_RPA,TrOp)

    return
end
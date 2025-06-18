function HF_ERPA(Int_Params::Vector{Any},NN_File::String,NNN_File::String,Params::Vector{Any},Orb::Vector{NOrb},Orb_NN::NNOrb,VNN::NNInt)
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

    # Transform transition operators to the HF basis ...
    @time TrOp = Transition_Operators_Transform(Params[4],Params[8],TrOp,Orb)

    # Count & pre-index all phonon states in JP subspaces ...
    @time N_nu, Orb_Phonon = HF_ERPA_Phonon_Count(Params,N_Phonon,Phonon,Particle,Hole)

    # Start ERPA iteration ...
    @time E_RPA, X_RPA, Y_RPA, pRho, nRho = HF_ERPA_Iteration(Int_Params,NN_File,NNN_File,Params,N_nu,Orb_Phonon,Phonon,N_Particle,Particle,N_Hole,Hole,Orb,Orb_NN,VNN,TrOp)

    # Elmag. reduced multipole operators ...
    @time rM = HF_ERPA_rM(Params,N_nu,Orb_Phonon,Phonon,Particle,Hole,X_RPA,Y_RPA,TrOp,pRho,nRho)

    # Elmag. reduced transition intensities ...
    @time rB = HF_ERPA_rB(Params,N_nu,rM)

    # Collectiviy of 1-phonon ERPA states ...
    @time s = HF_ERPA_Collectivity(Params,N_nu,X_RPA,Y_RPA)

    # ERPA ground state energy ...
    @time E_0_RPA = HF_ERPA_Energy(Params,N_nu,N_Particle,Particle,N_Hole,Hole,E_RPA,Y_RPA)

    # ERPA solutions export ...
    @time HF_ERPA_Export(Params,Orb,N_nu,E_0_RPA,E_RPA,X_RPA,Y_RPA,s,rB,pRho,nRho)

    # ERPA transition radial densities export ...
    # Requires manual control in corresponding function in HF_ERPA_Transitions.jl file
    #   !!! One has to manually choose the phonon transition to export !!!
    #@time HF_ERPA_Transition_Densities_Export(Params,Orb,N_nu,Orb_Phonon,Phonon,Particle,Hole,X_RPA,Y_RPA,TrOp)

    return
end
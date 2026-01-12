function HF_solver(Params::Parameters)
    # Make single-particle orbitals - NuHamil ordering ...
    Orb = orbitals_make(Params)

    # 1-body kinetic operator ...
    @time T = T1b(Params,Orb)

    # 2-body bare NN interaction & Orbitals ...
    @time V_NN, Orb_NN = V2b_read(Params,Orb)

    # 3-body bare NNN interaction & Orbitals ...
    @time V_NNN, Orb_NNN = V3b_no2b_read(Params,Orb)

    # Solve HF equations ...
    @time h, C, Rho, Iteration = HF_solve(Params,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN)

    # HF energy calculation ...
    @time E_HF = HF_energy(Params,Rho,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN)

    # Calculate the total HF kinetic energy ...
    @time T_HF = T1b_energy(Params,Rho,Orb,T)

    # Calculation summary ...
    @time HF_summary(Params,E_HF,T_HF,Iteration)

    # Calculate & export radial HF densities & radii ...
    Summary_File = "IO/" * Params.Calc.Path * "/HF/HF_Summary.dat"
    Densities_File = "IO/" * Params.Calc.Path * "/HF/Densities/HF_Radial_Densities.dat"
    @time OBDM_export(Params,Orb,Summary_File,Densities_File,Rho,C)

    # Calculate & export radial HF potential ...
        # To be refined ... export non-local V^HF_lj 
    #HF_Radial_Potential(Params,Orb,C,O1B((h.p .- t.p),(h.n .- t.n)))

    # Single-Particle States export ...
    HF_SPS_summary(Params,h,Orb)

    # Export of single-particle HF Hamiltonian & HF single-particle states ... in matrix C ...
    HF_export(Params,C,h,Orb)

    # Perform Beyond-mean-field HF-MBPT calculations & export the residual 2-body NN interation ...
    if Params.Calc.HF.BMF == true
        # Make density-dependent residual NN interaction ... NO2B approximation ...
        @time V_NN = V2b_residual_no2b(Params,Orb,Orb_NN,Orb_NNN,Rho,V_NN,V_NNN)

        # Transform the density-dependent NN interaction to the target HF basis ...
        @time V_NN_Res = O2b_transformation(Params,Orb,Orb_NN,V_NN,C)

        # Export residual 2-body interaction in the binary format ...
        V_NN_Export_Path = "IO/" * Params.Calc.Path * "/Bin/V2B_HF.bin"
        @time O2b_export(Params,Orb_NN,V_NN_Res,V_NN_Export_Path)

        # Perform HF-MBPT calculation ...
        @time HF_MBPT(Params,Orb,Orb_NN,V_NN_Res)

        # Export HF solutions ... orbitals & NN interaction in the Human-Readable Format (HRF) ...
        #   Maybe some refinement needed here???
        if Params.Calc.HF.HRF == true
            println("Exporting the HF single-particle orbitals & residual 2-body NN interaction in the Human-Readable-Format (HRF) ...")
            orbitals_export(Params,Orb)
            O2b_export_HRF(Params,Orb_NN,O_NN,"IO/" * Params.Calc.Path * "/Bin/V2B_HF.dat")
        end

        # Deallocate the residual NN interaction ...
        V_NN_Res = nothing

        # Perform the Garbace Collection ...
        GC.gc()

    end

    # Deallocate the NN & NNN interaction ...
    V_NN, V_NNN = nothing, nothing

    # Perform the Garbace Collection ...
    GC.gc()

    return
end

function HF_solve(Params::Parameters,Orb::Vector{Orb1B},Orb_NN::Orb2B,Orb_NNN::Orb3B,T::O1B,V_NN::O2B,V_NNN::Array{Vector{Vector{Float32}},4})
    # Read calculation parameters ...
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1) * (N_max + 2), 2)
    epsilon = Params.Calc.HF.Tol
    Iteration_max = Params.Calc.HF.IMax

    # Initialize local iteration parameters ...
    dE, Iteration = 1.0, 0

    # Preallocate arrays ...
        # Vectors for single-particle energies ...
    SPE = pnVector(zeros(Float64,a_max),zeros(Float64,a_max))
    SPE_old = pnVector(zeros(Float64,a_max),zeros(Float64,a_max))

        # HF single-particle Hamiltonian matrix ...
    h = O1B(zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max))

        # Initial guess on single-particle states & 1-body density matrix ... LHO orbitals ...
    C = O1B(diagm(ones(Float64,a_max)), diagm(ones(Float64,a_max)))
    Rho = HF_density_operator(a_max,C,Orb)

    # Iteration of spherical HF equations ...
    @time while (Iteration < Iteration_max) && (dE > epsilon)

        # Perform iteration of HF equations ... fills the HF Hamiltonian ...
        h = HF_allocate(Params,Rho,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN)

        # Diagonalize the HF Hamiltonians ...
        pSPE, pC = eigen(Symmetric(h.p),sortby=+)
        nSPE, nC = eigen(Symmetric(h.n),sortby=+)

        # Allocate single-particle energies & HF orbitals ...
        SPE = pnVector(pSPE,nSPE)
        C = O1B(pC,nC)

        # Reorder HF orbitals & SPEs ...
        C, SPE = HF_orbital_ordering(Orb,a_max,C,SPE)

        # Generate new HF density matrix ...
        Rho = HF_density_operator(a_max,C,Orb)
    
        # Check on convergence of HF SPEs ...
        dE = (sum(abs.(SPE.p .- SPE_old.p )) + sum(abs.(SPE.n .- SPE_old.n))) / Float64(2 * a_max)

        SPE_old  = pnVector(deepcopy(SPE.p), deepcopy(SPE.n))
        Iteration += 1

        @printf("\tHF iteration number: %4d   dE = %12.9f MeV\n", Iteration, dE)
    end

    println("\nIteration of HF eqs. with NO2B NN+NNN interaction has finished ...")

    return O1B(diagm(SPE.p),diagm(SPE.n)), C, Rho, Iteration
end
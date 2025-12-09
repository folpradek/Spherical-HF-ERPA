function BCS_Solver(Params::Parameters,Params_ref::Parameters)
    # Convergence precision epsilon ...
    epsilon = 1e-7

    # Make single-particle orbitals - NuHamil ordering ...
    Orb = Make_Orbitals(Params.Calc.A,Params.Calc.Z,Params.Int.Nmax)

    # Load 1-body kinetic operator ...
    T = T1B(Params.Int.Nmax,Orb,Params.Int.hw)

    # 2-body NN interaction & Orbitals ...
    @time VNN, Orb_NN = V2B_Read(Params,Orb)

    # 3-body NNN interaction & Orbitals ...
    @time VNNN, Orb_NNN = V3B_NO2B_Read(Params,Orb)

    # Solve HF-BCS equations ...
    @time E_MF, E_BCS, Lambda, SPE, SQE, C, U, V, Rho, Kappa, h, Delta = HF_BCS_Solve(Params,Params_ref,Orb,Orb_NN,Orb_NNN,T,VNN,VNNN,epsilon)

    # Calculate the total BCS ground-state kinetic energy ...
    T_BCS = Kinetic_Energy(Params,Rho,Orb,T)

    # Particle number fluctuation calculation ...
    dA = BCS_Particle_Number_Dispersion(Params,U,V,Orb)

    # Calculation summary ...
    BCS_Summary(Params,Params_ref,E_MF,E_BCS,T_BCS,Lambda,dA,epsilon)

    # Evaluate BCS charge radii & radial densities ...
    Summary_File = "IO/" * Params.Calc.Path * "/BCS/BCS_Summary.dat"
    Densities_File = "IO/" * Params.Calc.Path * "/BCS/Densities/BCS_Radial_Densities.dat"
    OBDM_Export(Params,Summary_File,Densities_File,pnMatrix(C.p * Rho.p * C.p', C.n * Rho.n * C.n'),C,Orb)

    # Export of single-quasiparticle energies, amplitudes U & V & also possibly radial densities ...
    BCS_SQS_Summary(Params,SPE,SQE,U,V,Orb)

    return
end

function HF_BCS_Solve(Params::Parameters,Params_ref::Parameters,Orb::Vector{NOrb},Orb_NN::NNOrb,Orb_NNN::NNNOrb,T::Matrix{Float64},VNN::NNInt,VNNN::Array{Vector{Vector{Float32}},4},epsilon::Float64)
    # Read calculation parameters ...
    Z_target, N_target = Params.Calc.Z, Params.Calc.A - Params.Calc.Z
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1) * (N_max + 2), 2)

    # Setup local iteration variables ...
    Iteration_BCS, Iteration_max = 0, 500
    
    # Preallocate some arrays ...
        # BCS amplitudes U & V vectors ...
    V = pnVector(zeros(Float64,a_max), zeros(Float64,a_max))
    U = pnVector(zeros(Float64,a_max), zeros(Float64,a_max))
        # BCS SQEs ...
    SQE = pnVector(zeros(Float64,a_max), zeros(Float64,a_max))
        # The pairing gap vector ...
    Delta = pnVector(0.5 .* ones(Float64,a_max), 0.5 .* ones(Float64,a_max))

    # Solve the HF-BCS approximation  ...
    println("\nStarting iteration of HF-BCS with NO2B NN+NNN interaction ...\n")

    # Define local function for reference HF calculation ...
    function HF_BCS_HF_Solve()
        # Setup single-particle orbitals for HF - NuHamil ordering ...
        Orb_HF = Make_Orbitals(Params_ref.Calc.A,Params_ref.Calc.Z,Params_ref.Int.Nmax)

        # Read 2-body NN interaction with CM correction for reference nucleus ...
            # Note that further self-consistent HF-BCS iterations are considered
            # with respect to the target nucleus ... the target values of A ...
        VNN_ref, Orb_NN_ref = V2B_Read(Params_ref,Orb_HF)

        # Call the HF Solver for reference closed-shell nucleus ...
        println("\nSolving the HF equations for reference closed-shell system ...")
        println("Reference nucleus:     A = " * string(Params_ref.Calc.A) * ",     Z = " * string(Params_ref.Calc.Z) * "\n")
        @time SPE, C, Rho, h, Iteration_HF = HF_Solve(Params_ref,Orb_HF,Orb_NN_ref,Orb_NNN,T,VNN_ref,VNNN,epsilon)

        # Calculate the HF mean-field energy ... for comparison
        println("\nCalculating the HF mean-field ground-state energy ... sanity check ...")
        @time E_HF = HF_Energy(Params_ref,Rho,Orb_HF,Orb_NN_ref,Orb_NNN,T,VNN_ref,VNNN)

        # Make residual density-dependent NN interaction in canonical HF basis... J = 0 - s-wave only ...
        println("\nMaking density-depenent residual NN interaction ... s-wave channel (J = 0) ...")
        @time VNN_res, Orb_NN_res = BCS_V2B_Res(Params_ref,Orb,Orb_NN_ref,Orb_NNN,VNN_ref,VNNN,C,Rho)

        # Drop VNN_ref and force Garbace Collection ...
        VNN_ref = nothing
        GC.gc()

        # Determine the initial value of chemical potential from the HF calculation ...
        Lambda = HF_BCS_Initialize_Chemical_Potential(Params,SPE,Orb_HF)

        return SPE, C, Rho, h, Lambda, VNN_res, Orb_NN_res
    end

    # Solve the HF equations for the reference closed-shell nucleus ...
    @time SPE, C, Rho, h, Lambda, VNN_res, Orb_NN_res = HF_BCS_HF_Solve()

    # Solve BCS equations ...
    println("\nInitializing the HF-BCS approximation ...")

    # Setup particle numbers Z & N ... exact from the HF iteration ...
    Z, N = Z_target, N_target

    # Setup particle number differences dZ & dN ...
    dZ, dN = abs(Z_target - Z) + 0.1, abs(N_target - N) + 0.1

    # Iteratively solve the BCS equations ...
    println("\nStarting iteration of BCS equations ...\n")
    while ((abs(dZ) > epsilon) || (abs(dN) > epsilon)) && (Iteration_BCS < Iteration_max)

        # Initial BCS iteration ... based on the initial pairing gap Delta ...
        if Iteration_BCS == 0
            # Evaluate the BCS SQEs ...
            SQE = BCS_Allocate_SQE(Params,SPE,Lambda,Delta)

            # Evaluate the BCS amplitudes U & V ...
            U, V = BCS_Allocate_Amplitudes(Params,Lambda,SPE,Delta,Orb)

            # Evaluate <Z> & <N> from BCS amplitudes V ...
            Z, N = BCS_Particle_Number(Params,V,Orb)

            # Evaluate the particle number deviations dZ & dN ...
            dZ, dN = Z_target - Z, N_target - N

            # Evaluate the chemical potential Lambda ...
            Lambda = BCS_Lambda(Params,SPE,Delta,Lambda,pnFloat(dZ,dN),Orb)

        # Regular BCS iteration ... starting with gap equation & input values of V & U amplitudes ...
        else
            # Calculate the pairing gap Delta ...
            Delta = BCS_Allocate_Delta(Params,U,V,Orb,Orb_NN_res,VNN_res)

            # Evaluate the BCS SQEs ...
            SQE = BCS_Allocate_SQE(Params,SPE,Lambda,Delta)

            # Evaluate the BCS amplitudes U & V ...
            U, V = BCS_Allocate_Amplitudes(Params,Lambda,SPE,Delta,Orb)

            # Evaluate <Z> & <N> from BCS amplitudes V ...
            Z, N = BCS_Particle_Number(Params,V,Orb)

            # Evaluate the chemical potential Lambda ...
            Lambda = BCS_Lambda(Params,SPE,Delta,Lambda,pnFloat(dZ,dN),Orb)
        
        end

        # Evaluate the BCS iteration ...
        Iteration_BCS += 1
        dZ, dN = Z_target - Z, N_target - N
        println("\nBCS iteration number:   " * string(Iteration_BCS) * "   dZ = " * string(round(dZ, sigdigits=8))* "   &   dN = " * string(round(dN, sigdigits=8)))
        println("\tCurrent particle numbers     ...     Z = " * string(Z) * ", N = " * string(N))
        println("\tCurrent chemical potentials        ...     pLambda = " * string(round(Lambda.p, digits = 6)) * " MeV, nLambda = " * string(round(Lambda.n, digits = 6)) * " MeV")

    end

    # Start self-consistent iteration of BCS equations ...
    if Params.Calc.Pairing.ScBCS == true
        println("\nStarting self-consistent iteration of BCS equations ...\n")

        # Initialize mean-field iteration variables ...
        Iteration_MF, dE = 0, 1.0

        # Initialize vector for old SPEs ...
        SPE_old = pnVector(zeros(Float64,a_max), zeros(Float64,a_max))

        # Start the HF-BCS self-consistent loop iteration ...
        while (dE > epsilon) && (Iteration_MF < Iteration_max)

            # Allocate density operators Rho & Kappa ...
            Rho, Kappa = BCS_Density_Operator(Params,U,V)

            # Transform Rho & Kappa to the reference LHO basis ...
            Rho = pnMatrix(C.p * Rho.p * C.p', C.n * Rho.n * C.n')
            Kappa = pnMatrix(C.p * Kappa.p * C.p', C.n * Kappa.n * C.n')

            # Evaluate new HF mean-field Hamiltonian h ...
            h = HF_BCS_Allocate(Params,Rho,Kappa,Orb,Orb_NN,Orb_NNN,T,VNN,VNNN)

            # Diagonalize the HF Hamiltonians ...
            pSPE, pC = eigen(Symmetric(h.p), sortby=+)
            nSPE, nC = eigen(Symmetric(h.n), sortby=+)

            # Reorder HF orbitals & SPEs ...
            C, SPE = HF_Orbital_Ordering(Orb,a_max,pnMatrix(pC,nC),pnVector(pSPE,nSPE))

            # Make residual density-dependent NN interaction in canonical HF basis... J = 0 - s-wave only ...
                # To avoid spam, Terminal Output is supressed for this call ...
            Out = "/dev/null"
            if Sys.iswindows()
                Out = "NUL"
            end
            open(Out, "w") do devnull_io
                redirect_stdout(devnull_io) do
                    redirect_stderr(devnull_io) do
                        VNN_res, Orb_NN_res = BCS_V2B_Res(Params,Orb,Orb_NN,Orb_NNN,VNN,VNNN,C,Rho)
                        return VNN_res, Orb_NN_res
                    end
                end
            end

            # Iterate the BCS equations ...
            Iteration_BCS = 0
            dZ, dN = 0.1, 0.1

            println("\n\tStarting iteration of BCS equations ...\n")
            while ((dZ > epsilon) || (dN > epsilon)) && (Iteration_BCS < Iteration_max)
                
                # Initial BCS iteration ... twofold for 2 initial guesses on Chemical Potentials ...
                if Iteration_BCS == 0
                    # Evaluate the BCS SQEs ...
                    SQE = BCS_Allocate_SQE(Params,SPE,Lambda,Delta)

                    # Evaluate the BCS amplitudes U & V ...
                    U, V = BCS_Allocate_Amplitudes(Params,Lambda,SPE,Delta,Orb)

                    # Evaluate <Z> & <N> from BCS amplitudes V ...
                    Z, N = BCS_Particle_Number(Params,V,Orb)

                    # Evaluate the particle number deviations dZ & dN ...
                    dZ, dN = abs(Z_target - Z), abs(N_target - N)

                    # Evaluate the chemical potential Lambda ...
                    Lambda = BCS_Lambda(Params,SPE,Delta,Lambda,pnFloat(dZ,dN),Orb)

                # Regular BCS iteration ... starting with gap equation & input values of V & U amplitudes ...
                elseif Iteration_BCS > 0
                    # Calculate the pairing gap Delta ...
                    Delta = BCS_Allocate_Delta(Params,U,V,Orb,Orb_NN_res,VNN_res)

                    # Evaluate the BCS SQEs ...
                    SQE = BCS_Allocate_SQE(Params,SPE,Lambda,Delta)

                    # Evaluate the BCS amplitudes U & V ...
                    U, V = BCS_Allocate_Amplitudes(Params,Lambda,SPE,Delta,Orb)

                    # Evaluate <Z> & <N> from BCS amplitudes V ...
                    Z, N = BCS_Particle_Number(Params,V,Orb)

                    # Evaluate the chemical potential Lambda ...
                    Lambda = BCS_Lambda(Params,SPE,Delta,Lambda,pnFloat(dZ,dN),Orb)

                end

                # Evaluate the BCS iteration ...
                Iteration_BCS += 1
                dZ, dN = abs(Params.Calc.Z - Z), abs(Params.Calc.A - Params.Calc.Z - N)
                println("\n\tBCS iteration number:   " * string(Iteration_BCS) * "   dZ = " * string(round(dZ, sigdigits=8))* "   &   dN = " * string(round(dN, sigdigits=8)))
                println("\t\tCurrent particle numbers     ...     Z = " * string(Z) * ", N = " * string(N))
                println("\t\tCurrent chemical potentials        ...     pLambda = " * string(round(Lambda.p, digits = 6)) * " MeV, nLambda = " * string(round(Lambda.n, digits = 6)) * " MeV")

            end

            # Evaluate mean-field iteration ...
            Iteration_MF += 1
            dE = (sum(abs.(SPE.p .- SPE_old.p )) + sum(abs.(SPE.n .- SPE_old.n))) / Float64(2 * a_max)
            println("\nMean-field iteration number:   " * string(Iteration_MF) * "   dE = " * string(round(dE, sigdigits=8)) * " MeV")
        
            # Save current SPEs ...
            SPE_old  = pnVector(deepcopy(SPE.p), deepcopy(SPE.n))
        
        end

    end

    if Iteration_BCS < Iteration_max
        println("\nBCS iteration with residual NN interaction has converged ...")
    elseif Iteration_BCS == Iteration_max
        println("\nBCS iteration with residual NN interaction has terminated ...")
        println("\n\t (!!!) CONVERGENCE WAS NOT REACHED (!!!) \n")
    end

    # Determine the BCS SQEs ...
    SQE = BCS_Allocate_SQE(Params,SPE,Lambda,Delta)

    # Determine the BCS amplitudes U & V ...
    U, V = BCS_Allocate_Amplitudes(Params,Lambda,SPE,Delta,Orb)

    # Determine the resulting density operators Rho & Kappa ...
    Rho, Kappa = BCS_Density_Operator(Params,U,V)

    # Express mean-field Hamiltonian h in the canonical HF basis ...
    h = pnMatrix(diagm(SPE.p),diagm(SPE.n))

    # Calculate the total HF mean-field + BCS pairing ground-state energy ...
    @time E_HF, E_BCS = BCS_Energy(Params,pnMatrix(C.p * Rho.p * C.p', C.n * Rho.n * C.n'),pnMatrix(C.p * Kappa.p * C.p', C.n * Kappa.n * C.n'),Orb,Orb_NN,Orb_NNN,T,VNN,VNNN)

    return E_HF, E_BCS, Lambda, SPE, SQE, C, U, V, Rho, Kappa, h, Delta
end
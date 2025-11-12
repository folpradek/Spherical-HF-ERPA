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
    @time E_MF, E_BCS, Lambda, SPE, SQE, C, U, V, Rho, Kappa, h, Delta, Iteration_BCS = HF_BCS_Solve(Params,Params_ref,Orb,Orb_NN,Orb_NNN,T,VNN,VNNN,epsilon)

    # Calculate the total BCS ground-state kinetic energy ...
    T_BCS = Kinetic_Energy(Params,Rho,Orb,T)

    # Particle number fluctuation calculation ...
    dA = BCS_Particle_Number_Dispersion(Params,U,V,Orb)

    # Calculation summary ...
    BCS_Summary(Params,Params_ref,E_MF,E_BCS,T_BCS,Lambda,dA,epsilon,Iteration_BCS)

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
    A, Z = Params.Calc.A, Params.Calc.Z
    hw = Params.Int.hw
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1) * (N_max + 2), 2)

    # Setup local iteration variables ...
    Iteration_max = 500
    Iteration_HF, Iteration_BCS = 0, 0
    
    # Preallocate arrays ...
        # Vectors for U & V BCS amplitudes ...
    pV, nV = zeros(Float64,a_max), zeros(Float64,a_max)
    pU, nU = zeros(Float64,a_max), zeros(Float64,a_max)
        # Vector for the pairing gap Delta ...
    pDelta, nDelta = 0.5 .* ones(Float64,a_max), 0.5 .* ones(Float64,a_max)

    # Solve the HF-BCS approximation  ...
    println("\nStarting iteration of HF-BCS with NO2B NN+NNN interaction ...\n")

    # Solve the HF equations for the reference closed-shell nucleus ...
    function HF_BCS_HF_Solve()
        # Setup single-particle orbitals for HF - NuHamil ordering ...
        Orb_HF = Make_Orbitals(Params_ref.Calc.A,Params_ref.Calc.Z,Params_ref.Int.Nmax)

        # Read 2-body NN interaction with CM correction for reference nucleus ...
            # Note that further self-consistent HF-BCS iterations are considered
            # with respect to the target nucleus ... the target values of A ...
        VNN_ref, Orb_NN_ref = V2B_Read(Params_ref,Orb)

        # Call the HF Solver for reference closed-shell nucleus ...
        println("\nSolving the HF equations for reference closed-shell system ...")
        println("Reference nucleus:     A = " * string(Params_ref.Calc.A) * ",     Z = " * string(Params_ref.Calc.Z) * "\n")
        @time SPE, C, Rho, h, Iteration_HF = HF_Solve(Params_ref,Orb_HF,Orb_NN_ref,Orb_NNN,T,VNN_ref,VNNN,epsilon)

        # Make residual density-dependent NN interaction in canonical HF basis... J = 0 - s-wave only ...
        println("\nMaking density-depenent residual NN interaction ... s-wave channel (J = 0) ...")
        @time VNN_res, Orb_NN_res = BCS_V2B_Res(Params,Orb,Orb_NN_ref,Orb_NNN,VNN_ref,VNNN,C,Rho)

        # Drop VNN_ref and force Garbace Collection ...
        VNN_ref = nothing
        GC.gc()

        return SPE, C, Rho, h, Iteration_HF, VNN_res, Orb_NN_res
    end

    @time SPE, C, Rho, h, Iteration_HF, VNN_res, Orb_NN_res = HF_BCS_HF_Solve()

    # Delete this ... old block ...
    #=
            # Read 2-body NN interaction with CM correction for reference nucleus ...
                # Note that further self-consistent HF-BCS iterations are considered
                # with respect to the target nucleus ... the target values of A ...
            @time VNN_ref, Orb_NN_ref = V2B_Read(Params_ref,Orb)

            # First solve the HF equations for reference closed-shell nucleus ...
            println("\nSolving the HF equations for reference closed-shell system ...")
            println("Reference nucleus:     A = " * string(Params_ref.Calc.A) * ",     Z = " * string(Params_ref.Calc.Z) * "\n")
            @time SPE, C, Rho, h, Iteration_HF = HF_Solve(Params_ref,Orb_HF,Orb_NN_ref,Orb_NNN,T,VNN_ref,VNNN,epsilon)

            # Make residual density-dependent NN interaction in canonical HF basis... J = 0 - s-wave only ...
            println("\nMaking density-depenent residual NN interaction ... s-wave channel (J = 0) ...")
            @time VNN_res, Orb_NN_res = BCS_V2B_Res(Params,Orb,Orb_NN_ref,Orb_NNN,VNN_ref,VNNN,C,Rho)
    =#

    # Solve BCS equations ...
    println("\nInitializing the HF-BCS approximation ...")

    # Initial guess on chemical potentials Lambda ...
    pLambda, nLambda = 0.2, 0.2

    # Setup particle numbers ...
    Z, Z_1, Z_2 = Z, 0, 0
    N, N_1, N_2 = (A - Z), 0, 0
    dZ, dN = abs(Params.Calc.Z - Params_ref.Calc.Z) + 0.01, abs(Params.Calc.A - Params.Calc.Z - Params_ref.Calc.A + Params_ref.Calc.Z) + 0.01

    # Iteratively solve the BCS equations ...
    println("\nStarting iteration of BCS equations ...\n")
    while ((dZ > epsilon) || (dN > epsilon)) && (Iteration_BCS < Iteration_max)
        # Initial BCS iteration ... twofold for 2 initial guesses on Chemical Potentials ...
        if Iteration_BCS == 0
            # Calculate U & V from initial guess ... from the initial guess
            @inbounds for a in 1:a_max
                pME = (SPE.p[a] - pLambda) / sqrt((SPE.p[a] - pLambda)^2 + pDelta[a]^2)
                pV[a] = sqrt(0.5 * (1.0 - pME))
                pU[a] = sqrt(0.5 * (1.0 + pME))

                nME = (SPE.n[a] - nLambda) / sqrt((SPE.n[a] - nLambda)^2 + nDelta[a]^2)
                nV[a] = sqrt(0.5 * (1.0 - nME))
                nU[a] = sqrt(0.5 * (1.0 + nME))
            end

            # Evaluate <Z> & <N> from BCS amplitudes V ...
            Z, N = BCS_Particle_Number(Params,pnVector(pV,nV),Orb)

            pLambda, nLambda = pLambda + 0.1 * (Params.Calc.Z - Z), nLambda + 0.1 * (Params.Calc.A - Params.Calc.Z - N)

        # Regular BCS iteration ... starting with gap equation & input values of V & U amplitudes ...
        elseif Iteration_BCS > 0
            # Calculate gap equation ...
            @inbounds for a in 1:a_max
                j_a = Orb[a].j
                pSum, nSum = 0.0, 0.0
                @inbounds for b in 1:a_max
                    j_b = Orb[b].j
                    ja_jb_hat = sqrt((Float64(j_b) + 1.0) / (Float64(j_a) + 1.0))
                    pME = - ja_jb_hat * V2B(a,a,b,b,0,1,VNN_res.pp,Orb,Orb_NN_res) * pU[b] * pV[b]
                    nME = - ja_jb_hat * V2B(a,a,b,b,0,1,VNN_res.nn,Orb,Orb_NN_res) * nU[b] * nV[b]
                    pSum += pME
                    nSum += nME
                end
                pDelta[a] = pSum
                nDelta[a] = nSum
            end

            # Recalculate U & V amplitudes ...
            @inbounds for a in 1:a_max
                pME = 0.5 * (SPE.p[a] - pLambda) / sqrt((SPE.p[a] - pLambda)^2 + pDelta[a]^2)
                pV[a] = sqrt((0.5 - pME))
                pU[a] = sqrt((0.5 + pME))

                nME = 0.5 * (SPE.n[a] - nLambda) / sqrt((SPE.n[a] - nLambda)^2 + nDelta[a]^2)
                nV[a] = sqrt((0.5 - nME))
                nU[a] = sqrt((0.5 + nME))
            end

            # Evaluate <Z> & <N> from BCS amplitudes V ...
            Z, N = BCS_Particle_Number(Params,pnVector(pV,nV),Orb)

            pLambda, nLambda = pLambda + 0.1 * (Params.Calc.Z - Z), nLambda + 0.1 * (Params.Calc.A - Params.Calc.Z - N)
        end

        Iteration_BCS += 1
        println("\nCurrent particle number values are ...")
        println("Z = " * string(Z))
        println("N = " * string(N))
        println("pLambda = " * string(pLambda))
        println("nLambda = " * string(nLambda))
        dZ, dN = abs(Params.Calc.Z - Z), abs(Params.Calc.A - Params.Calc.Z - N)

        println("BCS iteration number:   " * string(Iteration_BCS) * "   Proton number difference:   " * string(round(dZ, sigdigits=8))* "   &   Neutron number difference:   " * string(round(dN, sigdigits=8)))

    end

    # Start self-consistent iteration of BCS equations ...
    if Params.Calc.Pairing.ScBCS == true
        println("\nStarting self-consistent iteration of BCS equations ...\n")
        Iteration = 0
        dE = 1.0
        SPE_old = pnVector(zeros(Float64,a_max), zeros(Float64,a_max))
        while (dE > epsilon) && (Iteration < Iteration_max)

            # Allocate density operators Rho & Kappa ...
            Rho = BCS_Density_Operator(a_max,pnVector(pV,nV))
            Kappa = BCS_Pairing_Operator(a_max,pnVector(pU,nU),pnVector(pV,nV))

            # Transform Rho & Kappa to the reference LHO basis ...
            Rho = pnMatrix(C.p * Rho.p * C.p', C.n * Rho.n * C.n')
            Kappa = pnMatrix(C.p * Kappa.p * C.p', C.n * Kappa.n * C.n')

            # Evaluate new HF mean-field Hamiltonian
            h = HF_BCS_Allocate(Params,Rho,Kappa,Orb,Orb_NN,Orb_NNN,T,VNN,VNNN)

            # Diagonalize the HF Hamiltonians ...
            pSPE, pC = eigen(Symmetric(h.p), sortby=+)
            nSPE, nC = eigen(Symmetric(h.n), sortby=+)

            # Reorder HF orbitals & SPEs ...
            C, SPE = HF_Orbital_Ordering(Orb,a_max,pnMatrix(pC,nC),pnVector(pSPE,nSPE))

            # Make residual density-dependent NN interaction in canonical HF basis... J = 0 - s-wave only ...
                # For brevity, Terminal Output is supressed for this call ...
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

            println("\nStarting iteration of BCS equations ...\n")
            while ((dZ > epsilon) || (dN > epsilon)) && (Iteration_BCS < Iteration_max)
                # Initial BCS iteration ... twofold for 2 initial guesses on Chemical Potentials ...
                if Iteration_BCS == 0
                    # Calculate U & V from initial guess ... from the initial guess
                    @inbounds for a in 1:a_max
                        pME = (SPE.p[a] - pLambda) / sqrt((SPE.p[a] - pLambda)^2 + pDelta[a]^2)
                        pV[a] = sqrt(0.5 * (1.0 - pME))
                        pU[a] = sqrt(0.5 * (1.0 + pME))

                        nME = (SPE.n[a] - nLambda) / sqrt((SPE.n[a] - nLambda)^2 + nDelta[a]^2)
                        nV[a] = sqrt(0.5 * (1.0 - nME))
                        nU[a] = sqrt(0.5 * (1.0 + nME))
                    end

                    # Evaluate <Z> & <N> from BCS amplitudes V ...
                    Z, N = BCS_Particle_Number(Params,pnVector(pV,nV),Orb)

                    pLambda, nLambda = pLambda + 0.1 * (Params.Calc.Z - Z), nLambda + 0.1 * (Params.Calc.A - Params.Calc.Z - N)

                # Regular BCS iteration ... starting with gap equation & input values of V & U amplitudes ...
                elseif Iteration_BCS > 0
                    # Allocate the pairing gap Delta ...
                    pDelta, nDelta = BCS_Allocate_Delta(Params,pnVector(pU,nU),pnVector(pV,nV),Orb,Orb_NN_res,VNN_res)

                    # Recalculate U & V amplitudes ...
                    @inbounds for a in 1:a_max
                        pME = 0.5 * (SPE.p[a] - pLambda) / sqrt((SPE.p[a] - pLambda)^2 + pDelta[a]^2)
                        pV[a] = sqrt((0.5 - pME))
                        pU[a] = sqrt((0.5 + pME))

                        nME = 0.5 * (SPE.n[a] - nLambda) / sqrt((SPE.n[a] - nLambda)^2 + nDelta[a]^2)
                        nV[a] = sqrt((0.5 - nME))
                        nU[a] = sqrt((0.5 + nME))
                    end

                    # Evaluate <Z> & <N> from BCS amplitudes V ...
                    Z, N = BCS_Particle_Number(Params,pnVector(pV,nV),Orb)
                    
                    # Update the chemical potential ... simple quenching formula ... (e.g. see Suhonen Chapter 14)
                    pLambda, nLambda = pLambda + 0.1 * (Params.Calc.Z - Z), nLambda + 0.1 * (Params.Calc.A - Params.Calc.Z - N)
                end

                Iteration_BCS += 1
                println("\nCurrent particle number values are ...")
                println("Z = " * string(Z))
                println("N = " * string(N))
                println("pLambda = " * string(pLambda))
                println("nLambda = " * string(nLambda))
                dZ, dN = abs(Params.Calc.Z - Z), abs(Params.Calc.A - Params.Calc.Z - N)

                println("BCS iteration number:   " * string(Iteration_BCS) * "   Proton number difference:   " * string(round(dZ, sigdigits=8))* "   &   Neutron number difference:   " * string(round(dN, sigdigits=8)))

            end

            dE = (sum(abs.(SPE.p .- SPE_old.p )) + sum(abs.(SPE.n .- SPE_old.n))) / Float64(2 * a_max)
            SPE_old  = pnVector(deepcopy(SPE.p), deepcopy(SPE.n))

            Iteration += 1
            println("Iteration number:   " * string(Iteration) * "   Energy difference:   " * string(round(dE, sigdigits=8)) * " MeV")
        end

    end

    println("\nBCS iteration with residual NN interaction has converged ...")

    # Allocate resulting chemical potential Lambda ...
    Lambda = pnFloat(pLambda,nLambda)

    # Allocate BCS amplitudes U & V ...
    U, V = pnVector(pU,nU), pnVector(pV,nV)

    # Allocate the pairing gap Delta ...
    Delta = pnVector(pDelta,nDelta)

    # Determine the single-quasiparticle energies (SQE) ...
    SQE = BCS_SQE(a_max,SPE,Lambda,Delta)

    # Determine resulting density Rho ...
        # For Self-Consistent BCS ... determined from V amplitudes ...
    Rho = BCS_Density_Operator(a_max,V)

    # Determine resulting pairing tensor Kappa ...
    Kappa = BCS_Pairing_Operator(a_max,U,V)

    # Express mean-field Hamiltonian h in the canonical HF basis ...
    h = pnMatrix(diagm(SPE.p),diagm(SPE.n))

    # Calculate the HF mean-field energy ... Requires Rho expressed in the LHO basis ...
        # Note that E_HF != E_HF_ref ... due to the CMS correction! ...
    @time E_HF = HF_Energy(Params,pnMatrix(C.p * Rho.p * C.p', C.n * Rho.n * C.n'),Orb,Orb_NN,Orb_NNN,T,VNN,VNNN)

    # Calculate BCS ground-state pairing energy ...
    @time E_BCS = BCS_Energy(Params,Kappa,Orb,Orb_NN_res,VNN_res)

    return E_HF, E_BCS, Lambda, SPE, SQE, C, U, V, Rho, Kappa, h, Delta, Iteration_BCS
end
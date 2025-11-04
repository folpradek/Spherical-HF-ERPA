function BCS_Solver(Params::Parameters,Params_ref::Parameters)
    # Make single-particle orbitals - NuHamil ordering ...
    Orb = Make_Orbitals(Params.Calc.A,Params.Calc.Z,Params.Int.Nmax)

    # Load 1-body kinetic operator ...
    T = T1B(Params.Int.Nmax,Orb,Params.Int.hw)

    # 2-body NN interaction & Orbitals ...
    @time VNN, Orb_NN = V2B_Read(Params,Orb)

    # 3-body NNN interaction & Orbitals ...
    @time VNNN, Orb_NNN = V3B_NO2B_Read(Params,Orb)

    # Solve HF-BCS equations ...
    @time E_MF, E_BCS, lambda, SPE, SQE, C, U, V, Rho, Kappa, h, Delta, Iteration_BCS = HF_BCS_Solve(Params,Params_ref,Orb,Orb_NN,Orb_NNN,T,VNN,VNNN,1e-7)

    # Inspection print of amomal densities ...
    #display(C.p * Kappa.p * C.p')
    #display(C.n * Kappa.n * C.n')

    # Particle number fluctuation calculation ...
    dA = BCS_dN(Params,U,V,Orb)

    # Calculation summary ...
    BCS_Summary(Params,Params_ref,E_MF,E_BCS,lambda,dA,1e-7,Iteration_BCS)

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

    # Setup single-particle orbitals for HF - NuHamil ordering ...
    Orb_HF = Make_Orbitals(Params_ref.Calc.A,Params_ref.Calc.Z,Params_ref.Int.Nmax)

    # Setup local iteration variables ...
    Iteration_max = 500
    Iteration_HF, Iteration_BCS = 0, 0
    delta = 1.0
    
    # Preallocate arrays ...
        # Vectors for U & V BCS amplitudes ...
    pV, nV = zeros(Float64,a_max), zeros(Float64,a_max)
    pU, nU = zeros(Float64,a_max), zeros(Float64,a_max)
        # Vector for the pairing gap Delta ...
    pDelta, nDelta = 0.5 .* ones(Float64,a_max), 0.5 .* ones(Float64,a_max)

    # Solve the HF-BCS approximation  ...
    println("\nStarting iteration of HF-BCS with NO2B NN+NNN interaction ...\n")

    # First solve the HF equations for reference closed-shell nucleus ...
    println("\nSolving the HF equations for reference closed-shell system ...")
    println("Reference nucleus:     A = " * string(Params_ref.Calc.A) * ",     Z = " * string(Params_ref.Calc.Z) * "\n")
    @time SPE, C, Rho, h, Iteration_HF = HF_Solve(Params_ref,Orb_HF,Orb_NN,Orb_NNN,T,VNN,VNNN,epsilon)

    # Make residual density-dependent NN interaction in canonical HF basis... J = 0 - s-wave only ...
    println("\nMaking density-depenent residual NN interaction ... s-wave channel (J = 0) ...")
    @time VNN_Res, Orb_NN_Res = BCS_V2B_Res(Params,Orb,Orb_NN,Orb_NNN,VNN,VNNN,C,Rho)

    # Solve BCS equations ...
    println("\nInitializing the HF-BCS approximation ...")

    # Initial guess on chemical potentials lambda ...
    plambda, nlambda = 0.2, 0.2

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
                pME = (SPE.p[a] - plambda) / sqrt((SPE.p[a] - plambda)^2 + pDelta[a]^2)
                pV[a] = sqrt(0.5 * (1.0 - pME))
                pU[a] = sqrt(0.5 * (1.0 + pME))

                nME = (SPE.n[a] - nlambda) / sqrt((SPE.n[a] - nlambda)^2 + nDelta[a]^2)
                nV[a] = sqrt(0.5 * (1.0 - nME))
                nU[a] = sqrt(0.5 * (1.0 + nME))
            end

            # Calculate the average particle numbers for BCS amplitudes ...
            Z, N = 0, 0
            @inbounds for a in 1:a_max
                z = Float64(Orb[a].j + 1) * pV[a]^2
                Z += z

                n = Float64(Orb[a].j + 1) * nV[a]^2
                N += n
            end

            plambda, nlambda = plambda + 0.1 * (Params.Calc.Z - Z), nlambda + 0.1 * (Params.Calc.A - Params.Calc.Z - N)

        # Regular BCS iteration ... starting with gap equation & input values of V & U amplitudes ...
        elseif Iteration_BCS > 0
            # Calculate gap equation ...
            @inbounds for a in 1:a_max
                j_a = Orb[a].j
                pSum, nSum = 0.0, 0.0
                @inbounds for b in 1:a_max
                    j_b = Orb[b].j
                    ja_jb_hat = sqrt((Float64(j_b) + 1.0) / (Float64(j_a) + 1.0))
                    pME = - ja_jb_hat * V2B(a,a,b,b,0,1,VNN_Res.pp,Orb,Orb_NN_Res) * pU[b] * pV[b]
                    nME = - ja_jb_hat * V2B(a,a,b,b,0,1,VNN_Res.nn,Orb,Orb_NN_Res) * nU[b] * nV[b]
                    pSum += pME
                    nSum += nME
                end
                pDelta[a] = pSum
                nDelta[a] = nSum
            end

            # Recalculate U & V amplitudes ...
            @inbounds for a in 1:a_max
                pME = 0.5 * (SPE.p[a] - plambda) / sqrt((SPE.p[a] - plambda)^2 + pDelta[a]^2)
                pV[a] = sqrt((0.5 - pME))
                pU[a] = sqrt((0.5 + pME))

                nME = 0.5 * (SPE.n[a] - nlambda) / sqrt((SPE.n[a] - nlambda)^2 + nDelta[a]^2)
                nV[a] = sqrt((0.5 - nME))
                nU[a] = sqrt((0.5 + nME))
            end

            # Evaluate <Z> & <N> from BCS amplitudes U & V ...
            Z, N = 0, 0
            @inbounds for a in 1:a_max
                z = Float64(Orb[a].j + 1) * pV[a]^2
                Z += z
                n = Float64(Orb[a].j + 1) * nV[a]^2
                N += n
            end

            plambda, nlambda = plambda + 0.1 * (Params.Calc.Z - Z), nlambda + 0.1 * (Params.Calc.A - Params.Calc.Z - N)
        end

        Iteration_BCS += 1
        println("\nCurrently particle number values are ...")
        println("Z = " * string(Z))
        println("N = " * string(N))
        println("pLambda = " * string(plambda))
        println("nLambda = " * string(nlambda))
        dZ, dN = abs(Params.Calc.Z - Z), abs(Params.Calc.A - Params.Calc.Z - N)

        println("BCS iteration number:   " * string(Iteration_BCS) * "   Proton number difference:   " * string(round(dZ, sigdigits=8))* "   &   Neutron number difference:   " * string(round(dN, sigdigits=8)))

    end

    SCBCS = false

    # Start self-consistent iteration of BCS equations ...
    if SCBCS == true
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

            Out = "/dev/null"
            if Sys.iswindows()
                Out = "NUL"
            end

            # Make residual density-dependent NN interaction in canonical HF basis... J = 0 - s-wave only ...
            # For brevity, Terminal Output is supressed for this call ...
            open(Out, "w") do devnull_io
                redirect_stdout(devnull_io) do
                    redirect_stderr(devnull_io) do
                        VNN_Res, Orb_NN_Res = BCS_V2B_Res(Params,Orb,Orb_NN,Orb_NNN,VNN,VNNN,C,Rho)
                        return VNN_Res, Orb_NN_Res
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
                        pME = (SPE.p[a] - plambda) / sqrt((SPE.p[a] - plambda)^2 + pDelta[a]^2)
                        pV[a] = sqrt(0.5 * (1.0 - pME))
                        pU[a] = sqrt(0.5 * (1.0 + pME))

                        nME = (SPE.n[a] - nlambda) / sqrt((SPE.n[a] - nlambda)^2 + nDelta[a]^2)
                        nV[a] = sqrt(0.5 * (1.0 - nME))
                        nU[a] = sqrt(0.5 * (1.0 + nME))
                    end

                    # Calculate the average particle numbers for BCS amplitudes ...
                    Z, N = 0.0, 0.0
                    @inbounds for a in 1:a_max
                        z = Float64(Orb[a].j + 1) * pV[a]^2
                        Z += z

                        n = Float64(Orb[a].j + 1) * nV[a]^2
                        N += n
                    end

                    plambda, nlambda = plambda + 0.1 * (Params.Calc.Z - Z), nlambda + 0.1 * (Params.Calc.A - Params.Calc.Z - N)

                # Regular BCS iteration ... starting with gap equation & input values of V & U amplitudes ...
                elseif Iteration_BCS > 0
                    # Calculate gap equation ...
                    @inbounds for a in 1:a_max
                        j_a = Orb[a].j
                        pSum, nSum = 0.0, 0.0
                        @inbounds for b in 1:a_max
                            j_b = Orb[b].j
                            ja_jb_hat = sqrt((Float64(j_b) + 1.0) / (Float64(j_a) + 1.0))
                            pME = - ja_jb_hat * V2B(a,a,b,b,0,1,VNN_Res.pp,Orb,Orb_NN_Res) * pU[b] * pV[b]
                            nME = - ja_jb_hat * V2B(a,a,b,b,0,1,VNN_Res.nn,Orb,Orb_NN_Res) * nU[b] * nV[b]
                            pSum += pME
                            nSum += nME
                        end
                        pDelta[a] = pSum
                        nDelta[a] = nSum
                    end

                    # Recalculate U & V amplitudes ...
                    @inbounds for a in 1:a_max
                        pME = 0.5 * (SPE.p[a] - plambda) / sqrt((SPE.p[a] - plambda)^2 + pDelta[a]^2)
                        pV[a] = sqrt((0.5 - pME))
                        pU[a] = sqrt((0.5 + pME))

                        nME = 0.5 * (SPE.n[a] - nlambda) / sqrt((SPE.n[a] - nlambda)^2 + nDelta[a]^2)
                        nV[a] = sqrt((0.5 - nME))
                        nU[a] = sqrt((0.5 + nME))
                    end

                    # Evaluate <Z> & <N> from BCS amplitudes U & V ...
                    Z, N = 0.0, 0.0
                    @inbounds for a in 1:a_max
                        z = Float64(Orb[a].j + 1) * pV[a]^2
                        Z += z
                        n = Float64(Orb[a].j + 1) * nV[a]^2
                        N += n
                    end

                    plambda, nlambda = plambda + 0.1 * (Params.Calc.Z - Z), nlambda + 0.1 * (Params.Calc.A - Params.Calc.Z - N)
                end

                Iteration_BCS += 1
                println("\nCurrently particle number values are ...")
                println("Z = " * string(Z))
                println("N = " * string(N))
                println("pLambda = " * string(plambda))
                println("nLambda = " * string(nlambda))
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

    # Allocate resulting chemical potential lambda ...
    lambda = pnFloat(plambda,nlambda)

    # Allocate BCS amplitudes U & V ...
    U, V = pnVector(pU,nU), pnVector(pV,nV)

    # Allocate the pairing gap Delta ...
    Delta = pnVector(pDelta,nDelta)

    # Determine the single-quasiparticle energies (SQE) ...
    SQE = BCS_SQE(a_max,SPE,lambda,Delta)

    # Determine resulting density Rho ...
        # For Self-Consistent BCS ... determined from V amplitudes ...
    Rho = BCS_Density_Operator(a_max,V)

        #=
    if SCBCS == true
        Rho = BCS_Density_Operator(a_max,V)
        # For one-time BCS calculation ... the reference core is used ...
        # Rho transformed to the canonical basis ...
    elseif SCBCS == false
        Rho = pnMatrix(C.p' * Rho.p * C.p, C.n' * Rho.n * C.n)
    end
    =#

    # Determine resulting pairing tensor Kappa ...
    Kappa = BCS_Pairing_Operator(a_max,U,V)

    # Express mean-field Hamiltonian h in the canonical HF basis ...
    h = pnMatrix(diagm(SPE.p),diagm(SPE.n))

    # Calculate the HF mean-field energy ... Requires Rho expressed in the LHO basis ...
        # Note that E_HF != E_HF_ref ... due to the CMS correction! ...
    @time E_HF = HF_Energy(Params,pnMatrix(C.p * Rho.p * C.p', C.n * Rho.n * C.n'),Orb,Orb_NN,Orb_NNN,T,VNN,VNNN)

    # Calculate BCS ground-state pairing energy ...
    @time E_BCS = BCS_Energy(Params,Kappa,Orb,Orb_NN_Res,VNN_Res)

    return E_HF, E_BCS, lambda, SPE, SQE, C, U, V, Rho, Kappa, h, Delta, Iteration_BCS
end

function HF_BCS_Allocate(Params::Parameters,Rho::pnMatrix,Kappa::pnMatrix,Orb::Vector{NOrb},Orb_NN::NNOrb,Orb_NNN::NNNOrb,T::Matrix{Float64},VNN::NNInt,VNNN::Array{Vector{Vector{Float32}},4})
    # Read parameters ...
    A = Params.Calc.A
    N_max = Params.Calc.Nmax
    N_2max = Params.Calc.N2max
    N_3max = Params.Calc.N3max
    CMS = Params.Calc.CMS
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Allocate new HF Hamiltonian matrices ...
    pH, nH = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)

    # Allocate the single-particle field H ...
    @inbounds Threads.@threads for a in 1:a_max
        n_a = Orb[a].n
        l_a = Orb[a].l
        j_a = Orb[a].j
        ja_hat = 1.0 / (Float64(j_a) + 1.0)
        @inbounds for d in 1:a
            l_d = Orb[d].l
            j_d = Orb[d].j
            if l_a == l_d && j_a == j_d
                n_d = Orb[d].n
                pSum = 0.0
                nSum = 0.0
                @inbounds for b in 1:a_max
                    n_b = Orb[b].n
                    l_b = Orb[b].l
                    if (2*(n_a + n_b) + l_a + l_b) <= N_2max
                        j_b = Orb[b].j
                        @inbounds for e in 1:a_max
                            n_e = Orb[e].n
                            l_e = Orb[e].l
                            j_e = Orb[e].j

                            # Normal Density-dependent part ...
                            if l_b == l_e && j_b == j_e && (2*(n_d + n_e) + l_d + l_e) <= N_2max
                                pRho_be, nRho_be = Rho.p[b,e], Rho.n[b,e]
                                @inbounds for J = div(abs(j_d - j_e),2):div(j_d + j_e,2)

                                    # 2-body NN interaction part ...
                                    if rem(l_a + l_b, 2) == rem(l_d + l_e, 2)
                                        J_ja_hat = Float64(2*J + 1) * ja_hat
                                        @views pSum += J_ja_hat * V2B(a,b,d,e,J,1,VNN.pp,Orb,Orb_NN) * pRho_be
                                        @views pSum += J_ja_hat * V2B(a,b,d,e,J,0,VNN.pn,Orb,Orb_NN) * nRho_be
                                        @views nSum += J_ja_hat * V2B(a,b,d,e,J,1,VNN.nn,Orb,Orb_NN) * nRho_be
                                        @views nSum += J_ja_hat * V2B(b,a,e,d,J,0,VNN.pn,Orb,Orb_NN) * pRho_be

                                    end

                                    # 3-body NNN interaction part ...
                                    @inbounds for c in 1:a_max
                                        n_c = Orb[c].n
                                        l_c = Orb[c].l
                                        if (2*(n_a + n_b + n_c) + l_a + l_b + l_c) <= N_3max
                                            P = rem(l_a + l_b + l_c, 2) + 1
                                            j_c = Orb[c].j
                                            @inbounds for f in 1:a_max
                                                n_f = Orb[f].n
                                                l_f = Orb[f].l
                                                if (2*(n_d + n_e + n_f) + l_d + l_e + l_f) <= N_3max && l_c == l_f && P == (rem(l_d + l_e + l_f,2) + 1)
                                                    j_f = Orb[f].j
                                                    if j_c == j_f
                                                        pRho_cf = Rho.p[c,f]
                                                        nRho_cf = Rho.n[c,f]

                                                        ME001 = V3B_NO2B(a,b,c,0,d,e,f,0,J,1,P,VNNN,Orb,Orb_NNN)
                                                        ME101 = V3B_NO2B(a,b,c,1,d,e,f,0,J,1,P,VNNN,Orb,Orb_NNN)
                                                        ME011 = V3B_NO2B(a,b,c,0,d,e,f,1,J,1,P,VNNN,Orb,Orb_NNN)
                                                        ME111 = V3B_NO2B(a,b,c,1,d,e,f,1,J,1,P,VNNN,Orb,Orb_NNN)
                                                        ME113 = V3B_NO2B(a,b,c,1,d,e,f,1,J,3,P,VNNN,Orb,Orb_NNN)

                                                        pSum += ja_hat * (0.5*ME113*pRho_be*pRho_cf + 0.25 * (ME001 +
                                                                sqrt(1.0/3.0)*ME101 + sqrt(1.0/3.0)*ME011 + 1.0/3.0*ME111 + 2.0/3.0*ME113)*nRho_be*nRho_cf +
                                                                1.0/3.0 * (2*ME111 + ME113)*pRho_be*nRho_cf)

                                                        nSum += ja_hat * (0.5*ME113*nRho_be*nRho_cf + 0.25 * (ME001 +
                                                                sqrt(1.0/3.0)*ME101 + sqrt(1.0/3.0)*ME011 + 1.0/3.0*ME111 + 2.0/3.0*ME113)*pRho_be*pRho_cf +
                                                                1.0/3.0 * (2*ME111 + ME113)*nRho_be*pRho_cf)

                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end


                            # Anomal Density-dependent part ... Only 3-body NNN interaction part ...
                            if (2*(n_d + n_e) + l_d + l_e) <= N_2max
                                jb_je_hat = sqrt(Float64((j_b + 1) * (j_e + 1)))
                                @inbounds for c in 1:a_max
                                    n_c = Orb[c].n
                                    l_c = Orb[c].l
                                    if l_c == l_b && (2*(n_a + n_b + n_c) + l_a + l_b + l_c) <= N_3max
                                        P = rem(l_a + l_b + l_c, 2) + 1
                                        j_c = Orb[c].j
                                        if j_c == j_b
                                            pKappa_cb, nKappa_cb = Kappa.p[c,b], Kappa.n[c,b]
                                            @inbounds for f in 1:a_max
                                                n_f = Orb[f].n
                                                l_f = Orb[f].l
                                                if (2*(n_d + n_e + n_f) + l_d + l_e + l_f) <= N_3max && P == (rem(l_d + l_e + l_f,2) + 1) && l_e == l_f
                                                    j_f = Orb[f].j
                                                    if j_e == j_f
                                                        pKappa_ef, nKappa_ef = Kappa.p[e,f], Kappa.n[e,f]

                                                        ME111 = V3B_NO2B(b,c,a,1,e,f,d,1,0,1,P,VNNN,Orb,Orb_NNN)
                                                        ME113 = V3B_NO2B(b,c,a,1,e,f,d,1,0,3,P,VNNN,Orb,Orb_NNN)

                                                        pSum += -0.5 * jb_je_hat * ja_hat * (ME113 * pKappa_cb * pKappa_ef +
                                                                1.0 / 3.0 * (2.0 * ME111 + ME113) * nKappa_cb * nKappa_ef)
                                                        nSum += -0.5 * jb_je_hat * ja_hat * (ME113 * nKappa_cb * nKappa_ef +
                                                                1.0 / 3.0 * (2.0 * ME111 + ME113) * pKappa_cb * pKappa_ef)

                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end

                # 1-body kinetic energy & inclusion of Center-of-Mass motion (CM) correction ...
                if CMS == "CMS1+2B"
                    pH[a,d] = pSum + T[a,d] * (1.0 - 1.0 / A)
                    nH[a,d] = nSum + T[a,d] * (1.0 - 1.0 / A)
                    pH[d,a] = pH[a,d]
                    nH[d,a] = nH[a,d]
                elseif CMS == "CMS2B"
                    pH[a,d] = pSum
                    nH[a,d] = nSum
                    pH[d,a] = pH[a,d]
                    nH[d,a] = nH[a,d]
                else
                    pH[a,d] = pSum + T[a,d]
                    nH[a,d] = nSum + T[a,d]
                    pH[d,a] = pH[a,d]
                    nH[d,a] = nH[a,d]
                end

            end
        end
    end

    return pnMatrix(pH,nH)
end

function BCS_Density_Operator(a_max::Int64,V::pnVector)
    # Initialize density matrices ...
    pRho, nRho = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)

    # Compute the 1-body density matrix ... diagonal in canonical basis ...
    @inbounds for a in 1:a_max
        pRho[a,a] = V.p[a]^2
        nRho[a,a] = V.n[a]^2
    end

    return pnMatrix(pRho,nRho)
end

function BCS_Pairing_Operator(a_max::Int64,U::pnVector,V::pnVector)
    # Initialize pairing tensors ...
    pKappa, nKappa = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)

    # Compute the pairing tensor ... diagonal in canonical basis ...
    @inbounds for a in 1:a_max
        pKappa[a,a] = V.p[a] * U.p[a]
        nKappa[a,a] = V.n[a] * U.n[a]
    end

    return pnMatrix(pKappa,nKappa)
end

function BCS_Energy(Params::Parameters,Kappa::pnMatrix,Orb::Vector{NOrb},Orb_NN::NNOrb,VNN::NNInt)
   # Read parameters ...
    N_max = Params.Calc.Nmax
    N_2max = Params.Calc.N2max
    a_max = div((N_max + 1)*(N_max + 2),2)

    #Calculate the BCS ground-state energy ...
    println("\nCalculating total mean-field + BCS ground-state energy ...")

    E_BCS_partial = Threads.Atomic{Float64}[Threads.Atomic{Float64}(0.0) for _ in 1:Threads.nthreads()]

    @inbounds Threads.@threads for a = 1:a_max
        thread_id = Threads.threadid()
        n_a = Orb[a].n
        l_a = Orb[a].l
        j_a = Orb[a].j
        @inbounds for b = 1:a_max
            n_b = Orb[b].n
            l_b = Orb[b].l
            j_b = Orb[b].j
            if ((2*(n_a + n_b) + l_a + l_b) <= N_2max) && (l_a == l_b) && (j_a == j_b)
                @inbounds for d = 1:a_max
                    l_d = Orb[d].l
                    j_d = Orb[d].j
                    n_d = Orb[d].n
                    ja_jd_hat = sqrt(Float64((j_a + 1)*(j_d + 1)))
                    @inbounds for e = 1:a_max
                        n_e = Orb[e].n
                        l_e = Orb[e].l
                        j_e = Orb[e].j
                        if (l_d == l_e && j_d == j_e) && ((2*(n_d + n_e) + l_d + l_e) <= N_2max) && (rem(l_a + l_b,2) == rem(l_d + l_e,2))
                            @views E_BCS_partial[thread_id][] += 0.25 * ja_jd_hat * Kappa.p[a,b] * Kappa.p[d,e] * V2B(a,b,d,e,0,1,VNN.pp,Orb,Orb_NN)
                            @views E_BCS_partial[thread_id][] += 0.25 * ja_jd_hat * Kappa.n[a,b] * Kappa.n[d,e] * V2B(a,b,d,e,0,1,VNN.nn,Orb,Orb_NN)
                        end
                    end
                end
            end
        end
    end

    E_BCS = sum(x[] for x in E_BCS_partial)

    println("\nBCS ground-state pairing energy    ...   E_BCS = " * string(E_BCS) * " MeV")

    return E_BCS
end

function BCS_dN(Params::Parameters,U::pnVector,V::pnVector,Orb::Vector{NOrb})
    # Read parameters ...
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1)*(N_max + 2),2)
    dZ, dN = 0.0, 0.0

    # Calculate the particle number dispersion ...
    println("\nCalculating the BCS dispersion of proton & neutron particle numbers ...")

    @inbounds for a in 1:a_max
        j_a = Orb[a].j
        ja_hat = Float64(j_a + 1)
        dZ += 2.0 * ja_hat * U.p[a]^2 * V.p[a]^2
        dN += 2.0 * ja_hat * U.n[a]^2 * V.n[a]^2
    end

    # Calculate square roots of dispersion numbers ...
    dZ, dN = sqrt(dZ), sqrt(dN)

    println("\nBCS dispersion of nucleons numbers are ...")
    println("dZ = " * string(round(dZ,digits=5)))
    println("dN = " * string(round(dN,digits=5)))

    return pnFloat(dZ,dN)
end

function BCS_SQE(a_max::Int64,SPE::pnVector,lambda::pnFloat,Delta::pnVector)
    # Initialize vectors for SQEs ...
    pSQE, nSQE = zeros(Float64,a_max), zeros(Float64,a_max)

    # Calculate BCS single-quasiparticle energies (SQEs) ...
    println("\nEvalutiang BCS single-quasiparticle energies ...")
    @inbounds for a in 1:a_max
        pE_a = sqrt((SPE.p[a] - lambda.p)^2 + Delta.p[a]^2)
        nE_a = sqrt((SPE.n[a] - lambda.n)^2 + Delta.n[a]^2)
        pSQE[a], nSQE[a] = pE_a, nE_a
    end

    # Store SQEs ...
    SQE = pnVector(pSQE,nSQE)

    return SQE
end

function BCS_V2B_Res(Params::Parameters,Orb::Vector{NOrb},Orb_NN::NNOrb,Orb_NNN::NNNOrb,VNN::NNInt,VNNN::Array{Vector{Vector{Float32}},4},C::pnMatrix,Rho::pnMatrix)
    # Read parameters ...
    N_max = Params.Calc.Nmax
    N_2max = Params.Calc.N2max
    N_3max = Params.Calc.N3max
    a_max = div((N_max + 1)*(N_max + 2),2)
    J_max = 0

    c = Params.Calc.cV_res
    JP = JP_Ini(J_max)

    pRho, nRho = Rho.p, Rho.n

    VNN_Res = deepcopy(VNN)

    println("\nStarting calculation of residual NN interaction...\n")

    # Include NO2B NNN interaction to NN component ...
    println("\nMaking density dependent residual 2-body interaction...")

    if abs(1.0 - c) < 1e-5
        @time @inbounds for i in JP
            J = i[1]
            P = i[2]
            if P == 1
                println("Calculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = +")
            else
                println("Calculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = -")
            end
            @views N_t0 = Orb_NN.N[1,P,J+1]
            @views N_t1 = Orb_NN.N[2,P,J+1]
            @inbounds Threads.@threads for Bra in 1:max(N_t0, N_t1)

                if Bra <= N_t0
                    @views a = Orb_NN.Ind[1,P,J+1][Bra][1]
                    @views b = Orb_NN.Ind[1,P,J+1][Bra][2]
                    l_a = Orb[a].l
                    n_a = Orb[a].n
                    l_b = Orb[b].l
                    n_b = Orb[b].n
                    for Ket in 1:Bra
                        Ind = Bra + (Ket - 1) * N_t0 - div(Ket * (Ket - 1),2)
                        @views d = Orb_NN.Ind[1,P,J+1][Ket][1]
                        @views e = Orb_NN.Ind[1,P,J+1][Ket][2]
                        l_d = Orb[d].l
                        n_d = Orb[d].n
                        l_e = Orb[e].l
                        n_e = Orb[e].n
                        pnSum = 0.0
                        @inbounds for c in 1:a_max
                            n_c = Orb[c].n
                            l_c = Orb[c].l
                            j_c = Orb[c].j
                            if ((2*(n_a + n_b + n_c) + l_a + l_b + l_c) <= N_3max)
                                P3B = rem(l_a + l_b + l_c, 2) + 1
                                @inbounds for f in 1:a_max
                                    n_f = Orb[f].n
                                    l_f = Orb[f].l
                                    j_f = Orb[f].j
                                    if (l_c == l_f) && (j_c == j_f) && ((2*(n_d + n_e + n_f) + l_d + l_e + l_f) <= N_3max)

                                        Hat = 1.0/(Float64(2*J) + 1.0)
                                        ME001 = V3B_NO2B(a,b,c,0,d,e,f,0,J,1,P3B,VNNN,Orb,Orb_NNN)
                                        ME101 = V3B_NO2B(a,b,c,1,d,e,f,0,J,1,P3B,VNNN,Orb,Orb_NNN)
                                        ME011 = V3B_NO2B(a,b,c,0,d,e,f,1,J,1,P3B,VNNN,Orb,Orb_NNN)
                                        ME111 = V3B_NO2B(a,b,c,1,d,e,f,1,J,1,P3B,VNNN,Orb,Orb_NNN)
                                        ME113 = V3B_NO2B(a,b,c,1,d,e,f,1,J,3,P3B,VNNN,Orb,Orb_NNN)

                                        @views pnSum += Hat * ((1/2 * ME001 - 1/sqrt(12) * ME101 - 1/sqrt(12) * ME011 +
                                                    1/6 * ME111 + 1/3 * ME113) * pRho[c,f] + (1/2 * ME001 + 1/sqrt(12) *
                                                    ME101 + 1/sqrt(12) * ME011 + 1/6 * ME111 +  1/3 * ME113) * nRho[c,f])
                                    end
                                end
                            end
                        end
                        @views VNN_Res.pn[P,J+1][Ind] += pnSum
                    end
                end

                if Bra <= N_t1
                    @views a = Orb_NN.Ind[2,P,J+1][Bra][1]
                    @views b = Orb_NN.Ind[2,P,J+1][Bra][2]
                    l_a = Orb[a].l
                    n_a = Orb[a].n
                    l_b = Orb[b].l
                    n_b = Orb[b].n
                    for Ket in 1:Bra
                        Ind = Bra + (Ket - 1) * N_t1 - div(Ket * (Ket - 1),2)
                        @views d = Orb_NN.Ind[2,P,J+1][Ket][1]
                        @views e = Orb_NN.Ind[2,P,J+1][Ket][2]
                        l_d = Orb[d].l
                        n_d = Orb[d].n
                        l_e = Orb[e].l
                        n_e = Orb[e].n
                        ppSum = 0.0
                        nnSum = 0.0
                        @inbounds for c in 1:a_max
                            n_c = Orb[c].n
                            l_c = Orb[c].l
                            j_c = Orb[c].j
                            if ((2*(n_a + n_b + n_c) + l_a + l_b + l_c) <= N_3max)
                                P3B = rem(l_a + l_b + l_c, 2) + 1
                                @inbounds for f in 1:a_max
                                    n_f = Orb[f].n
                                    l_f = Orb[f].l
                                    j_f = Orb[f].j
                                    if (l_c == l_f) && (j_c == j_f) && ((2*(n_d + n_e + n_f) + l_d + l_e + l_f) <= N_3max)

                                        Hat = 1.0/(Float64(2*J) + 1.0)
                                        ME111 = V3B_NO2B(a,b,c,1,d,e,f,1,J,1,P3B,VNNN,Orb,Orb_NNN)
                                        ME113 = V3B_NO2B(a,b,c,1,d,e,f,1,J,3,P3B,VNNN,Orb,Orb_NNN)

                                        @views ppSum += Hat * (ME113 * pRho[c,f] + (2/3 * ME111 + 1/3 * ME113) * nRho[c,f])
                                        @views nnSum += Hat * (ME113 * nRho[c,f] + (2/3 * ME111 + 1/3 * ME113) * pRho[c,f])

                                    end
                                end
                            end
                        end
                        @views VNN_Res.pp[P,J+1][Ind] += ppSum
                        @views VNN_Res.nn[P,J+1][Ind] += nnSum
                    end
        
                end

            end
        end
    else
        @time @inbounds for i in JP
            J = i[1]
            P = i[2]
            if P == 1
                println("Calculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = +")
            else
                println("Calculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = -")
            end
            @views N_t0 = Orb_NN.N[1,P,J+1]
            @views N_t1 = Orb_NN.N[2,P,J+1]
            @inbounds Threads.@threads  for Bra in 1:max(N_t0, N_t1)

                if Bra <= N_t0
                    @views a = Orb_NN.Ind[1,P,J+1][Bra][1]
                    @views b = Orb_NN.Ind[1,P,J+1][Bra][2]
                    l_a = Orb[a].l
                    n_a = Orb[a].n
                    l_b = Orb[b].l
                    n_b = Orb[b].n
                    for Ket in 1:Bra
                        Ind = Bra + (Ket - 1) * N_t0 - div(Ket * (Ket - 1),2)
                        @views d = Orb_NN.Ind[1,P,J+1][Ket][1]
                        @views e = Orb_NN.Ind[1,P,J+1][Ket][2]
                        l_d = Orb[d].l
                        n_d = Orb[d].n
                        l_e = Orb[e].l
                        n_e = Orb[e].n
                        pnSum = 0.0
                        @inbounds for c in 1:a_max
                            n_c = Orb[c].n
                            l_c = Orb[c].l
                            j_c = Orb[c].j
                            if ((2*(n_a + n_b + n_c) + l_a + l_b + l_c) <= N_3max)
                                P3B = rem(l_a + l_b + l_c, 2) + 1
                                @inbounds for f in 1:a_max
                                    n_f = Orb[f].n
                                    l_f = Orb[f].l
                                    j_f = Orb[f].j
                                    if (l_c == l_f) && (j_c == j_f) && ((2*(n_d + n_e + n_f) + l_d + l_e + l_f) <= N_3max)

                                        Hat = 1.0/(Float64(2*J) + 1.0)
                                        ME001 = V3B_NO2B(a,b,c,0,d,e,f,0,J,1,P3B,VNNN,Orb,Orb_NNN)
                                        ME101 = V3B_NO2B(a,b,c,1,d,e,f,0,J,1,P3B,VNNN,Orb,Orb_NNN)
                                        ME011 = V3B_NO2B(a,b,c,0,d,e,f,1,J,1,P3B,VNNN,Orb,Orb_NNN)
                                        ME111 = V3B_NO2B(a,b,c,1,d,e,f,1,J,1,P3B,VNNN,Orb,Orb_NNN)
                                        ME113 = V3B_NO2B(a,b,c,1,d,e,f,1,J,3,P3B,VNNN,Orb,Orb_NNN)

                                        @views pnSum += Hat * ((1/2 * ME001 - 1/sqrt(12) * ME101 - 1/sqrt(12) * ME011 +
                                                    1/6 * ME111 + 1/3 * ME113) * pRho[c,f] + (1/2 * ME001 + 1/sqrt(12) *
                                                    ME101 + 1/sqrt(12) * ME011 + 1/6 * ME111 +  1/3 * ME113) * nRho[c,f])
                                    end
                                end
                            end
                        end
                        @views VNN.pn[P,J+1][Ind] += pnSum
                        @views VNN.pn[P,J+1][Ind] = c * VNN.pn[P,J+1][Ind]
                    end
                end

                if Bra <= N_t1
                    @views a = Orb_NN.Ind[2,P,J+1][Bra][1]
                    @views b = Orb_NN.Ind[2,P,J+1][Bra][2]
                    l_a = Orb[a].l
                    n_a = Orb[a].n
                    l_b = Orb[b].l
                    n_b = Orb[b].n
                    for Ket in 1:Bra
                        Ind = Bra + (Ket - 1) * N_t1 - div(Ket * (Ket - 1),2)
                        @views d = Orb_NN.Ind[2,P,J+1][Ket][1]
                        @views e = Orb_NN.Ind[2,P,J+1][Ket][2]
                        l_d = Orb[d].l
                        n_d = Orb[d].n
                        l_e = Orb[e].l
                        n_e = Orb[e].n
                        ppSum = 0.0
                        nnSum = 0.0
                        @inbounds for c in 1:a_max
                            n_c = Orb[c].n
                            l_c = Orb[c].l
                            j_c = Orb[c].j
                            if ((2*(n_a + n_b + n_c) + l_a + l_b + l_c) <= N_3max)
                                P3B = rem(l_a + l_b + l_c, 2) + 1
                                @inbounds for f in 1:a_max
                                    n_f = Orb[f].n
                                    l_f = Orb[f].l
                                    j_f = Orb[f].j
                                    if (l_c == l_f) && (j_c == j_f) && ((2*(n_d + n_e + n_f) + l_d + l_e + l_f) <= N_3max)

                                        Hat = 1.0/(Float64(2*J) + 1.0)
                                        ME111 = V3B_NO2B(a,b,c,1,d,e,f,1,J,1,P3B,VNNN,Orb,Orb_NNN)
                                        ME113 = V3B_NO2B(a,b,c,1,d,e,f,1,J,3,P3B,VNNN,Orb,Orb_NNN)

                                        @views ppSum += Hat * (ME113 * pRho[c,f] + (2/3 * ME111 + 1/3 * ME113) * nRho[c,f])
                                        @views nnSum += Hat * (ME113 * nRho[c,f] + (2/3 * ME111 + 1/3 * ME113) * pRho[c,f])

                                    end
                                end
                            end
                        end
                        @views VNN.pp[P,J+1][Ind] += ppSum
                        @views VNN.nn[P,J+1][Ind] += nnSum
                        @views VNN.pp[P,J+1][Ind] = c * VNN.pp[P,J+1][Ind]
                        @views VNN.nn[P,J+1][Ind] = c * VNN.nn[P,J+1][Ind]
                    end
        
                end

            end
        end
    end

    # Perform transformation of the residual NN interaction to the canonical mean-field basis ...
    println("\nTransforming residual 2-body interaction from the LHO to the target basis...")

    # Index 1
    VNN_Res_I, Orb_NN_Res = V2B_Res_Ind1(Params,JP,Orb,Orb_NN,VNN_Res,C)

    # Index 2
    VNN_Res_I, VNN_Res_II = V2B_Res_Ind2(Params,JP,Orb,Orb_NN_Res,VNN_Res_I,C)

    # Index 3
    VNN_Res_II = V2B_Res_Ind3(Params,JP,Orb,Orb_NN_Res,VNN_Res_I,VNN_Res_II,C)

    # Index 4
    VNN_Res, Orb_NN_Res = V2B_Res_Ind4(Params,JP,Orb,Orb_NN_Res,VNN_Res_II,C)
    
    println("\nResidual 2-body interaction ready...\n")

    return VNN_Res, Orb_NN_Res
end
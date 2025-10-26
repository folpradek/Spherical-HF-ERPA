function HFB_Solver(Params::Parameters)
    # Make single-particle orbitals - NuHamil ordering ...
    Orb = Make_Orbitals(Params.Calc.A,Params.Calc.Z,Params.Int.Nmax)

    # Load 1-body kinetic operator ...
    T = T1B(Params.Int.Nmax,Orb,Params.Int.hw)

    # 2-body NN interaction & Orbitals ...
    @time VNN, Orb_NN = V2B_Read(Params,Orb)

    # 3-body NNN interaction & Orbitals ...
    @time VNNN, Orb_NNN = V3B_NO2B_Read(Params,Orb)

    # Solve HFB equations ...
    @time Lambda, SQE, U, V, Rho, Kappa, H, Delta, Iteration = HFB_Solve(Params,Orb,Orb_NN,Orb_NNN,T,VNN,VNNN,1e-7)

    # Calculation of the total HFB mean-field ground-state energy ...
    @time E_HFB = HFB_Energy(Params,Rho,Kappa,Orb,Orb_NN,Orb_NNN,T,VNN,VNNN)

    # Calculate the total HFB kinetic energy ...


    # Particle number fluctuation calculation ...
    #dA = HFB_dN(Params,U,V,Orb)

    # Calculation summary ...
    #HFB_Summary(Params,E_HFB,lambda,dA,1e-7,Iteration)

    # Export of single-quasiparticle energies, amplitudes U & V & also possibly radial densities ...
    #HFB_SQS_Summary(Params,SPE,SQE,U,V,Orb)

    return
end

function HFB_Solve_BackUp(Params::Parameters,Orb::Vector{NOrb},Orb_NN::NNOrb,Orb_NNN::NNNOrb,T::Matrix{Float64},VNN::NNInt,VNNN::Array{Vector{Vector{Float32}},4},epsilon::Float64)
    # Read calculation parameters ...
    A_target, Z_target, N_target = Params.Calc.A, Params.Calc.Z, Params.Calc.A - Params.Calc.Z
    hw = Params.Int.hw
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1) * (N_max + 2), 2)

    # Setup local iteration variables ...
    Iteration, Iteration_max = 0, 150

    # Preallocate arrays ...
        # Matrices for U & V HFB ...
    V = pnMatrix(zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max))
    U = pnMatrix(zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max))
         # Matrices for densities Rho & Kappa ...
    Rho = pnMatrix(zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max))
    Kappa = pnMatrix(zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max))
        # Matrices for the single-particle field H & pairing field Delta ...
    H = pnMatrix(zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max))
    Delta = pnMatrix(zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max))
        # Vectors for single-(quasi)particle energies ...
    SQE = pnVector(zeros(Float64,a_max),zeros(Float64,a_max))

    # Setup particle numbers ...
    Z, N = 0, 0
    dZ, dN = 1.0, 1.0

    # Initial guess on chemical potentials lambda ...
    Lambda = pnFloat(0.1,0.1)

    # Initial guess on densities Rho & Kappa
    Rho, Kappa = HFB_Density_Operator_Initialize(Params,Orb)

    # Solve the spherical HFB equations ... by the means of self-consistent iteration ...
    println("\nStarting iteration of HFB equations ...\n")

    while ((dZ > epsilon) || (dN > epsilon)) && (Iteration < Iteration_max)

        # Allocate the single-particle field H and pairing field Delta ...
        H, Delta = HFB_Allocate(Params,Lambda,Rho,Kappa,Orb,Orb_NN,Orb_NNN,T,VNN,VNNN)

        # Diagonalize the HFB equations ... basis is reordered as needed ...
        SQE, U, V = HFB_Diagonalize(Params,H,Delta,Orb)

        # Generate new HF density matrix ...
        Rho_new, Kappa_new = HFB_Density_Operator(Params,U,V,Orb)
        Rho = pnMatrix(0.95 .* Rho.p + 0.05 * Rho_new.p, 0.95 .* Rho.n + 0.05 * Rho_new.n)
        Kappa = pnMatrix(0.95 .* Kappa.p + 0.05 * Kappa_new.p, 0.95 .* Kappa.n + 0.05 * Kappa_new.n)

        # Determine new average particle numbers ...
        Z, N = HFB_Particle_Number(Params,Rho,Orb)

        # Quench the chemical potential Lambda ...
        Lambda = pnFloat(Lambda.p + 0.5 * (Params.Calc.Z - Z), Lambda.n + 0.5 * (Params.Calc.A - Params.Calc.Z - N))
        #Lambda = HFB_Lambda(Params,Lambda,pnFloat(Params.Calc.Z - Z,Params.Calc.A - Params.Calc.Z - N),SQE,Delta,Orb)

        # Evaluate iteration ...
        Iteration += 1
        dZ, dN = abs(Z_target- Z), abs(N_target - N)

        println("\n\nHFB iteration number:   " * string(Iteration) * "   Proton number difference:   " * string(round(dZ, sigdigits=8))* "   &   Neutron number difference:   " * string(round(dN, sigdigits=8)))
        #println("Current values of particle numbers & chemical potentials are ...")
        println("Z = " * string(Z) * ",     pLambda = " * string(Lambda.p))
        println("N = " * string(N) * ",     nLambda = " * string(Lambda.n))

    end

    # Print arrays ...
    #display(Delta.n)
    #display(SQE.n)
    #display(U.n)
    #display(V.n)
    #display(Kappa.n)
    #display(Rho.n)
    #display(H.n)

    println("\nHFB iteration with residual NN interaction has converged ...")

    return Lambda, SQE, U, V, Rho, Kappa, H, Delta, Iteration
end

function HFB_Solve(Params::Parameters,Orb::Vector{NOrb},Orb_NN::NNOrb,Orb_NNN::NNNOrb,T::Matrix{Float64},VNN::NNInt,VNNN::Array{Vector{Vector{Float32}},4},epsilon::Float64)
    # Read calculation parameters ...
    A_target, Z_target, N_target = Params.Calc.A, Params.Calc.Z, Params.Calc.A - Params.Calc.Z
    hw = Params.Int.hw
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1) * (N_max + 2), 2)

    # Setup local iteration variables ...
    Iteration, Iteration_max = 0, 150
    dE, dZ, dN = 1.0, 1.0, 1.0

    # Preallocate arrays ...
        # Matrices for U & V HFB ...
    V = pnMatrix(zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max))
    U = pnMatrix(zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max))
         # Matrices for densities Rho & Kappa ...
    Rho = pnMatrix(zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max))
    Kappa = pnMatrix(zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max))
        # Matrices for the single-particle field H & pairing field Delta ...
    H = pnMatrix(zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max))
    Delta = pnMatrix(zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max))
        # Vectors for single-(quasi)particle energies ...
    SQE = pnVector(zeros(Float64,a_max),zeros(Float64,a_max))

    # Setup particle numbers ...
    Z, N = 0, 0

    # Initial guess on chemical potentials lambda ...
    Lambda = pnFloat(0.1,0.1)

    # Initial guess on densities Rho & Kappa
    Rho, Kappa = HFB_Density_Operator_Initialize(Params,Orb)

    # Solve the spherical HFB equations ... by the means of self-consistent iteration ...
    println("\nStarting iteration of HFB equations ...\n")

    while (dE > epsilon) && (Iteration < Iteration_max)

        Rho_old = pnMatrix(zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max))
        Kappa_old = pnMatrix(zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max))
        SQE_old = pnVector(deepcopy(SQE.p),deepcopy(SQE.n))


        Iteration_Lambda = 0
        Lambda_1, Lambda_2 = pnFloat(Lambda.p,Lambda.n), pnFloat(0.85 * Lambda.p, 0.85 *Lambda.n)
        Z_1, Z_2 = 0.0, 0.0
        N_1, N_2 = 0.0, 0.0
        while (dN > epsilon) || (dZ > epsilon) && (Iteration_Lambda < Iteration_max)

            if Iteration_Lambda == 0
                # Allocate the single-particle field H and pairing field Delta ...
                H_1, Delta_1 = HFB_Allocate(Params,Lambda_1,Rho,Kappa,Orb,Orb_NN,Orb_NNN,T,VNN,VNNN)

                # Diagonalize the HFB equations ... basis is reordered as needed ...
                SQE_1, U_1, V_1 = HFB_Diagonalize(Params,H_1,Delta_1,Orb)

                # Generate new HF density matrix ...
                Rho_1, Kappa_1 = HFB_Density_Operator(Params,U_1,V_1,Orb)

                # Determine new average particle numbers ...
                Z_1, N_1 = HFB_Particle_Number(Params,Rho_1,Orb)


                H_2, Delta_2 = HFB_Allocate(Params,Lambda_2,Rho,Kappa,Orb,Orb_NN,Orb_NNN,T,VNN,VNNN)

                # Diagonalize the HFB equations ... basis is reordered as needed ...
                SQE_2, U_2, V_2 = HFB_Diagonalize(Params,H_2,Delta_2,Orb)

                # Generate new HF density matrix ...
                Rho_2, Kappa_2 = HFB_Density_Operator(Params,U_2,V_2,Orb)

                # Determine new average particle numbers ...
                Z_2, N_2 = HFB_Particle_Number(Params,Rho_2,Orb)

                Z, N = Z_2, N_2
                Z_2, N_2 = Z_1, N_1

                Lambda_2 = pnFloat(Lambda_1.p, Lambda_1.n)

                Lambda = HFB_Lambda_Secant(pnFloat(Z_target - Z_1),pnFloat(Z_target - Z_2),Lambda_1,Lambda_2)

                SQE_old = pnVector(SQE_1.p, SQE_1.n)
                SQE = pnVector(SQE_2.p, SQE_2.n)

            else

                H, Delta = HFB_Allocate(Params,Lambda,Rho,Kappa,Orb,Orb_NN,Orb_NNN,T,VNN,VNNN)


            end

            # Diagonalize the HFB equations ... basis is reordered as needed ...
            SQE, U, V = HFB_Diagonalize(Params,H,Delta,Orb)

            # Generate new HF density matrix ...
            Rho_new, Kappa_new = HFB_Density_Operator(Params,U,V,Orb)

            # Determine new average particle numbers ...
            Z, N = HFB_Particle_Number(Params,Rho_new,Orb)

            # Quench the chemical potential Lambda ... (!!!) ... Include the secant method (!!!) ...
            #Lambda = pnFloat(Lambda.p + 0.5 * (Params.Calc.Z - Z), Lambda.n + 0.5 * (Params.Calc.A - Params.Calc.Z - N))
            #Lambda = HFB_Lambda(Params,Lambda,pnFloat(Params.Calc.Z - Z,Params.Calc.A - Params.Calc.Z - N),SQE,Delta,Orb)
            Lambda = HFB_Lambda_Secant(Lambda_1,Lambda_2,dA_1,dA_2)

            

            # Evaluate iteration ...
            Iteration_Lambda += 1
            dZ, dN = abs(Z_target- Z), abs(N_target - N)
            println("Z = " * string(Z) * ",     pLambda = " * string(Lambda.p))
            println("N = " * string(N) * ",     nLambda = " * string(Lambda.n))
            println("\nHFB lambda iteration number:     " * string(Iteration_Lambda))
        end

        #Rho = pnMatrix(0.8 .* Rho.p + 0.2 * Rho_new.p, 0.8 .* Rho.n + 0.2 * Rho_new.n)
        #Kappa = pnMatrix(0.8 .* Kappa.p + 0.2 * Kappa_new.p, 0.8 .* Kappa.n + 0.2 * Kappa_new.n)
        Rho, Kappa = pnMatrix(Rho_new.p,Rho_new.n), pnMatrix(Kappa_new.p,Kappa_new.n)

        Iteration += 1
        dE = (sum(abs.(SQE.p .- SQE_old.p )) + sum(abs.(SQE.n .- SQE_old.n))) / Float64(2 * a_max)
        println("\n\nHFB iteration number:   " * string(Iteration) * "   Single-quasiparticle energy difference:   " * string(round(dZ, sigdigits=8))* "   &   Neutron number difference:   " * string(round(dN, sigdigits=8)))
        
        #println("Current values of particle numbers & chemical potentials are ...")
        #println("Z = " * string(Z) * ",     pLambda = " * string(Lambda.p))
        #println("N = " * string(N) * ",     nLambda = " * string(Lambda.n))

    end

    # Print arrays ...
    #display(Delta.n)
    #display(SQE.n)
    #display(U.n)
    #display(V.n)
    #display(Kappa.n)
    #display(Rho.n)
    #display(H.n)

    println("\nHFB iteration with residual NN interaction has converged ...")

    return Lambda, SQE, U, V, Rho, Kappa, H, Delta, Iteration
end

function HFB_Lambda_Secant(dA_1::pnFloat,dA_2::pnFloat,Lambda_1::pnFloat,Lambda_2::pnFloat)
    pLambda = (Lambda_1.p * dA_2.p - Lambda_2.p * dA_1.p) / (dA_2.p - dA_1.p)
    nLambda = (Lambda_1.n * dA_2.n - Lambda_2.n * dA_1.n) / (dA_2.n - dA_1.n)
    return pnFloat(pLambda,nLambda)
end

function HFB_Density_Operator_Initialize(Params::Parameters,Orb::Vector{NOrb})
    # Read parameters ...
    Z_target, N_target = Params.Calc.Z, Params.Calc.A - Params.Calc.Z
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1) * (N_max + 2), 2)

    pRho, nRho = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)
    pKappa, nKappa = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)

    Z, N = 0, 0

    for a in 1:a_max
        if (Z - Z_target) < 0.0
            Z += Orb[a].j + 1
            pRho[a,a] = 1.0
        else
            break
        end
    end

    for a in 1:a_max
        if (N - N_target) < 0.0
            N += Orb[a].j + 1
            nRho[a,a] = 1.0
        else
            break
        end
    end

    for a in 1:a_max
        l_a, j_a = Orb[a].l, Orb[a].j
        for b in 1:a_max
            l_b, j_b = Orb[b].l, Orb[b].j
            if l_a == l_b && j_a == j_b && a == b
                pKappa[a,b], nKappa[a,b] = 0.01, 0.01
            end
        end
    end

    return pnMatrix(pRho,nRho), pnMatrix(pKappa,nKappa)
end

function HFB_Allocate(Params::Parameters,Lambda::pnFloat,Rho::pnMatrix,Kappa::pnMatrix,Orb::Vector{NOrb},Orb_NN::NNOrb,Orb_NNN::NNNOrb,T::Matrix{Float64},VNN::NNInt,VNNN::Array{Vector{Vector{Float32}},4})
    # Read parameters ...
    A = Params.Calc.A
    N_max = Params.Calc.Nmax
    N_2max = Params.Calc.N2max
    N_3max = Params.Calc.N3max
    CMS = Params.Calc.CMS
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Allocate new HF Hamiltonian matrices ...
    pH, nH = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)
    pDelta, nDelta = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)

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

                                                        pSum += -0.25 * jb_je_hat * ja_hat^2 * (ME113 * pKappa_cb * pKappa_ef +
                                                                1.0 / 3.0 * (2.0 * ME111 + ME113) * nKappa_cb * nKappa_ef)
                                                        nSum += -0.25 * jb_je_hat * ja_hat^2 * (ME113 * nKappa_cb * nKappa_ef +
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

                # Add 0-body chemical potential Lambda ...
                if a == d
                    pSum -= Lambda.p
                    nSum -= Lambda.n
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

    # Allocate the pairing field Delta ...
    @inbounds Threads.@threads for a in 1:a_max
        n_a = Orb[a].n
        l_a = Orb[a].l
        j_a = Orb[a].j
        ja_hat = sqrt(Float64(j_a + 1))
        @inbounds for b in 1:a_max
            n_b = Orb[b].n
            l_b = Orb[b].l
            if (l_a == l_b) && ((2*(n_a + n_b) + l_a + l_b) <= N_2max)
                j_b = Orb[b].j
                if (j_a == j_b)
                    pSum, nSum = 0.0, 0.0
                    @inbounds for d in 1:a_max
                        n_d = Orb[d].n
                        l_d = Orb[d].l
                        j_d = Orb[d].j
                        jd_hat = sqrt(Float64(j_d + 1))
                        @inbounds for e in 1:a_max
                            n_e = Orb[e].n
                            l_e = Orb[e].l
                            j_e = Orb[e].j
                            if (l_d == l_e) && (j_d == j_e) && ((2*(n_d + n_e) + l_d + l_e) <= N_2max)
                                pKappa_de, nKappa_de = Kappa.p[d,e], Kappa.n[d,e]
                                # 2-body NN interaction part ...
                                pSum += -0.25 * jd_hat / ja_hat * V2B(a,b,d,e,0,1,VNN.pp,Orb,Orb_NN) * pKappa_de
                                nSum += -0.25 * jd_hat / ja_hat * V2B(a,b,d,e,0,1,VNN.nn,Orb,Orb_NN) * nKappa_de

                                # 3-body NNN interaction part ...
                                @inbounds for c in 1:a_max
                                    n_c = Orb[c].n
                                    l_c = Orb[c].l
                                    if (2*(n_a + n_b + n_c) + l_a + l_b + l_c) <= N_3max
                                        P = rem(l_a + l_b + l_c, 2) + 1
                                        j_c = Orb[c].j
                                        jc_hat = Float64(j_c + 1)^2

                                        @inbounds for f in 1:a_max
                                            n_f = Orb[f].n
                                            l_f = Orb[f].l
                                            if (2*(n_d + n_e + n_f) + l_d + l_e + l_f) <= N_3max && l_c == l_f && P == (rem(l_d + l_e + l_f,2) + 1)
                                                j_f = Orb[f].j
                                                if j_c == j_f
                                                    pRho_cf, nRho_cf = Rho.p[c,f], Rho.n[c,f]

                                                    ME111 = V3B_NO2B(a,b,c,1,d,e,f,1,0,1,P,VNNN,Orb,Orb_NNN)
                                                    ME113 = V3B_NO2B(a,b,c,1,d,e,f,1,0,3,P,VNNN,Orb,Orb_NNN)

                                                    pSum -= 0.25 * ja_hat * jd_hat / jc_hat * (ME113 * pKappa_de * pRho_cf +
                                                            1.0 / 3.0 * (2.0 * ME111 + ME113) * pKappa_de * nRho_cf)
                                                    nSum -= 0.25 * ja_hat * jd_hat / jc_hat * (ME113 * nKappa_de * nRho_cf +
                                                            1.0 / 3.0 * (2.0 * ME111 + ME113) * nKappa_de * pRho_cf)
                                                    
                                                end
                                            end
                                        end

                                    end
                                end

                            end
                        end
                    end
                    pDelta[a,b] = pSum
                    nDelta[a,b] = nSum
                end
            end
        end
    end

    return pnMatrix(pH,nH), pnMatrix(pDelta,nDelta)
end

function HFB_Diagonalize(Params::Parameters,H::pnMatrix,Delta::pnMatrix,Orb::Vector{NOrb})
    # Read parameters ...
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1) * (N_max + 2), 2)

    # Preallocate arrays for solutions ...
    pU, pV, pSQE = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max), zeros(Float64,a_max)
    nU, nV, nSQE = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max), zeros(Float64,a_max)

    # Allocate the HFB eigenvalue system ...
    pHFB = [H.p Delta.p; -1.0 .* Delta.p -1.0 .* H.p]
    nHFB = [H.n Delta.n; -1.0 .* Delta.n -1.0 .* H.n]

    # Solve the HFB equations ...
    pE, pC = eigen(pHFB)
    nE, nC = eigen(nHFB)

    # Only positive energy solutions are extracted ...
    @inbounds for a in 1:a_max
        pSQE[a], nSQE[a] = pE[a+a_max], nE[a+a_max]
        pU[:,a], pV[:,a] = pC[1:a_max,a+a_max], pC[a_max+1:2*a_max,a+a_max]
        nU[:,a], nV[:,a] = nC[1:a_max,a+a_max], nC[a_max+1:2*a_max,a+a_max]
    end

    # Perform reordering - to match quantum numbers j & l & ascending in SQE ...
    SQE, U, V = HFB_Orbital_Ordering(Params,pnVector(pSQE,nSQE),pnMatrix(pU,nU),pnMatrix(pV,nV),Orb)

    return SQE, U, V
end

function HFB_Orbital_Ordering(Params::Parameters,SQE::pnVector,U::pnMatrix,V::pnMatrix,Orb::Vector{NOrb})
    # Read parameters ...
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1) * (N_max + 2), 2)

    # Read input arrays ...
    pU, pV, pSQE =  U.p, V.p, SQE.p
    nU, nV, nSQE =  U.n, V.n, SQE.n

    # Preallocate temporary arrays ...
    pOrb_order, nOrb_order = Vector{Int64}(undef,a_max), Vector{Int64}(undef,a_max)
    pOrb_mask, nOrb_mask = falses(a_max), falses(a_max)

    pl_values, pj_values = zeros(Float64,a_max), zeros(Float64,a_max)
    nl_values, nj_values = zeros(Float64,a_max), zeros(Float64,a_max)

    # Evaluate values of j & l for single-quasiparticle orbitals ...
    for a in 1:a_max
        pjSum, plSum = 0.0, 0.0
        njSum, nlSum = 0.0, 0.0
        for b in 1:a_max
            pjSum, plSum = pjSum + Float64(Orb[b].j) * (pU[b,a]^2 + pV[b,a]^2), plSum + Float64(Orb[b].l) * (pU[b,a]^2 + pV[b,a]^2)
            njSum, nlSum = njSum + Float64(Orb[b].j) * (nU[b,a]^2 + nV[b,a]^2), nlSum + Float64(Orb[b].l) * (nU[b,a]^2 + nV[b,a]^2)
        end
        pl_values[a], pj_values[a] = plSum, pjSum
        nl_values[a], nj_values[a] = nlSum, njSum
    end

    # Find the ordering for single-quasiparticle orbitals ...
    for a in 1:a_max
        pl, pj = pl_values[a], pj_values[a]
        nl, nj = nl_values[a], nj_values[a]
        for b in 1:a_max
            if pOrb_mask[b] == false && abs(Float64(Orb[b].l) - pl) < 1e-7 && abs(Float64(Orb[b].j) - pj) < 1e-7
                pOrb_order[b] = a
                pOrb_mask[b] = true
                break
            end
        end

        for b in 1:a_max
            if nOrb_mask[b] == false && abs(Float64(Orb[b].l) - nl) < 1e-7 && abs(Float64(Orb[b].j) - nj) < 1e-7
                nOrb_order[b] = a
                nOrb_mask[b] = true
                break
            end
        end
    end

    # Perform reordering of single-quasiparticle orbitals ...
    pSQE, pU, pV = pSQE[pOrb_order], pU[:,pOrb_order], pV[:,pOrb_order]
    nSQE, nU, nV = nSQE[nOrb_order], nU[:,nOrb_order], nV[:,nOrb_order]

    return pnVector(pSQE,nSQE), pnMatrix(pU,nU), pnMatrix(pV,nV)
end

function HFB_Density_Operator(Params::Parameters,U::pnMatrix,V::pnMatrix,Orb::Vector{NOrb})
    # Read parameters ...
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1) * (N_max + 2), 2)

    # Preallocate matrices for Rho & Kappa density operators ...
    pRho, nRho = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)
    pKappa, nKappa = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)

    # Allocate density operators Rho & Kappa ...
    for a in 1:a_max
        for b in 1:a_max
            if (Orb[a].j == Orb[b].j) && (Orb[a].l == Orb[b].l)
                pRhoSum, pKappaSum = 0.0, 0.0
                nRhoSum, nKappaSum = 0.0, 0.0
                for c in 1:a_max
                    pMERho, pMEKappa = V.p[a,c] * V.p[b,c], V.p[a,c] * U.p[b,c]
                    nMERho, nMEKappa = V.n[a,c] * V.n[b,c], V.n[a,c] * U.n[b,c]
                    pRhoSum, pKappaSum = pRhoSum + pMERho, pKappaSum + pMEKappa
                    nRhoSum, nKappaSum = nRhoSum + nMERho, nKappaSum + nMEKappa
                end
                pRho[a,b], pKappa[a,b] = pRhoSum, pKappaSum
                nRho[a,b], nKappa[a,b] = nRhoSum, nKappaSum
            end
        end
    end

    return pnMatrix(pRho,nRho), pnMatrix(pKappa,nKappa)
end

function HFB_Particle_Number(Params::Parameters,Rho::pnMatrix,Orb::Vector{NOrb})
    # Read parameters ...
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1) * (N_max + 2), 2)

    # Initialize particle numbers ...
    Z, N = 0,0, 0,0

    # Evaluate the angular-momentum weighted trace of density matrices ...
    for a in 1:a_max
        Z += Rho.p[a,a] * Float64(Orb[a].j + 1)
        N += Rho.n[a,a] * Float64(Orb[a].j + 1)
    end

    return Z, N
end

function HFB_Lambda(Params::Parameters,dA::pnFloat,Lambda::pnFloat,SQE::pnVector,Delta::pnMatrix,Orb::Vector{NOrb})
    # Read parameters ...
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1) * (N_max + 2), 2)

    # Read chemical potentials ...
    pLambda, nLambda = Lambda.p, Lambda.n

    # Quench chemical potentials ...
    pQ, nQ = 0.0, 0.0
    for a in 1:a_max
        j_a = Orb[a].j
        pQ += 0.5 * Float64(j_a + 1) * Delta.p[a,a]^2 / SQE.p[a]^3
        nQ += 0.5 * Float64(j_a + 1) * Delta.n[a,a]^2 / SQE.n[a]^3
    end
    #println("\n\npQ = " * string(pQ) * " MeV,     nQ = " * string(nQ) * " MeV\n\n")
    #throw("Stop here")
    if pQ < 0.01
        pQ = 1.0
    end
    if nQ < 0.01
        nQ = 1.0
    end
    pLambda += dA.p / pQ
    nLambda += dA.n / nQ
    
    return pnFloat(pLambda,nLambda)
end

function HFB_Energy(Params::Parameters,Rho::pnMatrix,Kappa::pnMatrix,Orb::Vector{NOrb},Orb_NN::NNOrb,Orb_NNN::NNNOrb,T::Matrix{Float64},VNN::NNInt,VNNN::Array{Vector{Vector{Float32}},4})
    # Read parameters ...
    A = Params.Calc.A
    N_max = Params.Calc.Nmax
    N_2max = Params.Calc.N2max
    N_3max = Params.Calc.N3max
    CMS = Params.Calc.CMS
    a_max = div((N_max + 1)*(N_max + 2),2)

    #Calculate the HF energy ...

    E_HFB = 0.0

    println("\nCalculating total HFB ground-state energy ...")

    E_h_partial = Threads.Atomic{Float64}[Threads.Atomic{Float64}(0.0) for _ in 1:Threads.nthreads()]
    E_Delta_partial = Threads.Atomic{Float64}[Threads.Atomic{Float64}(0.0) for _ in 1:Threads.nthreads()]

    # Single-particle mean-field energy ...
    @inbounds Threads.@threads for a = 1:a_max
        thread_id = Threads.threadid()
        n_a = Orb[a].n
        l_a = Orb[a].l
        j_a = Orb[a].j
        @inbounds for b = 1:a_max
            n_b = Orb[b].n
            l_b = Orb[b].l
            if (2*(n_a + n_b) + l_a + l_b) <= N_2max
                j_b = Orb[b].j
                @inbounds for d = 1:a_max
                    l_d = Orb[d].l
                    j_d = Orb[d].j
                    if (l_a == l_d) && (j_a == j_d)
                        n_d = Orb[d].n
                        pRho_ad = Rho.p[a,d]
                        nRho_ad = Rho.n[a,d]
                        @inbounds for e = 1:a_max
                            n_e = Orb[e].n
                            l_e = Orb[e].l
                            j_e = Orb[e].j
                            if (l_b == l_e && j_b == j_e) && ((2*(n_d + n_e) + l_d + l_e) <= N_2max)
                                pRho_be = Rho.p[b,e]
                                nRho_be = Rho.n[b,e]
                                @inbounds for J = div(abs(j_a - j_b),2):div((j_a + j_b),2)

                                    # 2-body NN interaction ...
                                    if (rem(l_a + l_b, 2) == rem(l_d + l_e, 2))

                                        Hat =  Float64(2*J + 1)
                                        @views E_h_partial[thread_id][] += 0.5 * Hat * pRho_ad * pRho_be * V2B(a,b,d,e,J,1,VNN.pp,Orb,Orb_NN)
                                        @views E_h_partial[thread_id][] += 0.5 * Hat * nRho_ad * nRho_be * V2B(a,b,d,e,J,1,VNN.nn,Orb,Orb_NN)
                                        @views E_h_partial[thread_id][] += Hat * pRho_ad * nRho_be * V2B(a,b,d,e,J,0,VNN.pn,Orb,Orb_NN)

                                    end
                                    
                                    # 3-body NNN interaction ...
                                    @inbounds for c = 1:a_max
                                        n_c = Orb[c].n
                                        l_c = Orb[c].l
                                        if ((2*(n_a + n_b + n_c) + l_a + l_b + l_c) <= N_3max)
                                            j_c = Orb[c].j
                                            P = rem(l_a + l_b + l_c, 2) + 1
                                            @inbounds for f = 1:a_max
                                                n_f = Orb[f].n
                                                l_f = Orb[f].l
                                                j_f = Orb[f].j
                                                if (2*(n_d + n_e + n_f) + l_d + l_e + l_f) <= N_3max && l_c == l_f && j_c == j_f && P == (rem(l_d + l_e + l_f, 2)+1)
                                                    pRho_cf, nRho_cf = Rho.p[c,f], Rho.n[c,f]

                                                    me1 = V3B_NO2B(a,b,c,0,d,e,f,0,J,1,P,VNNN,Orb,Orb_NNN)
                                                    me2 = V3B_NO2B(a,b,c,1,d,e,f,0,J,1,P,VNNN,Orb,Orb_NNN)
                                                    me3 = V3B_NO2B(a,b,c,0,d,e,f,1,J,1,P,VNNN,Orb,Orb_NNN)
                                                    me4 = V3B_NO2B(a,b,c,1,d,e,f,1,J,1,P,VNNN,Orb,Orb_NNN)
                                                    me5 = V3B_NO2B(a,b,c,1,d,e,f,1,J,3,P,VNNN,Orb,Orb_NNN)

                                                    @views E_h_partial[thread_id][] += 1.0/6.0 * (me5 * pRho_ad * pRho_be * pRho_cf + (2.0 * me4 + me5) *
                                                    pRho_ad * pRho_be * nRho_cf + (1.5 * me1 + sqrt(3.0/4.0) * me2 +
                                                    sqrt(3.0/4.0) * me3 + 0.5 * me4 + me5)  * pRho_ad * nRho_be * nRho_cf +
                                                    me5 * nRho_ad * nRho_be * nRho_cf)
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
            if CMS == "CMS1+2B"
                @views E_h_partial[thread_id][] += (Rho.p[a,b] + Rho.n[a,b]) * T[a,b] * Float64(j_a + 1) * (1.0 - 1.0 / Float64(A))
            elseif CMS == "CMS2B"
                @views E_h_partial[thread_id][] += 0
            else
                @views E_h_partial[thread_id][] += (Rho.p[a,b] + Rho.n[a,b]) * T[a,b] * Float64(j_a + 1)
            end
        end
    end

    # Pairing field energy ...
    @inbounds Threads.@threads for a = 1:a_max
        thread_id = Threads.threadid()
        n_a = Orb[a].n
        l_a = Orb[a].l
        j_a = Orb[a].j
        j_a_hat = sqrt(Float64(Orb[a].j + 1)) 
        @inbounds for b = 1:a_max
            n_b = Orb[b].n
            l_b = Orb[b].l
            if (2*(n_a + n_b) + l_a + l_b) <= N_2max && l_a == l_b
                j_b = Orb[b].j
                if j_a == j_b
                    pKappa_ba, nKappa_ba = Kappa.p[b,a], Kappa.n[b,a]
                    @inbounds for d = 1:a_max
                        n_d = Orb[d].n
                        l_d = Orb[d].l
                        j_d = Orb[d].j
                        j_d_hat = sqrt(Float64(Orb[d].j + 1)) 
                        @inbounds for e = 1:a_max
                            n_e = Orb[e].n
                            l_e = Orb[e].l
                            if ((2*(n_d + n_e) + l_d + l_e) <= N_2max) && l_d == l_e
                                j_e = Orb[e].j
                                if j_d == j_e
                                    pKappa_de, nKappa_de = Kappa.p[d,e], Kappa.n[d,e]

                                    # 2-body NN interaction ...
                                    @views E_Delta_partial[thread_id][] += -0.25 * j_a_hat * j_d_hat * pKappa_ba * pKappa_de * V2B(a,b,d,e,0,1,VNN.pp,Orb,Orb_NN)
                                    @views E_Delta_partial[thread_id][] += -0.25 * j_a_hat * j_d_hat * nKappa_ba * nKappa_de * V2B(a,b,d,e,0,1,VNN.nn,Orb,Orb_NN)
                                
                                    # 3-body NNN interaction ...
                                    @inbounds for c = 1:a_max
                                        n_c = Orb[c].n
                                        l_c = Orb[c].l
                                        if ((2*(n_a + n_b + n_c) + l_a + l_b + l_c) <= N_3max)
                                            j_c = Orb[c].j
                                            j_c_hat = Float64(Orb[c].j + 1)
                                            P = rem(l_a + l_b + l_c, 2) + 1
                                            @inbounds for f = 1:a_max
                                                n_f = Orb[f].n
                                                l_f = Orb[f].l
                                                j_f = Orb[f].j
                                                if (2*(n_d + n_e + n_f) + l_d + l_e + l_f) <= N_3max && l_c == l_f && j_c == j_f && P == (rem(l_d + l_e + l_f, 2)+1)
                                                    pRho_cf, nRho_cf = Rho.p[c,f], Rho.n[c,f]

                                                    ME111 = V3B_NO2B(a,b,c,1,d,e,f,1,0,1,P,VNNN,Orb,Orb_NNN)
                                                    ME113 = V3B_NO2B(a,b,c,1,d,e,f,1,0,3,P,VNNN,Orb,Orb_NNN)

                                                    @views E_Delta_partial[thread_id][] += -0.25 * j_a_hat * j_d_hat / j_c_hat * (pKappa_ba * pKappa_de *
                                                                                            pRho_cf * ME113 + pKappa_ba * pKappa_de * nRho_cf * 1.0 / 3.0 *
                                                                                            (2.0 * ME111 + ME113) + nKappa_ba * nKappa_de * nRho_cf * ME113 +
                                                                                            nKappa_ba * nKappa_de * pRho_cf * 1.0 / 3.0 *(2.0 * ME111 + ME113))
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
        end
    end

    E_HFB += sum(x[] for x in E_h_partial)
    E_HFB += sum(x[] for x in E_Delta_partial)

    println("\nHFB energy    ...   E_HFB = " * string(E_HFB) * " MeV")

    return E_HFB
end

function HFB_Kinetic_Energy(Params::Parameters,Rho::pnMatrix,Orb::Vector{NOrb},T::Matrix{Float64})
    # Read parameters ...
    A = Params.Calc.A
    hw = Params.Calc.hw
    N_max = Params.Calc.Nmax
    N_2max = Params.Calc.N2max
    CMS = Params.Calc.CMS
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Initialize matrix for kinetic operator ...
    pT, nT = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)

    # Allocate effective 1-body kinetic operator ...
    @inbounds Threads.@threads for a = 1:a_max
        n_a = Orb[a].n
        l_a = Orb[a].l
        j_a = Orb[a].j
        @inbounds for d = 1:a_max
            n_d = Orb[d].n
            l_d = Orb[d].l
            j_d = Orb[d].j
            if l_a == l_d && j_a == j_d
                pSum = 0.0
                nSum = 0.0
                @inbounds for b = 1:a_max
                    n_b = Orb[b].n
                    l_b = Orb[b].l
                    j_b = Orb[b].j
                    if (2*(n_a + n_b) + l_a + l_b) <= N_2max
                        @inbounds for e = 1:a_max
                            n_e = Orb[e].n
                            l_e = Orb[e].l
                            j_e = Orb[e].j
                            if l_b == l_e && j_b == j_e && (2*(n_d + n_e) + l_d + l_e) <= N_2max

                                # 2-body CMS correction
                                if rem(l_a + l_b, 2) == rem(l_d + l_e, 2)
                                    @inbounds for j = div(abs(j_a - j_b),2):div((j_a + j_b),2)

                                        if CMS == "CMS1+2B"
                                            TNN_sym =  T2B(Orb,a,b,d,e,j) * hw
    
                                            TNN_antisym = 1.0 / sqrt(Float64((1 + KroneckerDelta(a,b))*(1 + KroneckerDelta(d,e)))) * (T2B(Orb, a, b, d, e, j) -
                                                          Float64((-1)^(round(div(j_d + j_e,2) - j))) * T2B(Orb, a, b, e, d, j)) * hw
                                        elseif CMS == "CMS2B"
                                            Amp = hw / sqrt(Float64(1 + KroneckerDelta(a,b)) * Float64(1 + KroneckerDelta(d,e)))
                                            Amp_2 = 1.0 / sqrt(Float64(1 + KroneckerDelta(a,b)) * Float64(1 + KroneckerDelta(d,e)))
    
                                            TNN_sym = hw * T2B(Orb, a, b, d, e, j) + T[a,d] * KroneckerDelta(b,e) + T[b,e] * KroneckerDelta(a,d)

                                            TNN_antisym = (Amp * T2B(Orb,a,b,d,e,j) + Amp_2 * (KroneckerDelta(b,e) * T[a,d] + KroneckerDelta(a,d) * T[b,e])
                                                        - Float64((-1)^(div(j_d + j_e,2) - j)) * (Amp * T2B(Orb,a,b,e,d,j) + Amp_2 * (Float64(KroneckerDelta(b,d)) *
                                                        T[a,e] + Float64(KroneckerDelta(a,e)) * T[b,d]))) / Float64(A)
                                        else
                                            TNN_sym = 0.0
                                            TNN_antisym = 0.0
                                        end

                                        pSum += 1.0/Float64(A) * TNN_antisym * Rho.p[b,e] * Float64(2*j + 1) / Float64(j_a + 1)
                                        pSum += 1.0/Float64(A) * TNN_sym * Rho.n[b,e] * Float64(2*j + 1) / Float64(j_a + 1)
                                        nSum += 1.0/Float64(A) * TNN_antisym * Rho.n[b,e] * Float64(2*j + 1) / Float64(j_a + 1)
                                        nSum += 1.0/Float64(A) * TNN_sym * Rho.p[b,e] * Float64(2*j + 1) / Float64(j_a + 1)

                                    end
                                end

                            end
                        end
                    end
                end

                # 1-body kinetic operator & CMS correction
                if CMS == "CMS1+2B"
                    pT[a,d] = pSum + T[a,d] * (1.0 - 1.0/Float64(A))
                    nT[a,d] = nSum + T[a,d] * (1.0 - 1.0/Float64(A))
                elseif CMS == "CMS2B"
                    pT[a,d] = pSum
                    nT[a,d] = nSum
                else
                    pT[a,d] = pSum + T[a,d]
                    nT[a,d] = nSum + T[a,d]
                end

            end
        end
    end

    # Calculate the total ground-state kinetic energy ...
    T_HF = 0.0

    @inbounds for a = 1:a_max
        T_HF += (pT[a,a] * Rho.p[a,a] + nT[a,a] * Rho.n[a,a]) * Float64(Orb[a].j + 1)
    end

    println("\nTotal mean-field kinetic energy reads:    <T> = " * string(T_HF) * " MeV")

    return T_HF, pnMatrix(pT,nT)
end
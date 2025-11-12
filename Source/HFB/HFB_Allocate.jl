function HFB_Allocate_Indices(Params::Parameters,Orb::Vector{NOrb})
    # Read # initialize parameters ...
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1)*(N_max + 2),2)
    ab_count = 0

    # Count how many pairs of (ab) are there ...
    @inbounds for a in 1:a_max
        l_a = Orb[a].l
        j_a = Orb[a].j
        @inbounds for b in 1:a
            l_b = Orb[b].l
            j_b = Orb[b].j
            if j_a == j_b && l_a == l_b
                ab_count += 1
            end
        end
    end

    # Initialite the array ab ...
    ab, ab_count = zeros(Int64,3,3*ab_count), 0

    # Allocate the array ab for single-particle contributions to the field H ...
    @inbounds for a in 1:a_max
        l_a = Orb[a].l
        j_a = Orb[a].j
        @inbounds for b in 1:a
            l_b = Orb[b].l
            j_b = Orb[b].j
            if j_a == j_b && l_a == l_b
                ab_count += 1
                ab[1,ab_count] = a
                ab[2,ab_count] = b
                ab[3,ab_count] = 1
            end
        end
    end

    # Allocate the array ab for pairing contributions to the field H ...
    @inbounds for a in 1:a_max
        l_a = Orb[a].l
        j_a = Orb[a].j
        @inbounds for b in 1:a
            l_b = Orb[b].l
            j_b = Orb[b].j
            if j_a == j_b && l_a == l_b
                ab_count += 1
                ab[1,ab_count] = a
                ab[2,ab_count] = b
                ab[3,ab_count] = 2
            end
        end
    end

    # Allocate the array ab for contributions to the field Delta ...
    @inbounds for a in 1:a_max
        l_a = Orb[a].l
        j_a = Orb[a].j
        @inbounds for b in 1:a
            l_b = Orb[b].l
            j_b = Orb[b].j
            if j_a == j_b && l_a == l_b
                ab_count += 1
                ab[1,ab_count] = a
                ab[2,ab_count] = b
                ab[3,ab_count] = 3
            end
        end
    end

    return ab, ab_count
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
    pH, pH_A = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)
    nH, nH_A = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)
    pDelta, nDelta = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)

    # Allocate  indices for HFB iteration ...
    ad, ad_count = HFB_Allocate_Indices(Params,Orb)



    # Allocate the fields H & Delta ...
    @inbounds Threads.@threads for ad_i in 1:ad_count
        a, d, Type = ad[1,ad_i], ad[2,ad_i], ad[3,ad_i]
        n_a, l_a, j_a = Orb[a].n, Orb[a].l, Orb[a].j
        n_d, l_d, j_d = Orb[d].n, Orb[d].l, Orb[d].j

        j_a_hat = sqrt(Float64(j_a + 1))
        s_j_a_hat = Float64(j_a + 1)
        is_j_a_hat = 1.0 / (Float64(j_a) + 1.0)

        # Single-particle field H ... normal density contributions ...
        if Type == 1
            pHSum, nHSum = 0.0, 0.0

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
                                    J_j_a_hat = Float64(2*J + 1) * is_j_a_hat
                
                                    pHSum += J_j_a_hat * V2B(a,b,d,e,J,1,VNN.pp,Orb,Orb_NN) * pRho_be
                                    pHSum += J_j_a_hat * V2B(a,b,d,e,J,0,VNN.pn,Orb,Orb_NN) * nRho_be
                                    nHSum += J_j_a_hat * V2B(a,b,d,e,J,1,VNN.nn,Orb,Orb_NN) * nRho_be
                                    nHSum += J_j_a_hat * V2B(b,a,e,d,J,0,VNN.pn,Orb,Orb_NN) * pRho_be

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
                                                    pRho_cf, nRho_cf = Rho.p[c,f], Rho.n[c,f]

                                                    ME001 = V3B_NO2B(a,b,c,0,d,e,f,0,J,1,P,VNNN,Orb,Orb_NNN)
                                                    ME101 = V3B_NO2B(a,b,c,1,d,e,f,0,J,1,P,VNNN,Orb,Orb_NNN)
                                                    ME011 = V3B_NO2B(a,b,c,0,d,e,f,1,J,1,P,VNNN,Orb,Orb_NNN)
                                                    ME111 = V3B_NO2B(a,b,c,1,d,e,f,1,J,1,P,VNNN,Orb,Orb_NNN)
                                                    ME113 = V3B_NO2B(a,b,c,1,d,e,f,1,J,3,P,VNNN,Orb,Orb_NNN)

                                                    pHSum += is_j_a_hat * (0.5*ME113*pRho_be*pRho_cf + 0.25 * (ME001 +
                                                            sqrt(1.0/3.0)*ME101 + sqrt(1.0/3.0)*ME011 + ME111/3.0 + 2.0/3.0*ME113)*nRho_be*nRho_cf +
                                                            1.0/3.0 * (2.0*ME111 + ME113)*pRho_be*nRho_cf)

                                                    nHSum += is_j_a_hat * (0.5*ME113*nRho_be*nRho_cf + 0.25 * (ME001 +
                                                            sqrt(1.0/3.0)*ME101 + sqrt(1.0/3.0)*ME011 + ME111/3.0 + 2.0/3.0*ME113)*pRho_be*pRho_cf +
                                                            1.0/3.0 * (2.0*ME111 + ME113)*nRho_be*pRho_cf)

                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end

                        #=
                        # Anomal Density-dependent part ... Only 3-body NNN interaction part ...
                        if (2*(n_d + n_e) + l_d + l_e) <= N_2max
                            A_NNN_Amp = 0.25 * sqrt(Float64((j_b + 1) * (j_e + 1))) * is_j_a_hat^2
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

                                                    pHSum += A_NNN_Amp * (ME113 * pKappa_cb * pKappa_ef +
                                                            1.0 / 3.0 * (2.0 * ME111 + ME113) * nKappa_cb * nKappa_ef)
                                                    nHSum += A_NNN_Amp * (ME113 * nKappa_cb * nKappa_ef +
                                                            1.0 / 3.0 * (2.0 * ME111 + ME113) * pKappa_cb * pKappa_ef)

                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                        =#

                    end
                end
            end

            # Include the 1-body kinetic energy & inclusion of Center-of-Mass motion (CM) correction ...
                # Combined 1- + 2-body kinetic operator with CM correction ...
            if CMS == "CMS1+2B"
                pHSum += T[a,d] * (1.0 - 1.0 / A)
                nHSum += T[a,d] * (1.0 - 1.0 / A)
                # Pure 1-body kinetic operator with no CM correction ...
            elseif CMS != "CMS2B"
                pHSum += T[a,d]
                nHSum += T[a,d]
            end
                # No contribution for pure 2-body kinetic operator with CM correction ...


            # Add 0-body chemical potential Lambda ...
            if a == d
                pHSum -= Lambda.p
                nHSum -= Lambda.n
            end

            # Allocate pH & nH ...
            if a != d
                pH[a,d], pH[d,a] = pHSum, pHSum
                nH[a,d], nH[d,a] = nHSum, nHSum
            elseif a == d
                pH[a,a], nH[a,a] = pHSum, nHSum
            end
        end

        # Single-particle field H ... pairing density contributions ...
        ###=
        if Type == 2
            pH_ASum, nH_ASum = 0.0, 0.0

            @inbounds for b in 1:a_max
                n_b = Orb[b].n
                l_b = Orb[b].l
                if (2*(n_a + n_b) + l_a + l_b) <= N_2max
                    j_b = Orb[b].j
                    @inbounds for e in 1:a_max
                        n_e = Orb[e].n
                        l_e = Orb[e].l
                        j_e = Orb[e].j

                        # Anomal Density-dependent part ... Only 3-body NNN interaction part ...
                        if (2*(n_d + n_e) + l_d + l_e) <= N_2max
                            A_NNN_Amp = 0.25 * sqrt(Float64((j_b + 1) * (j_e + 1))) * is_j_a_hat^2
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

                                                    pH_ASum += A_NNN_Amp * (ME113 * pKappa_cb * pKappa_ef +
                                                            1.0 / 3.0 * (2.0 * ME111 + ME113) * nKappa_cb * nKappa_ef)
                                                    nH_ASum += A_NNN_Amp * (ME113 * nKappa_cb * nKappa_ef +
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

            # Allocate pH & nH ...
            if a != d
                pH_A[a,d], pH_A[d,a] = pH_ASum, pH_ASum
                nH_A[a,d], nH_A[d,a] = nH_ASum, nH_ASum
            elseif a == d
                pH_A[a,a], nH_A[a,a] = pH_ASum, nH_ASum
            end
        end
        ##=#

        # Single-quasiparticle pairing field Delta ...
        if Type == 3
            if ((2*(n_a + n_d) + l_a + l_d) <= N_2max)
                pDeltaSum, nDeltaSum = 0.0, 0.0

                @inbounds for b in 1:a_max
                    n_b = Orb[b].n
                    l_b = Orb[b].l
                    j_b = Orb[b].j
                    j_b_hat = sqrt(Float64(j_b + 1))
                    NN_Amp = 0.5 * j_b_hat / j_a_hat
                    @inbounds for e in 1:a_max
                        n_e = Orb[e].n
                        l_e = Orb[e].l
                        j_e = Orb[e].j
                        if (l_b == l_e) && (j_b == j_e) && ((2*(n_b + n_e) + l_b + l_e) <= N_2max)
                            pKappa_be, nKappa_be = Kappa.p[b,e], Kappa.n[b,e]

                            # 2-body NN interaction part ...
                            pDeltaSum += NN_Amp * V2B(a,d,b,e,0,1,VNN.pp,Orb,Orb_NN) * pKappa_be
                            nDeltaSum += NN_Amp * V2B(a,d,b,e,0,1,VNN.nn,Orb,Orb_NN) * nKappa_be

                            # 3-body NNN interaction part ...
                            @inbounds for c in 1:a_max
                                n_c = Orb[c].n
                                l_c = Orb[c].l
                                if (2*(n_a + n_d + n_c) + l_a + l_d + l_c) <= N_3max
                                    P = rem(l_a + l_d + l_c, 2) + 1
                                    j_c = Orb[c].j

                                    NNN_Amp = 0.5 * j_b_hat * j_a_hat / Float64(j_c + 1)^2

                                    @inbounds for f in 1:a_max
                                        n_f = Orb[f].n
                                        l_f = Orb[f].l
                                        if (2*(n_b + n_e + n_f) + l_b + l_e + l_f) <= N_3max && l_c == l_f && P == (rem(l_b + l_e + l_f,2) + 1)
                                            j_f = Orb[f].j
                                            if j_e == j_f
                                                pRho_cf, nRho_cf = Rho.p[c,f], Rho.n[c,f]

                                                ME111 = V3B_NO2B(a,d,c,1,b,e,f,1,0,1,P,VNNN,Orb,Orb_NNN)
                                                ME113 = V3B_NO2B(a,d,c,1,b,e,f,1,0,3,P,VNNN,Orb,Orb_NNN)

                                                pDeltaSum += NNN_Amp * (ME113 * pKappa_be * pRho_cf +
                                                            (2.0 * ME111 + ME113) / 3.0 * pKappa_be * nRho_cf)
                                                nDeltaSum += NNN_Amp * (ME113 * nKappa_be * nRho_cf +
                                                            (2.0 * ME111 + ME113) / 3.0 * nKappa_be * pRho_cf)
            
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end

                if a != d
                    pDelta[a,d], pDelta[d,a] = pDeltaSum, pDeltaSum
                    nDelta[a,d], nDelta[d,a] = nDeltaSum, nDeltaSum
                elseif a == d
                    pDelta[a,d] = pDeltaSum
                    nDelta[a,d] = nDeltaSum
                end

            end
        end

    end

    # Try these ... Co-pilot proposal ...
    # Better dynamic paralelization of the inner loop,
    # also separete anomal density iteartion ... Type 3 ??
    #
    # maybe needed ??? using Base.Threads
    #=
    @inbounds Threads.@threads for ad_i in 1:ad_count
        a, d, Type = ad[1,ad_i], ad[2,ad_i], ad[3,ad_i]
        n_a, l_a, j_a = Orb[a].n, Orb[a].l, Orb[a].j
        n_d, l_d, j_d = Orb[d].n, Orb[d].l, Orb[d].j

        j_a_hat = sqrt(Float64(j_a + 1))
        s_j_a_hat = Float64(j_a + 1)
        is_j_a_hat = 1.0 / (Float64(j_a) + 1.0)

        # Single-particle field H ...
        if Type == 1
            pHSum, nHSum = 0.0, 0.0

            pHSum_local = 0.0
            nHSum_local = 0.0

            # Thread-local accumulators for each thread
            local_pH = zeros(Float64, nthreads())
            local_nH = zeros(Float64, nthreads())

            @threads for b in 1:a_max
                tid = threadid()
                pHSum_local, nHSum_local = 0.0, 0.0
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
                                    J_j_a_hat = Float64(2*J + 1) * is_j_a_hat
                
                                    pHSum_local += J_j_a_hat * V2B(a,b,d,e,J,1,VNN.pp,Orb,Orb_NN) * pRho_be
                                    pHSum_local += J_j_a_hat * V2B(a,b,d,e,J,0,VNN.pn,Orb,Orb_NN) * nRho_be
                                    nHSum_local += J_j_a_hat * V2B(a,b,d,e,J,1,VNN.nn,Orb,Orb_NN) * nRho_be
                                    nHSum_local += J_j_a_hat * V2B(b,a,e,d,J,0,VNN.pn,Orb,Orb_NN) * pRho_be

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
                                                    pRho_cf, nRho_cf = Rho.p[c,f], Rho.n[c,f]

                                                    ME001 = V3B_NO2B(a,b,c,0,d,e,f,0,J,1,P,VNNN,Orb,Orb_NNN)
                                                    ME101 = V3B_NO2B(a,b,c,1,d,e,f,0,J,1,P,VNNN,Orb,Orb_NNN)
                                                    ME011 = V3B_NO2B(a,b,c,0,d,e,f,1,J,1,P,VNNN,Orb,Orb_NNN)
                                                    ME111 = V3B_NO2B(a,b,c,1,d,e,f,1,J,1,P,VNNN,Orb,Orb_NNN)
                                                    ME113 = V3B_NO2B(a,b,c,1,d,e,f,1,J,3,P,VNNN,Orb,Orb_NNN)

                                                    pHSum_local += is_j_a_hat * (0.5*ME113*pRho_be*pRho_cf + 0.25 * (ME001 +
                                                            sqrt(1.0/3.0)*ME101 + sqrt(1.0/3.0)*ME011 + ME111/3.0 + 2.0/3.0*ME113)*nRho_be*nRho_cf +
                                                            1.0/3.0 * (2.0*ME111 + ME113)*pRho_be*nRho_cf)

                                                    nHSum_local += is_j_a_hat * (0.5*ME113*nRho_be*nRho_cf + 0.25 * (ME001 +
                                                            sqrt(1.0/3.0)*ME101 + sqrt(1.0/3.0)*ME011 + ME111/3.0 + 2.0/3.0*ME113)*pRho_be*pRho_cf +
                                                            1.0/3.0 * (2.0*ME111 + ME113)*nRho_be*pRho_cf)

                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end

                        # Anomal Density-dependent part ... Only 3-body NNN interaction part ...
                        if (2*(n_d + n_e) + l_d + l_e) <= N_2max
                            A_NNN_Amp = 0.25 * sqrt(Float64((j_b + 1) * (j_e + 1))) * is_j_a_hat^2
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

                                                    pHSum_local += A_NNN_Amp * (ME113 * pKappa_cb * pKappa_ef +
                                                            1.0 / 3.0 * (2.0 * ME111 + ME113) * nKappa_cb * nKappa_ef)
                                                    nHSum_local += A_NNN_Amp * (ME113 * nKappa_cb * nKappa_ef +
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

                local_pH[tid] += pHSum_local
                local_nH[tid] += nHSum_local

            end

            pHSum = sum(local_pH)
            nHSum = sum(local_nH)

            # Include the 1-body kinetic energy & inclusion of Center-of-Mass motion (CM) correction ...
                # Combined 1- + 2-body kinetic operator with CM correction ...
            if CMS == "CMS1+2B"
                pHSum += T[a,d] * (1.0 - 1.0 / A)
                nHSum += T[a,d] * (1.0 - 1.0 / A)
                # Pure 1-body kinetic operator with no CM correction ...
            elseif CMS != "CMS2B"
                pHSum += T[a,d]
                nHSum += T[a,d]
            end
                # No contribution for pure 2-body kinetic operator with CM correction ...


            # Add 0-body chemical potential Lambda ...
            if a == d
                pHSum -= Lambda.p
                nHSum -= Lambda.n
            end

            # Allocate pH & nH ...
            if a != d
                pH[a,d], pH[d,a] = pHSum, pHSum
                nH[a,d], nH[d,a] = nHSum, nHSum
            elseif a == d
                pH[a,a], nH[a,a] = pHSum, nHSum
            end
        end

        # Single-quasiparticle pairing field Delta ...
        if Type == 2
            if ((2*(n_a + n_d) + l_a + l_d) <= N_2max)
                pDeltaSum, nDeltaSum = 0.0, 0.0

                @inbounds for b in 1:a_max
                    n_b = Orb[b].n
                    l_b = Orb[b].l
                    j_b = Orb[b].j
                    j_b_hat = sqrt(Float64(j_b + 1))
                    NN_Amp = 0.5 * j_b_hat / j_a_hat
                    @inbounds for e in 1:a_max
                        n_e = Orb[e].n
                        l_e = Orb[e].l
                        j_e = Orb[e].j
                        if (l_b == l_e) && (j_b == j_e) && ((2*(n_b + n_e) + l_b + l_e) <= N_2max)
                            pKappa_be, nKappa_be = Kappa.p[b,e], Kappa.n[b,e]

                            # 2-body NN interaction part ...
                            pDeltaSum += NN_Amp * V2B(a,d,b,e,0,1,VNN.pp,Orb,Orb_NN) * pKappa_be
                            nDeltaSum += NN_Amp * V2B(a,d,b,e,0,1,VNN.nn,Orb,Orb_NN) * nKappa_be

                            # 3-body NNN interaction part ...
                            @inbounds for c in 1:a_max
                                n_c = Orb[c].n
                                l_c = Orb[c].l
                                if (2*(n_a + n_d + n_c) + l_a + l_d + l_c) <= N_3max
                                    P = rem(l_a + l_d + l_c, 2) + 1
                                    j_c = Orb[c].j

                                    NNN_Amp = 0.5 * j_b_hat * j_a_hat / Float64(j_c + 1)^2

                                    @inbounds for f in 1:a_max
                                        n_f = Orb[f].n
                                        l_f = Orb[f].l
                                        if (2*(n_b + n_e + n_f) + l_b + l_e + l_f) <= N_3max && l_c == l_f && P == (rem(l_b + l_e + l_f,2) + 1)
                                            j_f = Orb[f].j
                                            if j_e == j_f
                                                pRho_cf, nRho_cf = Rho.p[c,f], Rho.n[c,f]

                                                ME111 = V3B_NO2B(a,d,c,1,b,e,f,1,0,1,P,VNNN,Orb,Orb_NNN)
                                                ME113 = V3B_NO2B(a,d,c,1,b,e,f,1,0,3,P,VNNN,Orb,Orb_NNN)

                                                pDeltaSum += NNN_Amp * (ME113 * pKappa_be * pRho_cf +
                                                            (2.0 * ME111 + ME113) / 3.0 * pKappa_be * nRho_cf)
                                                nDeltaSum += NNN_Amp * (ME113 * nKappa_be * nRho_cf +
                                                            (2.0 * ME111 + ME113) / 3.0 * nKappa_be * pRho_cf)
            
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end

                if a != d
                    pDelta[a,d], pDelta[d,a] = pDeltaSum, pDeltaSum
                    nDelta[a,d], nDelta[d,a] = nDeltaSum, nDeltaSum
                elseif a == d
                    pDelta[a,d] = pDeltaSum
                    nDelta[a,d] = nDeltaSum
                end

            end
        end

    end
    =#

    return pnMatrix(pH .+ pH_A, nH .+ nH_A), pnMatrix(pDelta, nDelta)
end
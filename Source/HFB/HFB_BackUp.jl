# Temporary script to store alternative versions of some HFB functions ...
function HFB_Allocate_Indices1(Params::Parameters,Orb::Vector{NOrb})
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
    ab, ab_count = zeros(Int64,2,ab_count), 0

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
            end
        end
    end

    return ab, ab_count
end

function HFB_Allocate_Indices2(Params::Parameters,Orb::Vector{NOrb})
    # Read # initialize parameters ...
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1)*(N_max + 2),2)
    ab_count, de_count = 0, 0

    # Count how many ab pairs are there ...
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

    # Count how many de pairs are there ...
    @inbounds for d in 1:a_max
        l_d = Orb[d].l
        j_d = Orb[d].j
        @inbounds for e in 1:a_max
            l_e = Orb[e].l
            j_e = Orb[e].j
            de_count += 1
        end
    end

    # Initialite the array ab ...
    ab, ab_count = zeros(Int64,2,ab_count), 0
    de, de_count = zeros(Int64,2,de_count), 0

    # Allocate the array ab ...
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
            end
        end
    end

    # Allocate the array de ...
    @inbounds for d in 1:a_max
        l_d = Orb[d].l
        j_d = Orb[d].j
        @inbounds for e in 1:a_max
            l_e = Orb[e].l
            j_e = Orb[e].j
            de_count += 1
            de[1,de_count] = d
            de[2,de_count] = e
        end
    end

    return ab, ab_count, de, de_count
end

function HFB_Allocate_Indices0(Params::Parameters,Orb::Vector{NOrb})
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

    # Allocate the array ab for the field H ... normal contributions ...
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

    # Allocate the array ab for the field H ... anomal contributions ...
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

    # Allocate the array ab for the field Delta ...
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

function HFB_Allocate1(Params::Parameters,Lambda::pnFloat,Rho::pnMatrix,Kappa::pnMatrix,Orb::Vector{NOrb},Orb_NN::NNOrb,Orb_NNN::NNNOrb,T::Matrix{Float64},VNN::NNInt,VNNN::Array{Vector{Vector{Float32}},4})
    # Read parameters ...
    A = Params.Calc.A
    N_max = Params.Calc.Nmax
    N_2max = Params.Calc.N2max
    N_3max = Params.Calc.N3max
    CMS = Params.Calc.CMS
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Allocate new HF Hamiltonian matrices ...
    pH, pDelta = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)
    nH, nDelta = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)

    # Allocate  indices for HFB iteration ...
    ad, ad_count = HFB_Allocate_Indices(Params,Orb)

    #=
    # Og imeplemetation ...
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
        #=
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
        =#

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
    =#

        # Updated implementation ... so far the fastest?
    # Allocate the fields H & Delta ...
    @inbounds for ad_i in 1:ad_count
        a, d = ad[1,ad_i], ad[2,ad_i]
        n_a, l_a, j_a = Orb[a].n, Orb[a].l, Orb[a].j
        n_d, l_d, j_d = Orb[d].n, Orb[d].l, Orb[d].j

        j_a_hat = sqrt(Float64(j_a + 1))
        s_j_a_hat = Float64(j_a + 1)
        is_j_a_hat = 1.0 / (Float64(j_a) + 1.0)


        # Thread-local accumulators for each thread ...
        pH_local = zeros(Float64, Threads.nthreads())
        nH_local = zeros(Float64, Threads.nthreads())

        # Single-particle field H ..
            # H field thread-local accumulators for each thread ...
        pH_local = zeros(Float64, Threads.nthreads())
        nH_local = zeros(Float64, Threads.nthreads())
        @inbounds Threads.@threads for b in 1:a_max
            tid = Threads.threadid()
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

            pH_local[tid] += pHSum_local
            nH_local[tid] += nHSum_local

        end

        # Sum over the Thread-local accumulators for H ...
        pHSum = sum(pH_local)
        nHSum = sum(nH_local)

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

        # Single-quasiparticle pairing field Delta ...
        if ((2*(n_a + n_d) + l_a + l_d) <= N_2max)
                # Delta field thread-local accumulators for each thread ...
            pDelta_local = zeros(Float64, Threads.nthreads())
            nDelta_local = zeros(Float64, Threads.nthreads())

            @inbounds Threads.@threads for b in 1:a_max
                tid = Threads.threadid()
                pDeltaSum_local, nDeltaSum_local = 0.0, 0.0
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
                        pDeltaSum_local += NN_Amp * V2B(a,d,b,e,0,1,VNN.pp,Orb,Orb_NN) * pKappa_be
                        nDeltaSum_local += NN_Amp * V2B(a,d,b,e,0,1,VNN.nn,Orb,Orb_NN) * nKappa_be

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

                                            pDeltaSum_local += NNN_Amp * (ME113 * pKappa_be * pRho_cf +
                                                        (2.0 * ME111 + ME113) / 3.0 * pKappa_be * nRho_cf)
                                            nDeltaSum_local += NNN_Amp * (ME113 * nKappa_be * nRho_cf +
                                                        (2.0 * ME111 + ME113) / 3.0 * nKappa_be * pRho_cf)
        
                                        end
                                    end
                                end
                            end
                        end
                    end
                end

                pDelta_local[tid] += pDeltaSum_local
                nDelta_local[tid] += nDeltaSum_local
            end

            # Sum over the Thread-local accumulators for Delta ...
            pDeltaSum = sum(pDelta_local)
            nDeltaSum = sum(nDelta_local)

            # Allocate Delta ...
            if a != d
                pDelta[a,d], pDelta[d,a] = pDeltaSum, pDeltaSum
                nDelta[a,d], nDelta[d,a] = nDeltaSum, nDeltaSum
            elseif a == d
                pDelta[a,d] = pDeltaSum
                nDelta[a,d] = nDeltaSum
            end

        end

    end


    return pnMatrix(pH, nH), pnMatrix(pDelta, nDelta)
end

function HFB_Allocate2(Params::Parameters,Lambda::pnFloat,Rho::pnMatrix,Kappa::pnMatrix,Orb::Vector{NOrb},Orb_NN::NNOrb,Orb_NNN::NNNOrb,T::Matrix{Float64},VNN::NNInt,VNNN::Array{Vector{Vector{Float32}},4})
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
    ad, ad_count, be, be_count = HFB_Allocate_Indices(Params,Orb)

        # New flattened inner loop threaded implementation ...
    # Allocate the fields H & Delta ...
    @inbounds for ad_i in 1:ad_count
        a, d = ad[1,ad_i], ad[2,ad_i]
        n_a, l_a, j_a = Orb[a].n, Orb[a].l, Orb[a].j
        n_d, l_d, j_d = Orb[d].n, Orb[d].l, Orb[d].j

        j_a_hat = sqrt(Float64(j_a + 1))
        s_j_a_hat = Float64(j_a + 1)
        is_j_a_hat = 1.0 / (Float64(j_a) + 1.0)


        # Thread-local accumulators for each thread ...
        pH_local = zeros(Float64, Threads.nthreads())
        nH_local = zeros(Float64, Threads.nthreads())

        # Single-particle field H ..
            # H field thread-local accumulators for each thread ...
        pH_local = zeros(Float64, Threads.nthreads())
        nH_local = zeros(Float64, Threads.nthreads())

        @inbounds Threads.@threads for be_i in 1:be_count
            b, e = be[1,be_i], be[2,be_i]
            n_b, l_b, j_b = Orb[b].n, Orb[b].l, Orb[b].j
            n_e, l_e, j_e = Orb[e].n, Orb[e].l, Orb[e].j

            tid = Threads.threadid()
            pHSum_local, nHSum_local = 0.0, 0.0

            if (2*(n_a + n_b) + l_a + l_b) <= N_2max

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

            pH_local[tid] += pHSum_local
            nH_local[tid] += nHSum_local

        end

        # Sum over the Thread-local accumulators for H ...
        pHSum = sum(pH_local)
        nHSum = sum(nH_local)

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

        # Single-quasiparticle pairing field Delta ...
        if ((2*(n_a + n_d) + l_a + l_d) <= N_2max)
                # Delta field thread-local accumulators for each thread ...
            pDelta_local = zeros(Float64, Threads.nthreads())
            nDelta_local = zeros(Float64, Threads.nthreads())

            @inbounds Threads.@threads for be_i in 1:be_count
                b, e = be[1,be_i], be[2,be_i]
                n_b, l_b, j_b = Orb[b].n, Orb[b].l, Orb[b].j
                n_e, l_e, j_e = Orb[e].n, Orb[e].l, Orb[e].j

                j_b_hat = sqrt(Float64(j_b + 1))
                NN_Amp = 0.5 * j_b_hat / j_a_hat

                tid = Threads.threadid()
                pDeltaSum_local, nDeltaSum_local = 0.0, 0.0

                if (l_b == l_e) && (j_b == j_e) && ((2*(n_b + n_e) + l_b + l_e) <= N_2max)
                    pKappa_be, nKappa_be = Kappa.p[b,e], Kappa.n[b,e]

                    # 2-body NN interaction part ...
                    pDeltaSum_local += NN_Amp * V2B(a,d,b,e,0,1,VNN.pp,Orb,Orb_NN) * pKappa_be
                    nDeltaSum_local += NN_Amp * V2B(a,d,b,e,0,1,VNN.nn,Orb,Orb_NN) * nKappa_be

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

                                        pDeltaSum_local += NNN_Amp * (ME113 * pKappa_be * pRho_cf +
                                                    (2.0 * ME111 + ME113) / 3.0 * pKappa_be * nRho_cf)
                                        nDeltaSum_local += NNN_Amp * (ME113 * nKappa_be * nRho_cf +
                                                    (2.0 * ME111 + ME113) / 3.0 * nKappa_be * pRho_cf)
    
                                    end
                                end
                            end
                        end
                    end
                end

                pDelta_local[tid] += pDeltaSum_local
                nDelta_local[tid] += nDeltaSum_local
            end

            # Sum over the Thread-local accumulators for Delta ...
            pDeltaSum = sum(pDelta_local)
            nDeltaSum = sum(nDelta_local)

            # Allocate Delta ...
            if a != d
                pDelta[a,d], pDelta[d,a] = pDeltaSum, pDeltaSum
                nDelta[a,d], nDelta[d,a] = nDeltaSum, nDeltaSum
            elseif a == d
                pDelta[a,d] = pDeltaSum
                nDelta[a,d] = nDeltaSum
            end

        end

    end

    return pnMatrix(pH, nH), pnMatrix(pDelta, nDelta)
end

function HFB_Allocate0(Params::Parameters,Lambda::pnFloat,Rho::pnMatrix,Kappa::pnMatrix,Orb::Vector{NOrb},Orb_NN::NNOrb,Orb_NNN::NNNOrb,T::Matrix{Float64},VNN::NNInt,VNNN::Array{Vector{Vector{Float32}},4})
    # Read parameters ...
    A = Params.Calc.A
    N_max = Params.Calc.Nmax
    N_2max = Params.Calc.N2max
    N_3max = Params.Calc.N3max
    CMS = Params.Calc.CMS
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Allocate new HF Hamiltonian matrices ...
    pH, pH_A, pDelta = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)
    nH, nH_A, nDelta = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)

    # Allocate  indices for HFB iteration ...
    ad, ad_count = HFB_Allocate_Indices(Params,Orb)

    # Og imeplemetation ...
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

    return pnMatrix(pH .+ pH_A, nH .+ nH_A), pnMatrix(pDelta, nDelta)
end

# Trivial density mixing update ... obsolete ... no longer used in the code ...
function HFB_Mixing_Update(a_max::Int64,Rho::pnMatrix,Rho_old::pnMatrix,Kappa::pnMatrix,Kappa_old::pnMatrix)
    # Basic parameters ...
    q, q_min = 0.9, 1e-8

    # Allocate finite differences for densities ...
    dRho = pnMatrix(Rho.p .- Rho_old.p, Rho.n .- Rho_old.n)
    dKappa = pnMatrix(Kappa.p .- Kappa_old.p, Kappa.n .- Kappa_old.n)

    # Evaluate current finite difference D ...
    D = (sum(abs.(dRho.p)) + sum(abs.(dKappa.p)) + sum(abs.(dRho.n)) + sum(abs.(dKappa.n))) / Float64(4*a_max^2)

    # Iterate the quenching of proton densities ...
    while q > q_min
        # Perform trial step ...
        pRho_trial, pKappa_trial = (1.0 - q) * Rho_old.p .+ q * Rho.p, (1.0 - q) * Kappa_old.p .+ q * Kappa.p
        nRho_trial, nKappa_trial = (1.0 - q) * Rho_old.n .+ q * Rho.n, (1.0 - q) * Kappa_old.n .+ q * Kappa.n

        # Symmetrize trial densities ...
        pRho_trial .= 0.5 * (pRho_trial .+ pRho_trial')
        pKappa_trial .= 0.5 * (pKappa_trial .+ pKappa_trial')
        nRho_trial .= 0.5 * (nRho_trial .+ nRho_trial')
        nKappa_trial .= 0.5 * (nKappa_trial .+ nKappa_trial')

        D_trial = (sum(abs.(pRho_trial .- Rho_old.p)) + sum(abs.(pKappa_trial .- Kappa_old.p)) + sum(abs.(nRho_trial .- Rho_old.n)) + sum(abs.(nKappa_trial .- Kappa_old.n))) / Float64(4*a_max^2)

        # Condition on accepting the current trial step ...
        if D_trial / D > 1.1
            q = 0.5 * q
        else
            return pnMatrix(pRho_trial,nRho_trial), pnMatrix(pKappa_trial,nKappa_trial)
        end
    end

    # No improvement due to quenching ... return old densities ...
    return pnMatrix(0.975 * Rho.p, 0.975 * Rho.n), pnMatrix(0.975 * Kappa.p, 0.975 * Kappa.n)
end
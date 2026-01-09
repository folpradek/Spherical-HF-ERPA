function HF_allocate(Params::Parameters,Rho::O1B,Orb::Vector{Orb1B},Orb_NN::Orb2B,Orb_NNN::Orb3B,T::O1B,V_NN::O2B,V_NNN::Array{Vector{Vector{Float32}},4})
    # Read parameters ...
    A = Float64(Params.Calc.A)
    N_max = Params.Calc.Nmax
    N_2max = Params.Calc.N2max
    N_3max = Params.Calc.N3max
    CMS = Params.Calc.CMS
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Allocate  indices for HF iteration ...
    ad, ad_count = HF_allocate_indices(Params,Orb)

    # Allocate new HF Hamiltonian matrices ...
    pH, nH = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)

    # Allocate the single-particle field H ...
    @inbounds Threads.@threads for ad_i in 1:ad_count
        a, d = ad[1,ad_i], ad[2,ad_i]
        n_a, l_a, j_a = Orb[a].n, Orb[a].l, Orb[a].j
        n_d, l_d, j_d = Orb[d].n, Orb[d].l, Orb[d].j
        is_j_a_hat = 1.0 / Float64(j_a + 1)
        pHSum, nHSum = 0.0, 0.0
        @inbounds for b in 1:a_max
            n_b = Orb[b].n
            l_b = Orb[b].l
            P = rem(l_a + l_b, 2) + 1
            if (2*(n_a + n_b) + l_a + l_b) <= N_2max
                j_b = Orb[b].j
                @inbounds for e in 1:a_max
                    n_e = Orb[e].n
                    l_e = Orb[e].l
                    j_e = Orb[e].j

                    if l_b == l_e && j_b == j_e && (2*(n_d + n_e) + l_d + l_e) <= N_2max
                        pRho_be, nRho_be = Rho.p[b,e], Rho.n[b,e]

                        @inbounds for J = div(abs(j_d - j_e),2):div(j_d + j_e,2)
                            
                            # 2-body NN interaction part ...
                            if rem(l_a + l_b, 2) == rem(l_d + l_e, 2)
                                J_j_a_hat = Float64(2*J + 1) * is_j_a_hat

                                pHSum += J_j_a_hat * pRho_be * O2b_pp(a,b,d,e,J,P,V_NN,Orb,Orb_NN)
                                pHSum += J_j_a_hat * nRho_be * O2b_pn(a,b,d,e,J,P,V_NN,Orb_NN)
                                nHSum += J_j_a_hat * nRho_be * O2b_nn(a,b,d,e,J,P,V_NN,Orb,Orb_NN)
                                nHSum += J_j_a_hat * pRho_be * O2b_pn(b,a,e,d,J,P,V_NN,Orb_NN)
                            end

                            # 3-body NNN interaction part ...
                            @inbounds for c in 1:a_max
                                n_c = Orb[c].n
                                l_c = Orb[c].l
                                if (2*(n_a + n_b + n_c) + l_a + l_b + l_c) <= N_3max
                                    P3B = rem(l_a + l_b + l_c, 2) + 1
                                    j_c = Orb[c].j
                                    @inbounds for f in 1:a_max
                                        n_f = Orb[f].n
                                        l_f = Orb[f].l
                                        if (2*(n_d + n_e + n_f) + l_d + l_e + l_f) <= N_3max && l_c == l_f && P3B == (rem(l_d + l_e + l_f,2) + 1)
                                            j_f = Orb[f].j
                                            if j_c == j_f
                                                pRho_cf, nRho_cf = Rho.p[c,f], Rho.n[c,f]

                                                ME001 = V3b_no2b(a,b,c,0,d,e,f,0,J,1,P3B,V_NNN,Orb,Orb_NNN)
                                                ME101 = V3b_no2b(a,b,c,1,d,e,f,0,J,1,P3B,V_NNN,Orb,Orb_NNN)
                                                ME011 = V3b_no2b(a,b,c,0,d,e,f,1,J,1,P3B,V_NNN,Orb,Orb_NNN)
                                                ME111 = V3b_no2b(a,b,c,1,d,e,f,1,J,1,P3B,V_NNN,Orb,Orb_NNN)
                                                ME113 = V3b_no2b(a,b,c,1,d,e,f,1,J,3,P3B,V_NNN,Orb,Orb_NNN)

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
            pHSum += T.p[a,d] * (1.0 - 1.0 / A)
            nHSum += T.n[a,d] * (1.0 - 1.0 / A)
            # Pure 1-body kinetic operator with no CM correction ...
        elseif CMS != "CMS2B"
            pHSum += T.p[a,d]
            nHSum += T.n[a,d]
        end
            # No contribution for pure 2-body kinetic operator with CM correction ...

        # Allocate pH & nH ...
        if a != d
            pH[a,d], pH[d,a] = pHSum, pHSum
            nH[a,d], nH[d,a] = nHSum, nHSum
        elseif a == d
            pH[a,a], nH[a,a] = pHSum, nHSum
        end

    end

    return O1B(pH,nH)
end

function HF_allocate_indices(Params::Parameters,Orb::Vector{Orb1B})
    # Read # initialize parameters ...
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Count how many pairs of (ab) are there ...
    ab_count = 0
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

    # Initialite the counting array ab ...
    ab, ab_count = zeros(Int64,2,ab_count), 0

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

    return ab, ab_count
end
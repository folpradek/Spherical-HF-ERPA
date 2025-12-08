function HFB_Energy(Params::Parameters,Rho::pnMatrix,Kappa::pnMatrix,Orb::Vector{NOrb},Orb_NN::NNOrb,Orb_NNN::NNNOrb,T::Matrix{Float64},VNN::NNInt,VNNN::Array{Vector{Vector{Float32}},4})
    # Read parameters ...
    A = Params.Calc.A
    N_max = Params.Calc.Nmax
    N_2max = Params.Calc.N2max
    N_3max = Params.Calc.N3max
    CMS = Params.Calc.CMS
    a_max = div((N_max + 1)*(N_max + 2),2)

    #Calculate the HF energy ...

    E_HFB = [0.0, 0.0]

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
                                    @views E_Delta_partial[thread_id][] += 0.25 * j_a_hat * j_d_hat * pKappa_ba * pKappa_de * V2B(a,b,d,e,0,1,VNN.pp,Orb,Orb_NN)
                                    @views E_Delta_partial[thread_id][] += 0.25 * j_a_hat * j_d_hat * nKappa_ba * nKappa_de * V2B(a,b,d,e,0,1,VNN.nn,Orb,Orb_NN)
                                
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

                                                    @views E_Delta_partial[thread_id][] += 0.25 * j_a_hat * j_d_hat * (pKappa_ba * pKappa_de *
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

    E_HFB[1] += sum(x[] for x in E_h_partial)
    E_HFB[2] += sum(x[] for x in E_Delta_partial)

    println("\nHFB energy    ...   E_HFB = " * string(sum(E_HFB)) * " MeV")
    println("\tMean-field energy    ...   E_MF  = " * string(E_HFB[1]) * " MeV")
    println("\tPairing energy       ...   E_Par = " * string(E_HFB[2]) * " MeV")

    return E_HFB
end
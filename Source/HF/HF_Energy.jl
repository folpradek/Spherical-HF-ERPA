function HF_Energy(Params::Vector{Any},pRho::Matrix{Float64},nRho::Matrix{Float64},Orb::Vector{NOrb},Orb_NN::NNOrb, Orb_NNN::NNNOrb, T::Matrix{Float64}, VNN::NNInt, VNNN::Array{Vector{Vector{Float32}},4})
    A = Params[5]
    N_max = Params[7]
    N_2max = Params[8]
    N_3max = Params[9]
    CMS = Params[10]

    a_max = div((N_max + 1)*(N_max + 2),2)

    E_HF = 0.0
    println("\nCalculating total HF energy ...")
    E_HF_partial = Threads.Atomic{Float64}[Threads.Atomic{Float64}(0.0) for _ in 1:Threads.nthreads()]
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
                        pRho_ad = pRho[a,d]
                        nRho_ad = nRho[a,d]
                        @inbounds for e = 1:a_max
                            n_e = Orb[e].n
                            l_e = Orb[e].l
                            j_e = Orb[e].j
                            if (l_b == l_e && j_b == j_e) && ((2*(n_d + n_e) + l_d + l_e) <= N_2max)
                                pRho_be = pRho[b,e]
                                nRho_be = nRho[b,e]
                                @inbounds for J = div(abs(j_a - j_b),2):div((j_a + j_b),2)

                                    # 2-body NN interaction
                                    if (rem(l_a + l_b, 2) == rem(l_d + l_e, 2))

                                        Hat =  Float64(2*J + 1)
                                        @views E_HF_partial[thread_id][] += 0.5 * Hat * pRho_ad * pRho_be * V2B(a,b,d,e,J,1,VNN.pp,Orb,Orb_NN)
                                        @views E_HF_partial[thread_id][] += 0.5 * Hat * nRho_ad * nRho_be * V2B(a,b,d,e,J,1,VNN.nn,Orb,Orb_NN)
                                        @views E_HF_partial[thread_id][] += Hat * pRho_ad * nRho_be * V2B(a,b,d,e,J,0,VNN.pn,Orb,Orb_NN)

                                    end
                                    
                                    # 3-body NNN interaction
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
                                                    pRho_cf = pRho[c,f]
                                                    nRho_cf = nRho[c,f]

                                                    me1 = V3B_NO2B(a,b,c,0,d,e,f,0,J,1,P,VNNN,Orb,Orb_NNN)
                                                    me2 = V3B_NO2B(a,b,c,1,d,e,f,0,J,1,P,VNNN,Orb,Orb_NNN)
                                                    me3 = V3B_NO2B(a,b,c,0,d,e,f,1,J,1,P,VNNN,Orb,Orb_NNN)
                                                    me4 = V3B_NO2B(a,b,c,1,d,e,f,1,J,1,P,VNNN,Orb,Orb_NNN)
                                                    me5 = V3B_NO2B(a,b,c,1,d,e,f,1,J,3,P,VNNN,Orb,Orb_NNN)

                                                    @views E_HF_partial[thread_id][] += 1.0/6.0 * (me5 * pRho_ad * pRho_be * pRho_cf + (2.0 * me4 + me5) *
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
                @views E_HF_partial[thread_id][] += (pRho[a,b] + nRho[a,b]) * T[a,b] * Float64(j_a + 1) * (1.0 - 1.0 / Float64(A))
            elseif CMS == "CMS2B"
                @views E_HF_partial[thread_id][] += 0
            else
                @views E_HF_partial[thread_id][] += (pRho[a,b] + nRho[a,b]) * T[a,b] * Float64(j_a + 1)
            end
        end
    end

    E_HF = sum(x[] for x in E_HF_partial)

    println("\nHF energy    ...   E_HF = " * string(E_HF) * " MeV")

    return E_HF
end

function HF_Kinetic_Energy(Params::Vector{Any},pRho::Matrix{Float64},nRho::Matrix{Float64},Orb::Vector{NOrb},T::Matrix{Float64})
    HbarOmega = Params[1]
    A = Params[5]
    N_max = Params[7]
    N_2max = Params[8]
    CMS = Params[10]

    a_max = div((N_max + 1)*(N_max + 2),2)
    pT = zeros(Float64,a_max,a_max)
    nT = zeros(Float64,a_max,a_max)

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
                                            TNN_sym =  T2B(Orb,a,b,d,e,j) * HbarOmega
    
                                            TNN_antisym = 1.0 / sqrt(Float64((1 + KroneckerDelta(a,b))*(1 + KroneckerDelta(d,e)))) * (T2B(Orb, a, b, d, e, j) -
                                                          Float64((-1)^(round(div(j_d + j_e,2) - j))) * T2B(Orb, a, b, e, d, j)) * HbarOmega
                                        elseif CMS == "CMS2B"
                                            Amp = HbarOmega / sqrt(Float64(1 + KroneckerDelta(a,b)) * Float64(1 + KroneckerDelta(d,e)))
                                            Amp_2 = 1.0 / sqrt(Float64(1 + KroneckerDelta(a,b)) * Float64(1 + KroneckerDelta(d,e)))
    
                                            TNN_sym =  HbarOmega * T2B(Orb, a, b, d, e, j) + T[a,d] * KroneckerDelta(b,e) + T[b,e] * KroneckerDelta(a,d)

                                            TNN_antisym = (Amp * T2B(Orb,a,b,d,e,j) + Amp_2 * (KroneckerDelta(b,e) * T[a,d] + KroneckerDelta(a,d) * T[b,e])
                                                        - Float64((-1)^(div(j_d + j_e,2) - j)) * (Amp * T2B(Orb,a,b,e,d,j) + Amp_2 * (Float64(KroneckerDelta(b,d)) *
                                                        T[a,e] + Float64(KroneckerDelta(a,e)) * T[b,d]))) / Float64(A)
                                        else
                                            TNN_sym = 0.0
                                            TNN_antisym = 0.0
                                        end

                                        pSum += 1.0/Float64(A) * TNN_antisym * pRho[b,e] * Float64(2*j + 1) / Float64(j_a + 1)
                                        pSum += 1.0/Float64(A) * TNN_sym * nRho[b,e] * Float64(2*j + 1) / Float64(j_a + 1)
                                        nSum += 1.0/Float64(A) * TNN_antisym * nRho[b,e] * Float64(2*j + 1) / Float64(j_a + 1)
                                        nSum += 1.0/Float64(A) * TNN_sym * pRho[b,e] * Float64(2*j + 1) / Float64(j_a + 1)

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

    T_HF = 0.0

    @inbounds for a = 1:a_max
        T_HF += (pT[a,a] * pRho[a,a] + nT[a,a] * nRho[a,a]) * Float64(Orb[a].j + 1)
    end

    println("\nTotal mean-field kinetic energy reads:    <T> = " * string(T_HF) * " MeV")

    return T_HF, pT, nT
end
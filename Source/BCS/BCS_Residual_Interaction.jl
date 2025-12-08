function BCS_V2B_Res(Params::Parameters,Orb::Vector{NOrb},Orb_NN::NNOrb,Orb_NNN::NNNOrb,VNN::NNInt,VNNN::Array{Vector{Vector{Float32}},4},C::pnMatrix,Rho::pnMatrix)
    # Read parameters ...
    N_max = Params.Calc.Nmax
    N_2max = Params.Calc.N2max
    N_3max = Params.Calc.N3max
    a_max = div((N_max + 1)*(N_max + 2),2)
    J_max = 0

    # Read the rescalling coefficients c, s2 & s3 ...
    c, s2, s3 = Params.Calc.cV_res, Params.Calc.Pairing.s2, Params.Calc.Pairing.s3

    # Initialize the array of values J & P ...
    JP = JP_Ini(J_max)

    # Read density matrix Rho ...
    pRho, nRho = Rho.p, Rho.n

    # Prepare array for residual interaction ...
    VNN_res = deepcopy(VNN)

    println("\nStarting calculation of residual NN interaction...\n")

    # Include NO2B NNN interaction to NN component ...
    println("\nMaking density dependent residual 2-body interaction...")

    if abs(1.0 - c) < 1e-5 && abs(1.0 - s2) < 1e-5 && abs(1.0 - s3) < 1e-5
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
                        @views VNN_res.pn[P,J+1][Ind] += pnSum
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
                        @views VNN_res.pp[P,J+1][Ind] += ppSum
                        @views VNN_res.nn[P,J+1][Ind] += nnSum
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
                        @views VNN_res.pn[P,J+1][Ind] = c * s2 * VNN_res.pn[P,J+1][Ind] + c * s3 * pnSum
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
                        @views VNN_res.pp[P,J+1][Ind] = c * s2 * VNN_res.pp[P,J+1][Ind] + c * s3 * ppSum
                        @views VNN_res.nn[P,J+1][Ind] = c * s2 * VNN_res.nn[P,J+1][Ind] + c * s3 * nnSum
                    end
        
                end

            end
        end
    end

    # Perform transformation of the residual NN interaction to the canonical mean-field basis ...
    println("\nTransforming residual 2-body interaction from the LHO to the target basis...")

    # Index 1
    VNN_res_I, Orb_NN_res = V2B_Res_Ind1(Params,JP,Orb,Orb_NN,VNN_res,C)

    # Index 2
    VNN_res_I, VNN_res_II = V2B_Res_Ind2(Params,JP,Orb,Orb_NN_res,VNN_res_I,C)

    # Index 3
    VNN_res_II = V2B_Res_Ind3(Params,JP,Orb,Orb_NN_res,VNN_res_I,VNN_res_II,C)

    # Index 4
    VNN_res, Orb_NN_res = V2B_Res_Ind4(Params,JP,Orb,Orb_NN_res,VNN_res_II,C)
    
    println("\nResidual 2-body interaction ready...\n")

    return VNN_res, Orb_NN_res
end
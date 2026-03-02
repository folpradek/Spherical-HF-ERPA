function HF_MBPT2_TBDM(Params::Parameters,h_N::O1B,V_NN::O2B,Orb::Vector{Orb1B},Orb_NN::Orb2B)
    # Read parameters ...
    N_2max = Params.Calc.N2max
    J_max = N_2max + 1

    # Iteraction list for J & P ...
    JP = JP_initialize(J_max)

    # Initialize 2-body operator O ...
    Rho_NN = O2b_initialize(Params,Orb;Make_Orb_NN=false)

    # Allocate the components of PT(2) 2-body density matrix Rho2B ...
    println("\nAllocating the components of MBPT(2) 2-body density matrix Rho_NN ...")
    @time @inbounds for i in JP
        J, P = i[1], i[2]
        if P == 1
            println("\tCalculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = +")
        else
            println("\tCalculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = -")
        end

        N_T0, N_T1 = Orb_NN.N[1,P,J+1], Orb_NN.N[2,P,J+1]

        @inbounds Threads.@threads for Bra in 1:max(N_T0, N_T1)
            @inbounds for Ket in 1:Bra

                # Case of pn interaction ... T = 0
                if Bra <= N_T0
                    Ind = Bra + (Ket - 1) * N_T0 - div(Ket * (Ket - 1),2)
                    a, b = Orb_NN.Ind[1,P,J+1][Bra][1], Orb_NN.Ind[1,P,J+1][Bra][2]
                    c, d = Orb_NN.Ind[1,P,J+1][Ket][1], Orb_NN.Ind[1,P,J+1][Ket][2]
                    #j_a, j_b, j_c, j_d = Orb[a].j, Orb[b].j, Orb[c].j, Orb[d].j
                    #l_a, l_b, l_c, l_d = Orb[a].l, Orb[b].l, Orb[c].l, Orb[d].l

                    # Case of Rho_pn ...
                    if (Orb[a].pO == 1 && Orb[b].nO == 1) && (Orb[c].pO == 0 && Orb[d].nO == 0)
                        pE_a, nE_b, pE_c ,nE_d = h_N.p[a,a], h_N.n[b,b], h_N.p[c,c], h_N.n[d,d]
                        ME = O2b_pn(a,b,c,d,J,P,V_NN,Orb_NN) / (pE_a + nE_b - pE_c - nE_d)
                        @views Rho_NN.pn[P,J+1][Ind] = ME
                    elseif (Orb[a].pO == 0 && Orb[b].nO == 0) && (Orb[c].pO == 1 && Orb[d].nO == 1)
                        pE_a, nE_b, pE_c ,nE_d = h_N.p[a,a], h_N.n[b,b], h_N.p[c,c], h_N.n[d,d]
                        ME = O2b_pn(a,b,c,d,J,P,V_NN,Orb_NN) / (pE_c + nE_d - pE_a - nE_b)
                        @views Rho_NN.pn[P,J+1][Ind] = ME
                    end

                end

                # Case of pp & nn interaction ... T = 1
                if Bra <= N_T1
                    Ind = Bra + (Ket - 1) * N_T1 - div(Ket * (Ket - 1),2)
                    a, b = Orb_NN.Ind[2,P,J+1][Bra][1], Orb_NN.Ind[2,P,J+1][Bra][2]
                    c, d = Orb_NN.Ind[2,P,J+1][Ket][1], Orb_NN.Ind[2,P,J+1][Ket][2]
                    #j_a, j_b, j_c, j_d = Orb[a].j, Orb[b].j, Orb[c].j, Orb[d].j
                    #l_a, l_b, l_c, l_d = Orb[a].l, Orb[b].l, Orb[c].l, Orb[d].l

                    # Case of Rho_pp ...
                     if (Orb[a].pO == 1 && Orb[b].pO == 1) && (Orb[c].pO == 0 && Orb[d].pO == 0)
                        pE_a, pE_b, pE_c ,pE_d = h_N.p[a,a], h_N.p[b,b], h_N.p[c,c], h_N.p[d,d]
                        ME = O2b_pp(a,b,c,d,J,P,V_NN,Orb,Orb_NN) / (pE_a + pE_b - pE_c - pE_d)
                        @views Rho_NN.pp[P,J+1][Ind] = ME
                    elseif (Orb[a].pO == 0 && Orb[b].pO == 0) && (Orb[c].pO == 1 && Orb[d].pO == 1)
                        pE_a, pE_b, pE_c ,pE_d = h_N.p[a,a], h_N.p[b,b], h_N.p[c,c], h_N.p[d,d]
                        ME = O2b_pp(a,b,c,d,J,P,V_NN,Orb,Orb_NN) / (pE_c + pE_d - pE_a - pE_b)
                        @views Rho_NN.pp[P,J+1][Ind] = ME
                    end

                    # Case of Rho_nn ...
                    if (Orb[a].nO == 1 && Orb[b].nO == 1) && (Orb[c].nO == 0 && Orb[d].nO == 0)
                        nE_a, nE_b, nE_c ,nE_d = h_N.n[a,a], h_N.n[b,b], h_N.n[c,c], h_N.n[d,d]
                        ME = O2b_nn(a,b,c,d,J,P,V_NN,Orb,Orb_NN) / (nE_a + nE_b - nE_c - nE_d)
                        @views Rho_NN.nn[P,J+1][Ind] = ME
                    elseif (Orb[a].nO == 0 && Orb[b].nO == 0) && (Orb[c].nO == 1 && Orb[d].nO == 1)
                        nE_a, nE_b, nE_c ,nE_d = h_N.n[a,a], h_N.n[b,b], h_N.n[c,c], h_N.n[d,d]
                        ME = O2b_nn(a,b,c,d,J,P,V_NN,Orb,Orb_NN) / (nE_c + nE_d - nE_a - nE_b)
                        @views Rho_NN.nn[P,J+1][Ind] = ME
                    end

                end

            end
        end

    end

    return Rho_NN
end

function h1b_correlations(Params::Parameters,h_N::O1B,Rho_NN::O2B,V_NNN::Array{Vector{Vector{Float32}},4},Orb::Vector{Orb1B},Orb_NN::Orb2B,Orb_NNN::Orb3B)
    # Read parameters ...
    N_max, N_2max, N_3max = Params.Calc.Nmax, Params.Calc.N2max, Params.Calc.N3max
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Import LHO -> HF transformation matrices ...
    @time C_HF = O1b_import(Params,Orb,"IO/" * Params.Calc.Path * "/Bin/C_HF.bin")

    Rho_NN_LHO = O2b_transformation(Params,Orb,Orb_NN,Rho_NN,O1B(C_HF.p',C_HF.n'))

    cf, cf_count = HF_allocate_indices(Params,Orb)

    h_p_LHO, h_n_LHO = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)

    function h1b_correlation_indiced(Params::Parameters)
        # Read parameters ...
        N_max = Params.Calc.Nmax
        a_max = div((N_max + 1)*(N_max + 2),2)

        ab_list = Vector{Tuple{Int64,Int64}}(undef,a_max*a_max)
        ab_count = 0
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                ab_count += 1
                ab_list[ab_count] = (a,b)
            end
        end

        return ab_list
    end

    ab_list = h1b_correlation_indiced(Params)

    println("\nCalculating correlated 1-body Hamiltonian h1b_correlations ...")
    @inbounds for cf_i in 1:cf_count
        c, f = cf[1,cf_i], cf[2,cf_i]
        n_c, l_c, j_c = Orb[c].n, Orb[c].l, Orb[c].j
        n_f, l_f, j_f = Orb[f].n, Orb[f].l, Orb[f].j

        j_c_hat = Float64(j_c + 1)

        ph_local = zeros(Float64,Threads.maxthreadid())
        nh_local = zeros(Float64,Threads.maxthreadid())

        @inbounds Threads.@threads :static for ab in 1:a_max^2
            a, b = ab_list[ab][1], ab_list[ab][2]
            n_a, l_a, j_a = Orb[a].n, Orb[a].l, Orb[a].j
            n_b, l_b, j_b = Orb[b].n, Orb[b].l, Orb[b].j

            Tid = Threads.threadid()
            pSum, nSum = 0.0, 0.0

            if (2*(n_a + n_b) + l_a + l_b) <= N_2max && (2*(n_a + n_b + n_c) + l_a + l_b + l_c) <= N_3max

                P2B = rem(l_a + l_b, 2) + 1
                P3B = rem(l_a + l_b + l_c, 2) + 1
                @inbounds for J in div(abs(j_a - j_b),2):div((j_a + j_b),2)

                    @inbounds for d in 1:a_max
                        n_d, l_d, j_d = Orb[d].n, Orb[d].l, Orb[d].j

                        @inbounds @simd for e in 1:a_max
                            n_e, l_e, j_e = Orb[e].n, Orb[e].l, Orb[e].j

                            if (2*(n_d + n_e) + l_d + l_e) <= N_2max && (2*(n_d + n_e + n_f) + l_d + l_e + l_f) <= N_3max &&
                            rem(l_d + l_e, 2) + 1 == P2B && rem(l_d + l_e + l_f, 2) + 1 == P3B && abs(j_d - j_e) <= 2*J && 2*J <= (j_d + j_e)

                                ppRho_abde = O2b_pp(a,b,d,e,J,P2B,Rho_NN_LHO,Orb,Orb_NN)
                                pnRho_abde = O2b_pn(a,b,d,e,J,P2B,Rho_NN_LHO,Orb_NN)
                                pnRho_baed = O2b_pn(b,a,e,d,J,P2B,Rho_NN_LHO,Orb_NN)
                                nnRho_abde = O2b_nn(a,b,d,e,J,P2B,Rho_NN_LHO,Orb,Orb_NN)

                                ME001 = V3b_no2b(a,b,c,0,d,e,f,0,J,1,P3B,V_NNN,Orb,Orb_NNN)
                                ME101 = V3b_no2b(a,b,c,1,d,e,f,0,J,1,P3B,V_NNN,Orb,Orb_NNN)
                                ME011 = V3b_no2b(a,b,c,0,d,e,f,1,J,1,P3B,V_NNN,Orb,Orb_NNN)
                                ME111 = V3b_no2b(a,b,c,1,d,e,f,1,J,1,P3B,V_NNN,Orb,Orb_NNN)
                                ME113 = V3b_no2b(a,b,c,1,d,e,f,1,J,3,P3B,V_NNN,Orb,Orb_NNN)

                                npME001 = V3b_no2b(b,a,c,0,e,d,f,0,J,1,P3B,V_NNN,Orb,Orb_NNN)
                                npME101 = V3b_no2b(b,a,c,1,e,d,f,0,J,1,P3B,V_NNN,Orb,Orb_NNN)
                                npME011 = V3b_no2b(b,a,c,0,e,d,f,1,J,1,P3B,V_NNN,Orb,Orb_NNN)
                                npME111 = V3b_no2b(b,a,c,1,e,d,f,1,J,1,P3B,V_NNN,Orb,Orb_NNN)
                                npME113 = V3b_no2b(b,a,c,1,e,d,f,1,J,3,P3B,V_NNN,Orb,Orb_NNN)

                                @inbounds pSum += 0.25 / j_c_hat * (ME113 * ppRho_abde + (2.0 * ME111 + ME113) / 3.0 * nnRho_abde +
                                        2.0 * (ME001 - ME101 / sqrt(3.0) - ME011 / sqrt(3.0) + ME111 / 3.0 + 2.0 / 3.0 * ME113) * pnRho_abde)
                                @inbounds nSum += 0.25 / j_c_hat * (ME113 * nnRho_abde + (2.0 * ME111 + ME113) / 3.0 * ppRho_abde +
                                        2.0 * (npME001 - npME101 / sqrt(3.0) - npME011 / sqrt(3.0) + npME111 / 3.0 + 2.0 / 3.0 * npME113) * pnRho_baed)

                            end

                        end

                    end

                end

            end


            ph_local[Tid] += pSum
            nh_local[Tid] += nSum

        end

        @inbounds h_p_LHO[c,f] += sum(ph_local)
        @inbounds h_n_LHO[c,f] += sum(nh_local)

        if c != f
            @inbounds h_p_LHO[f,c] += sum(ph_local)
            @inbounds h_n_LHO[f,c] += sum(nh_local)
        end


    end

    h_N.p .+= C_HF.p' * h_p_LHO * C_HF.p
    h_N.n .+= C_HF.n' * h_n_LHO * C_HF.n

    return h_N
end

function HF_ScRPA_TBDM_allocate(Params::Parameters,N_nu::Matrix{Int64},Orb_Phonon::Matrix{Vector{Int64}},Phonon::Vector{PhState},N_Particle::pnInteger,Particle::pnSVector,N_Hole::pnInteger,Hole::pnSVector,Orb::Vector{Orb1B},Orb_NN::Orb2B,X_RPA::Matrix{Matrix{ComplexF64}},Y_RPA::Matrix{Matrix{ComplexF64}})
    # Read parameters ...
    N_max = Params.Calc.Nmax
    N_2max = 2*N_max
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Preallocate arrays of allowed values of J & P ...
    JP_list = JP_initialize(Params.Calc.N2max + 1)
    
    # Initialite the 2-body NN correlation matrix Sigma_NN ...
    Sigma_NN = O2b_initialize(Params,Orb)

    function orbitals_ph_mappping(Params::Parameters,N_Particle::pnInteger,Particle::pnSVector,N_Hole::pnInteger,Hole::pnSVector)
        # Read parameters ...
        N_max = Params.Calc.Nmax
        a_max = div((N_max + 1)*(N_max + 2),2)

        orbitals_ph_list = Vector{Vector{Int64}}(undef,a_max)

        @inbounds for a in 1:a_max
            orbitals_ph_list[a] = Vector{Int64}(undef,2)
        end

        @inbounds for p in 1:N_Particle.p
            a = Particle.p[p].a
            orbitals_ph_list[a][1] = p
        end

        @inbounds for h in 1:N_Hole.p
            a = Hole.p[h].a
            orbitals_ph_list[a][1] = h
        end

        @inbounds for p in 1:N_Particle.n
            a = Particle.n[p].a
            orbitals_ph_list[a][2] = p
        end

        @inbounds for h in 1:N_Hole.n
            a = Hole.n[h].a
            orbitals_ph_list[a][2] = h
        end

        return orbitals_ph_list
    end

    orbitals_ph_list = orbitals_ph_mappping(Params,N_Particle,Particle,N_Hole,Hole)

    println("\nAllocating 2-body RPA correlation matrix Sigma_NN ...")


    @inbounds Threads.@threads for JP in JP_list;
        J, P = JP[1], JP[2]
        if P == 1
            println("\tCalculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = +")
        else
            println("\tCalculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = -")
        end
        N_T0, N_T1 = Orb_NN.N[1,P,J+1], Orb_NN.N[2,P,J+1]
        N_ph = N_nu[J+1,P]
        @inbounds for Bra in 1:max(N_T0, N_T1)
            @inbounds for Ket in 1:Bra

                # Case of pn interaction ... T = 0
                if Bra <= N_T0
                    Ind = Bra + (Ket - 1) * N_T0 - div(Ket * (Ket - 1),2)
                    a, b = Orb_NN.Ind[1,P,J+1][Bra][1], Orb_NN.Ind[1,P,J+1][Bra][2]
                    c, d = Orb_NN.Ind[1,P,J+1][Ket][1], Orb_NN.Ind[1,P,J+1][Ket][2]
                    l_a, j_a = Orb[a].l, Orb[a].j
                    l_b, j_b = Orb[b].l, Orb[b].j
                    l_c, j_c = Orb[c].l, Orb[c].j
                    l_d, j_d = Orb[d].l, Orb[d].j

                    pO_a, nO_b = Orb[a].pO, Orb[b].nO
                    pO_c, nO_d = Orb[c].pO, Orb[d].nO

                    Q = rem(l_a + l_c,2) + 1

                    pnSigmaSum = 0.0

                    # 2p-2h terms ...
                        if ((pO_a + nO_b) == 0 && (pO_c + nO_d) == 2) || ((pO_a + nO_b) == 2 && (pO_c + nO_d) == 0)

                            if (pO_a + nO_b) == 0 && (pO_c + nO_d) == 2
                                p, l_p, j_p = a, l_a, j_a
                                q, l_q, j_q = b, l_b, j_b
                                h, l_g, j_h = c, l_c, j_c
                                g, l_h, j_g = d, l_d, j_d
                            else
                                p, l_p, j_p = c, l_c, j_c
                                q, l_q, j_q = d, l_d, j_d
                                h, l_g, j_h = a, l_a, j_a
                                g, l_h, j_g = b, l_b, j_b
                            end

                            p_p, p_q = orbitals_ph_list[p][1], orbitals_ph_list[q][2]
                            h_h, h_g = orbitals_ph_list[h][1], orbitals_ph_list[g][2]

                            Q = rem(l_p + l_h,2) + 1
                            if Q == (rem(l_q + l_g,2) + 1)
                                @inbounds for I in max(div(abs(j_p - j_h),2),div(abs(j_q - j_g),2)):min(div((j_p + j_h),2),div((j_q + j_g),2))
                                    Amp = 0.5 * Float64(2*I + 1) * (-1)^(J + I + div(j_q + j_h,2)) * f6j(j_p,j_q,2*J,j_g,j_h,2*I)
                                    @inbounds for nu in 1:N_ph
                                        ME = Amp * (real(Y_RPA[I+1,Q][nu][p_p,h_h] * conj(X_RPA[I+1,Q][nu][p_q,h_g])) + real(X_RPA[I+1,Q][nu][p_p,h_h] * conj(Y_RPA[I+1,Q][nu][p_q,h_g])))
                                        #pnSigmaSum += ME
                                    end
                                end
                            end

                        end

                    #=
                        if Q == (rem(l_b + l_d,2) + 1)
                            if (pO_a + nO_b) == 0 && (pO_c + nO_d) == 2

                                @inbounds for I in max(div(abs(j_a - j_c),2),div(abs(j_b - j_d),2)):min(div((j_a + j_c),2),div((j_b + j_d),2))
                                    p_a, p_b = orbitals_ph_list[a][1], orbitals_ph_list[b][2]
                                    h_c, h_d = orbitals_ph_list[c][1], orbitals_ph_list[d][2]
                                    Phase = 0.5 * Float64(2*I + 1) * (-1)^(J + I + j_c)
                                    Amp1 = Phase * (-1)^(j_b) * f6j(j_a,j_b,2*J,j_d,j_c,2*I)
                                    Amp2 = Phase * (-1)^(j_a) * f6j(j_b,j_a,2*J,j_d,j_c,2*I)
                                    @inbounds for nu in 1:N_ph
                                        ME = Amp1 * real(Y_RPA[I+1,Q][nu][p_a,h_c] * conj(X_RPA[I+1,Q][nu][p_b,h_d])) + Amp2 * real(X_RPA[I+1,Q][nu][p_a,h_c] * conj(Y_RPA[I+1,Q][nu][p_b,h_d]))
                                        pnSigmaSum += ME
                                    end
                                end

                            elseif (pO_a + nO_b) == 2 && (pO_c + nO_d) == 0

                                    @inbounds for I in max(div(abs(j_a - j_c),2),div(abs(j_b - j_d),2)):min(div((j_a + j_c),2),div((j_b + j_d),2))
                                        h_a, h_b = orbitals_ph_list[a][1], orbitals_ph_list[b][2]
                                        p_c, p_d = orbitals_ph_list[c][1], orbitals_ph_list[d][2]
                                        Phase = 0.5 * Float64(2*I + 1) * (-1)^(J + I + j_a)
                                        Amp1 = Phase * (-1)^(j_d) * f6j(j_c,j_d,2*J,j_b,j_a,2*I)
                                        Amp2 = Phase * (-1)^(j_c) * f6j(j_d,j_c,2*J,j_b,j_a,2*I)
                                        @inbounds for nu in 1:N_ph
                                            ME = Amp1 * real(Y_RPA[I+1,Q][nu][p_c,h_a] * conj(X_RPA[I+1,Q][nu][p_d,h_b])) + Amp2 * real(X_RPA[I+1,Q][nu][p_c,h_a] * conj(Y_RPA[I+1,Q][nu][p_d,h_b]))
                                            pnSigmaSum += ME
                                        end
                                    end

                            end
                        end
                    =#

                    # 1p1-h terms ... no proton-neutron 1p-1h 2-body correlation function terms ...

                    @inbounds Sigma_NN.pn[P,J+1][Ind] = pnSigmaSum

                end

                # Case of pp & nn interaction ... T = 1
                if Bra <= N_T1
                    Ind = Bra + (Ket - 1) * N_T1 - div(Ket * (Ket - 1),2)
                    a, b = Orb_NN.Ind[2,P,J+1][Bra][1], Orb_NN.Ind[2,P,J+1][Bra][2]
                    c, d = Orb_NN.Ind[2,P,J+1][Ket][1], Orb_NN.Ind[2,P,J+1][Ket][2]
                    l_a, j_a = Orb[a].l, Orb[a].j
                    l_b, j_b = Orb[b].l, Orb[b].j
                    l_c, j_c = Orb[c].l, Orb[c].j
                    l_d, j_d = Orb[d].l, Orb[d].j

                    pO_a, pO_b = Orb[a].pO, Orb[b].pO
                    pO_c, pO_d = Orb[c].pO, Orb[d].pO

                    nO_a, nO_b = Orb[a].nO, Orb[b].nO
                    nO_c, nO_d = Orb[c].nO, Orb[d].nO

                    ppSigmaSum, nnSigmaSum = 0.0, 0.0

                    # 2p-2h channel ...
                        # proton-proton terms ...
                        if ((pO_a + pO_b) == 0 && (pO_c + pO_d) == 2) || ((pO_a + pO_b) == 2 && (pO_c + pO_d) == 0)

                            if (pO_a + pO_b) == 0 && (pO_c + pO_d) == 2
                                p, l_p, j_p = a, l_a, j_a
                                q, l_q, j_q = b, l_b, j_b
                                h, l_g, j_h = c, l_c, j_c
                                g, l_h, j_g = d, l_d, j_d
                            else
                                p, l_p, j_p = c, l_c, j_c
                                q, l_q, j_q = d, l_d, j_d
                                h, l_g, j_h = a, l_a, j_a
                                g, l_h, j_g = b, l_b, j_b
                            end

                            p_p, p_q = orbitals_ph_list[p][1], orbitals_ph_list[q][1]
                            h_h, h_g = orbitals_ph_list[h][1], orbitals_ph_list[g][1]

                            # Direct term ...
                            Q = rem(l_p + l_h,2) + 1
                            if Q == (rem(l_q + l_g,2) + 1)
                                @inbounds for I in max(div(abs(j_p - j_h),2),div(abs(j_q - j_g),2)):min(div((j_p + j_h),2),div((j_q + j_g),2))
                                    Amp = 0.5 * Float64(2*I + 1) * (-1)^(J + I + div(j_q + j_h,2)) * f6j(j_p,j_q,2*J,j_g,j_h,2*I)
                                    @inbounds for nu in 1:N_ph
                                        ME = Amp * real(Y_RPA[I+1,Q][nu][p_p,h_h] * conj(X_RPA[I+1,Q][nu][p_q,h_g]))
                                        ppSigmaSum += ME
                                    end
                                end
                            end

                            # Exchange term ...
                            Q = rem(l_p + l_g,2) + 1
                            if Q == (rem(l_q + l_h,2) + 1)
                                @inbounds for I in max(div(abs(j_p - j_g),2),div(abs(j_q - j_h),2)):min(div((j_p + j_g),2),div((j_q + j_h),2))
                                    Amp = 0.5 * Float64(2*I + 1) * (-1)^(I + div(j_q + j_g,2) + 1) * f6j(j_p,j_q,2*J,j_h,j_g,2*I)
                                    @inbounds for nu in 1:N_ph
                                        ME = Amp * real(Y_RPA[I+1,Q][nu][p_p,h_g] * conj(X_RPA[I+1,Q][nu][p_q,h_h]))
                                        ppSigmaSum += ME
                                    end
                                end
                            end

                        end

                        # Neutron-neutron terms
                        if ((nO_a + nO_b) == 0 && (nO_c + nO_d) == 2) || ((nO_a + nO_b) == 2 && (nO_c + nO_d) == 0)

                            if (nO_a + nO_b) == 0 && (nO_c + nO_d) == 2
                                p, l_p, j_p = a, l_a, j_a
                                q, l_q, j_q = b, l_b, j_b
                                h, l_g, j_h = c, l_c, j_c
                                g, l_h, j_g = d, l_d, j_d
                            else
                                p, l_p, j_p = c, l_c, j_c
                                q, l_q, j_q = d, l_d, j_d
                                h, l_g, j_h = a, l_a, j_a
                                g, l_h, j_g = b, l_b, j_b
                            end

                            p_p, p_q = orbitals_ph_list[p][2], orbitals_ph_list[q][2]
                            h_h, h_g = orbitals_ph_list[h][2], orbitals_ph_list[g][2]

                            # Direct term ...
                            Q = rem(l_p + l_h,2) + 1
                            if Q == (rem(l_q + l_g,2) + 1)
                                @inbounds for I in max(div(abs(j_p - j_h),2),div(abs(j_q - j_g),2)):min(div((j_p + j_h),2),div((j_q + j_g),2))
                                    Amp = 0.5 * Float64(2*I + 1) * (-1)^(J + I + div(j_q + j_h,2)) * f6j(j_p,j_q,2*J,j_g,j_h,2*I)
                                    @inbounds for nu in 1:N_ph
                                        ME = Amp * real(Y_RPA[I+1,Q][nu][p_p,h_h] * conj(X_RPA[I+1,Q][nu][p_q,h_g]))
                                        nnSigmaSum += ME
                                    end
                                end
                            end

                            # Exchange term ...
                            Q = rem(l_p + l_g,2) + 1
                            if Q == (rem(l_q + l_h,2) + 1)
                                @inbounds for I in max(div(abs(j_p - j_g),2),div(abs(j_q - j_h),2)):min(div((j_p + j_g),2),div((j_q + j_h),2))
                                    Amp = 0.5 * Float64(2*I + 1) * (-1)^(J + div(j_q + j_g,2) + 1) * f6j(j_p,j_q,2*J,j_h,j_g,2*I)
                                    @inbounds for nu in 1:N_ph
                                        ME = Amp * real(Y_RPA[I+1,Q][nu][p_p,h_g] * conj(X_RPA[I+1,Q][nu][p_q,h_h]))
                                        nnSigmaSum += ME
                                    end
                                end
                            end

                        end


                    # 2p-2h channel ...
                    #=
                        # proton-proton terms ...
                        if (pO_a + pO_b) == 0 && (pO_c + pO_d) == 2
                            # Direct term ...
                            Q = rem(l_a + l_c,2) + 1
                            if Q == (rem(l_b + l_d,2) + 1)
                                @inbounds for I in max(div(abs(j_a - j_c),2),div(abs(j_b - j_d),2)):min(div((j_a + j_c),2),div((j_b + j_d),2))

                                    p_a, p_b = orbitals_ph_list[a][1], orbitals_ph_list[b][1]
                                    h_c, h_d = orbitals_ph_list[c][1], orbitals_ph_list[d][1]
                                    Amp = 0.5 * Float64(2*I + 1) * (-1)^(J + I + div(j_b + j_c,2)) * f6j(j_a,j_b,2*J,j_d,j_c,2*I) / (sqrt((1 + kronecker_delta(a,b))*(1 + kronecker_delta(c,d))))
                                    
                                    @inbounds for nu in 1:N_ph
                                        ME = Amp * real(Y_RPA[I+1,Q][nu][p_a,h_c] * conj(X_RPA[I+1,Q][nu][p_b,h_d]))
                                        ppSigmaSum += ME
                                    
                                    end
                                end
                            end

                            # Exchange term ...
                            Q = rem(l_a + l_d,2) + 1
                            if Q == (rem(l_b + l_c,2) + 1)
                                @inbounds for I in max(div(abs(j_a - j_d),2),div(abs(j_b - j_c),2)):min(div((j_a + j_d),2),div((j_b + j_c),2))

                                    p_a, p_b = orbitals_ph_list[a][1], orbitals_ph_list[b][1]
                                    h_c, h_d = orbitals_ph_list[c][1], orbitals_ph_list[d][1]
                                    Amp = 0.5 * Float64(2*I + 1) * (-1)^(I + div(j_b + j_d,2) + 1) * f6j(j_a,j_b,2*J,j_c,j_d,2*I) / (sqrt((1 + kronecker_delta(a,b))*(1 + kronecker_delta(c,d))))
                                    @inbounds for nu in 1:N_ph
                                        ME = Amp *  real(Y_RPA[I+1,Q][nu][p_a,h_d] * conj(X_RPA[I+1,Q][nu][p_b,h_c]))
                                        ppSigmaSum += ME
                                    end

                                end
                            end
                            
                        elseif (pO_a + pO_b) == 2 && (pO_c + pO_d) == 0

                            # Direct term ...
                            Q = rem(l_a + l_c,2) + 1
                            if Q == (rem(l_b + l_d,2) + 1)
                                @inbounds for I in max(div(abs(j_a - j_c),2),div(abs(j_b - j_d),2)):min(div((j_a + j_c),2),div((j_b + j_d),2))

                                    h_a, h_b = orbitals_ph_list[a][1], orbitals_ph_list[b][1]
                                    p_c, p_d = orbitals_ph_list[c][1], orbitals_ph_list[d][1]
                                    Amp = 0.5 * Float64(2*I + 1) * (-1)^(J + I + div(j_a + j_d,2)) * f6j(j_c,j_d,2*J,j_b,j_a,2*I) / (sqrt((1 + kronecker_delta(a,b))*(1 + kronecker_delta(c,d))))
                                    @inbounds for nu in 1:N_ph
                                        ME = Amp * real(Y_RPA[I+1,Q][nu][p_c,h_a] * conj(X_RPA[I+1,Q][nu][p_d,h_b]))
                                        ppSigmaSum += ME
                                    end

                                end
                            end

                            # Exchange term ...
                            Q = rem(l_a + l_d,2) + 1
                            if Q == (rem(l_b + l_c,2) + 1)
                                @inbounds for I in max(div(abs(j_a - j_d),2),div(abs(j_b - j_c),2)):min(div((j_a + j_d),2),div((j_b + j_c),2))

                                    h_a, h_b = orbitals_ph_list[a][1], orbitals_ph_list[b][1]
                                    p_c, p_d = orbitals_ph_list[c][1], orbitals_ph_list[d][1]
                                    Amp = 0.5 * Float64(2*I + 1) * (-1)^(I + div(j_b + j_d,2) + 1) * f6j(j_c,j_d,2*J,j_a,j_b,2*I) / (sqrt((1 + kronecker_delta(a,b))*(1 + kronecker_delta(c,d))))
                                    @inbounds for nu in 1:N_ph
                                        ME = Amp *  real(Y_RPA[I+1,Q][nu][p_d,h_a] * conj(X_RPA[I+1,Q][nu][p_c,h_b]))
                                        ppSigmaSum += ME
                                    end

                                end
                            end

                        end

                        # neutron-neutron terms ...
                        if (nO_a + nO_b) == 0 && (nO_c + nO_d) == 2
                            # Direct term ...
                            Q = rem(l_a + l_c,2) + 1
                            if Q == (rem(l_b + l_d,2) + 1)
                                @inbounds for I in max(div(abs(j_a - j_c),2),div(abs(j_b - j_d),2)):min(div((j_a + j_c),2),div((j_b + j_d),2))

                                    p_a, p_b = orbitals_ph_list[a][2], orbitals_ph_list[b][2]
                                    h_c, h_d = orbitals_ph_list[c][2], orbitals_ph_list[d][2]
                                    Amp = 0.5 * Float64(2*I + 1) * (-1)^(J + I + div(j_b + j_c,2)) * f6j(j_a,j_b,2*J,j_d,j_c,2*I) / (sqrt((1 + kronecker_delta(a,b))*(1 + kronecker_delta(c,d))))
                                    
                                    @inbounds for nu in 1:N_ph
                                        ME = Amp * real(Y_RPA[I+1,Q][nu][p_a,h_c] * conj(X_RPA[I+1,Q][nu][p_b,h_d]))
                                        nnSigmaSum += ME
                                    
                                    end
                                end
                            end

                            # Exchange term ...
                            Q = rem(l_a + l_d,2) + 1
                            if Q == (rem(l_b + l_c,2) + 1)
                                @inbounds for I in max(div(abs(j_a - j_d),2),div(abs(j_b - j_c),2)):min(div((j_a + j_d),2),div((j_b + j_c),2))

                                    p_a, p_b = orbitals_ph_list[a][2], orbitals_ph_list[b][2]
                                    h_c, h_d = orbitals_ph_list[c][2], orbitals_ph_list[d][2]
                                    Amp = 0.5 * Float64(2*I + 1) * (-1)^(I + div(j_b + j_d,2) + 1) * f6j(j_a,j_b,2*J,j_c,j_d,2*I) / (sqrt((1 + kronecker_delta(a,b))*(1 + kronecker_delta(c,d))))
                                    @inbounds for nu in 1:N_ph
                                        ME = Amp *  real(Y_RPA[I+1,Q][nu][p_b,h_c] * conj(X_RPA[I+1,Q][nu][p_a,h_d]))
                                        nnSigmaSum += ME
                                    end

                                end
                            end
                            
                        elseif (nO_a + nO_b) == 2 && (nO_c + nO_d) == 0

                            # Direct term ...
                            Q = rem(l_a + l_c,2) + 1
                            if Q == (rem(l_b + l_d,2) + 1)
                                @inbounds for I in max(div(abs(j_a - j_c),2),div(abs(j_b - j_d),2)):min(div((j_a + j_c),2),div((j_b + j_d),2))

                                    h_a, h_b = orbitals_ph_list[a][2], orbitals_ph_list[b][2]
                                    p_c, p_d = orbitals_ph_list[c][2], orbitals_ph_list[d][2;]
                                    Amp = 0.5 * Float64(2*I + 1) * (-1)^(J + I + div(j_a + j_d,2)) * f6j(j_c,j_d,2*J,j_b,j_a,2*I) / (sqrt((1 + kronecker_delta(a,b))*(1 + kronecker_delta(c,d))))
                                    @inbounds for nu in 1:N_ph
                                        ME = Amp * real(Y_RPA[I+1,Q][nu][p_c,h_a] * conj(X_RPA[I+1,Q][nu][p_d,h_b]))
                                        nnSigmaSum += ME
                                    end

                                end
                            end

                            # Exchange term ...
                            Q = rem(l_a + l_d,2) + 1
                            if Q == (rem(l_b + l_c,2) + 1)
                                @inbounds for I in max(div(abs(j_a - j_d),2),div(abs(j_b - j_c),2)):min(div((j_a + j_d),2),div((j_b + j_c),2))

                                    h_a, h_b = orbitals_ph_list[a][2], orbitals_ph_list[b][2]
                                    p_c, p_d = orbitals_ph_list[c][2], orbitals_ph_list[d][2]
                                    Amp = 0.5 * Float64(2*I + 1) * (-1)^(I + div(j_b + j_d,2) + 1) * f6j(j_c,j_d,2*J,j_a,j_b,2*I) / (sqrt((1 + kronecker_delta(a,b))*(1 + kronecker_delta(c,d))))
                                    @inbounds for nu in 1:N_ph
                                        ME = Amp *  real(Y_RPA[I+1,Q][nu][p_d,h_a] * conj(X_RPA[I+1,Q][nu][p_c,h_b]))
                                        nnSigmaSum += ME
                                    end

                                end
                            end

                        end
                        =#

                    
                        #=
                    # 1p-1h channel ... does not contribute ???
                        # proton-proton terms ...
                        if (pO_a + pO_b) == 1 && (pO_c + pO_d) == 1

                            Phase = 1.0

                            if pO_a == 0
                                p, l_p, j_p = a, l_a, j_a
                                h, l_h, j_h = b, l_b, j_b
                            else
                                p, l_p, j_p = b, l_b, j_b
                                h, l_h, j_h = a, l_a, j_a
                                Phase = (-1)^(div(j_a + j_b,2) - J + 1) * Phase
                            end

                            if pO_c == 0
                                q, l_q, j_q = c, l_c, j_c
                                g, l_g, j_g = d, l_d, j_d
                            else
                                q, l_q, j_q = d, l_d, j_d
                                g, l_g, j_g = c, l_c, j_c
                                Phase = (-1)^(div(j_c + j_d,2) - J + 1) * Phase
                            end

                            Q = rem(l_p + l_g,2) + 1
                            if Q == (rem(l_q + l_h,2) + 1)
                                @inbounds for I in max(div(abs(j_p - j_g),2),div(abs(j_q - j_h),2)):min(div((j_p + j_g),2),div((j_q + j_h),2))
                                    if p != q && h != g
                                        Amp = (-1)^(div(j_h + j_q,2) + I + 1) * Phase * Float64(2*I + 1) * f6j(j_p,j_h,2*J,j_q,j_g,2*I)
                                        p_p, p_q = orbitals_ph_list[p][1], orbitals_ph_list[q][1]
                                        h_h, h_g = orbitals_ph_list[h][1], orbitals_ph_list[g][1]
                                        @inbounds for nu in 1:N_ph
                                            ME = Amp * real(Y_RPA[I+1,Q][nu][p_p,h_g] * conj(Y_RPA[I+1,Q][nu][p_q,h_h]))
                                            ppSigmaSum += ME
                                        end
                                    end
                                end
                            end


                        end

                        # neutron-neutron terms ...
                        if (nO_a + nO_b) == 1 && (nO_c + nO_d) == 1

                            Phase = 1.0

                            if nO_a == 0
                                p, l_p, j_p = a, l_a, j_a
                                h, l_h, j_h = b, l_b, j_b
                            else
                                p, l_p, j_p = b, l_b, j_b
                                h, l_h, j_h = a, l_a, j_a
                                Phase = (-1)^(div(j_a + j_b,2) - J + 1) * Phase
                            end

                            if nO_c == 0
                                q, l_q, j_q = c, l_c, j_c
                                g, l_g, j_g = d, l_d, j_d
                            else
                                q, l_q, j_q = d, l_d, j_d
                                g, l_g, j_g = c, l_c, j_c
                                Phase = (-1)^(div(j_c + j_d,2) - J + 1) * Phase
                            end

                            Q = rem(l_p + l_g,2) + 1
                            if Q == (rem(l_q + l_h,2) + 1)
                                @inbounds for I in max(div(abs(j_p - j_g),2),div(abs(j_q - j_h),2)):min(div((j_p + j_g),2),div((j_q + j_h),2))
                                    if p != q && h != g
                                        Amp = (-1)^(div(j_h + j_q,2) + I + 1) * Phase * Float64(2*I + 1) * f6j(j_p,j_h,2*J,j_q,j_g,2*I)
                                        p_p, p_q = orbitals_ph_list[p][2], orbitals_ph_list[q][2]
                                        h_h, h_g = orbitals_ph_list[h][2], orbitals_ph_list[g][2]
                                        @inbounds for nu in 1:N_ph
                                            ME = Amp * real(Y_RPA[I+1,Q][nu][p_p,h_g] * conj(Y_RPA[I+1,Q][nu][p_q,h_h]))
                                            nnSigmaSum += ME
                                        end
                                    end
                                end
                            end


                        end
                    =#

                    @inbounds Sigma_NN.pp[P,J+1][Ind] = ppSigmaSum
                    @inbounds Sigma_NN.nn[P,J+1][Ind] = nnSigmaSum

                end


            end
        end

    end


    println("\nSuccesfully allocate 2-body RPA correlation matrix Sigma_NN ...")
    
    return Sigma_NN
end

function V2b_correlation_energy(Params::Parameters,Orb::Vector{Orb1B},Orb_NN::Orb2B,Sigma::O2B,V_NN::O2B)
    # Read parameters ...
    N_max = Params.Calc.Nmax
    N_2max = Params.Calc.N2max
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Calculate the HF energy ...
    println("\nCalculating the correlation ground-state energy ...")

    E_corr, E_corr_threads = 0.0, zeros(Float64,Threads.maxthreadid())

    @inbounds Threads.@threads :static for a = 1:a_max
        Sum, Tid = 0.0, Threads.threadid()
        n_a = Orb[a].n
        l_a = Orb[a].l
        j_a = Orb[a].j
        @inbounds for b = 1:a_max
            n_b = Orb[b].n
            l_b = Orb[b].l
            P = rem(l_a + l_b, 2) + 1
            if (2*(n_a + n_b) + l_a + l_b) <= N_2max
                j_b = Orb[b].j
                @inbounds for c = 1:a_max
                    l_c = Orb[c].l
                    j_c = Orb[c].j
                    n_c = Orb[c].n
                    @inbounds for d = 1:a_max
                        n_d = Orb[d].n
                        l_d = Orb[d].l
                        j_d = Orb[d].j
                        if (2*(n_c + n_d) + l_c + l_d) <= N_2max
                            @inbounds for J in max(div(abs(j_a - j_b),2),div(abs(j_c - j_d),2)):min(div((j_a + j_b),2),div((j_c + j_d),2))
                                if (rem(l_a + l_b, 2) == rem(l_c + l_d, 2))
                                    Hat =  Float64(2*J + 1)
                                    Sum += 0.25 * Hat * O2b_pp(a,b,c,d,J,P,V_NN,Orb,Orb_NN) * O2b_pp(a,b,c,d,J,P,Sigma,Orb,Orb_NN)
                                    Sum += 0.25 * Hat * O2b_nn(a,b,c,d,J,P,V_NN,Orb,Orb_NN) * O2b_nn(a,b,c,d,J,P,Sigma,Orb,Orb_NN)
                                    Sum += Hat * O2b_pn(a,b,c,d,J,P,V_NN,Orb_NN) * O2b_pn(a,b,c,d,J,P,Sigma,Orb_NN)
                                end
                            end
                        end
                    end

                end
            end
        end
        E_corr_threads[Tid] += Sum
    end

    # Sum-up the total HF energy ...
    E_corr = sum(E_corr_threads)

    println("\tCorrelation energy ...   E_corr = " * string(E_corr) * " MeV")

    return E_corr
end
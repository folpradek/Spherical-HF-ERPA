function HF_BCS_Allocate(Params::Parameters,Rho::pnMatrix,Kappa::pnMatrix,Orb::Vector{NOrb},Orb_NN::NNOrb,Orb_NNN::NNNOrb,T::Matrix{Float64},VNN::NNInt,VNNN::Array{Vector{Vector{Float32}},4})
    # Read parameters ...
    A = Params.Calc.A
    N_max = Params.Calc.Nmax
    N_2max = Params.Calc.N2max
    N_3max = Params.Calc.N3max
    CMS = Params.Calc.CMS
    a_max = div((N_max + 1)*(N_max + 2),2)

    # BCS pairing rescale parameter s3 ...
    s3 = Params.Calc.Pairing.s3

    # Allocate  indices for allocation of H  ...
    ad, ad_count = HF_Allocate_Indices(Params,Orb)

    # Allocate new HF Hamiltonian matrices ...
    pH, nH = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)

    # Allocate the single-particle field H ...
    @inbounds Threads.@threads for ad_i in 1:ad_count
        a, d = ad[1,ad_i], ad[2,ad_i]
        n_a, l_a, j_a = Orb[a].n, Orb[a].l, Orb[a].j
        n_d, l_d, j_d = Orb[d].n, Orb[d].l, Orb[d].j

        is_j_a_hat = 1.0 / (Float64(j_a) + 1.0)

        # Single-particle field H ...
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

                    # Anomal Density-dependent part ... Only 3-body NNN interaction part ...
                    if (2*(n_d + n_e) + l_d + l_e) <= N_2max
                        A_NNN_Amp = 0.25 * sqrt(Float64((j_b + 1) * (j_e + 1))) * is_j_a_hat^2 * s3
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

        # Allocate pH & nH ...
        if a != d
            pH[a,d], pH[d,a] = pHSum, pHSum
            nH[a,d], nH[d,a] = nHSum, nHSum
        elseif a == d
            pH[a,a], nH[a,a] = pHSum, nHSum
        end

    end

    return pnMatrix(pH,nH)
end

function BCS_Allocate_Delta(Params::Parameters,U::pnVector,V::pnVector,Orb::Vector{NOrb},Orb_NN_res::NNOrb,VNN_res::NNInt)
    # Read parameters ...
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Initialize the pairing gap Delta ...
    pDelta, nDelta = zeros(Float64,a_max), zeros(Float64,a_max)

    # Allocate the pairing gap Delta ...
    @inbounds for a in 1:a_max
        j_a = Orb[a].j
        pSum, nSum = 0.0, 0.0
        @inbounds for b in 1:a_max
            j_b = Orb[b].j
            ja_jb_hat = sqrt((Float64(j_b) + 1.0) / (Float64(j_a) + 1.0))
            pME = - ja_jb_hat * V2B(a,a,b,b,0,1,VNN_res.pp,Orb,Orb_NN_res) * U.p[b] * V.p[b]
            nME = - ja_jb_hat * V2B(a,a,b,b,0,1,VNN_res.nn,Orb,Orb_NN_res) * U.n[b] * V.n[b]
            pSum += pME
            nSum += nME
        end
        pDelta[a] = pSum
        nDelta[a] = nSum
    end    

    return pDelta, nDelta
end

function BCS_SQE(a_max::Int64,SPE::pnVector,Lambda::pnFloat,Delta::pnVector)
    # Initialize vectors for SQEs ...
    pSQE, nSQE = zeros(Float64,a_max), zeros(Float64,a_max)

    # Calculate BCS single-quasiparticle energies (SQEs) ...
    println("\nEvalutiang BCS single-quasiparticle energies ...")
    @inbounds for a in 1:a_max
        pE_a = sqrt((SPE.p[a] - Lambda.p)^2 + Delta.p[a]^2)
        nE_a = sqrt((SPE.n[a] - Lambda.n)^2 + Delta.n[a]^2)
        pSQE[a], nSQE[a] = pE_a, nE_a
    end

    # Store SQEs ...
    SQE = pnVector(pSQE,nSQE)

    return SQE
end
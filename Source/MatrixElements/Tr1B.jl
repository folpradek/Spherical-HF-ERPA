function Tr1b_initialize(Params::Parameters,Orb::Vector{Orb1B})
    # Read parameters
    hw = Params.Calc.hw
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1)*(N_max+2),2)

    # Basic constants ...
    hc = 197.326980
    m_n = 939.565346
    m_p = 938.272013
    b_osc = sqrt(0.5 * (m_p + m_n) * hw) / hc

    println("\nConstructing electromagnetic transition operators E^lambda, M^lambda in the LHO basis ...")

    # Allocate E_lambda transition operators ...
    pE0 = E_lambda_initialize(0,b_osc,a_max,Orb)
    nE0 = E_lambda_initialize(0,b_osc,a_max,Orb)
    pE1 = E_lambda_initialize(1,b_osc,a_max,Orb)
    nE1 = E_lambda_initialize(1,b_osc,a_max,Orb)
    pE2 = E_lambda_initialize(2,b_osc,a_max,Orb)
    nE2 = E_lambda_initialize(2,b_osc,a_max,Orb)
    pE3 = E_lambda_initialize(3,b_osc,a_max,Orb)
    nE3 = E_lambda_initialize(3,b_osc,a_max,Orb)

    # Allocate M_lambda transition operators ...
    pM1 = M_lambda_initialize(1,-1,b_osc,a_max,Orb)
    nM1 = M_lambda_initialize(1,1,b_osc,a_max,Orb)
    pM2 = M_lambda_initialize(2,-1,b_osc,a_max,Orb)
    nM2 = M_lambda_initialize(2,1,b_osc,a_max,Orb)
    pM3 = M_lambda_initialize(3,-1,b_osc,a_max,Orb)
    nM3 = M_lambda_initialize(3,1,b_osc,a_max,Orb)
    
    TrE0 = O1B(pE0,nE0)
    TrE1 = O1B(pE1,nE1)
    TrE2 = O1B(pE2,nE2)
    TrE3 = O1B(pE3,nE3)
    TrM1 = O1B(pM1,nM1)
    TrM2 = O1B(pM2,nM2)
    TrM3 = O1B(pM3,nM3)

    TrOp = Tr1B(TrE0,TrE1,TrE2,TrE3,TrM1,TrM2,TrM3)

    println("\nTransition operators ready ...")

    return TrOp
end

function E_lambda_initialize(lambda::Int64,b_osc::Float64,a_max::Int64,Orb::Vector{Orb1B})
    E_lambda = zeros(Float64,a_max,a_max)
    @inbounds for a in 1:a_max
        l_a = Orb[a].l
        j_a = Orb[a].j
        n_a = Orb[a].n
        @inbounds for b in 1:a_max
            l_b = Orb[b].l
            j_b = Orb[b].j
            n_b = Orb[b].n
            ME = 0.0
            if rem(l_a + l_b + lambda,2) == 0
                Amp = Float64((-1)^(lambda + div(j_a - 1,2))) * sqrt(Float64((j_a + 1) * (j_b + 1))) * fCG(j_a,j_b,2*lambda,1,-1,0) / sqrt(4.0 * pi)
                ME = Amp * radial_moment_LHO(lambda,n_a,l_a,n_b,l_b,b_osc)
                if lambda == 0
                    Amp = Float64((-1)^(lambda + div(j_a - 1,2))) * sqrt(Float64((j_a + 1) * (j_b + 1))) * fCG(j_a,j_b,2*lambda,1,-1,0) / sqrt(4.0 * pi)
                    ME = Amp * radial_moment_LHO(lambda+2,n_a,l_a,n_b,l_b,b_osc)
                end
            end
            E_lambda[a,b] = ME
        end
    end
    return E_lambda
end

function M_lambda_initialize(lambda::Int64,t::Int64,b_osc::Float64,a_max::Int64,Orb::Vector{Orb1B})
    mu_N = 0.10515
    g_p = 5.586
    g_n = -3.826
    M_lambda = zeros(Float64,a_max,a_max)
    @inbounds for a in 1:a_max
        l_a = Orb[a].l
        j_a = Orb[a].j
        n_a = Orb[a].n
        @inbounds for b in 1:a
            l_b = Orb[b].l
            j_b = Orb[b].j
            n_b = Orb[b].n
            if rem(l_a + l_b + lambda + 1,2) == 0
                kappa = 0.5 * Float64((-1)^(div(j_a+1,2) + l_a) * (j_a + 1) + (-1)^(div(j_b+1,2) + l_b) * (j_b + 1))
                if t == 1
                    Amp = Float64(mu_N * (-1)^(lambda + div(j_a - 1,2))) * sqrt(Float64((j_a + 1) * (j_b + 1))) * fCG(j_a,j_b,2*lambda,1,-1,0) *
                            (lambda  - kappa) * (-0.5 * g_n) / sqrt(4.0 * pi)
                elseif t == -1
                    Amp = Float64(mu_N * (-1)^(lambda + div(j_a - 1,2))) * sqrt(Float64((j_a + 1) * (j_b + 1))) * fCG(j_a,j_b,2*lambda,1,-1,0) *
                            (lambda  - kappa) * (1.0 + kappa / (lambda + 1.0) - 0.5 * g_p) / sqrt(4.0 * pi)
                end
                ME = Amp * radial_moment_LHO(lambda,n_a,l_a,n_b,l_b,b_osc)
                M_lambda[a,b] = ME
                M_lambda[b,a] = ME
            end
        end
    end
    return M_lambda
end

function Tr1b_transformation(Params::Parameters,TrOp::Tr1B,Orb::Vector{Orb1B},C::O1B)
    # Parameter initialization...
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1)*(N_max + 2),2)

    pC, nC = C.p, C.n

    println("\nTransforming transition operators from the reference basis to the given target basis ...")

    pE0_new = zeros(Float64,a_max,a_max)
    nE0_new = zeros(Float64,a_max,a_max)
    pE1_new = zeros(Float64,a_max,a_max)
    nE1_new = zeros(Float64,a_max,a_max)
    pE2_new = zeros(Float64,a_max,a_max)
    nE2_new = zeros(Float64,a_max,a_max)
    pE3_new = zeros(Float64,a_max,a_max)
    nE3_new = zeros(Float64,a_max,a_max)

    pM1_new = zeros(Float64,a_max,a_max)
    nM1_new = zeros(Float64,a_max,a_max)
    pM2_new = zeros(Float64,a_max,a_max)
    nM2_new = zeros(Float64,a_max,a_max)
    pM3_new = zeros(Float64,a_max,a_max)
    nM3_new = zeros(Float64,a_max,a_max)

    @inbounds for a in 1:a_max
        l_a = Orb[a].l
        j_a = Orb[a].j
        @inbounds for b in 1:a_max
            l_b = Orb[b].l
            j_b = Orb[b].j

            pE0Sum = 0.0
            nE0Sum = 0.0
            pE1Sum = 0.0
            nE1Sum = 0.0
            pE2Sum = 0.0
            nE2Sum = 0.0
            pE3Sum = 0.0
            nE3Sum = 0.0

            pM1Sum = 0.0
            nM1Sum = 0.0
            pM2Sum = 0.0
            nM2Sum = 0.0
            pM3Sum = 0.0
            nM3Sum = 0.0

            @inbounds for k in 1:a_max
                l_k = Orb[k].l
                j_k = Orb[k].j
                if l_a == l_k && j_a == j_k
                    @inbounds for l in 1:a_max
                        l_l = Orb[l].l
                        j_l = Orb[l].j
                        if  l_b == l_l && j_b == j_l

                            # E0
                            if rem(l_a + l_b,2) == 0
                                pE0ME = TrOp.E0.p[k,l] * pC[k,a] * pC[l,b]
                                nE0ME = TrOp.E0.n[k,l] * nC[k,a] * nC[l,b]
                                pE0Sum += pE0ME
                                nE0Sum += nE0ME
                            end

                            # E1
                            if rem(l_a + l_b + 1,2) == 0
                                pE1ME = TrOp.E1.p[k,l] * pC[k,a] * pC[l,b]
                                nE1ME = TrOp.E1.n[k,l] * nC[k,a] * nC[l,b]
                                pE1Sum += pE1ME
                                nE1Sum += nE1ME
                            end

                            # E2
                            if rem(l_a + l_b + 2,2) == 0
                                pE2ME = TrOp.E2.p[k,l] * pC[k,a] * pC[l,b]
                                nE2ME = TrOp.E2.n[k,l] * nC[k,a] * nC[l,b]
                                pE2Sum += pE2ME
                                nE2Sum += nE2ME
                            end

                            # E3
                            if rem(l_a + l_b + 3,2) == 0
                                pE3ME = TrOp.E3.p[k,l] * pC[k,a] * pC[l,b]
                                nE3ME = TrOp.E3.n[k,l] * nC[k,a] * nC[l,b]
                                pE3Sum += pE3ME
                                nE3Sum += nE3ME
                            end

                            # M1
                            if rem(l_a + l_b + 2,2) == 0
                                pM1ME = TrOp.M1.p[k,l] * pC[k,a] * pC[l,b]
                                nM1ME = TrOp.M1.n[k,l] * nC[k,a] * nC[l,b]
                                pM1Sum += pM1ME
                                nM1Sum += nM1ME
                            end

                            # M2
                            if rem(l_a + l_b + 3,2) == 0
                                pM2ME = TrOp.M2.p[k,l] * pC[k,a] * pC[l,b]
                                nM2ME = TrOp.M2.n[k,l] * nC[k,a] * nC[l,b]
                                pM2Sum += pM2ME
                                nM2Sum += nM2ME
                            end

                            # M3
                            if rem(l_a + l_b + 4,2) == 0
                                pM3ME = TrOp.M3.p[k,l] * pC[k,a] * pC[l,b]
                                nM3ME = TrOp.M3.n[k,l] * nC[k,a] * nC[l,b]
                                pM3Sum += pM3ME
                                nM3Sum += nM3ME
                            end

                        end
                    end
                end
            end

            pE0_new[a,b] = pE0Sum
            nE0_new[a,b] = nE0Sum
            pE1_new[a,b] = pE1Sum
            nE1_new[a,b] = nE1Sum
            pE2_new[a,b] = pE2Sum
            nE2_new[a,b] = nE2Sum
            pE3_new[a,b] = pE3Sum
            nE3_new[a,b] = nE3Sum

            pM1_new[a,b] = pM1Sum
            nM1_new[a,b] = nM1Sum
            pM2_new[a,b] = pM2Sum
            nM2_new[a,b] = nM2Sum
            pM3_new[a,b] = pM3Sum
            nM3_new[a,b] = nM3Sum

        end
    end

    TrE0_new = O1B(pE0_new,nE0_new)
    TrE1_new = O1B(pE1_new,nE1_new)
    TrE2_new = O1B(pE2_new,nE2_new)
    TrE3_new = O1B(pE3_new,nE3_new)
    TrM1_new = O1B(pM1_new,nM1_new)
    TrM2_new = O1B(pM2_new,nM2_new)
    TrM3_new = O1B(pM3_new,nM3_new)

    TrOp_new = Tr1B(TrE0_new,TrE1_new,TrE2_new,TrE3_new,TrM1_new,TrM2_new,TrM3_new)

    println("\nTransition 1-body operators transformed to the given target basis ...")
    
    return TrOp_new
end



# Work in progress ...

struct qpTr1B
    E0::qpO1B
    E1::qpO1B
    E2::qpO1B
    E3::qpO1B
    M1::qpO1B
    M2::qpO1B
    M3::qpO1B
end

function qpTr1b_initialize(Params::Parameters,Orb::Vector{Orb1B},C::O1B,U::O1B,V::O1B)
    println("\nPreparing 1-body electromagnetic transition operators in quasiparticle representation ...")
    # Make the 1-body electromagnetic transitions operator in particle representation ...
    TrOp = Tr1b_initialize(Params,Orb)

    # Transform 1-body electromagnetic transitions operator in particle representation to the target particle basis ...
    TrOp = Tr1b_transformation(Params,TrOp,Orb,C)

    # Evaluate the quasiparticle representation of the 1-body electromagnetic transition operator components ...
        # E_lambda
    qpTrE0 = qpTr1b_E_lambda_initialize(Params,Orb,TrOp,U,V,0)
    qpTrE1 = qpTr1b_E_lambda_initialize(Params,Orb,TrOp,U,V,1)
    qpTrE2 = qpTr1b_E_lambda_initialize(Params,Orb,TrOp,U,V,2)
    qpTrE3 = qpTr1b_E_lambda_initialize(Params,Orb,TrOp,U,V,3)
        # M_lambda
    qpTrM1 = qpTr1b_M_lambda_initialize(Params,Orb,TrOp,U,V,1)
    qpTrM2 = qpTr1b_M_lambda_initialize(Params,Orb,TrOp,U,V,2)
    qpTrM3 = qpTr1b_M_lambda_initialize(Params,Orb,TrOp,U,V,3)


    qpTrOp = qpTr1B(qpTrE0,qpTrE1,qpTrE2,qpTrE3,qpTrM1,qpTrM2,qpTrM3)

    println("\nElectromagnetic 1-body transition operators in quasiparticle representation successfully prepared ...")

    return qpTrOp
end

function qpTr1b_E_lambda_initialize(Params::Parameters,Orb::Vector{Orb1B},TrOp::Tr1B,U::O1B,V::O1B,lambda::Int64)
    # Read parameters
        #hw = Params.Calc.hw
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1)*(N_max+2),2)

        # Basic constants ...
        #hc = 197.326980
        #m_n = 939.565346
        #m_p = 938.272013
        #b_osc = sqrt(0.5 * (m_p + m_n) * hw) / hc

    # Initialite the transition operators ...
    pE_lambda_11 = zeros(Float64,a_max,a_max)
    nE_lambda_11 = zeros(Float64,a_max,a_max)
    pE_lambda_20 = zeros(Float64,a_max,a_max)
    nE_lambda_20 = zeros(Float64,a_max,a_max)

    @inbounds for a in 1:a_max
        l_a , j_a = Orb[a].l, Orb[a].j
        @inbounds for b in 1:a_max
            l_b, j_b = Orb[b].l, Orb[b].j
            if div(abs(j_a - j_b),2) <= lambda && lambda <= div(j_a + j_b,2) && rem(l_a + l_b + lambda,2) == 0
                pM11Sum, pM20Sum = 0.0, 0.0
                nM11Sum, nM20Sum = 0.0, 0.0
                @inbounds for k in 1:a_max
                    l_k, j_k = Orb[k].l, Orb[k].j

                    if l_a == l_k && j_a == j_k

                        @inbounds for l in 1:a_max
                            l_l, j_l = Orb[l].l, Orb[l].j

                            if l_l == l_b && j_l == j_b
                                Phase = Float64((-1)^(div(j_a + j_b,2) + lambda))

                                pME11, nME11 = 0.0, 0.0
                                pME20, nME20 = 0.0, 0.0

                                # E0
                                if lambda == 0
                                    # 11 part ...
                                    pME11 = TrOp.E0.p[k,l] * (U.p[k,a] * U.p[l,b] - Phase * V.p[k,a] * V.p[l,b] + sqrt(Float64(j_a + 1)) * kronecker_delta(a,b) * V.p[k,a] * V.p[l,b])
                                    nME11 = TrOp.E0.n[k,l] * (U.n[k,a] * U.n[l,b] - Phase * V.n[k,a] * V.n[l,b] + sqrt(Float64(j_a + 1)) * kronecker_delta(a,b) * V.n[k,a] * V.n[l,b])

                                    # 20 part ...
                                    pME20 = TrOp.E0.p[k,l] * V.p[k,a] * U.p[l,b]
                                    nME20 = TrOp.E0.n[k,l] * V.n[k,a] * U.n[l,b]
                                # E1
                                elseif lambda == 1
                                    # 11 part ...
                                    pME11 = TrOp.E1.p[k,l] * (U.p[k,a] * U.p[l,b] - Phase * V.p[k,a] * V.p[l,b]) * sqrt(3.0)
                                    nME11 = TrOp.E1.n[k,l] * (U.n[k,a] * U.n[l,b] - Phase * V.n[k,a] * V.n[l,b]) * sqrt(3.0)

                                    # 20 part ...
                                    pME20 = TrOp.E1.p[k,l] * V.p[k,a] * U.p[l,b] * sqrt(3.0)
                                    nME20 = TrOp.E1.n[k,l] * V.n[k,a] * U.n[l,b] * sqrt(3.0)
                                # E2
                                elseif lambda == 2
                                    # 11 part ...
                                    pME11 = TrOp.E2.p[k,l] * (U.p[k,a] * U.p[l,b] - Phase * V.p[k,a] * V.p[l,b]) * sqrt(5.0)
                                    nME11 = TrOp.E2.n[k,l] * (U.n[k,a] * U.n[l,b] - Phase * V.n[k,a] * V.n[l,b]) * sqrt(5.0)

                                    # 20 part ...
                                    pME20 = TrOp.E2.p[k,l] * V.p[k,a] * U.p[l,b] * sqrt(5.0)
                                    nME20 = TrOp.E2.n[k,l] * V.n[k,a] * U.n[l,b] * sqrt(5.0)
                                # E3
                                elseif lambda == 3
                                    # 11 part ...
                                    pME11 = TrOp.E3.p[k,l] * (U.p[k,a] * U.p[l,b] - Phase * V.p[k,a] * V.p[l,b])  * sqrt(7.0)
                                    nME11 = TrOp.E3.n[k,l] * (U.n[k,a] * U.n[l,b] - Phase * V.n[k,a] * V.n[l,b])  * sqrt(7.0)

                                    # 20 part ...
                                    pME20 = TrOp.E3.p[k,l] * V.p[k,a] * U.p[l,b]  * sqrt(7.0)
                                    nME20 = TrOp.E3.n[k,l] * V.n[k,a] * U.n[l,b]  * sqrt(7.0)
                                end

                                # Update the sums ...
                                pM11Sum += pME11
                                nM11Sum += nME11
                                pM20Sum += pME20
                                nM20Sum += nME20
                            end
                        end
                    end
                end
                
                pE_lambda_11[a,b] = pM11Sum
                nE_lambda_11[a,b] = nM11Sum
                pE_lambda_20[a,b] = pM20Sum
                nE_lambda_20[a,b] = nM20Sum

            end
        end
    end

    return qpO1B(O1B(pE_lambda_11,pE_lambda_20),O1B(nE_lambda_11,nE_lambda_20))
end

function qpTr1b_M_lambda_initialize(Params::Parameters,Orb::Vector{Orb1B},TrOp::Tr1B,U::O1B,V::O1B,lambda::Int64)
    # Read parameters
        #hw = Params.Calc.hw
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1)*(N_max+2),2)

        # Basic constants ...
        #hc = 197.326980
        #m_n = 939.565346
        #m_p = 938.272013
        #b_osc = sqrt(0.5 * (m_p + m_n) * hw) / hc

    # Initialite the transition operators ...
    pM_lambda_11 = zeros(Float64,a_max,a_max)
    nM_lambda_11 = zeros(Float64,a_max,a_max)
    pM_lambda_20 = zeros(Float64,a_max,a_max)
    nM_lambda_20 = zeros(Float64,a_max,a_max)

    @inbounds for a in 1:a_max
        l_a , j_a = Orb[a].l, Orb[a].j
        @inbounds for b in 1:a_max
            l_b, j_b = Orb[b].l, Orb[b].j
            if div(abs(j_a - j_b),2) <= lambda && lambda <= div(j_a + j_b,2) && rem(l_a + l_b + lambda,2) == 1
                pM11Sum, pM20Sum = 0.0, 0.0
                nM11Sum, nM20Sum = 0.0, 0.0
                @inbounds for k in 1:a_max
                    l_k, j_k = Orb[k].l, Orb[k].j

                    if l_a == l_k && j_a == j_k

                        @inbounds for l in 1:a_max
                            l_l, j_l = Orb[l].l, Orb[l].j

                            if l_l == l_b && j_l == j_b
                                Phase = Float64((-1)^(div(j_a + j_b,2) + lambda))

                                pME11, nME11 = 0.0, 0.0
                                pME20, nME20 = 0.0, 0.0

                                # E0
                                if lambda == 0
                                    # 11 part ...
                                    pME11 = TrOp.E0.p[k,l] * (U.p[k,a] * U.p[l,b] - Phase * V.p[k,a] * V.p[l,b] - sqrt(Float64(j_a + 1)) * kronecker_delta(a,b) * V.p[k,a] * V.p[l,b])
                                    nME11 = TrOp.E0.n[k,l] * (U.n[k,a] * U.n[l,b] - Phase * V.n[k,a] * V.n[l,b] - sqrt(Float64(j_a + 1)) * kronecker_delta(a,b) * V.n[k,a] * V.n[l,b])

                                    # 20 part ...
                                    pME20 = TrOp.E0.p[k,l] * (V.p[k,a] * U.p[l,b] - Phase * U.p[k,a] * V.p[l,b])
                                    nME20 = TrOp.E0.n[k,l] * (V.n[k,a] * U.n[l,b] - Phase * U.n[k,a] * V.n[l,b])
                                # E1
                                elseif lambda == 1
                                    # 11 part ...
                                    pME11 = TrOp.E1.p[k,l] * (U.p[k,a] * U.p[l,b] - Phase * V.p[k,a] * V.p[l,b])
                                    nME11 = TrOp.E1.n[k,l] * (U.n[k,a] * U.n[l,b] - Phase * V.n[k,a] * V.n[l,b])

                                    # 20 part ...
                                    pME20 = TrOp.E1.p[k,l] * (V.p[k,a] * U.p[l,b] - Phase * U.p[k,a] * V.p[l,b])
                                    nME20 = TrOp.E1.n[k,l] * (V.n[k,a] * U.n[l,b] - Phase * U.n[k,a] * V.n[l,b])
                                # E2
                                elseif lambda == 2
                                    # 11 part ...
                                    pME11 = TrOp.E2.p[k,l] * (U.p[k,a] * U.p[l,b] - Phase * V.p[k,a] * V.p[l,b])
                                    nME11 = TrOp.E2.n[k,l] * (U.n[k,a] * U.n[l,b] - Phase * V.n[k,a] * V.n[l,b])

                                    pM11Sum += pME11
                                    nM11Sum += nME11

                                    # 20 part ...
                                    pME20 = TrOp.E2.p[k,l] * (V.p[k,a] * U.p[l,b] - Phase * U.p[k,a] * V.p[l,b])
                                    nME20 = TrOp.E2.n[k,l] * (V.n[k,a] * U.n[l,b] - Phase * U.n[k,a] * V.n[l,b])
                                # E3
                                elseif lambda == 3
                                    # 11 part ...
                                    pME11 = TrOp.E3.p[k,l] * (U.p[k,a] * U.p[l,b] - Phase * V.p[k,a] * V.p[l,b])
                                    nME11 = TrOp.E3.n[k,l] * (U.n[k,a] * U.n[l,b] - Phase * V.n[k,a] * V.n[l,b])

                                    # 20 part ...
                                    pME20 = TrOp.E3.p[k,l] * (V.p[k,a] * U.p[l,b] - Phase * U.p[k,a] * V.p[l,b])
                                    nME20 = TrOp.E3.n[k,l] * (V.n[k,a] * U.n[l,b] - Phase * U.n[k,a] * V.n[l,b])
                                end

                                # Update the sums ...
                                pM11Sum += pME11
                                nM11Sum += nME11
                                pM20Sum += pME20
                                nM20Sum += nME20

                            end
                        end
                    end
                end
                pM_lambda_11[a,b] = pM11Sum
                nM_lambda_11[a,b] = nM11Sum
                pM_lambda_20[a,b] = pM20Sum
                nM_lambda_20[a,b] = nM20Sum
            end
        end
    end

    return qpO1B(O1B(pM_lambda_11, pM_lambda_20),O1B(nM_lambda_11, nM_lambda_20))
end
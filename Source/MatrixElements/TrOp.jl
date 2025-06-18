function Transition_Operators_Ini(HbarOmega::Float64,N_max::Int64,Orb::Vector{NOrb})
    # Basic constants ...
    HbarC = 197.326980
    m_n = 939.565346
    m_p = 938.272013
    b_osc = sqrt(0.5 * (m_p + m_n) * HbarOmega) / HbarC
    a_max = div((N_max + 1)*(N_max+2),2)

    println("\nConstructing electromagnetic transition operators E^lambda, M^lambda ...")

    # Allocate Elambda transition operators ...
    pE0 = Elambda_Ini(0,b_osc,a_max,Orb)
    nE0 = Elambda_Ini(0,b_osc,a_max,Orb)
    pE1 = Elambda_Ini(1,b_osc,a_max,Orb)
    nE1 = Elambda_Ini(1,b_osc,a_max,Orb)
    pE2 = Elambda_Ini(2,b_osc,a_max,Orb)
    nE2 = Elambda_Ini(2,b_osc,a_max,Orb)
    pE3 = Elambda_Ini(3,b_osc,a_max,Orb)
    nE3 = Elambda_Ini(3,b_osc,a_max,Orb)

    # Allocate Mlambda transition operators ...
    pM1 = Mlambda_Ini(1,-1,b_osc,a_max,Orb)
    nM1 = Mlambda_Ini(1,1,b_osc,a_max,Orb)
    pM2 = Mlambda_Ini(2,-1,b_osc,a_max,Orb)
    nM2 = Mlambda_Ini(2,1,b_osc,a_max,Orb)
    pM3 = Mlambda_Ini(3,-1,b_osc,a_max,Orb)
    nM3 = Mlambda_Ini(3,1,b_osc,a_max,Orb)
    
    TrE0 = pnMatrix(pE0,nE0)
    TrE1 = pnMatrix(pE1,nE1)
    TrE2 = pnMatrix(pE2,nE2)
    TrE3 = pnMatrix(pE3,nE3)
    TrM1 = pnMatrix(pM1,nM1)
    TrM2 = pnMatrix(pM2,nM2)
    TrM3 = pnMatrix(pM3,nM3)

    TrOp = TranOper(TrE0,TrE1,TrE2,TrE3,TrM1,TrM2,TrM3)

    println("\nTransition operators ready ...")

    return TrOp
end

function Elambda_Ini(lambda::Int64,b_osc::Float64,a_max::Int64,Orb::Vector{NOrb})
    Elambda = zeros(Float64,a_max,a_max)
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
                ME = Amp * Radial_Moment_LHO(lambda,n_a,l_a,n_b,l_b,b_osc)
                if lambda == 0
                    Amp = Float64((-1)^(lambda + div(j_a - 1,2))) * sqrt(Float64((j_a + 1) * (j_b + 1))) * fCG(j_a,j_b,2*lambda,1,-1,0) / sqrt(4.0 * pi)
                    ME = Amp * Radial_Moment_LHO(lambda+2,n_a,l_a,n_b,l_b,b_osc)
                end
            end
            Elambda[a,b] = ME
        end
    end
    return Elambda
end

function Mlambda_Ini(lambda::Int64,t::Int64,b_osc::Float64,a_max::Int64,Orb::Vector{NOrb})
    mu_N = 0.10515
    g_p = 5.586
    g_n = -3.826
    Mlambda = zeros(Float64,a_max,a_max)
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
                ME = Amp * Radial_Moment_LHO(lambda,n_a,l_a,n_b,l_b,b_osc)
                Mlambda[a,b] = ME
                Mlambda[b,a] = ME
            end
        end
    end
    return Mlambda
end

function Transition_Operators_Transform(N_max::Int64,Input_File::String,TrOp::TranOper,Orb::Vector{NOrb})
    # Parameter initialization...
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Import HF basis transformation matrices ...
    pU_Import = Input_File * "/Bin/pU.bin"
    pU = Matrix{Float64}(undef,a_max,a_max)
    open(pU_Import, "r") do Read_File
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                ME = read(Read_File, Float64)
                pU[a,b] = ME
            end
        end
    end

    nU_Import = Input_File * "/Bin/nU.bin"
    nU = Matrix{Float64}(undef,a_max,a_max)
    open(nU_Import, "r") do Read_File
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                ME = read(Read_File, Float64)
                nU[a,b] = ME
            end
        end
    end

    println("\nTransforming transition operators to the HF basis ...")

    pE0_HF = zeros(Float64,a_max,a_max)
    nE0_HF = zeros(Float64,a_max,a_max)
    pE1_HF = zeros(Float64,a_max,a_max)
    nE1_HF = zeros(Float64,a_max,a_max)
    pE2_HF = zeros(Float64,a_max,a_max)
    nE2_HF = zeros(Float64,a_max,a_max)
    pE3_HF = zeros(Float64,a_max,a_max)
    nE3_HF = zeros(Float64,a_max,a_max)

    pM1_HF = zeros(Float64,a_max,a_max)
    nM1_HF = zeros(Float64,a_max,a_max)
    pM2_HF = zeros(Float64,a_max,a_max)
    nM2_HF = zeros(Float64,a_max,a_max)
    pM3_HF = zeros(Float64,a_max,a_max)
    nM3_HF = zeros(Float64,a_max,a_max)

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
                                pE0ME = TrOp.E0.p[k,l] * pU[k,a] * pU[l,b]
                                nE0ME = TrOp.E0.n[k,l] * nU[k,a] * nU[l,b]
                                pE0Sum += pE0ME
                                nE0Sum += nE0ME
                            end

                            # E1
                            if rem(l_a + l_b + 1,2) == 0
                                pE1ME = TrOp.E1.p[k,l] * pU[k,a] * pU[l,b]
                                nE1ME = TrOp.E1.n[k,l] * nU[k,a] * nU[l,b]
                                pE1Sum += pE1ME
                                nE1Sum += nE1ME
                            end

                            # E2
                            if rem(l_a + l_b + 2,2) == 0
                                pE2ME = TrOp.E2.p[k,l] * pU[k,a] * pU[l,b]
                                nE2ME = TrOp.E2.n[k,l] * nU[k,a] * nU[l,b]
                                pE2Sum += pE2ME
                                nE2Sum += nE2ME
                            end

                            # E3
                            if rem(l_a + l_b + 3,2) == 0
                                pE3ME = TrOp.E3.p[k,l] * pU[k,a] * pU[l,b]
                                nE3ME = TrOp.E3.n[k,l] * nU[k,a] * nU[l,b]
                                pE3Sum += pE3ME
                                nE3Sum += nE3ME
                            end

                            # M1
                            if rem(l_a + l_b + 2,2) == 0
                                pM1ME = TrOp.M1.p[k,l] * pU[k,a] * pU[l,b]
                                nM1ME = TrOp.M1.n[k,l] * nU[k,a] * nU[l,b]
                                pM1Sum += pM1ME
                                nM1Sum += nM1ME
                            end

                            # M2
                            if rem(l_a + l_b + 3,2) == 0
                                pM2ME = TrOp.M2.p[k,l] * pU[k,a] * pU[l,b]
                                nM2ME = TrOp.M2.n[k,l] * nU[k,a] * nU[l,b]
                                pM2Sum += pM2ME
                                nM2Sum += nM2ME
                            end

                            # M3
                            if rem(l_a + l_b + 4,2) == 0
                                pM3ME = TrOp.M3.p[k,l] * pU[k,a] * pU[l,b]
                                nM3ME = TrOp.M3.n[k,l] * nU[k,a] * nU[l,b]
                                pM3Sum += pM3ME
                                nM3Sum += nM3ME
                            end

                        end
                    end
                end
            end

            pE0_HF[a,b] = pE0Sum
            nE0_HF[a,b] = nE0Sum
            pE1_HF[a,b] = pE1Sum
            nE1_HF[a,b] = nE1Sum
            pE2_HF[a,b] = pE2Sum
            nE2_HF[a,b] = nE2Sum
            pE3_HF[a,b] = pE3Sum
            nE3_HF[a,b] = nE3Sum

            pM1_HF[a,b] = pM1Sum
            nM1_HF[a,b] = nM1Sum
            pM2_HF[a,b] = pM2Sum
            nM2_HF[a,b] = nM2Sum
            pM3_HF[a,b] = pM3Sum
            nM3_HF[a,b] = nM3Sum

        end
    end

    TrE0_HF = pnMatrix(pE0_HF,nE0_HF)
    TrE1_HF = pnMatrix(pE1_HF,nE1_HF)
    TrE2_HF = pnMatrix(pE2_HF,nE2_HF)
    TrE3_HF = pnMatrix(pE3_HF,nE3_HF)
    TrM1_HF = pnMatrix(pM1_HF,nM1_HF)
    TrM2_HF = pnMatrix(pM2_HF,nM2_HF)
    TrM3_HF = pnMatrix(pM3_HF,nM3_HF)

    TrOp_HF = TranOper(TrE0_HF,TrE1_HF,TrE2_HF,TrE3_HF,TrM1_HF,TrM2_HF,TrM3_HF)

    println("\nTransition operators transformed to the HF basis ...")
    
    return TrOp_HF
end

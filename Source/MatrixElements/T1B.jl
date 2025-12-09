function T1B(N_max::Int64,Orb::Vector{NOrb},HbarOmega::Float64)
    a_max = div((N_max + 1)*(N_max + 2),2)
    T = zeros(Float64,a_max,a_max)
    for k = 1:a_max
        n_k = Orb[k].n
        l_k = Orb[k].l
        j_k = Orb[k].j
        T[k,k] += 0.5 * HbarOmega * (2 * n_k + l_k + 1.5)
        for l = 1:a_max
            n_l = Orb[l].n
            l_l = Orb[l].l
            j_l = Orb[l].j
            if l_k == l_l
                if j_k == j_l
                    if n_k == (n_l + 1)
                        T[k,l] += 0.5 * HbarOmega * sqrt(n_k * (n_k + l_k + 0.5))
                    end
                    if (n_k + 1) == n_l
                        T[k,l] += 0.5 * HbarOmega * sqrt(n_l * (n_l + l_l + 0.5))
                    end
                end
            end
        end
    end
    return T
end

function Kinetic_Energy(Params::Parameters,Rho::pnMatrix,Orb::Vector{NOrb},T::Matrix{Float64})
    # Read parameters ...
    A = Params.Calc.A
    hw = Params.Calc.hw
    N_max = Params.Calc.Nmax
    N_2max = Params.Calc.N2max
    CMS = Params.Calc.CMS
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Calculate the total kinetic energy ...
    t_partial = Threads.Atomic{Float64}[Threads.Atomic{Float64}(0.0) for _ in 1:Threads.nthreads()]

    # Include the 2-body CM correction to the kinetic energy ...
    @inbounds Threads.@threads for a = 1:a_max
        thread_id = Threads.threadid()
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

                # 2-body CM correction to the kinetic energy ...
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

                                if rem(l_a + l_b, 2) == rem(l_d + l_e, 2)
                                    @inbounds for J = div(abs(j_a - j_b),2):div((j_a + j_b),2)
                                        T_NN_Amp = Float64(2*J + 1) / Float64(A)

                                        if CMS == "CMS1+2B"
                                            TNN_sym =  T2B(Orb,a,b,d,e,J) * hw
    
                                            TNN_antisym = 1.0 / sqrt(Float64((1 + KroneckerDelta(a,b))*(1 + KroneckerDelta(d,e)))) * (T2B(Orb,a,b,d,e,J) -
                                                          Float64((-1)^(round(div(j_d + j_e,2) - J))) * T2B(Orb,a,b,e,d,J)) * hw
                                        elseif CMS == "CMS2B"
                                            Amp = hw / sqrt(Float64(1 + KroneckerDelta(a,b)) * Float64(1 + KroneckerDelta(d,e)))
                                            Amp_2 = 1.0 / sqrt(Float64(1 + KroneckerDelta(a,b)) * Float64(1 + KroneckerDelta(d,e)))
    
                                            TNN_sym = hw * T2B(Orb,a,b,d,e,J) + T[a,d] * KroneckerDelta(b,e) + T[b,e] * KroneckerDelta(a,d)

                                            TNN_antisym = (Amp * T2B(Orb,a,b,d,e,J) + Amp_2 * (KroneckerDelta(b,e) * T[a,d] + KroneckerDelta(a,d) * T[b,e])
                                                        - Float64((-1)^(div(j_d + j_e,2) - J)) * (Amp * T2B(Orb,a,b,e,d,J) + Amp_2 * (Float64(KroneckerDelta(b,d)) *
                                                        T[a,e] + Float64(KroneckerDelta(a,e)) * T[b,d]))) / Float64(A)
                                        else
                                            TNN_sym = 0.0
                                            TNN_antisym = 0.0
                                        end

                                        @views t_partial[thread_id][] += T_NN_Amp * TNN_antisym * Rho.p[b,e]
                                        @views t_partial[thread_id][] += T_NN_Amp * TNN_sym * Rho.n[b,e]
                                        @views t_partial[thread_id][] += T_NN_Amp * TNN_antisym * Rho.n[b,e]
                                        @views t_partial[thread_id][] += T_NN_Amp * TNN_sym * Rho.p[b,e]

                                    end
                                end

                            end
                        end
                    end
                end

                # 1-body kinetic operator & CMS correction ...
                if CMS == "CMS1+2B"
                    @views t_partial[thread_id][] += T[a,d] * (1.0 - 1.0 / Float64(A)) * Rho.p[a,d] * Float64(Orb[a].j + 1)
                    @views t_partial[thread_id][] += T[a,d] * (1.0 - 1.0 / Float64(A)) * Rho.n[a,d] * Float64(Orb[a].j + 1)
                elseif CMS != "CMS2B"
                    @views t_partial[thread_id][] += T[a,d] * Rho.p[a,d]
                    @views t_partial[thread_id][] += T[a,d] * Rho.n[a,d]
                end

            end
        end
    end

    t = sum(x[] for x in t_partial)

    println("\nTotal kinetic energy reads:    <T> = " * string(t) * " MeV")

    return t
end

function T1B_Transform(HbarOmega::Float64,N_max::Int64,Input_File::String,Orb::Vector{NOrb})
    # Parameter initialization...
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Import HF basis transformation matrices ...
    pU_Import = Input_File * "/pU.bin"
    pU = Matrix{Float64}(undef,a_max,a_max)
    open(pU_Import, "r") do Read_File
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                ME = read(Read_File, Float64)
                pU[a,b] = ME
            end
        end
    end

    nU_Import = Input_File * "/nU.bin"
    nU = Matrix{Float64}(undef,a_max,a_max)
    open(nU_Import, "r") do Read_File
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                ME = read(Read_File, Float64)
                nU[a,b] = ME
            end
        end
    end

    println("\nTransforming 1-body kinetic operators to the HF basis ...")

    T = T1B(N_max,Orb,HbarOmega)

    pT_HF = zeros(Float64,a_max,a_max)
    nT_HF = zeros(Float64,a_max,a_max)

    @inbounds for a in 1:a_max
        l_a = Orb[a].l
        j_a = Orb[a].j
        @inbounds for b in 1:a_max
            l_b = Orb[b].l
            j_b = Orb[b].j

            pTSum = 0.0
            nTSum = 0.0

            @inbounds for k in 1:a_max
                l_k = Orb[k].l
                j_k = Orb[k].j
                if l_a == l_k && j_a == j_k
                    @inbounds for l in 1:a_max
                        l_l = Orb[l].l
                        j_l = Orb[l].j
                        if  l_b == l_l && j_b == j_l
                            pTME = T[k,l] * pU[k,a] * pU[l,b]
                            nTME = T[k,l] * nU[k,a] * nU[l,b]
                            pTSum += pTME
                            nTSum += nTME
                        end
                    end
                end
            end
            pT_HF[a,b] = pTSum
            nT_HF[a,b] = nTSum
        end
    end

    println("\n1-body kinetic operators transformed to the HF basis ...")

    return pT_HF, nT_HF
end
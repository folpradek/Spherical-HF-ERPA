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
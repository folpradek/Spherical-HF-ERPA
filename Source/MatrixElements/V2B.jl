function V2B_Ini(Orb::Vector{NOrb},N_max::Int64)
    N_2max = 2 * N_max
    a_max = div((N_max + 1)*(N_max + 2),2)
    J_max = Int64(2*N_max + 1)

    N_Orb_NN = zeros(Int64,2,2,J_max+1)

    Orb_NN_Dic = Dict{Tuple{Int8, Int8, Int8, Int16, Int16}, Int32}()

    for a in 1:a_max
        n_a = Orb[a].n
        l_a = Orb[a].l
        j_a = Orb[a].j
        if (2*n_a + l_a) <= N_max
            for b = 1:a_max
                n_b = Orb[b].n
                l_b = Orb[b].l
                j_b = Orb[b].j
                if ((2*(n_a + n_b) + l_a + l_b) <= N_2max)
                    P = rem(l_a + l_b, 2) + 1
                    for J = div(abs(j_a - j_b),2):div((j_a + j_b),2)
                        if b <= a
                            N_Orb_NN[2,P,J+1] += 1
                        end
                        N_Orb_NN[1,P,J+1] += 1
                    end
                end
            end
        end
    end

    Ind_Orb_NN = Array{Vector{Vector{Int64}}}(undef,2,2,J_max+1)
    for P = 1:2
        for J = 0:J_max
            Ind_Orb_NN[1,P,J+1] = Vector{Vector{Int64}}(undef,N_Orb_NN[1,P,J+1])
            Ind_Orb_NN[2,P,J+1] = Vector{Vector{Int64}}(undef,N_Orb_NN[2,P,J+1])
        end
    end

    N_Orb_NN = zeros(Int64,2,2,J_max+1)

    for a in 1:a_max
        n_a = Orb[a].n
        l_a = Orb[a].l
        j_a = Orb[a].j
        if (2*n_a + l_a) <= N_max
            for b = 1:a_max
                n_b = Orb[b].n
                l_b = Orb[b].l
                j_b = Orb[b].j
                if ((2*(n_a + n_b) + l_a + l_b) <= N_2max)
                    P = rem(l_a + l_b, 2) + 1
                    for J = Int64(abs(j_a - j_b)/2):Int64((j_a + j_b)/2)
                        if b <= a
                            key_1 = (Int8(1),Int8(P),Int8(J),Int16(a),Int16(b))
                            N_Orb_NN[2,P,J+1] += 1
                            Orb_NN_Dic[key_1] = Int32(N_Orb_NN[2,P,J+1])
                            Ind_Orb_NN[2,P,J+1][N_Orb_NN[2,P,J+1]] = zeros(Int64, 2)
                            Ind_Orb_NN[2,P,J+1][N_Orb_NN[2,P,J+1]][1] = a
                            Ind_Orb_NN[2,P,J+1][N_Orb_NN[2,P,J+1]][2] = b
                        end
                        key_0 = (Int8(0),Int8(P),Int8(J),Int16(a),Int16(b))
                        N_Orb_NN[1,P,J+1] += 1
                        Orb_NN_Dic[key_0] = Int32(N_Orb_NN[1,P,J+1])
                        Ind_Orb_NN[1,P,J+1][N_Orb_NN[1,P,J+1]] = zeros(Int64, 2)
                        Ind_Orb_NN[1,P,J+1][N_Orb_NN[1,P,J+1]][1] = a
                        Ind_Orb_NN[1,P,J+1][N_Orb_NN[1,P,J+1]][2] = b
                    end
                end
            end
        end
    end

    Vpp = Matrix{Vector{Float64}}(undef,2,J_max+1)
    Vpn = Matrix{Vector{Float64}}(undef,2,J_max+1)
    Vnn = Matrix{Vector{Float64}}(undef,2,J_max+1)
    
    for P = 1:2
        for J=0:J_max
            Length1 = div(N_Orb_NN[1,P,J+1]*(N_Orb_NN[1,P,J+1]+1),2)
            Length2 = div(N_Orb_NN[2,P,J+1]*(N_Orb_NN[2,P,J+1]+1),2)
            Vpp[P,J+1] = zeros(Float64, Length2)
            Vpn[P,J+1] = zeros(Float64, Length1)
            Vnn[P,J+1] = zeros(Float64, Length2)
        end
    end

    VNN = NNInt(Vpp, Vpn, Vnn)

    Orb_NN = NNOrb(Orb_NN_Dic, N_Orb_NN, Ind_Orb_NN)

    return VNN, Orb_NN

end

function V2B_Index(a::Int64,b::Int64,c::Int64,d::Int64,J::Int64,P::Int64,t::Int64,Orb_NN::NNOrb)
    @views N = Orb_NN.N[t+1,P,J+1]
    t = Int8(t)
    P = Int8(P)
    J = Int8(J)
    a = Int16(a)
    b = Int16(b)
    c = Int16(c)
    d = Int16(d)
    Bra = Int64(Orb_NN.Dic[(t,P,J,a,b)])
    Ket = Int64(Orb_NN.Dic[(t,P,J,c,d)])
    if Bra < Ket
        Ket, Bra = Bra, Ket
    end
    Ind = Bra + (Ket - 1) * N - div(Ket * (Ket - 1),2)
    return Ind
end

function V2B(a::Int64,b::Int64,c::Int64,d::Int64,J::Int64,t::Int64,V::Matrix{Vector{Float64}},Orb::Vector{NOrb},Orb_NN::NNOrb)
    P = rem(Orb[a].l + Orb[b].l, 2) + 1

    if t == 1

        if a >= b && c >= d
            @inbounds Ind = V2B_Index(a,b,c,d,J,P,1,Orb_NN)
            @views v = V[P,J+1][Ind]
            return v
        elseif a < b && c >= d
            Amp = Float64((-1)^(J + 1 + div(Orb[a].j + Orb[b].j,2)))
            @inbounds Ind = V2B_Index(b,a,c,d,J,P,1,Orb_NN)
            @views v = Amp * V[P,J+1][Ind]
            return v
        elseif a >= b && c < d
            Amp = Float64((-1)^(J + 1 + div(Orb[c].j + Orb[d].j,2)))
            @inbounds Ind = V2B_Index(a,b,d,c,J,P,1,Orb_NN)
            @views v = Amp * V[P,J+1][Ind]
            return v
        elseif a < b && c < d
            Amp = Float64((-1)^(div(Orb[a].j + Orb[b].j + Orb[c].j + Orb[d].j,2)))
            @inbounds Ind = V2B_Index(b,a,d,c,J,P,1,Orb_NN)
            @views v = Amp * V[P,J+1][Ind]
            return v
        end

    elseif t == 0
        @inbounds Ind = V2B_Index(a,b,c,d,J,P,0,Orb_NN)
        @views v = V[P,J+1][Ind]
        return v
    end

end

function V2B_Read(Int_File::String,Params::Vector{Any},Orb::Vector{NOrb})
    # Read interactions params ...
    HbarOmega = Params[1]
    N_max = Params[2]
    N_2max = Params[3]
    A = Params[5]
    CMS = Params[10]

    a_max = div((N_max+1)*(N_max+2),2)
    J_max = N_2max + 1

    # Initialize 1-body kinetic operator ...
    T = T1B(N_max,Orb,HbarOmega)

    # Initialize NN interaction array ...
    println("\nPreparing 2-body NN interaction array...")
    VNN, Orb_NN = V2B_Ini(Orb,N_max)

    # Precalculate chunks of bytes to skip ...
    N_Chunk_skip = V2B_Count(N_max,N_2max,Orb)

    # Read NN interaction matrix elements ...
    println("\nReading 2-body NN interaction file...")
    @inbounds Threads.@threads for a in 1:a_max
        n_a = Orb[a].n
        l_a = Orb[a].l
        j_a = Orb[a].j
        open(Int_File, "r") do Bin_Read
            seek(Bin_Read, N_Chunk_skip[a])
            @inbounds for b in 1:a
                n_b = Orb[b].n
                l_b = Orb[b].l
                j_b = Orb[b].j
                P = rem(l_a + l_b,2) + 1
                if (2*(n_a + n_b) + l_a + l_b) <= N_2max
                    @inbounds for c in 1:a
                        n_c = Orb[c].n
                        l_c = Orb[c].l
                        j_c = Orb[c].j
                        d_max = c

                        if a == c
                            d_max = b
                        end

                        @inbounds for d in 1:d_max
                            n_d = Orb[d].n
                            l_d = Orb[d].l
                            j_d = Orb[d].j
                            if ((2*(n_c + n_d) + l_c + l_d) <= N_2max) && (rem(l_a + l_b,2) == rem(l_c + l_d,2))
                                @inbounds for J in div(max(abs(j_a - j_b), abs(j_c - j_d)),2):div(min(j_a + j_b, j_c + j_d),2)
                                    ME_00 = read(Bin_Read,Float64)
                                    ME_nn = read(Bin_Read,Float64)
                                    ME_10 = read(Bin_Read,Float64)
                                    ME_pp = read(Bin_Read,Float64)

                                    Amp = 1.0
                                    if a == b
                                        Amp = Amp / sqrt(2)
                                    end
                                    if c == d
                                        Amp = Amp / sqrt(2)
                                    end

                                    if (abs(Amp - 1.0) < 1e-4) || (rem(J,2) == 0)
                                        Ind = V2B_Index(a,b,c,d,J,P,1,Orb_NN)
                                        VNN.pp[P,J+1][Ind] += ME_pp
                                        VNN.nn[P,J+1][Ind] += ME_nn
                                    end

                                    if (abs(Amp - 1.0) < 1e-4) || (rem(J,2) == 0)
                                        Ind = V2B_Index(a,b,c,d,J,P,0,Orb_NN)
                                        VNN.pn[P,J+1][Ind] += 0.5 * ME_10

                                        if c != d
                                            Phase = (-1)^(div(j_c + j_d,2) + J + 1)
                                            Ind = V2B_Index(a,b,d,c,J,P,0,Orb_NN)
                                            VNN.pn[P,J+1][Ind] += 0.5 * Phase * ME_10
                                        end

                                        if (a != b) && ( c!= d)
                                            Phase = (-1)^(div(j_a + j_b + j_c + j_d,2))
                                            Ind = V2B_Index(b,a,d,c,J,P,0,Orb_NN)
                                            VNN.pn[P,J+1][Ind] += 0.5 * Phase * ME_10
                                        end

                                        if (a != b) && ((a != c) || (b != d))
                                            Phase = (-1)^(div(j_a + j_b,2) + J + 1)
                                            Ind = V2B_Index(b,a,c,d,J,P,0,Orb_NN)
                                            VNN.pn[P,J+1][Ind] += 0.5 * Phase * ME_10
                                        end


                                    end

                                    if (abs(Amp - 1.0) < 1e-4) || (rem(J,2) == 1)
                                        Ind = V2B_Index(a,b,c,d,J,P,0,Orb_NN)
                                        VNN.pn[P,J+1][Ind] += 0.5 * ME_00

                                        if c != d
                                            Phase = (-1)^(div(j_c + j_d,2) + J)
                                            Ind = V2B_Index(a,b,d,c,J,P,0,Orb_NN)
                                            VNN.pn[P,J+1][Ind] += 0.5 * Phase * ME_00
                                        end

                                        if (a != b) && (c != d)
                                            Phase = (-1)^(div(j_a + j_b + j_c + j_d,2))
                                            Ind = V2B_Index(b,a,d,c,J,P,0,Orb_NN)
                                            VNN.pn[P,J+1][Ind] += 0.5 * Phase * ME_00
                                        end

                                        if (a != b) && ((a != c) || (b != d))
                                            Phase = (-1)^(div(j_a + j_b,2) + J)
                                            Ind = V2B_Index(b,a,c,d,J,P,0,Orb_NN)
                                            VNN.pn[P,J+1][Ind] += 0.5 * Phase * ME_00
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

    println("\nNN 2-body interaction loaded...")

    if CMS == "CMS1+2B" || CMS == "CMS2B"
        println("\nLoading 2-body CMS correction...")
        JP_List = JP_Ini(J_max)
        @inbounds Threads.@threads for JP in JP_List
            J = JP[1]
            P = JP[2]
            N_t0 = Orb_NN.N[1,P,J+1]
            N_t1 = Orb_NN.N[2,P,J+1]
            @inbounds for Bra in 1:max(N_t0,N_t1)
                @inbounds for Ket in 1:Bra

                    if Bra <= N_t0
                        Ind = Bra + (Ket - 1) * N_t0 - div(Ket * (Ket - 1),2)
                        a = Orb_NN.Ind[1,P,J+1][Bra][1]
                        b = Orb_NN.Ind[1,P,J+1][Bra][2]
                        c = Orb_NN.Ind[1,P,J+1][Ket][1]
                        d = Orb_NN.Ind[1,P,J+1][Ket][2]
                        if CMS == "CMS1+2B"
                            T_sym = T2B(Orb,a,b,c,d,J) * HbarOmega / Float64(A)
                        elseif CMS == "CMS2B"
                            T_sym = (HbarOmega * T2B(Orb,a,b,c,d,J) + T[a,c] * Float64(KroneckerDelta(b,d)) + T[b,d] * Float64(KroneckerDelta(a,c))) / Float64(A)
                        end
                        VNN.pn[P,J+1][Ind] += T_sym
                    end

                    if Bra <= N_t1
                        Ind = Bra + (Ket - 1) * N_t1 - div(Ket * (Ket - 1),2)
                        a = Orb_NN.Ind[2,P,J+1][Bra][1]
                        b = Orb_NN.Ind[2,P,J+1][Bra][2]
                        c = Orb_NN.Ind[2,P,J+1][Ket][1]
                        d = Orb_NN.Ind[2,P,J+1][Ket][2]
                        j_c = Orb[c].j
                        j_d = Orb[d].j
                        if CMS == "CMS1+2B"
                            j_c = Orb[c].j
                            j_d = Orb[d].j
                            T_antisym = 1.0 / sqrt(Float64(1 + KroneckerDelta(a,b)) * Float64(1 + KroneckerDelta(c,d))) * HbarOmega * 
                                        (T2B(Orb,a,b,c,d,J) - Float64((-1)^(div(j_c + j_d,2) - J)) * T2B(Orb,a,b,d,c,J)) / Float64(A)
                        elseif CMS == "CMS2B"
                            j_c = Orb[c].j
                            j_d = Orb[d].j
                            Amp = HbarOmega / sqrt(Float64(1 + KroneckerDelta(a,b)) * Float64(1 + KroneckerDelta(c,d)))
                            Amp_2 = 1.0 / sqrt(Float64(1 + KroneckerDelta(a,b)) * Float64(1 + KroneckerDelta(c,d)))
                            T_antisym = (Amp * T2B(Orb,a,b,c,d,J) + Amp_2 * (KroneckerDelta(b,d) * T[a,c] + KroneckerDelta(a,c) * T[b,d])
                                        - Float64((-1)^(div(j_c + j_d,2) - J)) * (Amp * T2B(Orb,a,b,d,c,J) + Amp_2 * (Float64(KroneckerDelta(b,c)) *
                                        T[a,d] + Float64(KroneckerDelta(a,d)) * T[b,c]))) / Float64(A)
                        end
                        ME =  T_antisym * sqrt(Float64(1 + KroneckerDelta(a,b)) * Float64(1 + KroneckerDelta(c,d)))
                        VNN.pp[P,J+1][Ind] += ME
                        VNN.nn[P,J+1][Ind] += ME
                    end

                end
            end
        end
        println("\nCMS 2-body correction added...")
    end

    return VNN, Orb_NN
end

function V2B_Count(N_max::Int64,N_2max::Int64,Orb::Vector{NOrb})
    a_max = div((N_max+1)*(N_max+2),2)
    N_Chunk = Vector{Int64}(undef,a_max)
    N_Chunk_skip = Vector{Int64}(undef,a_max)

    @inbounds Threads.@threads for a in 1:a_max
        n_a = Orb[a].n
        l_a = Orb[a].l
        j_a = Orb[a].j
        N = 0
        @inbounds for b in 1:a
            n_b = Orb[b].n
            l_b = Orb[b].l
            j_b = Orb[b].j
            if (2*(n_a + n_b) + l_a + l_b) <= N_2max
                @inbounds for c in 1:a
                    n_c = Orb[c].n
                    l_c = Orb[c].l
                    j_c = Orb[c].j
                    d_max = c

                    if a == c
                        d_max = b
                    end

                    @inbounds for d in 1:d_max
                        n_d = Orb[d].n
                        l_d = Orb[d].l
                        j_d = Orb[d].j
                        if ((2*(n_c + n_d) + l_c + l_d) <= N_2max) && (rem(l_a + l_b,2) == rem(l_c + l_d,2))
                            @inbounds for J in div(max(abs(j_a - j_b), abs(j_c - j_d)),2):div(min(j_a + j_b, j_c + j_d),2)
                                N += 4
                            end
                        end
                    end
                end
            end
        end
        N_Chunk[a] = N
    end

    @inbounds for a in 1:a_max
        Sum = 0
        @inbounds for b in 1:(a-1)
            Sum += 8 * N_Chunk[b]
        end
        N_Chunk_skip[a] = Sum
    end

    return N_Chunk_skip
end
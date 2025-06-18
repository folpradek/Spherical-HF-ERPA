function V3B_NO2B_Ini(N_max::Int64,N_2max::Int64,N_3max::Int64,Orb::Vector{NOrb})
    a_max = div((N_max + 1)*(N_max + 2),2)
    l_max = N_max
    J_max = N_2max + 1

    N_Orb_NNN = Array{Vector}(undef,2,l_max+1,J_max+1,2)
    for P = 1:2
        for l = 0:l_max
            for J = 0:J_max
                for T = 1:2
                    j_count = 0
                    for j in Int64(abs(2*l-1)):2:Int64(2*l+1)
                        j_count += 1
                    end
                    N_Orb_NNN[P,l+1,J+1,T] = zeros(Int64,j_count)
                end
            end
        end
    end

    Orb_NNN_Dic = Dict()

    for a = 1:a_max
        n_a = Orb[a].n
        l_a = Orb[a].l
        j_a = Orb[a].j
        if (2*n_a + l_a) <= N_max
            for b = 1:a
                n_b = Orb[b].n
                l_b = Orb[b].l
                j_b = Orb[b].j
                if (2*(n_a + n_b) + l_a + l_b) <= N_2max
                    for c = 1:a_max
                        n_c = Orb[c].n
                        l_c = Orb[c].l
                        j_c = Orb[c].j
                        P = rem(l_a + l_b + l_c,2) + 1
                        if (2*(n_a + n_b + n_c) + l_a + l_b + l_c) <= N_3max &&
                           (2*(n_a + n_c) + l_a + l_c) <= N_2max &&
                           (2*(n_b + n_c) + l_b + l_c) <= N_2max
                           for J = Int64(abs(j_a - j_b)/2):Int64((j_a + j_b)/2) 
                                for T in 1:2
                                    for T_ab in (T-1):1
                                        key = (Int8(P),Int8(J),Int8(2*T-1),Int8(T_ab),Int16(a),Int16(b),Int8(l_c),Int8(j_c),Int8(n_c))
                                        if j_c == Int64(abs(2*l_c-1))
                                            N_Orb_NNN[P,l_c+1,J+1,T][1] += 1
                                            Orb_NNN_Dic[key] = Int32(N_Orb_NNN[P,l_c+1,J+1,T][1])
                                        elseif j_c == Int64(2*l_c+1)
                                            N_Orb_NNN[P,l_c+1,J+1,T][2] += 1
                                            Orb_NNN_Dic[key] = Int32(N_Orb_NNN[P,l_c+1,J+1,T][2])
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

    VNNN = Array{Vector{Vector{Float32}}}(undef,2,l_max+1,J_max+1,2)

    for P = 1:2
        for l = 0:l_max
            for J = 0:J_max
                for T = 1:2

                    j_count = 0
                    for j in Int64(abs(2*l-1)):2:Int64(2*l+1)
                        j_count += 1
                    end

                    VNNN[P,l+1,J+1,T] = Vector{Vector{Float32}}(undef,j_count)

                    j_count = 0
                    for j in Int64(abs(2*l-1)):2:Int64(2*l+1)
                        j_count += 1
                        
                        N_count = N_Orb_NNN[P,l+1,J+1,T][j_count]
                        Length = div(N_count*(N_count+1),2)
                        VNNN[P,l+1,J+1,T][j_count] = zeros(Float32, Length)
                    end

                end
            end
        end
    end

    Orb_NNN = NNNOrb(Orb_NNN_Dic, N_Orb_NNN)

    return VNNN, Orb_NNN

end

function V3B_NO2B_Index(a::Int64,b::Int64,c::Int64,T_ab::Int64,d::Int64,e::Int64,f::Int64,T_de::Int64,J::Int64,T::Int64,P::Int64,Orb::Vector{NOrb},Orb_NNN::NNNOrb)
    n_c = Int8(Orb[c].n)
    l_c = Int8(Orb[c].l)
    j_c = Int8(Orb[c].j)
    n_f = Int8(Orb[f].n) 
    Bra = Int64(Orb_NNN.Dic[(P,J,T,T_ab,a,b,l_c,j_c,n_c)])
    Ket = Int64(Orb_NNN.Dic[(P,J,T,T_de,d,e,l_c,j_c,n_f)])
    if Bra < Ket
        Bra, Ket = Ket, Bra
    end
    T_ind = div(T+1,2)
    if j_c == abs(2*l_c-1)
        @views N = Orb_NNN.N[P,l_c+1,J+1,T_ind][1]
    elseif j_c == (2*l_c+1)
        @views N = Orb_NNN.N[P,l_c+1,J+1,T_ind][2]
    end
    Ind = Bra + (Ket - 1)*N - div(Ket * (Ket - 1),2)
    return Ind
end

function V3B_NO2B(a::Int64,b::Int64,c::Int64,T_ab::Int64,d::Int64,e::Int64,f::Int64,T_de::Int64,J_de::Int64,T::Int64,P::Int64,VNNN::Array{Vector{Vector{Float32}},4},Orb::Vector{NOrb},Orb_NNN::NNNOrb)
    if a >= b
        if d >= e
            @inbounds Ind = V3B_NO2B_Index(a,b,c,T_ab,d,e,f,T_de,J_de,T,P,Orb,Orb_NNN)
            phase = 1.0
        else
            j_d = Orb[d].j
            j_e = Orb[e].j
            @inbounds Ind = V3B_NO2B_Index(a,b,c,T_ab,e,d,f,T_de,J_de,T,P,Orb,Orb_NNN)
            phase = Float64((-1)^(div(j_d + j_e,2) - J_de - T_de))
        end
    else
        if d >= e
            j_a = Orb[a].j
            j_b = Orb[b].j
            @inbounds Ind = V3B_NO2B_Index(b,a,c,T_ab,d,e,f,T_de,J_de,T,P,Orb,Orb_NNN)
            phase = Float64((-1)^(div(j_a + j_b,2) - J_de - T_ab))
        else
            j_a = Orb[a].j
            j_b = Orb[b].j
            j_d = Orb[d].j
            j_e = Orb[e].j
            @inbounds Ind = V3B_NO2B_Index(b,a,c,T_ab,e,d,f,T_de,J_de,T,P,Orb,Orb_NNN)
            phase = Float64((-1)^(div(j_a + j_b + j_d + j_e,2) - T_de - T_ab))
        end
    end
    T_ind = div((T+1),2)
    l_c = Orb[c].l
    j_c = Orb[c].j
    if j_c == abs(2*l_c-1)
        @views V = phase * Float64(VNNN[P,l_c+1,J_de+1,T_ind][1][Ind])
        return V
    elseif j_c == (2*l_c+1)
        @views V = phase * Float64(VNNN[P,l_c+1,J_de+1,T_ind][2][Ind])
        return V
    end
end

function V3B_NO2B_Read(Int_File::String,Params::Vector{Any},Orb::Vector{NOrb})
    N_max = Params[2]
    N_2max = Params[3]
    N_3max = Params[4]
    a_max = div((N_max+1)*(N_max+2),2)

    println("\nPreparing 3-body NNN interaction array ...")
    VNNN, Orb_NNN = V3B_NO2B_Ini(N_max,N_2max,N_3max,Orb)

    println("\nPrecounting # of 3-body NNN interaction matrix elements ...")
    N_Chunk_skip = V3B_NO2B_Count(N_max,N_2max,N_3max,Orb)

    println("\nReading 3-body NO2B NNN interaction file. ..")
    @inbounds Threads.@threads for a in 1:a_max
        n_a = Orb[a].n
        l_a = Orb[a].l
        j_a = Orb[a].j
        N_a = 2*n_a + l_a
        if (N_a <= N_max)
            open(Int_File, "r") do Bin_Read
                seek(Bin_Read, N_Chunk_skip[a])
                @inbounds for b in 1:a
                    n_b = Orb[b].n
                    l_b = Orb[b].l
                    j_b = Orb[b].j
                    N_b = 2*n_b + l_b
                    if ((N_a + N_b) <= N_2max)
                        @inbounds for c in 1:a_max
                            n_c = Orb[c].n
                            l_c = Orb[c].l
                            j_c = Orb[c].j
                            N_c = 2*n_c + l_c
                            if ((N_a + N_c) <= N_2max) && ((N_b + N_c) <= N_2max) && ((N_a + N_b + N_c) <= N_3max)
                                P_abc = rem(l_a + l_b + l_c,2) + 1
                                @inbounds for d in 1:a 
                                    n_d = Orb[d].n
                                    l_d = Orb[d].l
                                    j_d = Orb[d].j
                                    N_d = 2*n_d + l_d
                                    @inbounds for e in 1:d
                                        n_e = Orb[e].n
                                        l_e = Orb[e].l
                                        j_e = Orb[e].j
                                        N_e = 2*n_e + l_e
                                        if ((N_d + N_e) <= N_2max)
                                            @inbounds for f in 1:a_max
                                                n_f = Orb[f].n
                                                l_f = Orb[f].l
                                                j_f = Orb[f].j
                                                N_f = 2*n_f + l_f
                                                if (j_c == j_f) && (l_c == l_f) && ((N_d + N_f) <= N_2max) && ((N_e + N_f) <= N_2max) && ((N_d + N_e + N_f) <= N_3max)
                                                    P_def = rem(l_d + l_e + l_f,2) + 1
                                                    if (P_abc == P_def)
                                                        @inbounds for J in div(max(abs(j_a - j_b),abs(j_d - j_e)),2):div(min((j_a + j_b),(j_d + j_e)),2)
                                                            @inbounds for T_ab in 0:1
                                                                @inbounds for T_de in 0:1
                                                                    @inbounds for T in max(abs(2*T_ab-1),abs(2*T_de-1)):2:min((2*T_ab+1),(2*T_de+1))

                                                                        V = read(Bin_Read,Float32)

                                                                        T_ind = div((T+1),2)
                                                                        P = rem(l_a + l_b + l_c,2) + 1
                                                                        Ind = V3B_NO2B_Index(a,b,c,T_ab,d,e,f,T_de,J,T,P,Orb,Orb_NNN)
                                                                        if j_c == abs(2*l_c-1)
                                                                            VNNN[P,l_c+1,J+1,T_ind][1][Ind] = V
                                                                        elseif j_c == (2*l_c+1)
                                                                            VNNN[P,l_c+1,J+1,T_ind][2][Ind] = V
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
                        end
                    end
                end
            end
        end
    end

    println("\nFinished loading of 3-body NO2B NNN interaction...")

    return VNNN, Orb_NNN
end

function V3B_NO2B_Count(N_max::Int64,N_2max::Int64,N_3max::Int64,Orb::Vector{NOrb})
    a_max = div((N_max+1)*(N_max+2),2)
    N_Chunk = Vector{Int64}(undef,a_max)
    N_Chunk_skip = Vector{Int64}(undef,a_max)
    
    @inbounds Threads.@threads for a in 1:a_max
        n_a = Orb[a].n
        l_a = Orb[a].l
        j_a = Orb[a].j
        N_a = 2*n_a + l_a
        N = 0
        if (N_a <= N_max)
            @inbounds for b in 1:a
                n_b = Orb[b].n
                l_b = Orb[b].l
                j_b = Orb[b].j
                N_b = 2*n_b + l_b
                if ((N_a + N_b) <= N_2max)
                    @inbounds for c in 1:a_max
                        n_c = Orb[c].n
                        l_c = Orb[c].l
                        j_c = Orb[c].j
                        N_c = 2*n_c + l_c
                        if ((N_a + N_c) <= N_2max) && ((N_b + N_c) <= N_2max) && ((N_a + N_b + N_c) <= N_3max)
                            P_abc = rem(l_a + l_b + l_c,2) + 1
                            @inbounds for d in 1:a 
                                n_d = Orb[d].n
                                l_d = Orb[d].l
                                j_d = Orb[d].j
                                N_d = 2*n_d + l_d
                                @inbounds for e in 1:d
                                    n_e = Orb[e].n
                                    l_e = Orb[e].l
                                    j_e = Orb[e].j
                                    N_e = 2*n_e + l_e
                                    if ((N_d + N_e) <= N_2max)
                                        @inbounds for f in 1:a_max
                                            n_f = Orb[f].n
                                            l_f = Orb[f].l
                                            j_f = Orb[f].j
                                            N_f = 2*n_f + l_f
                                            if (j_c == j_f) && (l_c == l_f) && ((N_d + N_f) <= N_2max) && ((N_e + N_f) <= N_2max) && ((N_d + N_e + N_f) <= N_3max)
                                                P_def = rem(l_d + l_e + l_f,2) + 1
                                                if (P_abc == P_def)
                                                    @inbounds for J in div(max(abs(j_a - j_b),abs(j_d - j_e)),2):div(min((j_a + j_b),(j_d + j_e)),2)
                                                        @inbounds for T_ab in 0:1
                                                            @inbounds for T_de in 0:1
                                                                @inbounds for T in max(abs(2*T_ab-1),abs(2*T_de-1)):2:min((2*T_ab+1),(2*T_de+1))
                                                                    N += 1
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
                end
            end
        end
        N_Chunk[a] = N
    end

    @inbounds for a in 1:a_max
        Sum = 0
        @inbounds for b in 1:(a-1)
            Sum += 4 * N_Chunk[b]
        end
        N_Chunk_skip[a] = Sum
    end

    return N_Chunk_skip
end
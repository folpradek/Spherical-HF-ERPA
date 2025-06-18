function V2B_Res_Density(Params::Vector{Any},Orb::Vector{NOrb},Orb_NN::NNOrb,Orb_NNN::NNNOrb,VNN::NNInt,VNNN::Array{Vector{Vector{Float32}},4},pRho::Matrix{Float64},nRho::Matrix{Float64})
    println("\nStarting calculation of residual NN interaction...\n")
    N_max = Params[7]
    N_2max = Params[8]
    N_3max = Params[9]
    a_max = div((N_max + 1)*(N_max + 2),2)
    J_max = N_2max + 1

    JP = JP_Ini(J_max)

    # Include NO2B NNN interaction to VNN...
    println("\nMaking density dependent residual 2-body interaction...")

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
        @inbounds Threads.@threads  for Bra in 1:max(N_t0, N_t1)

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
                    @views VNN.pn[P,J+1][Ind] += pnSum
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
                    @views VNN.pp[P,J+1][Ind] += ppSum
                    @views VNN.nn[P,J+1][Ind] += nnSum
                end

    
            end

        end
    end

    # Perhaps faster version of previous loop ... I yet need to test its functionality ...
    #=
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
        @inbounds Threads.@threads for Bra in 1:N_t0
            t1 = true
            if Bra > N_t1
                t1 = false
            end
            @views a_0 = Orb_NN.Ind[1,P,J+1][Bra][1]
            @views b_0 = Orb_NN.Ind[1,P,J+1][Bra][2]
            l_a_0 = Orb[a_0].l
            n_a_0 = Orb[a_0].n
            l_b_0 = Orb[b_0].l
            n_b_0 = Orb[b_0].n
            if t1
                @views a_1 = Orb_NN.Ind[2,P,J+1][Bra][1]
                @views b_1 = Orb_NN.Ind[2,P,J+1][Bra][2]
                l_a_1 = Orb[a_1].l
                n_a_1 = Orb[a_1].n
                l_b_1 = Orb[b_1].l
                n_b_1 = Orb[b_1].n
            end

            @inbounds for Ket in 1:Bra
                Ind_0 = Bra + (Ket - 1) * N_t0 - div(Ket * (Ket - 1),2)
                @views d_0 = Orb_NN.Ind[1,P,J+1][Ket][1]
                @views e_0 = Orb_NN.Ind[1,P,J+1][Ket][2]
                l_d_0 = Orb[d_0].l
                n_d_0 = Orb[d_0].n
                l_e_0 = Orb[e_0].l
                n_e_0 = Orb[e_0].n
                pnSum = 0.0
                if t1
                    Ind_1 = Bra + (Ket - 1) * N_t1 - div(Ket * (Ket - 1),2)
                    @views d_1 = Orb_NN.Ind[2,P,J+1][Ket][1]
                    @views e_1 = Orb_NN.Ind[2,P,J+1][Ket][2]
                    l_d_1 = Orb[d_1].l
                    n_d_1 = Orb[d_1].n
                    l_e_1 = Orb[e_1].l
                    n_e_1 = Orb[e_1].n
                    ppSum = 0.0
                    nnSum = 0.0
                end
                @inbounds for c in 1:a_max
                    n_c = Orb[c].n
                    l_c = Orb[c].l
                    j_c = Orb[c].j
                    if ((2*(n_a_0 + n_b_0 + n_c) + l_a_0 + l_b_0 + l_c) <= N_3max)
                        P3B_0 = rem(l_a_0 + l_b_0 + l_c, 2) + 1
                        if t1
                            if ((2*(n_a_1 + n_b_1 + n_c) + l_a_1 + l_b_1 + l_c) <= N_3max)
                                P3B_1 = rem(l_a_1 + l_b_1 + l_c, 2) + 1
                            else
                                t1 = false
                            end
                        end
                        @inbounds for f in 1:a_max
                            n_f = Orb[f].n
                            l_f = Orb[f].l
                            j_f = Orb[f].j
                            if (l_c == l_f) && (j_c == j_f) && ((2*(n_d_0 + n_e_0 + n_f) + l_d_0 + l_e_0 + l_f) <= N_3max)

                                Hat = 1.0/(Float64(2*J) + 1.0)
                                ME001 = V3B_NO2B(a_0,b_0,c,0,d_0,e_0,f,0,J,1,P3B_0,VNNN,Orb,Orb_NNN)
                                ME101 = V3B_NO2B(a_0,b_0,c,1,d_0,e_0,f,0,J,1,P3B_0,VNNN,Orb,Orb_NNN)
                                ME011 = V3B_NO2B(a_0,b_0,c,0,d_0,e_0,f,1,J,1,P3B_0,VNNN,Orb,Orb_NNN)
                                ME111 = V3B_NO2B(a_0,b_0,c,1,d_0,e_0,f,1,J,1,P3B_0,VNNN,Orb,Orb_NNN)
                                ME113 = V3B_NO2B(a_0,b_0,c,1,d_0,e_0,f,1,J,3,P3B_0,VNNN,Orb,Orb_NNN)

                                @views pnSum += Hat * ((1/2 * ME001 - 1/sqrt(12) * ME101 - 1/sqrt(12) * ME011 +
                                            1/6 * ME111 + 1/3 * ME113) * pRho[c,f] + (1/2 * ME001 + 1/sqrt(12) *
                                            ME101 + 1/sqrt(12) * ME011 + 1/6 * ME111 +  1/3 * ME113) * nRho[c,f])

                                if t1
                                    if ((2*(n_d_1 + n_e_1 + n_f) + l_d_1 + l_e_1 + l_f) <= N_3max)
                                        ME111 = V3B_NO2B(a_1,b_1,c,1,d_1,e_1,f,1,J,1,P3B_1,VNNN,Orb,Orb_NNN)
                                        ME113 = V3B_NO2B(a_1,b_1,c,1,d_1,e_1,f,1,J,3,P3B_1,VNNN,Orb,Orb_NNN)
    
                                        @views ppSum += Hat * (ME113 * pRho[c,f] + (2/3 * ME111 + 1/3 * ME113) * nRho[c,f])
                                        @views nnSum += Hat * (ME113 * nRho[c,f] + (2/3 * ME111 + 1/3 * ME113) * pRho[c,f])
                                    else
                                        t1 = false
                                    end
                                end

                            end
                        end
                    end
                end
                @views VNN.pn[P,J+1][Ind_0] += pnSum
                if t1
                    @views VNN.pp[P,J+1][Ind_1] += ppSum
                    @views VNN.nn[P,J+1][Ind_1] += nnSum
                end
            end

        end
    end
    =#

    return VNN
end

function V2B_Res_Ini(Orb::Vector{NOrb},N_max::Int64)
    N_2max = 2 * N_max
    a_max = div((N_max + 1)*(N_max + 2),2)
    J_max = Int64(N_2max + 1)

    N_Orb_NN = zeros(Int64,2,J_max+1)
    Orb_NN_Dic = Dict{Tuple{Int8,Int8,Int16,Int16},Int32}()

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
                        N_Orb_NN[P,J+1] += 1
                    end
                end
            end
        end
    end

    Ind_Orb_NN = Matrix{Vector{Vector{Int64}}}(undef,2,J_max+1)
    for P = 1:2
        for J = 0:J_max
            Ind_Orb_NN[P,J+1] = Vector{Vector{Int64}}(undef,N_Orb_NN[P,J+1])
        end
    end

    N_Orb_NN = zeros(Int64,2,J_max+1)

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
                        key = (Int8(P),Int8(J),Int16(a),Int16(b))
                        N_Orb_NN[P,J+1] += 1
                        Orb_NN_Dic[key] = Int32(N_Orb_NN[P,J+1])
                        Ind_Orb_NN[P,J+1][N_Orb_NN[P,J+1]] = zeros(Int64, 2)
                        Ind_Orb_NN[P,J+1][N_Orb_NN[P,J+1]][1] = a
                        Ind_Orb_NN[P,J+1][N_Orb_NN[P,J+1]][2] = b
                    end
                end
            end
        end
    end

    Vpp = Matrix{Matrix{Float64}}(undef,2,J_max+1)
    Vpn = Matrix{Matrix{Float64}}(undef,2,J_max+1)
    Vnn = Matrix{Matrix{Float64}}(undef,2,J_max+1)
    
    for P = 1:2
        for J=0:J_max
            Length = N_Orb_NN[P,J+1]
            Vpp[P,J+1] = zeros(Float64, Length, Length)
            Vpn[P,J+1] = zeros(Float64, Length, Length)
            Vnn[P,J+1] = zeros(Float64, Length, Length)
        end
    end

    VNN = NNInt_Res(Vpp, Vpn, Vnn)

    Orb_NN = NNOrb_Res(Orb_NN_Dic, N_Orb_NN, Ind_Orb_NN)

    return VNN, Orb_NN

end

function V2B_Res_Index(a::Int64,b::Int64,J::Int64,P::Int64,Orb_NN_Res::NNOrb_Res)
    P = Int8(P)
    J = Int8(J)
    a = Int16(a)
    b = Int16(b)
    Ind = Int64(Orb_NN_Res.Dic[(P,J,a,b)])
    return Ind
end

function Orb_PreComp(a_max::Int64,j::Int64,l::Int64,Orb::Vector{NOrb})
    Orb_x = Vector{Int64}()
    @inbounds for x in 1:a_max
        orb = Orb[x]
        if (l == orb.l) && (j == orb.j)
            push!(Orb_x, orb.a)
        end
    end
    return Orb_x
end

function V2B_Res_Ind1(Params::Vector{Any},Orb::Vector{NOrb},Orb_NN::NNOrb,VNN::NNInt)
    # Parameter initialization...
    N_max = Params[7]
    N_2max = 2*N_max
    Output_File = Params[13]
    a_max = div((N_max + 1)*(N_max + 2),2)
    J_max = N_2max + 1

    # Import HF basis transformation matrices ...
    pU_Import = "IO/" * Output_File * "/Bin/pU.bin"
    pU = Matrix{Float64}(undef,a_max,a_max)
    open(pU_Import, "r") do Read_File
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                ME = read(Read_File, Float64)
                pU[a,b] = ME
            end
        end
    end

    nU_Import = "IO/" * Output_File * "/Bin/nU.bin"
    nU = Matrix{Float64}(undef,a_max,a_max)
    open(nU_Import, "r") do Read_File
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                ME = read(Read_File, Float64)
                nU[a,b] = ME
            end
        end
    end

    JP = JP_Ini(J_max)

    # Make temporary arrays for NN interaction...
    VNN_Res_I, Orb_NN_Res = V2B_Res_Ini(Orb,N_max)

    # Index 1
    println("\nTransformation in 1st index...")
    @time @inbounds Threads.@threads for i in JP
        J = i[1]
        P = i[2]
        if P == 1
            println("Calculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = +")
        else
            println("Calculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = -")
        end
        @views N = Orb_NN_Res.N[P,J+1]
        @inbounds for Bra in 1:N
            @views a = Orb_NN_Res.Ind[P,J+1][Bra][1]
            @views b = Orb_NN_Res.Ind[P,J+1][Bra][2]
            l_a = Orb[a].l
            j_a = Orb[a].j

            Orb_x = Orb_PreComp(a_max,j_a,l_a,Orb)
            @inbounds for Ket in 1:N
                c = Orb_NN_Res.Ind[P,J+1][Ket][1]
                d = Orb_NN_Res.Ind[P,J+1][Ket][2]

                ppSum = 0.0
                pnSum = 0.0
                nnSum = 0.0
                @inbounds for i in Orb_x
                    MEpp, MEpn, MEnn = V2B_Res_Ind1_MEs(a,b,c,d,i,J,pU,nU,VNN,Orb,Orb_NN)
                    ppSum += MEpp
                    pnSum += MEpn
                    nnSum += MEnn
                end
                @views VNN_Res_I.pp[P,J+1][Bra,Ket] = ppSum
                @views VNN_Res_I.pn[P,J+1][Bra,Ket] = pnSum
                @views VNN_Res_I.nn[P,J+1][Bra,Ket] = nnSum
            end
        end
    end
    return VNN_Res_I, Orb_NN_Res
end

function V2B_Res_Ind2(Params::Vector{Any},Orb::Vector{NOrb},Orb_NN_Res::NNOrb_Res,VNN_Res_I::NNInt_Res)
    # Parameter initialization...
    N_max = Params[7]
    N_2max = 2*N_max
    Output_File = Params[13]
    a_max = div((N_max + 1)*(N_max + 2),2)
    J_max = N_2max + 1

    # Import HF basis transformation matrices ...
    pU_Import = "IO/" * Output_File * "/Bin/pU.bin"
    pU = Matrix{Float64}(undef,a_max,a_max)
    open(pU_Import, "r") do Read_File
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                ME = read(Read_File, Float64)
                pU[a,b] = ME
            end
        end
    end

    nU_Import = "IO/" * Output_File * "/Bin/nU.bin"
    nU = Matrix{Float64}(undef,a_max,a_max)
    open(nU_Import, "r") do Read_File
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                ME = read(Read_File, Float64)
                nU[a,b] = ME
            end
        end
    end

    JP = JP_Ini(J_max)

    # Make 2nd temporary array for NN interaction...
    VNN_Res_II = NNInt_Res(deepcopy(VNN_Res_I.pp), deepcopy(VNN_Res_I.pn), deepcopy(VNN_Res_I.nn))

    # Index 2
    println("\nTransformation in 2nd index...")
    @time @inbounds Threads.@threads for i in JP
        J = i[1]
        P = i[2]
        if P == 1
            println("Calculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = +")
        else
            println("Calculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = -")
        end
        N = Orb_NN_Res.N[P,J+1]
        @inbounds for Bra in 1:N
            a = @views Orb_NN_Res.Ind[P,J+1][Bra][1]
            b = @views Orb_NN_Res.Ind[P,J+1][Bra][2]
            l_b = Orb[b].l
            j_b = Orb[b].j
            Orb_x = Orb_PreComp(a_max,j_b,l_b,Orb)
            @inbounds for Ket in 1:N
                ppSum, pnSum, nnSum = 0.0, 0.0, 0.0
                @inbounds for i in Orb_x
                    Bra_ai = V2B_Res_Index(a,i,J,P,Orb_NN_Res)
                    MEpp, MEpn, Mnn = V2B_Res_Ind2_MEs(P,J,Bra_ai,Ket,i,b,pU,nU,VNN_Res_I)
                    ppSum += MEpp
                    pnSum += MEpn
                    nnSum += Mnn
                end
                @views VNN_Res_II.pp[P,J+1][Bra,Ket] = ppSum
                @views VNN_Res_II.pn[P,J+1][Bra,Ket] = pnSum
                @views VNN_Res_II.nn[P,J+1][Bra,Ket] = nnSum
            end
        end
    end

    VNN_Res_I = NNInt_Res(deepcopy(VNN_Res_II.pp),deepcopy(VNN_Res_II.pn),deepcopy(VNN_Res_II.nn))

    return VNN_Res_I, VNN_Res_II
end

function V2B_Res_Ind3(Params::Vector{Any},Orb::Vector{NOrb},Orb_NN_Res::NNOrb_Res,VNN_Res_I::NNInt_Res,VNN_Res_II::NNInt_Res)
    # Parameter initialization...
    N_max = Params[7]
    N_2max = 2*N_max
    Output_File = Params[13]
    a_max = div((N_max + 1)*(N_max + 2),2)
    J_max = N_2max + 1

    # Import HF basis transformation matrices ...
    pU_Import = "IO/" * Output_File * "/Bin/pU.bin"
    pU = Matrix{Float64}(undef,a_max,a_max)
    open(pU_Import, "r") do Read_File
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                ME = read(Read_File, Float64)
                pU[a,b] = ME
            end
        end
    end

    nU_Import = "IO/" * Output_File * "/Bin/nU.bin"
    nU = Matrix{Float64}(undef,a_max,a_max)
    open(nU_Import, "r") do Read_File
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                ME = read(Read_File, Float64)
                nU[a,b] = ME
            end
        end
    end

    JP = JP_Ini(J_max)

    # Index 3
    println("\nTransformation in 3rd index...")
    @time @inbounds Threads.@threads for i in JP
        J = i[1]
        P = i[2]
        if P == 1
            println("Calculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = +")
        else
            println("Calculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = -")
        end
        @views N = Orb_NN_Res.N[P,J+1]
        @inbounds for Ket in 1:N
            @views c = Orb_NN_Res.Ind[P,J+1][Ket][1]
            @views d = Orb_NN_Res.Ind[P,J+1][Ket][2]
            l_c = Orb[c].l
            j_c = Orb[c].j
            Orb_x = Orb_PreComp(a_max,j_c,l_c,Orb)
            @inbounds for Bra in 1:N
                ppSum = 0.0
                pnSum = 0.0
                nnSum = 0.0
                @inbounds for i in Orb_x
                    Ket_id = V2B_Res_Index(i,d,J,P,Orb_NN_Res)
                    MEpp, MEpn, MEnn = V2B_Res_Ind3_MEs(P,J,Bra,Ket_id,i,c,pU,nU,VNN_Res_I)
                    ppSum += MEpp
                    pnSum += MEpn
                    nnSum += MEnn
                end
                @views VNN_Res_II.pp[P,J+1][Bra,Ket] = ppSum
                @views VNN_Res_II.pn[P,J+1][Bra,Ket] = pnSum
                @views VNN_Res_II.nn[P,J+1][Bra,Ket] = nnSum
            end
        end
    end

    return VNN_Res_II
end

function V2B_Res_Ind4(Params::Vector{Any},Orb::Vector{NOrb},Orb_NN_Res::NNOrb_Res,VNN_Res_II::NNInt_Res)
    # Parameter initialization...
    N_max = Params[7]
    N_2max = 2*N_max
    Output_File = Params[13]
    a_max = div((N_max + 1)*(N_max + 2),2)
    J_max = N_2max + 1
    
    # Import HF basis transformation matrices ...
    pU_Import = "IO/" * Output_File * "/Bin/pU.bin"
    pU = Matrix{Float64}(undef,a_max,a_max)
    open(pU_Import, "r") do Read_File
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                ME = read(Read_File, Float64)
                pU[a,b] = ME
            end
        end
    end

    nU_Import = "IO/" * Output_File * "/Bin/nU.bin"
    nU = Matrix{Float64}(undef,a_max,a_max)
    open(nU_Import, "r") do Read_File
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                ME = read(Read_File, Float64)
                nU[a,b] = ME
            end
        end
    end

    JP = JP_Ini(J_max)

    VNN, Orb_NN = V2B_Ini(Orb,N_max)

    # Index 4
    println("\nTransformation in 4th index...")
    @time @inbounds Threads.@threads for i in JP
        J = i[1]
        P = i[2]
        if P == 1
            println("Calculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = +")
        else
            println("Calculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = -")
        end
        @views N_t0 = Orb_NN.N[1,P,J+1]
        @views N_t1 = Orb_NN.N[2,P,J+1]
        @inbounds for Bra in 1:max(N_t0, N_t1)
            @inbounds for Ket in 1:Bra

                if Bra <= N_t0
                    Ind = Bra + (Ket - 1) * N_t0 - div(Ket * (Ket - 1),2)
                    @views a = Orb_NN.Ind[1,P,J+1][Bra][1]
                    @views b = Orb_NN.Ind[1,P,J+1][Bra][2]
                    @views c = Orb_NN.Ind[1,P,J+1][Ket][1]
                    @views d = Orb_NN.Ind[1,P,J+1][Ket][2]
                    l_d = Orb[d].l
                    j_d = Orb[d].j
                    Bra_ab = V2B_Res_Index(a,b,J,P,Orb_NN_Res)
                    Orb_x = Orb_PreComp(a_max,j_d,l_d,Orb)
                    pnSum = 0.0
                    @inbounds for i in Orb_x
                        Ket_ci = V2B_Res_Index(c,i,J,P,Orb_NN_Res)
                        MEpn = V2B_Res_Ind4_MEs_t0(P,J,Bra_ab,Ket_ci,i,d,nU,VNN_Res_II)
                        pnSum += MEpn
                    end
                    @views VNN.pn[P,J+1][Ind] = pnSum
                end

                if Bra <= N_t1
                    Ind = Bra + (Ket - 1) * N_t1 - div(Ket * (Ket - 1),2)
                    @views a = Orb_NN.Ind[2,P,J+1][Bra][1]
                    @views b = Orb_NN.Ind[2,P,J+1][Bra][2]
                    @views c = Orb_NN.Ind[2,P,J+1][Ket][1]
                    @views d = Orb_NN.Ind[2,P,J+1][Ket][2]
                    l_d = Orb[d].l
                    j_d = Orb[d].j
                    Bra_ab = V2B_Res_Index(a,b,J,P,Orb_NN_Res)
                    Orb_x = Orb_PreComp(a_max,j_d,l_d,Orb)
                    ppSum = 0.0
                    nnSum = 0.0
                    @inbounds for i in Orb_x
                        Ket_ci = V2B_Res_Index(c,i,J,P,Orb_NN_Res)
                        MEpp, MEnn = V2B_Res_Ind4_MEs_t1(P,J,Bra_ab,Ket_ci,i,d,pU,nU,VNN_Res_II)
                        ppSum += MEpp
                        nnSum += MEnn
                    end
                    @views VNN.pp[P,J+1][Ind] = ppSum
                    @views VNN.nn[P,J+1][Ind] = nnSum
                end

            end
        end

    end

    return VNN, Orb_NN
end

function V2B_Res(Params::Vector{Any},Orb::Vector{NOrb},Orb_NN::NNOrb,VNN::NNInt)
    # Transformation of VNN to the HF basis...
    println("\nTransforming residual 2-body interaction to the HF basis...")

    # Index 1
    VNN_Res_I, Orb_NN_Res = V2B_Res_Ind1(Params,Orb,Orb_NN,VNN)

    # Index 2
    VNN_Res_I, VNN_Res_II = V2B_Res_Ind2(Params,Orb,Orb_NN_Res,VNN_Res_I)

    # Index 3
    VNN_Res_II = V2B_Res_Ind3(Params,Orb,Orb_NN_Res,VNN_Res_I,VNN_Res_II)

    # Index 4
    VNN, Orb_NN = V2B_Res_Ind4(Params,Orb,Orb_NN_Res,VNN_Res_II)
    
    println("\nResidual 2-body interaction ready...\n")

    return VNN, Orb_NN

end

@inline function V2B_Res_Ind1_MEs(a::Int64,b::Int64,c::Int64,d::Int64,i::Int64,J::Int64,pU::Matrix{Float64},nU::Matrix{Float64},VNN::NNInt,Orb::Vector{NOrb},Orb_NN::NNOrb)
    return @views V2B(i,b,c,d,J,1,VNN.pp,Orb,Orb_NN) * pU[i,a], V2B(i,b,c,d,J,0,VNN.pn,Orb,Orb_NN) * pU[i,a], V2B(i,b,c,d,J,1,VNN.nn,Orb,Orb_NN) * nU[i,a]
end

@inline  function V2B_Res_Ind2_MEs(P::Int64,J::Int64,Bra::Int64,Ket::Int64,i::Int64,b::Int64,pU::Matrix{Float64},nU::Matrix{Float64},VNN_Res_I::NNInt_Res)
    return @views pU[i,b] * VNN_Res_I.pp[P,J+1][Bra,Ket], nU[i,b] * VNN_Res_I.pn[P,J+1][Bra,Ket], nU[i,b] * VNN_Res_I.nn[P,J+1][Bra,Ket]
end

@inline  function V2B_Res_Ind3_MEs(P::Int64,J::Int64,Bra::Int64,Ket::Int64,i::Int64,c::Int64,pU::Matrix{Float64},nU::Matrix{Float64},VNN_Res_I::NNInt_Res)
    return @views VNN_Res_I.pp[P,J+1][Bra,Ket] * pU[i,c], VNN_Res_I.pn[P,J+1][Bra,Ket] * pU[i,c], VNN_Res_I.nn[P,J+1][Bra,Ket] * nU[i,c]
end

@inline  function V2B_Res_Ind4_MEs_t0(P::Int64,J::Int64,Bra::Int64,Ket::Int64,i::Int64,d::Int64,nU::Matrix{Float64},VNN_Res_II::NNInt_Res)
    return @views VNN_Res_II.pn[P,J+1][Bra,Ket] * nU[i,d]
end

@inline  function V2B_Res_Ind4_MEs_t1(P::Int64,J::Int64,Bra::Int64,Ket::Int64,i::Int64,d::Int64,pU::Matrix{Float64},nU::Matrix{Float64},VNN_Res_II::NNInt_Res)
    return @views VNN_Res_II.pp[P,J+1][Bra,Ket] * pU[i,d], VNN_Res_II.nn[P,J+1][Bra,Ket] * nU[i,d]
end

function V2B_Res_Export(Params::Vector{Any},Orb_NN::NNOrb,VNN_Res::NNInt)
    Format = Params[12]

    if Format == "Bin"
        V2B_Res_Export_Bin(Params,Orb_NN,VNN_Res)
    elseif Format == "HRBin"
        V2B_Res_Export_HRBin(Params,Orb_NN,VNN_Res)
    elseif Format == "HR"
        V2B_Res_Export_HR(Params,Orb_NN,VNN_Res)
    end

    return
end

function V2B_Res_Export_Bin(Params::Vector{Any},Orb_NN::NNOrb,VNN_Res::NNInt)
    # Read parameters ...
    N_max = Params[7]
    N_2max = 2*N_max
    Output_File = Params[13]

    J_max = N_2max + 1
    V2B_Res_Export = "IO/" * Output_File * "/Bin/V2B_Res.bin"

    # Export residual 2-body interaction ...

    println("\nExporting residual 2-body interaction in internal binary format ...")

    open(V2B_Res_Export, "w") do Export_File
        # ppV
        @inbounds for J in 0:J_max
            @inbounds for P in 1:2
                N_t1 = Orb_NN.N[2,P,J+1]
                @inbounds for Bra in 1:N_t1
                    @inbounds for Ket in 1:Bra
                        Ind = Bra + (Ket - 1) * N_t1 - div(Ket * (Ket - 1),2)
                        ME = @views VNN_Res.pp[P,J+1][Ind]
                        write(Export_File, Float64(ME))
                    end
                end
            end
        end

        # pnV
        @inbounds for J in 0:J_max
            @inbounds for P in 1:2
                N_t0 = Orb_NN.N[1,P,J+1]
                @inbounds for Bra in 1:N_t0
                    @inbounds for Ket in 1:N_t0
                        Ind = Bra + (Ket - 1) * N_t0 - div(Ket * (Ket - 1),2)
                        ME = @views VNN_Res.pn[P,J+1][Ind]
                        write(Export_File, Float64(ME))
                    end
                end
            end
        end

        # nnV
        @inbounds for J in 0:J_max
            @inbounds for P in 1:2
                N_t1 = Orb_NN.N[2,P,J+1]
                @inbounds for Bra in 1:N_t1
                    @inbounds for Ket in 1:Bra
                        Ind = Bra + (Ket - 1) * N_t1 - div(Ket * (Ket - 1),2)
                        ME = @views VNN_Res.nn[P,J+1][Ind]
                        write(Export_File, Float64(ME))
                    end
                end
            end
        end

    end

    println("\nResidual 2-body interaction succesfully exported into file V2B_Res.bin ...")

    return
end

function V2B_Res_Export_HR(Params::Vector{Any},Orb_NN::NNOrb,VNN_Res::NNInt)
    # Read parameters ...
    N_max = Params[7]
    N_2max = 2*N_max
    Output_File = Params[13]

    J_max = N_2max + 1
    V2B_Res_Export = "IO/" * Output_File * "/V2B_Res_HR.dat"

    # Export residual 2-body interaction ...

    println("\nExporting residual 2-body interaction in human reeadable format ...")

    open(V2B_Res_Export, "w") do Export_File
        println(Export_File, "t\ta\tb\tc\td\tJ\tV")
        # ppV
        @inbounds for J in 0:J_max
            @inbounds for P in 1:2
                N_t1 = Orb_NN.N[2,P,J+1]
                @inbounds for Bra in 1:N_t1
                    a = 2*Orb_NN.Ind[2,P,J+1][Bra][1] - 1
                    b = 2*Orb_NN.Ind[2,P,J+1][Bra][2] - 1
                    @inbounds for Ket in 1:Bra
                        c = 2*Orb_NN.Ind[2,P,J+1][Ket][1] - 1
                        d = 2*Orb_NN.Ind[2,P,J+1][Ket][2] - 1
                        Ind = Bra + (Ket - 1) * N_t1 - div(Ket * (Ket - 1),2)
                        ME = VNN_Res.pp[P,J+1][Ind]
                        Row = "-1\t" * string(a) * "\t" * string(b) * "\t" *
                                string(c) * "\t" * string(d) * "\t" * string(J) *
                                "\t" * string(round(ME, digits = 10))
                        println(Export_File, Row)
                    end
                end
            end
        end


        # pnV
        @inbounds for J in 0:J_max
            @inbounds for P in 1:2
                N_t0 = Orb_NN.N[1,P,J+1]
                @inbounds for Bra in 1:N_t0
                    a = 2*Orb_NN.Ind[1,P,J+1][Bra][1] - 1
                    b = 2*Orb_NN.Ind[1,P,J+1][Bra][2]
                    @inbounds for Ket in 1:N_t0
                        c = 2*Orb_NN.Ind[1,P,J+1][Ket][1] - 1
                        d = 2*Orb_NN.Ind[1,P,J+1][Ket][2]
                        Ind = Bra + (Ket - 1) * N_t0 - div(Ket * (Ket - 1),2)
                        ME = VNN_Res.pn[P,J+1][Ind]
                        Row = "0\t" * string(a) * "\t" * string(b) * "\t" *
                                string(c) * "\t" * string(d) * "\t" * string(J) *
                                "\t" * string(round(ME, digits = 10))
                        println(Export_File, Row)
                    end
                end
            end
        end

        # nnV
        @inbounds for J in 0:J_max
            @inbounds for P in 1:2
                N_t1 = Orb_NN.N[2,P,J+1]
                @inbounds for Bra in 1:N_t1
                    a = 2*Orb_NN.Ind[2,P,J+1][Bra][1]
                    b = 2*Orb_NN.Ind[2,P,J+1][Bra][2]
                    @inbounds for Ket in 1:Bra
                        c = 2*Orb_NN.Ind[2,P,J+1][Ket][1]
                        d = 2*Orb_NN.Ind[2,P,J+1][Ket][2]
                        Ind = Bra + (Ket - 1) * N_t1 - div(Ket * (Ket - 1),2)
                        ME = VNN_Res.nn[P,J+1][Ind]
                        Row = "1\t" * string(a) * "\t" * string(b) * "\t" *
                                string(c) * "\t" * string(d) * "\t" * string(J) *
                                "\t" * string(round(ME, digits = 10))
                        println(Export_File, Row)
                    end
                end
            end
        end

    end

    println("\nResidual 2-body interaction succesfully exported into file V2B_Res.bin ...")

    return
end

function V2B_Res_Export_HRBin(Params::Vector{Any},Orb_NN::NNOrb,VNN_Res::NNInt)
    # Read parameters ...
    N_max = Params[7]
    N_2max = 2*N_max
    Output_File = Params[13]

    J_max = N_2max + 1
    V2B_Res_Export = "IO/" * Output_File * "/Bin/V2B_Res_HR.bin"

    # Export residual 2-body interaction ...

    println("\nExporting residual 2-body interaction in Human Readable binary format ...")

    open(V2B_Res_Export, "w") do Export_File
        # ppV

        @inbounds for J in 0:J_max
            @inbounds for P in 1:2
                N_t1 = Orb_NN.N[2,P,J+1]
                @inbounds for Bra in 1:N_t1
                    a = 2*Orb_NN.Ind[2,P,J+1][Bra][1] - 1
                    b = 2*Orb_NN.Ind[2,P,J+1][Bra][2] - 1
                    @inbounds for Ket in 1:Bra
                        c = 2*Orb_NN.Ind[2,P,J+1][Ket][1] - 1
                        d = 2*Orb_NN.Ind[2,P,J+1][Ket][2] - 1
                        Ind = Bra + (Ket - 1) * N_t1 - div(Ket * (Ket - 1),2)
                        ME = VNN_Res.pp[P,J+1][Ind]
                        write(Export_File, Int16(a))
                        write(Export_File, Int16(b))
                        write(Export_File, Int16(c))
                        write(Export_File, Int16(d))
                        write(Export_File, Int16(J))
                        write(Export_File, Float64(ME))
                    end
                end
            end
        end

        # pnV
        @inbounds for J in 0:J_max
            @inbounds for P in 1:2
                N_t0 = Orb_NN.N[1,P,J+1]
                @inbounds for Bra in 1:N_t0
                    a = 2*Orb_NN.Ind[1,P,J+1][Bra][1] - 1
                    b = 2*Orb_NN.Ind[1,P,J+1][Bra][2]
                    @inbounds for Ket in 1:N_t0
                        c = 2*Orb_NN.Ind[1,P,J+1][Ket][1] - 1
                        d = 2*Orb_NN.Ind[1,P,J+1][Ket][2]
                        Ind = Bra + (Ket - 1) * N_t0 - div(Ket * (Ket - 1),2)
                        ME = VNN_Res.pn[P,J+1][Ind]
                        write(Export_File, Int16(a))
                        write(Export_File, Int16(b))
                        write(Export_File, Int16(c))
                        write(Export_File, Int16(d))
                        write(Export_File, Int16(J))
                        write(Export_File, Float64(ME))
                    end
                end
            end
        end

        # nnV
        @inbounds for J in 0:J_max
            @inbounds for P in 1:2
                N_t1 = Orb_NN.N[2,P,J+1]
                @inbounds for Bra in 1:N_t1
                    a = 2*Orb_NN.Ind[2,P,J+1][Bra][1]
                    b = 2*Orb_NN.Ind[2,P,J+1][Bra][2]
                    @inbounds for Ket in 1:Bra
                        c = 2*Orb_NN.Ind[2,P,J+1][Ket][1]
                        d = 2*Orb_NN.Ind[2,P,J+1][Ket][2]
                        Ind = Bra + (Ket - 1) * N_t1 - div(Ket * (Ket - 1),2)
                        ME = VNN_Res.nn[P,J+1][Ind]
                        write(Export_File, Int16(a))
                        write(Export_File, Int16(b))
                        write(Export_File, Int16(c))
                        write(Export_File, Int16(d))
                        write(Export_File, Int16(J))
                        write(Export_File, Float64(ME))
                    end
                end
            end
        end


    end

    println("\nResidual 2-body interaction succesfully exported into file V2B_Res.bin ...")

    return
end

function V2B_Res_Import(N_max::Int64,Import_File::String,Orb::Vector{NOrb})
    # Read parameters ...
    N_2max = 2*N_max
    J_max = N_2max + 1
    V2B_Res_Import = Import_File * "/Bin/V2B_Res.bin"

    # Read Residual 2-body NN interaction ...
    println("\nReading Residual 2-body NN interaction ...")

    # Initialize NN residual interaction arrays ...
    VNN_Res, Orb_NN = V2B_Ini(Orb,N_max)

    # Iteraction list for J & P ...
    JP_list = JP_Ini(J_max)

    # Precalculate chunks of data for parallelized reading ...
    N_Chunk_skip = V2B_Res_Count(N_max,Orb_NN)

    # Load NN residual interaction ...

    # ppV
    @inbounds Threads.@threads for JP in JP_list
        J = JP[1]
        P = JP[2]
        N_t1 = Orb_NN.N[2,P,J+1]
        open(V2B_Res_Import, "r") do Bin_Read
            seek(Bin_Read, N_Chunk_skip[1,P,J+1])
            @inbounds for Bra in 1:N_t1
                @inbounds for Ket in 1:Bra
                    Ind = Bra + (Ket - 1) * N_t1 - div(Ket * (Ket - 1),2)
                    ME = read(Bin_Read, Float64)
                    VNN_Res.pp[P,J+1][Ind] = ME
                end
            end
        end
    end

    # pnV
    @inbounds Threads.@threads for JP in JP_list
        J = JP[1]
        P = JP[2]
        N_t0 = Orb_NN.N[1,P,J+1]
        open(V2B_Res_Import, "r") do Bin_Read
            seek(Bin_Read, N_Chunk_skip[2,P,J+1])
            @inbounds for Bra in 1:N_t0
                @inbounds for Ket in 1:N_t0
                    Ind = Bra + (Ket - 1) * N_t0 - div(Ket * (Ket - 1),2)
                    ME = read(Bin_Read, Float64)
                    VNN_Res.pn[P,J+1][Ind] = ME
                end
            end
        end
    end

    # nnV
    @inbounds Threads.@threads for JP in JP_list
        J = JP[1]
        P = JP[2]
        N_t1 = Orb_NN.N[2,P,J+1]
        open(V2B_Res_Import, "r") do Bin_Read
            seek(Bin_Read, N_Chunk_skip[3,P,J+1])
            @inbounds for Bra in 1:N_t1
                @inbounds for Ket in 1:Bra
                    Ind = Bra + (Ket - 1) * N_t1 - div(Ket * (Ket - 1),2)
                    ME = read(Bin_Read, Float64)
                    VNN_Res.nn[P,J+1][Ind] = ME
                end
            end
        end
    end

    println("\nResidual 2-body NN interaction succesfully loaded ...")

    return VNN_Res, Orb_NN
end

function V2B_Res_Count(N_max::Int64,Orb_NN::NNOrb)
    # Read parameters ...
    N_2max = 2*N_max
    J_max = N_2max + 1

    N_Chunk = Array{Int64}(undef,3,2,J_max+1)
    N_Chunk_skip = Array{Int64}(undef,3,2,J_max+1)

    JP_list = JP_Ini(J_max)

    # ppV
    @inbounds for JP in JP_list
        N = 0
        J = JP[1]
        P = JP[2]
        N_t1 = Orb_NN.N[2,P,J+1]
        @inbounds for Bra in 1:N_t1
            @inbounds for Ket in 1:Bra
                N += 1
            end
        end
        N_Chunk[1,P,J+1] = N
    end

    # pnV
    @inbounds for JP in JP_list
        N = 0
        J = JP[1]
        P = JP[2]
        N_t0 = Orb_NN.N[1,P,J+1]
        @inbounds for Bra in 1:N_t0
            @inbounds for Ket in 1:N_t0
                N += 1
            end
        end
        N_Chunk[2,P,J+1] = N
    end

    # nnV
    @inbounds for JP in JP_list
        N = 0
        J = JP[1]
        P = JP[2]
        N_t1 = Orb_NN.N[2,P,J+1]
        @inbounds for Bra in 1:N_t1
            @inbounds for Ket in 1:Bra
                N += 1
            end
        end
        N_Chunk[3,P,J+1] = N
    end

    @inbounds for t_1 in 1:3
        @inbounds for JP_1 in 1:(2*(J_max+1))
            J_1 = JP_list[JP_1][1]
            P_1 = JP_list[JP_1][2]
            Sum = 0
            @inbounds for t_2 in 1:t_1
                if t_2 < t_1
                    @inbounds for JP_2 in 1:(2*(J_max+1))
                        J_2 = JP_list[JP_2][1]
                        P_2 = JP_list[JP_2][2]
                        Sum += 8 * N_Chunk[t_2,P_2,J_2+1]
                    end
                elseif  t_2 == t_1
                    @inbounds for JP_2 in 1:(JP_1-1)
                        J_2 = JP_list[JP_2][1]
                        P_2 = JP_list[JP_2][2]
                        Sum += 8 * N_Chunk[t_2,P_2,J_2+1]
                    end
                end
            end
            N_Chunk_skip[t_1,P_1,J_1+1] = Sum
        end
    end

    return N_Chunk_skip
end
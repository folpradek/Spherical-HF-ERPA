struct W2B
    pp::Matrix{Vector{Float64}}
    ppn::Matrix{Vector{Float64}}
    npn::Matrix{Vector{Float64}}
    nn::Matrix{Vector{Float64}}
end

struct V2B_Temp
    pp::Matrix{Matrix{Float64}}
    pn::Matrix{Matrix{Float64}}
    nn::Matrix{Matrix{Float64}}
end

struct W2B_Temp
    pp::Matrix{Matrix{Float64}}
    ppn::Matrix{Matrix{Float64}}
    npn::Matrix{Matrix{Float64}}
    nn::Matrix{Matrix{Float64}}
end

struct NNOrb_Temp
    Dic::Dict{Tuple{Int8,Int8,Int16,Int16},Int32}
    N::Matrix{Int64}
    Ind::Matrix{Vector{Vector{Int64}}}
end

struct qpH40
    pp::Matrix{Vector{Float64}}
    pn::Matrix{Vector{Float64}}
    nn::Matrix{Vector{Float64}}
end

struct qpH31
    pp::Matrix{Vector{Float64}}
    pn2011::Matrix{Vector{Float64}}
    pn1120::Matrix{Vector{Float64}}
    nn::Matrix{Vector{Float64}}
end

struct qpH22
    pp::Matrix{Vector{Float64}}
    pn2002::Matrix{Vector{Float64}}
    pn1111::Matrix{Vector{Float64}}
    pn0220::Matrix{Vector{Float64}}
    nn::Matrix{Vector{Float64}}
end

struct qpH2B
    H40::qpH40
    H31::qpH31
    H22::qpH22
end

function W2b_initialize(Orb::Vector{NOrb},N_max::Int64)
    # Read parameters ...
    N_2max = 2 * N_max
    a_max = div((N_max + 1)*(N_max + 2),2)
    J_max = 0

    # Initiliaze array for counting of NN orbitals ...
    N_Orb_NN = zeros(Int64,2,2,J_max+1)

    # Initiliaze dictionary for NN orbitals ...
    Orb_NN_Dic = Dict{Tuple{Int8,Int8,Int8,Int16,Int16},Int32}()

    # Count the number of NN orbitals ...
    @inbounds for a in 1:a_max
        n_a = Orb[a].n
        l_a = Orb[a].l
        j_a = Orb[a].j
        if (2*n_a + l_a) <= N_max
            @inbounds for b = 1:a_max
                n_b = Orb[b].n
                l_b = Orb[b].l
                j_b = Orb[b].j
                if ((2*(n_a + n_b) + l_a + l_b) <= N_2max)
                    P = rem(l_a + l_b, 2) + 1
                    J = 0
                    if b <= a
                        N_Orb_NN[2,P,J+1] += 1
                    end
                    N_Orb_NN[1,P,J+1] += 1
                end
            end
        end
    end

    # Initialite the array for NN orbital indices ...
    Ind_Orb_NN = Array{Vector{Vector{Int64}}}(undef,2,2,J_max+1)
    @inbounds for P = 1:2
        @inbounds for J = 0:J_max
            Ind_Orb_NN[1,P,J+1] = Vector{Vector{Int64}}(undef,N_Orb_NN[1,P,J+1])
            Ind_Orb_NN[2,P,J+1] = Vector{Vector{Int64}}(undef,N_Orb_NN[2,P,J+1])
        end
    end

    # Reinitialize array for counting of NN orbitals ...
    N_Orb_NN = zeros(Int64,2,2,J_max+1)

    # Allocate the dictionary and inidices for NN orbitals ...
    @inbounds for a in 1:a_max
        n_a = Orb[a].n
        l_a = Orb[a].l
        j_a = Orb[a].j
        if (2*n_a + l_a) <= N_max
            @inbounds for b = 1:a_max
                n_b = Orb[b].n
                l_b = Orb[b].l
                j_b = Orb[b].j
                if ((2*(n_a + n_b) + l_a + l_b) <= N_2max)
                    if j_a == j_b
                        P = rem(l_a + l_b, 2) + 1
                        J = 0
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

    # Initialize the W interaction arrays ...
    W_pp = Matrix{Vector{Float64}}(undef,2,1)
    W_ppn = Matrix{Vector{Float64}}(undef,2,1)
    W_npn = Matrix{Vector{Float64}}(undef,2,1)
    W_nn = Matrix{Vector{Float64}}(undef,2,1)
    
    # Allocate entries of W ...
    @inbounds for P = 1:2
        LengthT0 = div(N_Orb_NN[1,P,1]*(N_Orb_NN[1,P,1]+1),2)
        LengthT1 = div(N_Orb_NN[2,P,1]*(N_Orb_NN[2,P,1]+1),2)
        W_pp[P,1] = zeros(Float64,LengthT1)
        W_ppn[P,1] = zeros(Float64,LengthT0)
        W_npn[P,1] = zeros(Float64,LengthT0)
        W_nn[P,1] = zeros(Float64,LengthT1)
    end

    # Allocate W ...
    W_NN = W2B(W_pp,W_ppn,W_npn,W_nn)

    return W_NN
end

function H2b_res_no2b(Params::Parameters,Orb::Vector{NOrb},Orb_NN::NNOrb,Orb_NNN::NNNOrb,V_NN::NNInt,V_NNN::Array{Vector{Vector{Float32}},4},Rho::pnMatrix,Kappa::pnMatrix)
    # Read parameters ...
    N_max = Params.Calc.Nmax
    N_2max = Params.Calc.N2max
    N_3max = Params.Calc.N3max
    a_max = div((N_max + 1)*(N_max + 2),2)
    J_max = N_2max + 1
    c = Params.Calc.cV_res

    # Initialize list of J,P values ...
    JP = JP_Ini(J_max)

    # Initialize the anomal density dependent NN interaction W ...
        # Yet to be written ...
    W_NN = W2b_initialize(Orb,N_max)

    println("\nInitializing the evaluation of the residual density dependent 2-body Hamiltonian ...\n")

    # Make normal density-dependent NN interaction V ...
    println("\nAllocating the normal density-dependent NN interaction V ...")
    @time @inbounds for i in JP
        J, P = i[1], i[2]
        if P == 1
            println("Calculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = +")
        else
            println("Calculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = -")
        end
        N_T0, N_T1 = Orb_NN.N[1,P,J+1], Orb_NN.N[2,P,J+1]

        @inbounds Threads.@threads for Bra in 1:max(N_T0, N_T1)

            if Bra <= N_T0
                a, b = Orb_NN.Ind[1,P,J+1][Bra][1], Orb_NN.Ind[1,P,J+1][Bra][2]
                n_a, l_a = Orb[a].n, Orb[a].l
                n_b, l_b = Orb[b].n, Orb[b].l

                @inbounds for Ket in 1:Bra
                    Ind = Bra + (Ket - 1) * N_T0 - div(Ket * (Ket - 1),2)
                    d, e = Orb_NN.Ind[1,P,J+1][Ket][1], Orb_NN.Ind[1,P,J+1][Ket][2]
                    n_d, l_d = Orb[d].n, Orb[d].l
                    n_e, l_e = Orb[e].n, Orb[e].l
                    pnSum = 0.0

                    @inbounds for c in 1:a_max
                        n_c, l_c, j_c = Orb[c].n, Orb[c].l, Orb[c].j

                        if ((2*(n_a + n_b + n_c) + l_a + l_b + l_c) <= N_3max)
                            P3B = rem(l_a + l_b + l_c, 2) + 1
                            @inbounds for f in 1:a_max
                                n_f, l_f, j_f = Orb[f].n, Orb[f].l, Orb[f].j

                                if (l_c == l_f) && (j_c == j_f) && ((2*(n_d + n_e + n_f) + l_d + l_e + l_f) <= N_3max)
                                    Hat = 1.0 / Float64(2*J + 1)
                                    ME001 = V3B_NO2B(a,b,c,0,d,e,f,0,J,1,P3B,V_NNN,Orb,Orb_NNN)
                                    ME101 = V3B_NO2B(a,b,c,1,d,e,f,0,J,1,P3B,V_NNN,Orb,Orb_NNN)
                                    ME011 = V3B_NO2B(a,b,c,0,d,e,f,1,J,1,P3B,V_NNN,Orb,Orb_NNN)
                                    ME111 = V3B_NO2B(a,b,c,1,d,e,f,1,J,1,P3B,V_NNN,Orb,Orb_NNN)
                                    ME113 = V3B_NO2B(a,b,c,1,d,e,f,1,J,3,P3B,V_NNN,Orb,Orb_NNN)

                                    pnSum += Hat * ((0.5 * ME001 - 1.0/sqrt(12.0) * ME101 - 1.0/sqrt(12.0) * ME011 +
                                                1.0/6.0 * ME111 + 1.0/3.0 * ME113) * Rho.p[c,f] + (0.5 * ME001 + 1.0/sqrt(12.0) *
                                                ME101 + 1.0/sqrt(12.0) * ME011 + 1.0/6.0 * ME111 +  1.0/3.0 * ME113) * Rho.n[c,f])
                                end
                            end
                        end
                    end
                    @views V_NN.pn[P,J+1][Ind] += pnSum
                end
            end

            if Bra <= N_T1
                a, b = Orb_NN.Ind[2,P,J+1][Bra][1], Orb_NN.Ind[2,P,J+1][Bra][2]
                n_a, l_a = Orb[a].n, Orb[a].l
                n_b, l_b = Orb[b].n, Orb[b].l

                @inbounds for Ket in 1:Bra
                    Ind = Bra + (Ket - 1) * N_T1 - div(Ket * (Ket - 1),2)
                    d, e = Orb_NN.Ind[2,P,J+1][Ket][1], Orb_NN.Ind[2,P,J+1][Ket][2]
                    n_d, l_d = Orb[d].n, Orb[d].l
                    n_e, l_e = Orb[e].n, Orb[e].l
                    ppSum, nnSum = 0.0, 0.0

                    @inbounds for c in 1:a_max
                        n_c, l_c, j_c = Orb[c].n, Orb[c].l, Orb[c].j

                        if ((2*(n_a + n_b + n_c) + l_a + l_b + l_c) <= N_3max)
                            P3B = rem(l_a + l_b + l_c, 2) + 1
                            @inbounds for f in 1:a_max
                                n_f, l_f, j_f = Orb[f].n, Orb[f].l, Orb[f].j

                                if (l_c == l_f) && (j_c == j_f) && ((2*(n_d + n_e + n_f) + l_d + l_e + l_f) <= N_3max)

                                    Hat = 1.0 / Float64(2*J + 1)
                                    ME111 = V3B_NO2B(a,b,c,1,d,e,f,1,J,1,P3B,V_NNN,Orb,Orb_NNN)
                                    ME113 = V3B_NO2B(a,b,c,1,d,e,f,1,J,3,P3B,V_NNN,Orb,Orb_NNN)

                                    ppSum += Hat * (ME113 * Rho.p[c,f] + (2.0/3.0 * ME111 + 1.0/3.0 * ME113) * Rho.n[c,f])
                                    nnSum += Hat * (ME113 * Rho.n[c,f] + (2.0/3.0 * ME111 + 1.0/3.0 * ME113) * Rho.p[c,f])

                                end
                            end
                        end
                    end
                    @views V_NN.pp[P,J+1][Ind] += ppSum
                    @views V_NN.nn[P,J+1][Ind] += nnSum
                end
    
            end

        end
    end

    # Make anomal density-dependent NN interaction W ...
    println("\nAllocating the anomal density-dependent NN interaction W ...")
    @time @inbounds for P in 1:2
        if P == 1
            println("Calculating   ...   J = " * string(0) * "/" * string(N_2max+1) * "\tP = +")
        else
            println("Calculating   ...   J = " * string(0) * "/" * string(N_2max+1) * "\tP = -")
        end
        N_T0, N_T1 = Orb_NN.N[1,P,1], Orb_NN.N[2,P,1]

        @inbounds Threads.@threads for Bra in 1:max(N_T0, N_T1)

            if Bra <= N_T0
                a, b = Orb_NN.Ind[1,P,1][Bra][1], Orb_NN.Ind[1,P,1][Bra][2]
                n_a, l_a = Orb[a].n, Orb[a].l
                n_b, l_b = Orb[b].n, Orb[b].l

                @inbounds for Ket in 1:Bra
                    Ind = Bra + (Ket - 1) * N_T0 - div(Ket * (Ket - 1),2)
                    d, e = Orb_NN.Ind[1,P,1][Ket][1], Orb_NN.Ind[1,P,1][Ket][2]
                    n_d, l_d = Orb[d].n, Orb[d].l
                    n_e, l_e = Orb[e].n, Orb[e].l

                    ppnSum, npnSum = 0.0, 0.0

                    @inbounds for c in 1:a_max
                        n_c, l_c, j_c = Orb[c].n, Orb[c].l, Orb[c].j

                        if ((2*(n_a + n_b + n_c) + l_a + l_b + l_c) <= N_3max)
                            @inbounds for f in 1:a_max
                                n_f, l_f, j_f = Orb[f].n, Orb[f].l, Orb[f].j

                                if (l_c == l_f) && (j_c == j_f) && ((2*(n_d + n_e + n_f) + l_d + l_e + l_f) <= N_3max)
                                    if Orb[b].j == Orb[e].j && l_b == l_e
                                        ppnP3B = rem(l_a + l_d + l_b, 2) + 1
                                        ppnME111 = V3B_NO2B(a,d,b,1,c,f,e,1,0,1,ppnP3B,V_NNN,Orb,Orb_NNN)
                                        ppnME113 = V3B_NO2B(a,d,b,1,c,f,e,1,0,3,ppnP3B,V_NNN,Orb,Orb_NNN)
                                        ppnSum +=  (2.0 * ppnME111 + 1.0 * ppnME113) / 3.0 * Kappa.p[c,f]
                                    end

                                    if Orb[a].j == Orb[d].j && l_a == l_d
                                        npnP3B = rem(l_b + l_e + l_a, 2) + 1
                                        npnME001 = V3B_NO2B(b,e,a,0,c,f,d,0,0,1,npnP3B,V_NNN,Orb,Orb_NNN)
                                        npnME101 = V3B_NO2B(b,e,a,1,c,f,d,0,0,1,npnP3B,V_NNN,Orb,Orb_NNN)
                                        npnME011 = V3B_NO2B(b,e,a,0,c,f,d,1,0,1,npnP3B,V_NNN,Orb,Orb_NNN)
                                        npnME111 = V3B_NO2B(b,e,a,1,c,f,d,1,0,1,npnP3B,V_NNN,Orb,Orb_NNN)
                                        npnME113 = V3B_NO2B(b,e,a,1,c,f,d,1,0,3,npnP3B,V_NNN,Orb,Orb_NNN)

                                        ppnSum +=  (2.0 * ppnME111 + 1.0 * ppnME113) / 3.0 * Kappa.p[c,f]
                                        npnSum +=  (0.5 * npnME001 - 1.0 / sqrt(12.0) * npnME101 - 1.0 / sqrt(12.0) * npnME011 +
                                                    1.0 / 6.0 * npnME111 + 1.0 / 3.0 * npnME113) * Kappa.n[c,f]
                                    end

                                end
                            end
                        end
                    end
                    @views W_NN.ppn[P,1][Ind] = ppnSum
                    @views W_NN.npn[P,1][Ind] = npnSum
                end

            end

            if Bra <= N_T1
                a, b = Orb_NN.Ind[2,P,1][Bra][1], Orb_NN.Ind[2,P,1][Bra][2]
                n_a, l_a = Orb[a].n, Orb[a].l
                n_b, l_b = Orb[b].n, Orb[b].l

                @inbounds for Ket in 1:Bra
                    Ind = Bra + (Ket - 1) * N_T1 - div(Ket * (Ket - 1),2)
                    d, e = Orb_NN.Ind[2,P,1][Ket][1], Orb_NN.Ind[2,P,1][Ket][2]
                    n_d, l_d = Orb[d].n, Orb[d].l
                    n_e, l_e = Orb[e].n, Orb[e].l
                    P3B = rem(l_a + l_b + l_d, 2) + 1
                    ppSum, nnSum = 0.0, 0.0

                    @inbounds for c in 1:a_max
                        n_c, l_c, j_c = Orb[c].n, Orb[c].l, Orb[c].j

                        if ((2*(n_a + n_b + n_c) + l_a + l_b + l_c) <= N_3max)
                            @inbounds for f in 1:a_max
                                n_f, l_f, j_f = Orb[f].n, Orb[f].l, Orb[f].j

                                if (l_c == l_f) && (j_c == j_f) && ((2*(n_d + n_e + n_f) + l_d + l_e + l_f) <= N_3max)
                                    if Orb[d].j == Orb[e].j && l_d == l_e
                                        ME113 = V3B_NO2B(a,b,d,1,c,f,e,1,0,3,P3B,V_NNN,Orb,Orb_NNN)

                                        ppSum += ME113 * Kappa.p[c,f]
                                        nnSum += ME113 * Kappa.n[c,f]
                                    end

                                end
                            end
                        end
                    end
                    @views W_NN.pp[P,1][Ind] = ppSum
                    @views W_NN.nn[P,1][Ind] = nnSum
                end
            end

        end
    end

    return V_NN, W_NN
end

function V2b_temp_initialize(Orb::Vector{NOrb},N_max::Int64;Make_Orb_NN::Bool=false)
    # Read parameters ...
    N_2max = 2 * N_max
    a_max = div((N_max + 1)*(N_max + 2),2)
    J_max = 2 * N_max + 1

    # Initiliaze array for counting of NN orbitals ...
    N_Orb_NN = zeros(Int64,2,J_max+1)

    # Initiliaze dictionary for NN orbitals ...
    Orb_NN_Dic = Dict{Tuple{Int8,Int8,Int16,Int16},Int32}()

    # Count the number of NN orbitals ...
    @inbounds for a in 1:a_max
        n_a = Orb[a].n
        l_a = Orb[a].l
        j_a = Orb[a].j
        if (2*n_a + l_a) <= N_max
            @inbounds for b = 1:a_max
                n_b = Orb[b].n
                l_b = Orb[b].l
                j_b = Orb[b].j
                if ((2*(n_a + n_b) + l_a + l_b) <= N_2max)
                    P = rem(l_a + l_b, 2) + 1
                    @inbounds for J = Int64(abs(j_a - j_b)/2):Int64((j_a + j_b)/2)
                        N_Orb_NN[P,J+1] += 1
                    end
                end
            end
        end
    end

    # Initialite the array for NN orbital indices ...
    Ind_Orb_NN = Matrix{Vector{Vector{Int64}}}(undef,2,J_max+1)
    @inbounds for P = 1:2
        @inbounds for J = 0:J_max
            Ind_Orb_NN[P,J+1] = Vector{Vector{Int64}}(undef,N_Orb_NN[P,J+1])
        end
    end

    # Reinitialize array for counting of NN orbitals ...
    N_Orb_NN = zeros(Int64,2,J_max+1)

    # Allocate the dictionary and inidices for NN orbitals ...
    @inbounds for a in 1:a_max
        n_a = Orb[a].n
        l_a = Orb[a].l
        j_a = Orb[a].j
        if (2*n_a + l_a) <= N_max
            @inbounds for b = 1:a_max
                n_b = Orb[b].n
                l_b = Orb[b].l
                j_b = Orb[b].j
                if ((2*(n_a + n_b) + l_a + l_b) <= N_2max)
                    P = rem(l_a + l_b, 2) + 1
                    @inbounds for J = Int64(abs(j_a - j_b)/2):Int64((j_a + j_b)/2)
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

    # Initialize the NN interaction arrays for V ...
    V_pp = Matrix{Matrix{Float64}}(undef,2,J_max+1)
    V_pn = Matrix{Matrix{Float64}}(undef,2,J_max+1)
    V_nn = Matrix{Matrix{Float64}}(undef,2,J_max+1)
    
    # Allocate entries of V ...
    @inbounds for P = 1:2
        @inbounds Threads.@threads for J=0:J_max
            Length = N_Orb_NN[P,J+1]
            V_pp[P,J+1] = zeros(Float64,Length,Length)
            V_pn[P,J+1] = zeros(Float64,Length,Length)
            V_nn[P,J+1] = zeros(Float64,Length,Length)
        end
    end

    # Allocate V & W ...
    V_NN_t = V2B_Temp(V_pp,V_pn,V_nn)

    if Make_Orb_NN == false
        return V_NN_t
    elseif Make_Orb_NN == true
        # Allocate the array for NN orbitals ...
        Orb_NN_t = NNOrb_Temp(Orb_NN_Dic,N_Orb_NN,Ind_Orb_NN)

        return V_NN_t, Orb_NN_t
    end

end

function W2b_temp_initialize(Orb::Vector{NOrb},N_max::Int64)
    # Read parameters ...
    N_2max = 2 * N_max
    a_max = div((N_max + 1)*(N_max + 2),2)
    J_max = 2*N_2max + 1

    # Initiliaze array for counting of NN orbitals ...
    N_Orb_NN = zeros(Int64,2,J_max+1)

    # Initiliaze dictionary for NN orbitals ...
    Orb_NN_Dic = Dict{Tuple{Int8,Int8,Int16,Int16},Int32}()

    # Count the number of NN orbitals ...
    @inbounds for a in 1:a_max
        n_a = Orb[a].n
        l_a = Orb[a].l
        j_a = Orb[a].j
        if (2*n_a + l_a) <= N_max
            @inbounds for b = 1:a_max
                n_b = Orb[b].n
                l_b = Orb[b].l
                j_b = Orb[b].j
                if ((2*(n_a + n_b) + l_a + l_b) <= N_2max)
                    P = rem(l_a + l_b, 2) + 1
                    @inbounds for J = Int64(abs(j_a - j_b)/2):Int64((j_a + j_b)/2)
                        N_Orb_NN[P,J+1] += 1
                    end
                end
            end
        end
    end

    # Initialite the array for NN orbital indices ...
    Ind_Orb_NN = Matrix{Vector{Vector{Int64}}}(undef,2,J_max+1)
    @inbounds for P = 1:2
        @inbounds for J = 0:J_max
            Ind_Orb_NN[P,J+1] = Vector{Vector{Int64}}(undef,N_Orb_NN[P,J+1])
        end
    end

    # Reinitialize array for counting of NN orbitals ...
    N_Orb_NN = zeros(Int64,2,J_max+1)

    # Allocate the dictionary and inidices for NN orbitals ...
    @inbounds for a in 1:a_max
        n_a = Orb[a].n
        l_a = Orb[a].l
        j_a = Orb[a].j
        if (2*n_a + l_a) <= N_max
            @inbounds for b = 1:a_max
                n_b = Orb[b].n
                l_b = Orb[b].l
                j_b = Orb[b].j
                if ((2*(n_a + n_b) + l_a + l_b) <= N_2max)
                    P = rem(l_a + l_b, 2) + 1
                    @inbounds for J = Int64(abs(j_a - j_b)/2):Int64((j_a + j_b)/2)
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

    # Initialize the NN interaction arrays for W ...
    W_pp = Matrix{Matrix{Float64}}(undef,2,1)
    W_ppn = Matrix{Matrix{Float64}}(undef,2,1)
    W_npn = Matrix{Matrix{Float64}}(undef,2,1)
    W_nn = Matrix{Matrix{Float64}}(undef,2,1)
    
    # Allocate entries of W ...
    @inbounds for P = 1:2
        W_pp[P,1] = zeros(Float64,N_Orb_NN[P,1],N_Orb_NN[P,1])
        W_ppn[P,1] = zeros(Float64,N_Orb_NN[P,1],N_Orb_NN[P,1])
        W_npn[P,1] = zeros(Float64,N_Orb_NN[P,1],N_Orb_NN[P,1])
        W_nn[P,1] = zeros(Float64,N_Orb_NN[P,1],N_Orb_NN[P,1])
    end

    # Allocate V & W ...
    W_NN_t = W2B_Temp(W_pp,W_ppn,W_npn,W_nn)

    return W_NN_t
end

function V2b_temp_index(a::Int64,b::Int64,J::Int64,P::Int64,Orb_NN_t::NNOrb_Temp)
    P, J, a, b = Int8(P), Int8(J), Int16(a), Int16(b)
    Ind = Int64(Orb_NN_t.Dic[(P,J,a,b)])
    return Ind
end

function qpH_initialize(Params::Parameters,Orb::Vector{NOrb})
    # Read parameters ...
    N_max = Params.Calc.Nmax
    N_2max = 2 * N_max
    a_max = div((N_max + 1)*(N_max + 2),2)
    J_max = Int64(2*N_max + 1)

    # Initiliaze array for counting of NN orbitals ...
    N_Orb_NN = zeros(Int64,2,2,J_max+1)

    # Initiliaze dictionary for NN orbitals ...
    Orb_NN_Dic = Dict{Tuple{Int8,Int8,Int8,Int16,Int16},Int32}()

    # Count the number of NN orbitals ...
    @inbounds for a in 1:a_max
        n_a = Orb[a].n
        l_a = Orb[a].l
        j_a = Orb[a].j
        if (2*n_a + l_a) <= N_max
            @inbounds for b = 1:a_max
                n_b = Orb[b].n
                l_b = Orb[b].l
                j_b = Orb[b].j
                if ((2*(n_a + n_b) + l_a + l_b) <= N_2max)
                    P = rem(l_a + l_b, 2) + 1
                    @inbounds for J = div(abs(j_a - j_b),2):div((j_a + j_b),2)
                        if b <= a
                            N_Orb_NN[2,P,J+1] += 1
                        end
                        N_Orb_NN[1,P,J+1] += 1
                    end
                end
            end
        end
    end

    # Initialite the array for NN orbital indices ...
    Ind_Orb_NN = Array{Vector{Vector{Int64}}}(undef,2,2,J_max+1)
    @inbounds for P = 1:2
        @inbounds for J = 0:J_max
            Ind_Orb_NN[1,P,J+1] = Vector{Vector{Int64}}(undef,N_Orb_NN[1,P,J+1])
            Ind_Orb_NN[2,P,J+1] = Vector{Vector{Int64}}(undef,N_Orb_NN[2,P,J+1])
        end
    end

    # Reinitialize array for counting of NN orbitals ...
    N_Orb_NN = zeros(Int64,2,2,J_max+1)

    # Allocate the dictionary and inidices for NN orbitals ...
    @inbounds for a in 1:a_max
        n_a = Orb[a].n
        l_a = Orb[a].l
        j_a = Orb[a].j
        if (2*n_a + l_a) <= N_max
            @inbounds for b = 1:a_max
                n_b = Orb[b].n
                l_b = Orb[b].l
                j_b = Orb[b].j
                if ((2*(n_a + n_b) + l_a + l_b) <= N_2max)
                    P = rem(l_a + l_b, 2) + 1
                    @inbounds for J = Int64(abs(j_a - j_b)/2):Int64((j_a + j_b)/2)
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

    # Initialize the NN interaction matrices ...
    H40_pp = Matrix{Vector{Float64}}(undef,2,J_max+1)
    H40_pn = Matrix{Vector{Float64}}(undef,2,J_max+1)
    H40_nn = Matrix{Vector{Float64}}(undef,2,J_max+1)
    H31_pp = Matrix{Vector{Float64}}(undef,2,J_max+1)
    H31_pn_2011 = Matrix{Vector{Float64}}(undef,2,J_max+1)
    H31_pn_1120 = Matrix{Vector{Float64}}(undef,2,J_max+1)
    H31_nn = Matrix{Vector{Float64}}(undef,2,J_max+1)
    H22_pp = Matrix{Vector{Float64}}(undef,2,J_max+1)
    H22_pn2002 = Matrix{Vector{Float64}}(undef,2,J_max+1)
    H22_pn1111 = Matrix{Vector{Float64}}(undef,2,J_max+1)
    H22_pn0220 = Matrix{Vector{Float64}}(undef,2,J_max+1)
    H22_nn = Matrix{Vector{Float64}}(undef,2,J_max+1)
    
    # Allocate entries of H ...
    @inbounds for P = 1:2
        @inbounds Threads.@threads for J=0:J_max
            Length_T0 = div(N_Orb_NN[1,P,J+1]*(N_Orb_NN[1,P,J+1]+1),2)
            Length_T1 = div(N_Orb_NN[2,P,J+1]*(N_Orb_NN[2,P,J+1]+1),2)

            H40_pp[P,J+1] = zeros(Float64,Length_T1)
            H40_pn[P,J+1] = zeros(Float64,Length_T0)
            H40_nn[P,J+1] = zeros(Float64,Length_T1)
            H31_pp[P,J+1] = zeros(Float64,Length_T1)
            H31_pn_2011[P,J+1] = zeros(Float64,Length_T0)
            H31_pn_1120[P,J+1] = zeros(Float64,Length_T0)
            H31_nn[P,J+1] = zeros(Float64,Length_T1)
            H22_pp[P,J+1] = zeros(Float64,Length_T1)
            H22_pn2002[P,J+1] = zeros(Float64,Length_T0)
            H22_pn1111[P,J+1] = zeros(Float64,Length_T0)
            H22_pn0220[P,J+1] = zeros(Float64,Length_T0)
            H22_nn[P,J+1] = zeros(Float64,Length_T1)
        end
    end

    # Allocate the NN interaction in quasiparticle picture ...
    H_NN = qpH2B(qpH40(H40_pp,H40_pn,H40_nn),
                 qpH31(H31_pp,H31_pn_2011,H31_pn_1120,H31_nn),
                 qpH22(H22_pp,H22_pn2002,H22_pn1111,H22_pn0220,H22_nn))

    return H_NN
end

function qpH(Params::Parameters,Orb::Vector{NOrb},Orb_NN::NNOrb,V_NN::NNInt,W_NN::W2B,C::pnMatrix,U::pnMatrix,V::pnMatrix)
    println("\nPreparing the Hamiltonian expressed in the J-scheme quasiparticle basis ...")

    # Perform the Canonical Transformation of V & W into the canonical basis ...
    V_NN, W_NN = qpH_canonical_transformation(Params,Orb,Orb_NN,V_NN,W_NN,C)

    # Initialize the 2-body quasiparticle Hamiltonian H^(kl) ...
    H_NN = qpH_initialize(Params,Orb)

    # Allocate H^(40) ...
    H_NN = qpH_allocate_H40(Params,Orb,Orb_NN,H_NN,V_NN,W_NN,U,V)

    # Allocate H^(31) ...
    H_NN = qpH_allocate_H31(Params,Orb,Orb_NN,H_NN,V_NN,U,V)

    # Allocate H^(22) ...
    H_NN = qpH_allocate_H22(Params,Orb,Orb_NN,H_NN,V_NN,U,V)
    
    println("\nQuasiparticle 2-body Hamiltonian H^(kl) has been fully initialized and allocated ...\n")

    return H_NN
end

function qpH_canonical_transformation(Params::Parameters,Orb::Vector{NOrb},Orb_NN::NNOrb,V_NN::NNInt,W_NN::W2B,C::pnMatrix)
    # Transformation of V_NN & W_NN to the canonical HFB basis...
    println("\nPerforming the Canonical Transformation of residual NN V & W interactions into the target canonical basis ...")

    # Preallocate arrays of allowed values of J & P ...
    JP = JP_initialize(Params.Calc.N2max + 1)

    # Transformation of the 1st index ...
    V_NN_t1, W_NN_t1, Orb_NN_t = qpH_canonical_transformation_ind1(Params,JP,Orb,Orb_NN,V_NN,W_NN,C)

    # Transformation of the 2nd index ...
    V_NN_t2, W_NN_t2 = qpH_canonical_transformation_ind2(Params,JP,Orb,Orb_NN_t,V_NN_t1,W_NN_t1,C)

    # Transformation of the 3rd index ...
    V_NN_t1, W_NN_t1 = qpH_canonical_transformation_ind3(Params,JP,Orb,Orb_NN_t,V_NN_t1,W_NN_t1,V_NN_t2,W_NN_t2,C)

    # Transformation of the 4th index ...
    V_NN, W_NN = qpH_canonical_transformation_ind4(Params,JP,Orb,Orb_NN,Orb_NN_t,V_NN,W_NN,V_NN_t1,W_NN_t1,C)

    # Deallocate temporary interaction arrays ...
    V_NN_t1, W_NN_t1 = nothing, nothing
    V_NN_t2, W_NN_t2 = nothing, nothing

    # Perform the Garbage collection ...
    GC.gc()
    
    println("\nResidual NN 2-body interaction V & W ready ...\n")

    return V_NN, W_NN
end

function qpH_canonical_transformation_ind1(Params::Parameters,JP::Vector{Vector{Int64}},Orb::Vector{NOrb},Orb_NN::NNOrb,V_NN::NNInt,W_NN::W2B,C::pnMatrix)
    # Parameter initialization...
    N_max = Params.Calc.Nmax
    N_2max = 2*N_max
    a_max = div((N_max + 1)*(N_max + 2),2)
    J_max = N_2max + 1

    # Make temporary arrays for NN interactions ...
        # Case of V ...
    V_NN_t1, Orb_NN_t = V2b_temp_initialize(Orb,N_max,Make_Orb_NN = true)
        # Case of W ...
    W_NN_t1 = W2b_temp_initialize(Orb,N_max)

    # Define function for the transformation MEs for the 1st index ...
        # Case of V ...
    @inline function V2b_temp_transformation_ind1_MEs(J::Int64,i::Int64,a::Int64,b::Int64,c::Int64,d::Int64,C::pnMatrix,V_NN::NNInt,Orb::Vector{NOrb},Orb_NN::NNOrb)
        return @views C.p[i,a] * V2B(i,b,c,d,J,1,V_NN.pp,Orb,Orb_NN), C.p[i,a] * V2B(i,b,c,d,J,0,V_NN.pn,Orb,Orb_NN),
            C.n[i,a] * V2B(i,b,c,d,J,1,V_NN.nn,Orb,Orb_NN)
    end
        # Case of W ...
    @inline function W2b_temp_transformation_ind1_MEs(i::Int64,a::Int64,b::Int64,c::Int64,d::Int64,C::pnMatrix,W_NN::W2B,Orb::Vector{NOrb},Orb_NN::NNOrb)
        return @views C.p[i,a] * V2B(i,b,c,d,0,1,W_NN.pp,Orb,Orb_NN), C.p[i,a] * V2B(i,b,c,d,0,0,W_NN.ppn,Orb,Orb_NN),
                      C.p[i,a] * V2B(i,b,c,d,0,0,W_NN.npn,Orb,Orb_NN), C.n[i,a] * V2B(i,b,c,d,0,1,W_NN.nn,Orb,Orb_NN)
    end

    # Transformation of the 1st index ...
    println("\nPerforming Canonical Transformation of the residual 2-body NN interaction in the 1st index...")
    @inbounds Threads.@threads for i in JP
        J, P = i[1], i[2]
        if P == 1
            println("Calculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = +")
        else
            println("Calculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = -")
        end
        N = Orb_NN_t.N[P,J+1]
        @inbounds for Bra in 1:N
            a = Orb_NN_t.Ind[P,J+1][Bra][1]
            b = Orb_NN_t.Ind[P,J+1][Bra][2]
            l_a, j_a = Orb[a].l, Orb[a].j

            Orb_x = Orb_PreComp(a_max,j_a,l_a,Orb)
            @inbounds for Ket in 1:N
                c = Orb_NN_t.Ind[P,J+1][Ket][1]
                d = Orb_NN_t.Ind[P,J+1][Ket][2]

                VppSum, VpnSum, VnnSum = 0.0, 0.0, 0.0
                if J == 0
                    WppSum, WppnSum, WnpnSum, WnnSum = 0.0, 0.0, 0.0, 0.0
                end
                
                @inbounds for i in Orb_x
                    VppME, VpnME, VnnME = V2b_temp_transformation_ind1_MEs(J,i,a,b,c,d,C,V_NN,Orb,Orb_NN)
                    VppSum += VppME
                    VpnSum += VpnME
                    VnnSum += VnnME

                    if J == 0
                        WppME, WppnME, WnpnME, WnnME = W2b_temp_transformation_ind1_MEs(i,a,b,c,d,C,W_NN,Orb,Orb_NN)
                        WppSum += WppME
                        WppnSum += WppnME
                        WnpnSum += WnpnME
                        WnnSum += WnnME
                    end

                end
                @views V_NN_t1.pp[P,J+1][Bra,Ket] = VppSum
                @views V_NN_t1.pn[P,J+1][Bra,Ket] = VpnSum
                @views V_NN_t1.nn[P,J+1][Bra,Ket] = VnnSum

                if J == 0
                    @views W_NN_t1.pp[P,1][Bra,Ket] = WppSum
                    @views W_NN_t1.ppn[P,1][Bra,Ket] = WppnSum
                    @views W_NN_t1.npn[P,1][Bra,Ket] = WnpnSum
                    @views W_NN_t1.nn[P,1][Bra,Ket] = WnnSum
                end

            end
        end
    end

    return V_NN_t1, W_NN_t1, Orb_NN_t
end

function qpH_canonical_transformation_ind2(Params::Parameters,JP::Vector{Vector{Int64}},Orb::Vector{NOrb},Orb_NN_t::NNOrb_Temp,V_NN_t1::V2B_Temp,W_NN_t1::W2B_Temp,C::pnMatrix)
    # Parameter initialization...
    N_max = Params.Calc.Nmax
    N_2max = 2*N_max
    a_max = div((N_max + 1)*(N_max + 2),2)
    J_max = N_2max + 1

    # Make another temporary arrays for NN interactions ...
        # Case of V ...
    V_NN_t2, Orb_NN_t2 = V2b_temp_initialize(Orb,N_max,Make_Orb_NN = true)
        # Case of W ...
    W_NN_t2 = W2b_temp_initialize(Orb,N_max)

    # Define function for the transformation MEs for the 2nd index ...
        # Case of V ...
    @inline function V2b_temp_transformation_ind2_MEs(P::Int64,J::Int64,Bra::Int64,Ket::Int64,i::Int64,b::Int64,C::pnMatrix,V_NN_t::V2B_Temp)
        return @views C.p[i,b] * V_NN_t.pp[P,J+1][Bra,Ket], C.n[i,b] * V_NN_t.pn[P,J+1][Bra,Ket], C.n[i,b] * V_NN_t.nn[P,J+1][Bra,Ket]
    end
        # Case of W ...
    @inline function W2b_temp_transformation_ind2_MEs(P::Int64,Bra::Int64,Ket::Int64,i::Int64,b::Int64,C::pnMatrix,W_NN_t::W2B_Temp)
        return C.p[i,b] * W_NN_t.pp[P,1][Bra,Ket], C.n[i,b] * W_NN_t.ppn[P,1][Bra,Ket], C.n[i,b] * W_NN_t.npn[P,1][Bra,Ket], C.n[i,b] * W_NN_t.nn[P,1][Bra,Ket]
    end

    # Transformation of the 2nd index ...
    println("\nPerforming Canonical Transformation of the residual 2-body NN interaction in the 2nd index...")

    @time @inbounds Threads.@threads for i in JP
        J, P = i[1], i[2]
        if P == 1
            println("Calculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = +")
        else
            println("Calculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = -")
        end
        N = Orb_NN_t.N[P,J+1]
        @inbounds for Bra in 1:N
            a = Orb_NN_t.Ind[P,J+1][Bra][1]
            b = Orb_NN_t.Ind[P,J+1][Bra][2]
            l_b, j_b = Orb[b].l, Orb[b].j

            Orb_x = Orb_PreComp(a_max,j_b,l_b,Orb)
            @inbounds for Ket in 1:N
                c = Orb_NN_t.Ind[P,J+1][Ket][1]
                d = Orb_NN_t.Ind[P,J+1][Ket][2]

                VppSum, VpnSum, VnnSum = 0.0, 0.0, 0.0
                if J == 0
                    WppSum, WppnSum, WnpnSum, WnnSum = 0.0, 0.0, 0.0, 0.0
                end
                
                @inbounds for i in Orb_x
                    Bra_ai = V2b_temp_index(a,i,J,P,Orb_NN_t)
                    VppME, VpnME, VnnME = V2b_temp_transformation_ind2_MEs(P,J,Bra_ai,Ket,i,b,C,V_NN_t1)
                    VppSum += VppME
                    VpnSum += VpnME
                    VnnSum += VnnME
                    if J == 0
                        WppME, WppnME, WnpnME, WnnME = W2b_temp_transformation_ind2_MEs(P,Bra_ai,Ket,i,b,C,W_NN_t1)
                        WppSum += WppME
                        WppnSum += WppnME
                        WnpnSum += WnpnME
                        WnnSum += WnnME
                    end
                end
                @views V_NN_t2.pp[P,J+1][Bra,Ket] = VppSum
                @views V_NN_t2.pn[P,J+1][Bra,Ket] = VpnSum
                @views V_NN_t2.nn[P,J+1][Bra,Ket] = VnnSum
                if J == 0
                    @views W_NN_t2.pp[P,1][Bra,Ket] = WppSum
                    @views W_NN_t2.ppn[P,1][Bra,Ket] = WppnSum
                    @views W_NN_t2.npn[P,1][Bra,Ket] = WnpnSum
                    @views W_NN_t2.nn[P,1][Bra,Ket] = WnnSum
                end
            end
        end
    end

    return V_NN_t2, W_NN_t2
end

function qpH_canonical_transformation_ind3(Params::Parameters,JP::Vector{Vector{Int64}},Orb::Vector{NOrb},Orb_NN_t::NNOrb_Temp,V_NN_t1::V2B_Temp,W_NN_t1::W2B_Temp,V_NN_t2::V2B_Temp,W_NN_t2::W2B_Temp,C::pnMatrix)
    # Parameter initialization...
    N_max = Params.Calc.Nmax
    N_2max = 2*N_max
    a_max = div((N_max + 1)*(N_max + 2),2)
    J_max = N_2max + 1

    # Define function for the transformation MEs for the 3rd index ...
        # Case of V ...
    @inline function V2b_temp_transformation_ind3_MEs(P::Int64,J::Int64,Bra::Int64,Ket::Int64,i::Int64,c::Int64,C::pnMatrix,V_NN_t::V2B_Temp)
        return @views C.p[i,c] * V_NN_t.pp[P,J+1][Bra,Ket], C.p[i,c] * V_NN_t.pn[P,J+1][Bra,Ket], C.n[i,c] * V_NN_t.nn[P,J+1][Bra,Ket]
    end
        # Case of W ...
    @inline function W2b_temp_transformation_ind3_MEs(P::Int64,Bra::Int64,Ket::Int64,i::Int64,c::Int64,C::pnMatrix,W_NN_t::W2B_Temp)
        return @views C.p[i,c] * W_NN_t.pp[P,1][Bra,Ket], C.p[i,c] * W_NN_t.ppn[P,1][Bra,Ket], C.p[i,c] * W_NN_t.npn[P,1][Bra,Ket], C.n[i,c] * W_NN_t.nn[P,1][Bra,Ket]
    end

    # Transformation of the 3rd index ...
    println("\nPerforming Canonical Transformation of the residual 2-body NN interaction in the 3rd index...")

    @time @inbounds Threads.@threads for i in JP
        J, P = i[1], i[2]
        if P == 1
            println("Calculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = +")
        else
            println("Calculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = -")
        end
        N = Orb_NN_t.N[P,J+1]
        @inbounds for Ket in 1:N
            c = Orb_NN_t.Ind[P,J+1][Ket][1]
            d = Orb_NN_t.Ind[P,J+1][Ket][2]
            l_c, j_c = Orb[c].l, Orb[c].j

            Orb_x = Orb_PreComp(a_max,j_c,l_c,Orb)

            @inbounds for Bra in 1:N
                VppSum, VpnSum, VnnSum = 0.0, 0.0, 0.0
                if J == 0
                    WppSum, WppnSum, WnpnSum, WnnSum = 0.0, 0.0, 0.0, 0.0
                end
                @inbounds for i in Orb_x
                    Ket_id = V2b_temp_index(i,d,J,P,Orb_NN_t)
                    VppME, VpnME, VnnME = V2b_temp_transformation_ind3_MEs(P,J,Bra,Ket_id,i,c,C,V_NN_t2)
                    VppSum += VppME
                    VpnSum += VpnME
                    VnnSum += VnnME
                    if J == 0
                        WppME, WppnME, WnpnME, WnnME = W2b_temp_transformation_ind3_MEs(P,Bra,Ket_id,i,c,C,W_NN_t2)
                        WppSum += WppME
                        WppnSum += WppnME
                        WnpnSum += WnpnME
                        WnnSum += WnnME
                    end
                end
                @views V_NN_t1.pp[P,J+1][Bra,Ket] = VppSum
                @views V_NN_t1.pn[P,J+1][Bra,Ket] = VpnSum
                @views V_NN_t1.nn[P,J+1][Bra,Ket] = VnnSum
                if J == 0
                    @views W_NN_t1.pp[P,1][Bra,Ket] = WppSum
                    @views W_NN_t1.ppn[P,1][Bra,Ket] = WppnSum
                    @views W_NN_t1.npn[P,1][Bra,Ket] = WnpnSum
                    @views W_NN_t1.nn[P,1][Bra,Ket] = WnnSum
                end
            end
        end
    end

    return V_NN_t1, W_NN_t1
end

function qpH_canonical_transformation_ind4(Params::Parameters,JP::Vector{Vector{Int64}},Orb::Vector{NOrb},Orb_NN::NNOrb,Orb_NN_t::NNOrb_Temp,V_NN::NNInt,W_NN::W2B,V_NN_t1::V2B_Temp,W_NN_t1::W2B_Temp,C::pnMatrix)
    # Parameter initialization...
    N_max = Params.Calc.Nmax
    N_2max = 2*N_max
    a_max = div((N_max + 1)*(N_max + 2),2)
    J_max = N_2max + 1

    # Define function for the transformation MEs for the 4th index ...
        # Case of V ... T = 0 ...
    @inline function V2b_temp_transformation_ind4_T0_MEs(P::Int64,J::Int64,Bra::Int64,Ket::Int64,i::Int64,d::Int64,C::pnMatrix,V_NN_t::V2B_Temp)
        return @views C.n[i,d] * V_NN_t.pn[P,J+1][Bra,Ket]
    end
        # Case of V ... T = 1 ...
    @inline function V2b_temp_transformation_ind4_T1_MEs(P::Int64,J::Int64,Bra::Int64,Ket::Int64,i::Int64,d::Int64,C::pnMatrix,V_NN_t::V2B_Temp)
        return @views C.p[i,d] * V_NN_t.pp[P,J+1][Bra,Ket], C.n[i,d] * V_NN_t.nn[P,J+1][Bra,Ket]
    end
        # Case of W ... T = 0 ...
    @inline function W2b_temp_transformation_ind4_T0_MEs(P::Int64,Bra::Int64,Ket::Int64,i::Int64,d::Int64,C::pnMatrix,W_NN_t::W2B_Temp)
        return @views C.n[i,d] * W_NN_t.ppn[P,1][Bra,Ket], C.n[i,d] * W_NN_t.npn[P,1][Bra,Ket]
    end
        # Case of W ... T = 1 ...
    @inline function W2b_temp_transformation_ind4_T1_MEs(P::Int64,Bra::Int64,Ket::Int64,i::Int64,d::Int64,C::pnMatrix,W_NN_t::W2B_Temp)
        return @views C.p[i,d] * W_NN_t.pp[P,1][Bra,Ket], C.n[i,d] * W_NN_t.nn[P,1][Bra,Ket]
    end

    # Index 4
    println("\nPerforming Canonical Transformation of the residual 2-body NN interaction in the 4th index...")

    @time @inbounds Threads.@threads for i in JP
        J, P = i[1], i[2]
        if P == 1
            println("Calculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = +")
        else
            println("Calculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = -")
        end
        N_T0, N_T1 = Orb_NN.N[1,P,J+1], Orb_NN.N[2,P,J+1]
        @inbounds for Bra in 1:max(N_T0, N_T1)
            @inbounds for Ket in 1:Bra

                # Case of pn interaction ... T = 0
                if Bra <= N_T0
                    Ind = Bra + (Ket - 1) * N_T0 - div(Ket * (Ket - 1),2)
                    a, b = Orb_NN.Ind[1,P,J+1][Bra][1], Orb_NN.Ind[1,P,J+1][Bra][2]
                    c, d = Orb_NN.Ind[1,P,J+1][Ket][1], Orb_NN.Ind[1,P,J+1][Ket][2]
                    l_d, j_d = Orb[d].l, Orb[d].j

                    Bra_ab = V2b_temp_index(a,b,J,P,Orb_NN_t)
                    Orb_x = Orb_PreComp(a_max,j_d,l_d,Orb)

                    VpnSum = 0.0

                    if J == 0
                        WppnSum, WnpnSum = 0.0, 0.0
                    end

                    @inbounds for i in Orb_x
                        Ket_ci = V2b_temp_index(c,i,J,P,Orb_NN_t)
                        VpnME = V2b_temp_transformation_ind4_T0_MEs(P,J,Bra_ab,Ket_ci,i,d,C,V_NN_t1)

                        VpnSum += VpnME

                        if J == 0
                            WppnME, WnpnME = W2b_temp_transformation_ind4_T0_MEs(P,Bra_ab,Ket_ci,i,d,C,W_NN_t1)

                            WppnSum += WppnME
                            WnpnSum += WnpnME
                        end
                    end

                    if J == 0
                        @views W_NN.ppn[P,1][Ind] = WppnSum
                        @views W_NN.npn[P,1][Ind] = WnpnSum
                    end

                    @views V_NN.pn[P,J+1][Ind] = VpnSum
                end

                # Case of pp & nn interaction ... T = 1
                if Bra <= N_T1
                    Ind = Bra + (Ket - 1) * N_T1 - div(Ket * (Ket - 1),2)
                    a, b = Orb_NN.Ind[2,P,J+1][Bra][1], Orb_NN.Ind[2,P,J+1][Bra][2]
                    c, d = Orb_NN.Ind[2,P,J+1][Ket][1], Orb_NN.Ind[2,P,J+1][Ket][2]
                    l_d, j_d = Orb[d].l, Orb[d].j

                    Bra_ab = V2b_temp_index(a,b,J,P,Orb_NN_t)
                    Orb_x = Orb_PreComp(a_max,j_d,l_d,Orb)
                    VppSum, VnnSum = 0.0, 0.0

                    if J == 0
                        WppSum, WnnSum = 0.0, 0.0
                    end

                    @inbounds for i in Orb_x
                        Ket_ci = V2b_temp_index(c,i,J,P,Orb_NN_t)
                        VppME, VnnME = V2b_temp_transformation_ind4_T1_MEs(P,J,Bra_ab,Ket_ci,i,d,C,V_NN_t1)

                        VppSum += VppME
                        VnnSum += VnnME

                        if J == 0
                            WppME, WnnME = W2b_temp_transformation_ind4_T1_MEs(P,Bra_ab,Ket_ci,i,d,C,W_NN_t1)

                            WppSum += WppME
                            WnnSum += WnnME
                        end

                    end

                    @views V_NN.pp[P,J+1][Ind] = VppSum
                    @views V_NN.nn[P,J+1][Ind] = VnnSum

                    if J == 0
                        @views W_NN.pp[P,1][Ind] = WppSum
                        @views W_NN.nn[P,1][Ind] = WnnSum
                    end

                end

            end
        end

    end

    return V_NN, W_NN
end

struct V2B_H40_Temp
    pp::Matrix{Matrix{Float64}}
    pn::Matrix{Matrix{Float64}}
    nn::Matrix{Matrix{Float64}}
end

struct W2B_H40_Temp
    pp::Vector{Matrix{Matrix{Float64}}}
    ppn::Vector{Matrix{Matrix{Float64}}}
    npn::Vector{Matrix{Matrix{Float64}}}
    nn::Vector{Matrix{Matrix{Float64}}}
end

function qpH_allocate_H40(Params::Parameters,Orb::Vector{NOrb},Orb_NN::NNOrb,H_NN::qpH2B,V_NN::NNInt,W_NN::W2B,U::pnMatrix,V::pnMatrix)
    println("\nAllocating H^(40) components of the quasiparticle Hamiltonian ...")

    # Preallocate arrays of allowed values of J & P ...
    JP = JP_initialize(Params.Calc.N2max + 1)

    # Initialize temporary V & W interaction arrays ...
    V_NN_t1, V_NN_t2, W_NN_t1, W_NN_t2, Orb_NN_t = qpH_initialize_H40(Params,Orb)


    # Transformation of the 1st index ...
    @time V_NN_t1, W_NN_t1 = qpH_allocate_H40_ind1(Params,JP,Orb,Orb_NN,Orb_NN_t,V_NN,W_NN,V_NN_t1,W_NN_t1,U,V)

    # Transformation of the 2nd index ...
    @time V_NN_t2, W_NN_t2 = qpH_allocate_H40_ind2(Params,JP,Orb,Orb_NN_t,V_NN_t1,W_NN_t1,V_NN_t2,W_NN_t2,U,V)

    # Transformation of the 3rd index ...
    @time V_NN_t1, W_NN_t1 = qpH_allocate_H40_ind3(Params,JP,Orb,Orb_NN_t,V_NN_t1,W_NN_t1,V_NN_t2,W_NN_t2,U,V)

    # Transformation of the 4th index ...
    @time H_NN = qpH_allocate_H40_ind4(Params,JP,Orb,Orb_NN,Orb_NN_t,H_NN,V_NN_t1,W_NN_t1,U,V)

    # Deallocate V_NN_t & W_NN_t ...
    V_NN_t1, W_NN_t1 = nothing, nothing
    V_NN_t2, W_NN_t2 = nothing, nothing

    # Perform the Garbage Collection ...
    GC.gc()

    println("\nAllocation of H^(40) components of the quasiparticle Hamiltonian done ...")

    return H_NN
end

function qpH_initialize_H40(Params::Parameters,Orb::Vector{NOrb})
    # Read parameters ...
    N_max = Params.Calc.Nmax
    N_2max = 2 * N_max
    a_max = div((N_max + 1)*(N_max + 2),2)
    J_max = Int64(2*N_max + 1)

    # Initiliaze array for counting of NN orbitals ...
    N_Orb_NN = zeros(Int64,2,J_max+1)

    # Initiliaze dictionary for NN orbitals ...
    Orb_NN_Dic = Dict{Tuple{Int8,Int8,Int16,Int16},Int32}()

    # Count the number of NN orbitals ...
    @inbounds for a in 1:a_max
        n_a = Orb[a].n
        l_a = Orb[a].l
        j_a = Orb[a].j
        if (2*n_a + l_a) <= N_max
            @inbounds for b = 1:a_max
                n_b = Orb[b].n
                l_b = Orb[b].l
                j_b = Orb[b].j
                if ((2*(n_a + n_b) + l_a + l_b) <= N_2max)
                    P = rem(l_a + l_b, 2) + 1
                    @inbounds for J = Int64(abs(j_a - j_b)/2):Int64((j_a + j_b)/2)
                        N_Orb_NN[P,J+1] += 1
                    end
                end
            end
        end
    end

    # Initialite the array for NN orbital indices ...
    Ind_Orb_NN = Matrix{Vector{Vector{Int64}}}(undef,2,J_max+1)
    @inbounds for P = 1:2
        @inbounds for J = 0:J_max
            Ind_Orb_NN[P,J+1] = Vector{Vector{Int64}}(undef,N_Orb_NN[P,J+1])
        end
    end

    # Reinitialize array for counting of NN orbitals ...
    N_Orb_NN = zeros(Int64,2,J_max+1)

    # Allocate the dictionary and inidices for NN orbitals ...
    @inbounds for a in 1:a_max
        n_a = Orb[a].n
        l_a = Orb[a].l
        j_a = Orb[a].j
        if (2*n_a + l_a) <= N_max
            @inbounds for b = 1:a_max
                n_b = Orb[b].n
                l_b = Orb[b].l
                j_b = Orb[b].j
                if ((2*(n_a + n_b) + l_a + l_b) <= N_2max)
                    P = rem(l_a + l_b, 2) + 1
                    @inbounds for J = Int64(abs(j_a - j_b)/2):Int64((j_a + j_b)/2)
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

    # Initialize the NN interaction matrices ...
        # Case of V ...
    V_pp_t1 = Matrix{Matrix{Float64}}(undef,2,J_max+1)
    V_pn_t1 = Matrix{Matrix{Float64}}(undef,2,J_max+1)
    V_nn_t1 = Matrix{Matrix{Float64}}(undef,2,J_max+1)
    V_pp_t2 = Matrix{Matrix{Float64}}(undef,2,J_max+1)
    V_pn_t2 = Matrix{Matrix{Float64}}(undef,2,J_max+1)
    V_nn_t2 = Matrix{Matrix{Float64}}(undef,2,J_max+1)
        # Case of W ...
    W_pp_t1 = Vector{Matrix{Matrix{Float64}}}(undef,2)
    W_ppn_t1 = Vector{Matrix{Matrix{Float64}}}(undef,2)
    W_npn_t1 = Vector{Matrix{Matrix{Float64}}}(undef,2)
    W_nn_t1 = Vector{Matrix{Matrix{Float64}}}(undef,2)
    W_pp_t2 = Vector{Matrix{Matrix{Float64}}}(undef,2)
    W_ppn_t2 = Vector{Matrix{Matrix{Float64}}}(undef,2)
    W_npn_t2 = Vector{Matrix{Matrix{Float64}}}(undef,2)
    W_nn_t2 = Vector{Matrix{Matrix{Float64}}}(undef,2)
    @inbounds for i in 1:2
        W_pp_t1[i] = Matrix{Matrix{Float64}}(undef,2,1)
        W_ppn_t1[i] = Matrix{Matrix{Float64}}(undef,2,1)
        W_npn_t1[i] = Matrix{Matrix{Float64}}(undef,2,1)
        W_nn_t1[i] = Matrix{Matrix{Float64}}(undef,2,1)
        W_pp_t2[i] = Matrix{Matrix{Float64}}(undef,2,1)
        W_ppn_t2[i] = Matrix{Matrix{Float64}}(undef,2,1)
        W_npn_t2[i] = Matrix{Matrix{Float64}}(undef,2,1)
        W_nn_t2[i] = Matrix{Matrix{Float64}}(undef,2,1)
    end
 
    # Allocate entries of V & W ...
    @inbounds for P = 1:2
        @inbounds Threads.@threads for J=0:J_max
            Length = N_Orb_NN[P,J+1]

            V_pp_t1[P,J+1] = zeros(Float64,Length,Length)
            V_pn_t1[P,J+1] = zeros(Float64,Length,Length)
            V_nn_t1[P,J+1] = zeros(Float64,Length,Length)
            V_pp_t2[P,J+1] = zeros(Float64,Length,Length)
            V_pn_t2[P,J+1] = zeros(Float64,Length,Length)
            V_nn_t2[P,J+1] = zeros(Float64,Length,Length)

            if J == 0
                @inbounds for i in 1:2
                    W_pp_t1[i][P,1] = zeros(Float64,Length,Length)
                    W_ppn_t1[i][P,1] = zeros(Float64,Length,Length)
                    W_npn_t1[i][P,1] = zeros(Float64,Length,Length)
                    W_nn_t1[i][P,1] = zeros(Float64,Length,Length)
                    W_pp_t2[i][P,1] = zeros(Float64,Length,Length)
                    W_ppn_t2[i][P,1] = zeros(Float64,Length,Length)
                    W_npn_t2[i][P,1] = zeros(Float64,Length,Length)
                    W_nn_t2[i][P,1] = zeros(Float64,Length,Length)
                end
            end

        end

    end

    # Allocate the NN interaction in quasiparticle picture ...
    V_NN_t1 = V2B_H40_Temp(V_pp_t1,V_pn_t1,V_nn_t1)
    V_NN_t2 = V2B_H40_Temp(V_pp_t2,V_pn_t2,V_nn_t2)
    W_NN_t1 = W2B_H40_Temp(W_pp_t1,W_ppn_t1,W_npn_t1,W_nn_t1)
    W_NN_t2 = W2B_H40_Temp(W_pp_t2,W_ppn_t2,W_npn_t2,W_nn_t2)

    # Allocate the array for NN orbitals ...
    Orb_NN_t = NNOrb_Temp(Orb_NN_Dic,N_Orb_NN,Ind_Orb_NN)

    return V_NN_t1, V_NN_t2, W_NN_t1, W_NN_t2, Orb_NN_t
end

function qpH_allocate_H40_ind1(Params::Parameters,JP::Vector{Vector{Int64}},Orb::Vector{NOrb},Orb_NN::NNOrb,Orb_NN_t::NNOrb_Temp,V_NN::NNInt,W_NN::W2B,V_NN_t1::V2B_H40_Temp,W_NN_t1::W2B_H40_Temp,U::pnMatrix,V::pnMatrix)
    # Parameter initialization...
    N_max = Params.Calc.Nmax
    N_2max = 2*N_max
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Define function for the transformation MEs for the 1st index ...
        # Case of V ...
    @inline function V2b_temp_transformation_ind1_MEs(J::Int64,i::Int64,a::Int64,b::Int64,c::Int64,d::Int64,T::pnMatrix,V_NN::NNInt,Orb::Vector{NOrb},Orb_NN::NNOrb)
        return @views T.p[i,a] * V2B(i,b,c,d,J,1,V_NN.pp,Orb,Orb_NN), T.p[i,a] * V2B(i,b,c,d,J,0,V_NN.pn,Orb,Orb_NN), T.n[i,a] * V2B(i,b,c,d,J,1,V_NN.nn,Orb,Orb_NN)
    end
        # Case of W ...
    @inline function W2b_temp_transformation_ind1_MEs(t::Int64,i::Int64,a::Int64,b::Int64,c::Int64,d::Int64,T::Matrix{Float64},W_NN::Matrix{Vector{Float64}},Orb::Vector{NOrb},Orb_NN::NNOrb)
        return @views T[i,a] * V2B(i,b,c,d,0,t,W_NN,Orb,Orb_NN)
    end
    # Transformation of the 1st index ...
    println("\nPerforming Quasiparticle Transformation of the residual 2-body NN interaction H40 in the 1st index...")
    @inbounds Threads.@threads for i in JP
        J, P = i[1], i[2]
        if P == 1
            println("Calculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = +")
        else
            println("Calculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = -")
        end
        N = Orb_NN_t.N[P,J+1]
        @inbounds for Bra in 1:N
            a = Orb_NN_t.Ind[P,J+1][Bra][1]
            b = Orb_NN_t.Ind[P,J+1][Bra][2]
            l_a, j_a = Orb[a].l, Orb[a].j

            Orb_x = Orb_PreComp(a_max,j_a,l_a,Orb)
            @inbounds for Ket in 1:N
                c = Orb_NN_t.Ind[P,J+1][Ket][1]
                d = Orb_NN_t.Ind[P,J+1][Ket][2]

                VppSum, VpnSum, VnnSum = 0.0, 0.0, 0.0
                if J == 0
                    WppSum, WppnSum1, WppnSum2, WnpnSum1, WnpnSum2, WnnSum = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
                end
                
                @inbounds for i in Orb_x
                    #Bra_ib = V2b_temp_index(i,b,J,P,Orb_NN_t)
                    #VppME, VpnME, VnnME = V2b_temp_transformation_ind1_MEs(P,J,Bra_ib,Ket,i,a,U,V_NN)
                    
                    VppME, VpnME, VnnME = V2b_temp_transformation_ind1_MEs(J,i,a,b,c,d,U,V_NN,Orb,Orb_NN)
                    VppSum += VppME
                    VpnSum += VpnME
                    VnnSum += VnnME

                    if J == 0
                        #=
                        WppME = W2b_temp_transformation_ind1_MEs(P,Bra_ib,Ket,i,a,U.p,W_NN.pp)
                        WppnME1 = W2b_temp_transformation_ind1_MEs(P,Bra_ib,Ket,i,a,U.p,W_NN.ppn)
                        WppnME2 = W2b_temp_transformation_ind1_MEs(P,Bra_ib,Ket,i,a,V.p,W_NN.ppn)
                        WnpnME1 = W2b_temp_transformation_ind1_MEs(P,Bra_ib,Ket,i,a,U.p,W_NN.npn)
                        WnpnME2 = W2b_temp_transformation_ind1_MEs(P,Bra_ib,Ket,i,a,U.p,W_NN.npn)
                        WnnME = W2b_temp_transformation_ind1_MEs(P,Bra_ib,Ket,i,a,U.n,W_NN.nn)
                        =#

                        WppME = W2b_temp_transformation_ind1_MEs(1,i,a,b,c,d,U.p,W_NN.pp,Orb,Orb_NN)
                        WppnME1 = W2b_temp_transformation_ind1_MEs(0,i,a,b,c,d,U.p,W_NN.ppn,Orb,Orb_NN)
                        WppnME2 = W2b_temp_transformation_ind1_MEs(0,i,a,b,c,d,V.p,W_NN.ppn,Orb,Orb_NN)
                        WnpnME1 = W2b_temp_transformation_ind1_MEs(0,i,a,b,c,d,U.p,W_NN.npn,Orb,Orb_NN)
                        WnpnME2 = W2b_temp_transformation_ind1_MEs(0,i,a,b,c,d,U.p,W_NN.npn,Orb,Orb_NN)
                        WnnME = W2b_temp_transformation_ind1_MEs(1,i,a,b,c,d,U.n,W_NN.nn,Orb,Orb_NN)


                        WppSum += WppME
                        WppnSum1 += WppnME1
                        WppnSum2 += WppnME2
                        WnpnSum1 += WnpnME1
                        WnpnSum2 += WnpnME2
                        WnnSum += WnnME
                    end

                end

                @views V_NN_t1.pp[P,J+1][Bra,Ket] = VppSum
                @views V_NN_t1.pn[P,J+1][Bra,Ket] = VpnSum
                @views V_NN_t1.nn[P,J+1][Bra,Ket] = VnnSum

                if J == 0
                    @views W_NN_t1.pp[1][P,1][Bra,Ket] = WppSum
                    @views W_NN_t1.ppn[1][P,1][Bra,Ket] = WppnSum1
                    @views W_NN_t1.npn[1][P,1][Bra,Ket] = WnpnSum1
                    @views W_NN_t1.nn[1][P,1][Bra,Ket] = WnnSum

                    @views W_NN_t1.pp[2][P,1][Bra,Ket] = WppSum
                    @views W_NN_t1.ppn[2][P,1][Bra,Ket] = WppnSum2
                    @views W_NN_t1.npn[2][P,1][Bra,Ket] = WnpnSum2
                    @views W_NN_t1.nn[2][P,1][Bra,Ket] = WnnSum
                end

            end
        end
    end

    return V_NN_t1, W_NN_t1
end

function qpH_allocate_H40_ind2(Params::Parameters,JP::Vector{Vector{Int64}},Orb::Vector{NOrb},Orb_NN_t::NNOrb_Temp,V_NN_t1::V2B_H40_Temp,W_NN_t1::W2B_H40_Temp,V_NN_t2::V2B_H40_Temp,W_NN_t2::W2B_H40_Temp,U::pnMatrix,V::pnMatrix)
    # Parameter initialization...
    N_max = Params.Calc.Nmax
    N_2max = 2*N_max
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Define function for the transformation MEs for the 2nd index ...
        # Case of V ...
    @inline function V2b_temp_transformation_ind2_MEs(P::Int64,J::Int64,Bra::Int64,Ket::Int64,i::Int64,b::Int64,T::pnMatrix,V_NN_t::V2B_H40_Temp)
        return @views T.p[i,b] * V_NN_t.pp[P,J+1][Bra,Ket], T.n[i,b] * V_NN_t.pn[P,J+1][Bra,Ket], T.n[i,b] * V_NN_t.nn[P,J+1][Bra,Ket]
    end
        # Case of W ...
    @inline function W2b_temp_transformation_ind2_MEs(Ind::Int64,P::Int64,Bra::Int64,Ket::Int64,i::Int64,b::Int64,T::Matrix{Float64},W_NN_t::Vector{Matrix{Matrix{Float64}}})
        return @views T[i,b] * W_NN_t[Ind][P,1][Bra,Ket]
    end

    # Transformation of the 2nd index ...
    println("\nPerforming Quasiparticle Transformation of the residual 2-body NN interaction H40 in the 2nd index...")

    @inbounds Threads.@threads for i in JP
        J, P = i[1], i[2]
        if P == 1
            println("Calculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = +")
        else
            println("Calculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = -")
        end
        N = Orb_NN_t.N[P,J+1]
        @inbounds for Bra in 1:N
            a = Orb_NN_t.Ind[P,J+1][Bra][1]
            b = Orb_NN_t.Ind[P,J+1][Bra][2]
            l_b, j_b = Orb[b].l, Orb[b].j

            Orb_x = Orb_PreComp(a_max,j_b,l_b,Orb)
            @inbounds for Ket in 1:N
                c = Orb_NN_t.Ind[P,J+1][Ket][1]
                d = Orb_NN_t.Ind[P,J+1][Ket][2]

                VppSum, VpnSum, VnnSum = 0.0, 0.0, 0.0

                if J == 0
                    WppSum1, WppnSum1, WnpnSum1, WnnSum1 = 0.0, 0.0, 0.0, 0.0
                    WppSum2, WppnSum2, WnpnSum2, WnnSum2 = 0.0, 0.0, 0.0, 0.0
                end
                
                @inbounds for i in Orb_x
                    Bra_ai = V2b_temp_index(a,i,J,P,Orb_NN_t)
                    VppME, VpnME, VnnME = V2b_temp_transformation_ind2_MEs(P,J,Bra_ai,Ket,i,b,U,V_NN_t1)
                    
                    VppSum += VppME
                    VpnSum += VpnME
                    VnnSum += VnnME

                    if J == 0
                        WppME1 = W2b_temp_transformation_ind2_MEs(1,P,Bra_ai,Ket,i,b,U.p,W_NN_t1.pp)
                        WppnME1 = W2b_temp_transformation_ind2_MEs(1,P,Bra_ai,Ket,i,b,U.n,W_NN_t1.ppn)
                        WnpnME1 = W2b_temp_transformation_ind2_MEs(1,P,Bra_ai,Ket,i,b,U.n,W_NN_t1.npn)
                        WnnME1 = W2b_temp_transformation_ind2_MEs(1,P,Bra_ai,Ket,i,b,U.n,W_NN_t1.nn)

                        WppME2 = W2b_temp_transformation_ind2_MEs(2,P,Bra_ai,Ket,i,b,V.p,W_NN_t1.pp)
                        WppnME2 = W2b_temp_transformation_ind2_MEs(2,P,Bra_ai,Ket,i,b,U.n,W_NN_t1.ppn)
                        WnpnME2 = W2b_temp_transformation_ind2_MEs(2,P,Bra_ai,Ket,i,b,V.n,W_NN_t1.npn)
                        WnnME2 = W2b_temp_transformation_ind2_MEs(2,P,Bra_ai,Ket,i,b,V.n,W_NN_t1.nn)

                        WppSum1 += WppME1
                        WppnSum1 += WppnME1
                        WnpnSum1 += WnpnME1
                        WnnSum1 += WnnME1

                        WppSum2 += WppME2
                        WppnSum2 += WppnME2
                        WnpnSum2 += WnpnME2
                        WnnSum2 += WnnME2
                    end

                end

                @views V_NN_t2.pp[P,J+1][Bra,Ket] = VppSum
                @views V_NN_t2.pn[P,J+1][Bra,Ket] = VpnSum
                @views V_NN_t2.nn[P,J+1][Bra,Ket] = VnnSum

                if J == 0
                    @views W_NN_t2.pp[1][P,1][Bra,Ket] = WppSum1
                    @views W_NN_t2.ppn[1][P,1][Bra,Ket] = WppnSum1
                    @views W_NN_t2.npn[1][P,1][Bra,Ket] = WnpnSum1
                    @views W_NN_t2.nn[1][P,1][Bra,Ket] = WnnSum1

                    @views W_NN_t2.pp[2][P,1][Bra,Ket] = WppSum2
                    @views W_NN_t2.ppn[2][P,1][Bra,Ket] = WppnSum2
                    @views W_NN_t2.npn[2][P,1][Bra,Ket] = WnpnSum2
                    @views W_NN_t2.nn[2][P,1][Bra,Ket] = WnnSum2
                end

            end
        end
    end

    return V_NN_t2, W_NN_t2
end

function qpH_allocate_H40_ind3(Params::Parameters,JP::Vector{Vector{Int64}},Orb::Vector{NOrb},Orb_NN_t::NNOrb_Temp,V_NN_t1::V2B_H40_Temp,W_NN_t1::W2B_H40_Temp,V_NN_t2::V2B_H40_Temp,W_NN_t2::W2B_H40_Temp,U::pnMatrix,V::pnMatrix)
    # Parameter initialization...
    N_max = Params.Calc.Nmax
    N_2max = 2*N_max
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Define function for the transformation MEs for the 3rd index ...
        # Case of V ...
    @inline function V2b_temp_transformation_ind3_MEs(P::Int64,J::Int64,Bra::Int64,Ket::Int64,i::Int64,c::Int64,T::pnMatrix,V_NN_t::V2B_H40_Temp)
        return @views T.p[i,c] * V_NN_t.pp[P,J+1][Bra,Ket], T.p[i,c] * V_NN_t.pn[P,J+1][Bra,Ket], T.n[i,c] * V_NN_t.nn[P,J+1][Bra,Ket]
    end
        # Case of W ...
    @inline function W2b_temp_transformation_ind3_MEs(Ind::Int64,P::Int64,Bra::Int64,Ket::Int64,i::Int64,c::Int64,T::Matrix{Float64},W_NN_t::Vector{Matrix{Matrix{Float64}}})
        return @views T[i,c] * W_NN_t[Ind][P,1][Bra,Ket]
    end

    # Transformation of the 3rd index ...
    println("\nPerforming Quasiparticle Transformation of the residual 2-body NN interaction H40 in the 3rd index...")

    @inbounds Threads.@threads for i in JP
        J, P = i[1], i[2]
        if P == 1
            println("Calculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = +")
        else
            println("Calculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = -")
        end
        N = Orb_NN_t.N[P,J+1]
        @inbounds for Ket in 1:N
            c = Orb_NN_t.Ind[P,J+1][Ket][1]
            d = Orb_NN_t.Ind[P,J+1][Ket][2]
            l_c, j_c = Orb[c].l, Orb[c].j

            Orb_x = Orb_PreComp(a_max,j_c,l_c,Orb)

            @inbounds for Bra in 1:N
                VppSum, VpnSum, VnnSum = 0.0, 0.0, 0.0

                if J == 0
                    WppSum1, WppnSum1, WnpnSum1, WnnSum1 = 0.0, 0.0, 0.0, 0.0
                    WppSum2, WppnSum2, WnpnSum2, WnnSum2 = 0.0, 0.0, 0.0, 0.0
                end

                @inbounds for i in Orb_x
                    Ket_id = V2b_temp_index(i,d,J,P,Orb_NN_t)
                    VppME, VpnME, VnnME = V2b_temp_transformation_ind3_MEs(P,J,Bra,Ket_id,i,c,V,V_NN_t2)

                    VppSum += VppME
                    VpnSum += VpnME
                    VnnSum += VnnME

                    if J == 0
                        WppME1 = W2b_temp_transformation_ind3_MEs(1,P,Bra,Ket_id,i,c,U.p,W_NN_t2.pp)
                        WppnME1 = W2b_temp_transformation_ind3_MEs(1,P,Bra,Ket_id,i,c,U.p,W_NN_t2.ppn)
                        WnpnME1 = W2b_temp_transformation_ind3_MEs(1,P,Bra,Ket_id,i,c,V.p,W_NN_t2.npn)
                        WnnME1 = W2b_temp_transformation_ind3_MEs(1,P,Bra,Ket_id,i,c,U.n,W_NN_t2.nn)

                        WppME2 = W2b_temp_transformation_ind3_MEs(2,P,Bra,Ket_id,i,c,V.p,W_NN_t2.pp)
                        WppnME2 = W2b_temp_transformation_ind3_MEs(2,P,Bra,Ket_id,i,c,V.p,W_NN_t2.ppn)
                        WnpnME2 = W2b_temp_transformation_ind3_MEs(2,P,Bra,Ket_id,i,c,V.p,W_NN_t2.npn)
                        WnnME2 = W2b_temp_transformation_ind3_MEs(2,P,Bra,Ket_id,i,c,V.n,W_NN_t2.nn)
                        
                        WppSum1 += WppME1
                        WppnSum1 += WppnME1
                        WnpnSum1 += WnpnME1
                        WnnSum1 += WnnME1
                        
                        WppSum2 += WppME2
                        WppnSum2 += WppnME2
                        WnpnSum2 += WnpnME2
                        WnnSum2 += WnnME2
                    end

                end
                @views V_NN_t1.pp[P,J+1][Bra,Ket] = VppSum
                @views V_NN_t1.pn[P,J+1][Bra,Ket] = VpnSum
                @views V_NN_t1.nn[P,J+1][Bra,Ket] = VnnSum

                if J == 0
                    @views W_NN_t1.pp[1][P,1][Bra,Ket] = WppSum1
                    @views W_NN_t1.ppn[1][P,1][Bra,Ket] = WppnSum1
                    @views W_NN_t1.npn[1][P,1][Bra,Ket] = WnpnSum1
                    @views W_NN_t1.nn[1][P,1][Bra,Ket] = WnnSum1

                    @views W_NN_t1.pp[2][P,1][Bra,Ket] = WppSum2
                    @views W_NN_t1.ppn[2][P,1][Bra,Ket] = WppnSum2
                    @views W_NN_t1.npn[2][P,1][Bra,Ket] = WnpnSum2
                    @views W_NN_t1.nn[2][P,1][Bra,Ket] = WnnSum2
                end

            end
        end
    end

    return V_NN_t1, W_NN_t1
end

function qpH_allocate_H40_ind4(Params::Parameters,JP::Vector{Vector{Int64}},Orb::Vector{NOrb},Orb_NN::NNOrb,Orb_NN_t::NNOrb_Temp,H_NN::qpH2B,V_NN_t1::V2B_H40_Temp,W_NN_t1::W2B_H40_Temp,U::pnMatrix,V::pnMatrix)
    # Parameter initialization...
    N_max = Params.Calc.Nmax
    N_2max = 2*N_max
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Define function for the transformation MEs for the 4th index ...
        # Case of V ... T = 0
    @inline function V2b_temp_transformation_ind4_T0_MEs(P::Int64,J::Int64,Bra::Int64,Ket::Int64,i::Int64,d::Int64,T::pnMatrix,V_NN_t::V2B_H40_Temp)
        return @views T.n[i,d] * V_NN_t.pn[P,J+1][Bra,Ket]
    end
        # Case of V ... T = 0
    @inline function V2b_temp_transformation_ind4_T1_MEs(P::Int64,J::Int64,Bra::Int64,Ket::Int64,i::Int64,d::Int64,T::pnMatrix,V_NN_t::V2B_H40_Temp)
        return @views T.p[i,d] * V_NN_t.pp[P,J+1][Bra,Ket], T.n[i,d] * V_NN_t.nn[P,J+1][Bra,Ket]
    end
        # Case of W ...
    @inline function W2b_temp_transformation_ind4_MEs(Ind::Int64,P::Int64,Bra::Int64,Ket::Int64,i::Int64,d::Int64,T::Matrix{Float64},W_NN_t::Vector{Matrix{Matrix{Float64}}})
        return @views T[i,d] * W_NN_t[Ind][P,1][Bra,Ket]
    end

    # Index 4
    println("\nPerforming Quasiparticle Transformation of the residual 2-body NN interaction H40 in the 4th index...")

    @inbounds Threads.@threads for i in JP
        J, P = i[1], i[2]
        if P == 1
            println("Calculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = +")
        else
            println("Calculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = -")
        end
        N_T0, N_T1 = Orb_NN.N[1,P,J+1], Orb_NN.N[2,P,J+1]

        @inbounds for Bra in 1:max(N_T0,N_T1)
            @inbounds for Ket in 1:Bra

                # Case of pn interaction ... T = 0
                if Bra <= N_T0
                    Ind = Bra + (Ket - 1) * N_T0 - div(Ket * (Ket - 1),2)
                    a, b = Orb_NN.Ind[1,P,J+1][Bra][1], Orb_NN.Ind[1,P,J+1][Bra][2]
                    c, d = Orb_NN.Ind[1,P,J+1][Ket][1], Orb_NN.Ind[1,P,J+1][Ket][2]
                    l_d, j_d = Orb[d].l, Orb[d].j

                    Bra_ab = V2b_temp_index(a,b,J,P,Orb_NN_t)
                    Orb_x = Orb_PreComp(a_max,j_d,l_d,Orb)

                    HpnSum = 0.0

                    @inbounds for i in Orb_x
                        Ket_ci = V2b_temp_index(c,i,J,P,Orb_NN_t)
                        VpnME = V2b_temp_transformation_ind4_T0_MEs(P,J,Bra_ab,Ket_ci,i,d,V,V_NN_t1)

                        HpnSum += VpnME

                        if J == 0
                            WppnME1 = W2b_temp_transformation_ind4_MEs(1,P,Bra_ab,Ket_ci,i,d,V.n,W_NN_t1.ppn)
                            WnpnME1 = W2b_temp_transformation_ind4_MEs(1,P,Bra_ab,Ket_ci,i,d,U.n,W_NN_t1.npn)

                            WppnME2 = W2b_temp_transformation_ind4_MEs(2,P,Bra_ab,Ket_ci,i,d,V.n,W_NN_t1.ppn)
                            WnpnME2 = W2b_temp_transformation_ind4_MEs(2,P,Bra_ab,Ket_ci,i,d,V.n,W_NN_t1.npn)

                            HpnSum += 0.25 * (WppnME1 - WppnME2 + WnpnME1 - WnpnME2)
                        end

                    end

                    @views H_NN.H40.pn[P,J+1][Ind] = HpnSum
                end

                # Case of pp & nn interaction ... T = 1
                if Bra <= N_T1
                    Ind = Bra + (Ket - 1) * N_T1 - div(Ket * (Ket - 1),2)
                    a, b = Orb_NN.Ind[2,P,J+1][Bra][1], Orb_NN.Ind[2,P,J+1][Bra][2]
                    c, d = Orb_NN.Ind[2,P,J+1][Ket][1], Orb_NN.Ind[2,P,J+1][Ket][2]
                    l_d, j_d = Orb[d].l, Orb[d].j

                    Bra_ab = V2b_temp_index(a,b,J,P,Orb_NN_t)
                    Orb_x = Orb_PreComp(a_max,j_d,l_d,Orb)

                    HppSum, HnnSum = 0.0, 0.0

                    @inbounds for i in Orb_x
                        Ket_ci = V2b_temp_index(c,i,J,P,Orb_NN_t)
                        VppME, VnnME = V2b_temp_transformation_ind4_T1_MEs(P,J,Bra_ab,Ket_ci,i,d,V,V_NN_t1)

                        HppSum += -0.25 * VppME
                        HnnSum += -0.25 * VnnME

                        if J == 0
                            WppME1 = W2b_temp_transformation_ind4_MEs(1,P,Bra_ab,Ket_ci,i,d,V.p,W_NN_t1.pp)
                            WnnME1 = W2b_temp_transformation_ind4_MEs(1,P,Bra_ab,Ket_ci,i,d,V.n,W_NN_t1.nn)

                            WppME2 = W2b_temp_transformation_ind4_MEs(2,P,Bra_ab,Ket_ci,i,d,V.p,W_NN_t1.pp)
                            WnnME2 = W2b_temp_transformation_ind4_MEs(2,P,Bra_ab,Ket_ci,i,d,V.n,W_NN_t1.nn)

                            HppSum += WppME1 / 3.0
                            HnnSum += WnnME1 / 3.0

                            HppSum += WppME2 / 3.0
                            HnnSum += WnnME2 / 3.0
                        end

                    end

                    @views H_NN.H40.pp[P,J+1][Ind] = HppSum
                    @views H_NN.H40.nn[P,J+1][Ind] = HnnSum

                end

            end
        end

    end

    return H_NN
end
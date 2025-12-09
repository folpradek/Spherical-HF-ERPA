function H2b_res_no2b(Params::Parameters,Orb::Vector{NOrb},Orb_NN::NNOrb,Orb_NNN::NNNOrb,V_NN::NNInt,V_NNN::Array{Vector{Vector{Float32}},4},Rho::pnMatrix)
    # Read parameters ...
    N_max = Params.Calc.Nmax
    N_2max = Params.Calc.N2max
    N_3max = Params.Calc.N3max
    a_max = div((N_max + 1)*(N_max + 2),2)
    J_max = N_2max + 1
    c = Params.Calc.cV_res

    # Initialize list of J,P values ...
    JP = JP_Ini(J_max)

    println("\nInitializing the evaluation of the residual density dependent 2-body Hamiltonian ...\n")

    # Make (normal) density-dependent NN interaction V ...
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

    return V_NN
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

function V2b_temp_index(a::Int64,b::Int64,J::Int64,P::Int64,Orb_NN_t::NNOrb_Temp)
    P, J, a, b = Int8(P), Int8(J), Int16(a), Int16(b)
    Ind = Int64(Orb_NN_t.Dic[(P,J,a,b)])
    return Ind
end

function qpH2b_initialize(Params::Parameters,Orb::Vector{NOrb})
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

function qpH2b(Params::Parameters,Orb::Vector{NOrb},Orb_NN::NNOrb,V_NN::NNInt,C::pnMatrix,U::pnMatrix,V::pnMatrix)
    println("\nPreparing the Hamiltonian expressed in the J-scheme quasiparticle basis ...")

    # Perform the Canonical Transformation of V & W into the canonical basis ...
    V_NN = qpH2b_canonical_transformation(Params,Orb,Orb_NN,V_NN,C)

    # Initialize the 2-body quasiparticle Hamiltonian H^(kl) ...
    H_NN = qpH2b_initialize(Params,Orb)

    # Allocate H^(40) ...
    H_NN = qpH2b_allocate_H40(Params,Orb,Orb_NN,H_NN,V_NN,U,V)

    # Allocate H^(31) ...
    H_NN = qpH2b_allocate_H31(Params,Orb,Orb_NN,H_NN,V_NN,U,V)

    # Allocate H^(22) ...
    H_NN = qpH2b_allocate_H22(Params,Orb,Orb_NN,H_NN,V_NN,U,V)
    
    println("\nQuasiparticle 2-body Hamiltonian H^(kl) has been fully initialized and allocated ...\n")

    return H_NN
end

function qpH2b_export(Params::Parameters,Orb_NN::NNOrb,H_NN::qpH2B)
    # Read parameters ...
    N_max = Params.Calc.Nmax
    N_2max = 2*N_max
    J_max = N_2max + 1
    Output_File = Params.Calc.Path

    # Set up the export paths ...
    qpH2B_40_Export = "IO/" * Output_File * "/Bin/qpH2B40.bin"
    qpH2B_31_Export = "IO/" * Output_File * "/Bin/qpH2B31.bin"
    qpH2B_22_Export = "IO/" * Output_File * "/Bin/qpH2B22.bin"

    # Export residual quasiparticle 2-body Hamiltonian ...
    println("\nExporting quasiparticle 2-body Hamiltonian in binary format ...")

    # Export H^(40) part of qp residual 2-body Hamiltonian ...
    open(qpH2B_40_Export, "w") do Export_File
        # ppH40
        @inbounds for J in 0:J_max
            @inbounds for P in 1:2
                N_T1 = Orb_NN.N[2,P,J+1]
                @inbounds for Bra in 1:N_T1
                    @inbounds for Ket in 1:Bra
                        Ind = Bra + (Ket - 1) * N_T1 - div(Ket * (Ket - 1),2)
                        ME = @views H_NN.H40.pp[P,J+1][Ind]
                        write(Export_File, Float64(ME))
                    end
                end
            end
        end

        # pnH40
        @inbounds for J in 0:J_max
            @inbounds for P in 1:2
                N_T0 = Orb_NN.N[1,P,J+1]
                @inbounds for Bra in 1:N_T0
                    @inbounds for Ket in 1:N_T0
                        Ind = Bra + (Ket - 1) * N_T0 - div(Ket * (Ket - 1),2)
                        ME = @views H_NN.H40.pn[P,J+1][Ind]
                        write(Export_File, Float64(ME))
                    end
                end
            end
        end

        # nnH40
        @inbounds for J in 0:J_max
            @inbounds for P in 1:2
                N_T1 = Orb_NN.N[2,P,J+1]
                @inbounds for Bra in 1:N_T1
                    @inbounds for Ket in 1:Bra
                        Ind = Bra + (Ket - 1) * N_T1 - div(Ket * (Ket - 1),2)
                        ME = @views H_NN.H40.nn[P,J+1][Ind]
                        write(Export_File, Float64(ME))
                    end
                end
            end
        end

    end

    # Export H^(31) part of qp residual 2-body Hamiltonian ...
    open(qpH2B_31_Export, "w") do Export_File
        # ppH31
        @inbounds for J in 0:J_max
            @inbounds for P in 1:2
                N_T1 = Orb_NN.N[2,P,J+1]
                @inbounds for Bra in 1:N_T1
                    @inbounds for Ket in 1:Bra
                        Ind = Bra + (Ket - 1) * N_T1 - div(Ket * (Ket - 1),2)
                        ME = @views H_NN.H31.pp[P,J+1][Ind]
                        write(Export_File, Float64(ME))
                    end
                end
            end
        end

        # pnH2011
        @inbounds for J in 0:J_max
            @inbounds for P in 1:2
                N_T0 = Orb_NN.N[1,P,J+1]
                @inbounds for Bra in 1:N_T0
                    @inbounds for Ket in 1:N_T0
                        Ind = Bra + (Ket - 1) * N_T0 - div(Ket * (Ket - 1),2)
                        ME = @views H_NN.H31.pn2011[P,J+1][Ind]
                        write(Export_File, Float64(ME))
                    end
                end
            end
        end

        # pnH1120
        @inbounds for J in 0:J_max
            @inbounds for P in 1:2
                N_T0 = Orb_NN.N[1,P,J+1]
                @inbounds for Bra in 1:N_T0
                    @inbounds for Ket in 1:N_T0
                        Ind = Bra + (Ket - 1) * N_T0 - div(Ket * (Ket - 1),2)
                        ME = @views H_NN.H31.pn1120[P,J+1][Ind]
                        write(Export_File, Float64(ME))
                    end
                end
            end
        end

        # nnH31
        @inbounds for J in 0:J_max
            @inbounds for P in 1:2
                N_T1 = Orb_NN.N[2,P,J+1]
                @inbounds for Bra in 1:N_T1
                    @inbounds for Ket in 1:Bra
                        Ind = Bra + (Ket - 1) * N_T1 - div(Ket * (Ket - 1),2)
                        ME = @views H_NN.H31.nn[P,J+1][Ind]
                        write(Export_File, Float64(ME))
                    end
                end
            end
        end

    end

    # Export H^(22) part of qp residual 2-body Hamiltonian ...
    open(qpH2B_22_Export, "w") do Export_File
        # ppH22
        @inbounds for J in 0:J_max
            @inbounds for P in 1:2
                N_T1 = Orb_NN.N[2,P,J+1]
                @inbounds for Bra in 1:N_T1
                    @inbounds for Ket in 1:Bra
                        Ind = Bra + (Ket - 1) * N_T1 - div(Ket * (Ket - 1),2)
                        ME = @views H_NN.H22.pp[P,J+1][Ind]
                        write(Export_File, Float64(ME))
                    end
                end
            end
        end

        # pnH2002
        @inbounds for J in 0:J_max
            @inbounds for P in 1:2
                N_T0 = Orb_NN.N[1,P,J+1]
                @inbounds for Bra in 1:N_T0
                    @inbounds for Ket in 1:N_T0
                        Ind = Bra + (Ket - 1) * N_T0 - div(Ket * (Ket - 1),2)
                        ME = @views H_NN.H22.pn2002[P,J+1][Ind]
                        write(Export_File, Float64(ME))
                    end
                end
            end
        end

        # pnH1111
        @inbounds for J in 0:J_max
            @inbounds for P in 1:2
                N_T0 = Orb_NN.N[1,P,J+1]
                @inbounds for Bra in 1:N_T0
                    @inbounds for Ket in 1:N_T0
                        Ind = Bra + (Ket - 1) * N_T0 - div(Ket * (Ket - 1),2)
                        ME = @views H_NN.H22.pn1111[P,J+1][Ind]
                        write(Export_File, Float64(ME))
                    end
                end
            end
        end

        # pnH0220
        @inbounds for J in 0:J_max
            @inbounds for P in 1:2
                N_T0 = Orb_NN.N[1,P,J+1]
                @inbounds for Bra in 1:N_T0
                    @inbounds for Ket in 1:N_T0
                        Ind = Bra + (Ket - 1) * N_T0 - div(Ket * (Ket - 1),2)
                        ME = @views H_NN.H22.pn0220[P,J+1][Ind]
                        write(Export_File, Float64(ME))
                    end
                end
            end
        end

        # nnH22
        @inbounds for J in 0:J_max
            @inbounds for P in 1:2
                N_T1 = Orb_NN.N[2,P,J+1]
                @inbounds for Bra in 1:N_T1
                    @inbounds for Ket in 1:Bra
                        Ind = Bra + (Ket - 1) * N_T1 - div(Ket * (Ket - 1),2)
                        ME = @views H_NN.H22.nn[P,J+1][Ind]
                        write(Export_File, Float64(ME))
                    end
                end
            end
        end

    end

    println("\nQuasiparticle 2-body Hamiltonian succesfully exported in binary format ...")

    return
end

# WORK IN PROGRESS ...
function qpH_import(Params::Parameters,Orb::Vector{NOrb})
    # Read parameters ...
    Import_File = Params.Calc.Path
    N_max = Params.Calc.Nmax
    N_2max = 2*N_max
    J_max = N_2max + 1
    V2B_res_Import = "IO/" * Import_File * "/Bin/V2B_res_HF.bin"

    # Read Residual 2-body NN interaction ...
    println("\nReading Residual 2-body NN interaction ...")

    # Initialize NN residual interaction arrays ...
    VNN_res, Orb_NN = V2B_Ini(Orb,N_max)

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
        open(V2B_res_Import, "r") do Bin_Read
            seek(Bin_Read, N_Chunk_skip[1,P,J+1])
            @inbounds for Bra in 1:N_t1
                @inbounds for Ket in 1:Bra
                    Ind = Bra + (Ket - 1) * N_t1 - div(Ket * (Ket - 1),2)
                    ME = read(Bin_Read, Float64)
                    VNN_res.pp[P,J+1][Ind] = ME
                end
            end
        end
    end

    # pnV
    @inbounds Threads.@threads for JP in JP_list
        J = JP[1]
        P = JP[2]
        N_t0 = Orb_NN.N[1,P,J+1]
        open(V2B_res_Import, "r") do Bin_Read
            seek(Bin_Read, N_Chunk_skip[2,P,J+1])
            @inbounds for Bra in 1:N_t0
                @inbounds for Ket in 1:N_t0
                    Ind = Bra + (Ket - 1) * N_t0 - div(Ket * (Ket - 1),2)
                    ME = read(Bin_Read, Float64)
                    VNN_res.pn[P,J+1][Ind] = ME
                end
            end
        end
    end

    # nnV
    @inbounds Threads.@threads for JP in JP_list
        J = JP[1]
        P = JP[2]
        N_t1 = Orb_NN.N[2,P,J+1]
        open(V2B_res_Import, "r") do Bin_Read
            seek(Bin_Read, N_Chunk_skip[3,P,J+1])
            @inbounds for Bra in 1:N_t1
                @inbounds for Ket in 1:Bra
                    Ind = Bra + (Ket - 1) * N_t1 - div(Ket * (Ket - 1),2)
                    ME = read(Bin_Read, Float64)
                    VNN_res.nn[P,J+1][Ind] = ME
                end
            end
        end
    end

    println("\nResidual 2-body NN interaction succesfully loaded ...")

    return VNN_res, Orb_NN
end

function qpH2b_count(N_max::Int64,Orb_NN::NNOrb)
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
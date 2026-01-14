function qpO2b_40_allocate(Params::Parameters,Orb::Vector{Orb1B},Orb_NN::Orb2B,Q_NN::qpO2B,O_NN::O2B,U::O1B,V::O1B)
    println("\nAllocating (40) components of the given operator in the quasiparticle representation ...")

    # Preallocate arrays of allowed values of J & P ...
    JP = JP_initialize(Params.Calc.N2max + 1)

    # Initialize temporary operator arrays ...
    Q_NN_t1, Q_NN_t2, Orb_NN_t = qpO2b_40_initialize(Params,Orb)

    # Transformation of the 1st index ...
    @time Q_NN_t1 = qpO2b_40_allocate_ind1(Params,JP,Orb,Orb_NN,Orb_NN_t,O_NN,Q_NN_t1,U,V)

    # Transformation of the 2nd index ...
    @time Q_NN_t2 = qpO2b_40_allocate_ind2(Params,JP,Orb,Orb_NN_t,Q_NN_t1,Q_NN_t2,U,V)

    # Transformation of the 3rd index ...
    @time Q_NN_t1 = qpO2b_40_allocate_ind3(Params,JP,Orb,Orb_NN_t,Q_NN_t1,Q_NN_t2,U,V)

    # Transformation of the 4th index ...
    @time Q_NN = qpO2b_40_allocate_ind4(Params,JP,Orb,Orb_NN,Orb_NN_t,Q_NN,Q_NN_t1,U,V)

    # Deallocate temporary arrays ...
    Q_NN_t1 = nothing
    Q_NN_t2 = nothing
    
    # Perform the Garbage Collection ...
    GC.gc()

    println("\nAllocation of (40) components of the given 2-body NN quasiparticle operator done ...")

    return Q_NN
end

function qpO2b_40_initialize(Params::Parameters,Orb::Vector{Orb1B})
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

    # Initialize the NN operator matrices ...
    O_pp_t1 = Matrix{Matrix{Float64}}(undef,2,J_max+1)
    O_pn_t1 = Matrix{Matrix{Float64}}(undef,2,J_max+1)
    O_nn_t1 = Matrix{Matrix{Float64}}(undef,2,J_max+1)
    O_pp_t2 = Matrix{Matrix{Float64}}(undef,2,J_max+1)
    O_pn_t2 = Matrix{Matrix{Float64}}(undef,2,J_max+1)
    O_nn_t2 = Matrix{Matrix{Float64}}(undef,2,J_max+1)
 
    # Allocate entries of O ...
    @inbounds for P = 1:2
        @inbounds Threads.@threads for J=0:J_max
            Length = N_Orb_NN[P,J+1]

            O_pp_t1[P,J+1] = zeros(Float64,Length,Length)
            O_pn_t1[P,J+1] = zeros(Float64,Length,Length)
            O_nn_t1[P,J+1] = zeros(Float64,Length,Length)
            O_pp_t2[P,J+1] = zeros(Float64,Length,Length)
            O_pn_t2[P,J+1] = zeros(Float64,Length,Length)
            O_nn_t2[P,J+1] = zeros(Float64,Length,Length)
        end
    end

    # Allocate the NN operator in the quasiparticle picture ...
    O_NN_t1 = O2B_40_Temp(O_pp_t1,O_pn_t1,O_nn_t1)
    O_NN_t2 = O2B_40_Temp(O_pp_t2,O_pn_t2,O_nn_t2)

    # Allocate the array for NN orbitals ...
    Orb_NN_t = Orb2B_Temp(Orb_NN_Dic,N_Orb_NN,Ind_Orb_NN)

    return O_NN_t1, O_NN_t2, Orb_NN_t
end

function qpO2b_40_allocate_ind1(Params::Parameters,JP::Vector{Vector{Int64}},Orb::Vector{Orb1B},Orb_NN::Orb2B,Orb_NN_t::Orb2B_Temp,O_NN::O2B,O_NN_t1::O2B_40_Temp,U::O1B,V::O1B)
    # Parameter initialization...
    N_max = Params.Calc.Nmax
    N_2max = 2*N_max
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Define function for the transformation MEs for the 1st index ...
    @inline function O2b_temp_transformation_ind1_MEs(J::Int64,P::Int64,i::Int64,a::Int64,b::Int64,c::Int64,d::Int64,T::O1B,O_NN::O2B,Orb::Vector{Orb1B},Orb_NN::Orb2B)
        return @views T.p[i,a] * O2b_pp(i,b,c,d,J,P,O_NN,Orb,Orb_NN), T.p[i,a] * O2b_pn(i,b,c,d,J,P,O_NN,Orb_NN), T.n[i,a] * O2b_nn(i,b,c,d,J,P,O_NN,Orb,Orb_NN)
    end

    # Transformation of the 1st index ...
    println("\nPerforming Quasiparticle Transformation of the residual 2-body NN interaction H40 in the 1st index...")
    @inbounds Threads.@threads for i in JP
        J, P = i[1], i[2]
        if P == 1
            println("\tCalculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = +")
        else
            println("\tCalculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = -")
        end
        N = Orb_NN_t.N[P,J+1]
        @inbounds for Bra in 1:N
            a = Orb_NN_t.Ind[P,J+1][Bra][1]
            b = Orb_NN_t.Ind[P,J+1][Bra][2]
            l_a, j_a = Orb[a].l, Orb[a].j
            l_b = Orb[b].l
            P_ab = rem(l_a + l_b,2) + 1

            Orb_x = Orb_PreComp(a_max,j_a,l_a,Orb)
            @inbounds for Ket in 1:N
                c = Orb_NN_t.Ind[P,J+1][Ket][1]
                d = Orb_NN_t.Ind[P,J+1][Ket][2]

                OppSum, OpnSum, OnnSum = 0.0, 0.0, 0.0

                @inbounds for i in Orb_x
                    OppME, OpnME, OnnME = O2b_temp_transformation_ind1_MEs(J,P_ab,i,a,b,c,d,U,O_NN,Orb,Orb_NN)

                    OppSum += OppME
                    OpnSum += OpnME
                    OnnSum += OnnME
                end

                @views O_NN_t1.pp[P,J+1][Bra,Ket] = OppSum
                @views O_NN_t1.pn[P,J+1][Bra,Ket] = OpnSum
                @views O_NN_t1.nn[P,J+1][Bra,Ket] = OnnSum
            end
        end
    end

    return O_NN_t1
end

function qpO2b_40_allocate_ind2(Params::Parameters,JP::Vector{Vector{Int64}},Orb::Vector{Orb1B},Orb_NN_t::Orb2B_Temp,O_NN_t1::O2B_40_Temp,O_NN_t2::O2B_40_Temp,U::O1B,V::O1B)
    # Parameter initialization...
    N_max = Params.Calc.Nmax
    N_2max = 2*N_max
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Define function for the transformation MEs for the 2nd index ...
    @inline function O2b_temp_transformation_ind2_MEs(P::Int64,J::Int64,Bra::Int64,Ket::Int64,i::Int64,b::Int64,T::O1B,O_NN_t::O2B_40_Temp)
        return @views T.p[i,b] * O_NN_t.pp[P,J+1][Bra,Ket], T.n[i,b] * O_NN_t.pn[P,J+1][Bra,Ket], T.n[i,b] * O_NN_t.nn[P,J+1][Bra,Ket]
    end

    # Transformation of the 2nd index ...
    println("\nPerforming Quasiparticle Transformation of the residual 2-body NN interaction H40 in the 2nd index...")

    @inbounds Threads.@threads for i in JP
        J, P = i[1], i[2]
        if P == 1
            println("\tCalculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = +")
        else
            println("\tCalculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = -")
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

                OppSum, OpnSum, OnnSum = 0.0, 0.0, 0.0
                
                @inbounds for i in Orb_x
                    Bra_ai = O2b_temp_index(a,i,J,P,Orb_NN_t)
                    OppME, OpnME, OnnME = O2b_temp_transformation_ind2_MEs(P,J,Bra_ai,Ket,i,b,U,O_NN_t1)
                    
                    OppSum += OppME
                    OpnSum += OpnME
                    OnnSum += OnnME
                end

                @views O_NN_t2.pp[P,J+1][Bra,Ket] = OppSum
                @views O_NN_t2.pn[P,J+1][Bra,Ket] = OpnSum
                @views O_NN_t2.nn[P,J+1][Bra,Ket] = OnnSum

            end
        end
    end

    return O_NN_t2
end

function qpO2b_40_allocate_ind3(Params::Parameters,JP::Vector{Vector{Int64}},Orb::Vector{Orb1B},Orb_NN_t::Orb2B_Temp,O_NN_t1::O2B_40_Temp,O_NN_t2::O2B_40_Temp,U::O1B,V::O1B)
    # Parameter initialization...
    N_max = Params.Calc.Nmax
    N_2max = 2*N_max
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Define function for the transformation MEs for the 3rd index ...
    @inline function O2b_temp_transformation_ind3_MEs(P::Int64,J::Int64,Bra::Int64,Ket::Int64,i::Int64,c::Int64,T::O1B,O_NN_t::O2B_40_Temp)
        return @views T.p[i,c] * O_NN_t.pp[P,J+1][Bra,Ket], T.p[i,c] * O_NN_t.pn[P,J+1][Bra,Ket], T.n[i,c] * O_NN_t.nn[P,J+1][Bra,Ket]
    end

    # Transformation of the 3rd index ...
    println("\nPerforming Quasiparticle Transformation of the residual 2-body NN interaction H40 in the 3rd index...")

    @inbounds Threads.@threads for i in JP
        J, P = i[1], i[2]
        if P == 1
            println("\tCalculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = +")
        else
            println("\tCalculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = -")
        end
        N = Orb_NN_t.N[P,J+1]
        @inbounds for Ket in 1:N
            c = Orb_NN_t.Ind[P,J+1][Ket][1]
            d = Orb_NN_t.Ind[P,J+1][Ket][2]
            l_c, j_c = Orb[c].l, Orb[c].j

            Orb_x = Orb_PreComp(a_max,j_c,l_c,Orb)

            @inbounds for Bra in 1:N
                OppSum, OpnSum, OnnSum = 0.0, 0.0, 0.0

                @inbounds for i in Orb_x
                    Ket_id = O2b_temp_index(i,d,J,P,Orb_NN_t)
                    OppME, OpnME, OnnME = O2b_temp_transformation_ind3_MEs(P,J,Bra,Ket_id,i,c,V,O_NN_t2)
                    OppSum += OppME
                    OpnSum += OpnME
                    OnnSum += OnnME

                end

                @views O_NN_t1.pp[P,J+1][Bra,Ket] = OppSum
                @views O_NN_t1.pn[P,J+1][Bra,Ket] = OpnSum
                @views O_NN_t1.nn[P,J+1][Bra,Ket] = OnnSum
            end
        end
    end

    return O_NN_t1
end

function qpO2b_40_allocate_ind4(Params::Parameters,JP::Vector{Vector{Int64}},Orb::Vector{Orb1B},Orb_NN::Orb2B,Orb_NN_t::Orb2B_Temp,O_NN::qpO2B,O_NN_t1::O2B_40_Temp,U::O1B,V::O1B)
    # Parameter initialization...
    N_max = Params.Calc.Nmax
    N_2max = 2*N_max
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Define function for the transformation MEs for the 4th index ...
        # Case of O ... T = 0
    @inline function O2b_temp_transformation_ind4_T0_MEs(P::Int64,J::Int64,Bra::Int64,Ket::Int64,i::Int64,d::Int64,T::O1B,O_NN_t::O2B_40_Temp)
        return @views T.n[i,d] * O_NN_t.pn[P,J+1][Bra,Ket]
    end
        # Case of O ... T = 0
    @inline function O2b_temp_transformation_ind4_T1_MEs(P::Int64,J::Int64,Bra::Int64,Ket::Int64,i::Int64,d::Int64,T::O1B,O_NN_t::O2B_40_Temp)
        return @views T.p[i,d] * O_NN_t.pp[P,J+1][Bra,Ket], T.n[i,d] * O_NN_t.nn[P,J+1][Bra,Ket]
    end

    # Index 4
    println("\nPerforming Quasiparticle Transformation of the residual 2-body NN interaction H40 in the 4th index...")

    @inbounds Threads.@threads for i in JP
        J, P = i[1], i[2]
        if P == 1
            println("\tCalculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = +")
        else
            println("\tCalculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = -")
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

                    Bra_ab = O2b_temp_index(a,b,J,P,Orb_NN_t)
                    Orb_x = Orb_PreComp(a_max,j_d,l_d,Orb)

                    OpnSum = 0.0

                    @inbounds for i in Orb_x
                        Ket_ci = O2b_temp_index(c,i,J,P,Orb_NN_t)
                        OpnME = O2b_temp_transformation_ind4_T0_MEs(P,J,Bra_ab,Ket_ci,i,d,V,O_NN_t1)

                        OpnSum += OpnME
                    end

                    Amp = Float64((-1)^J)
                    OpnSum = Amp * OpnSum

                    @views O_NN.qp40.pn[P,J+1][Ind] = OpnSum
                end

                # Case of pp & nn interaction ... T = 1
                if Bra <= N_T1
                    Ind = Bra + (Ket - 1) * N_T1 - div(Ket * (Ket - 1),2)
                    a, b = Orb_NN.Ind[2,P,J+1][Bra][1], Orb_NN.Ind[2,P,J+1][Bra][2]
                    c, d = Orb_NN.Ind[2,P,J+1][Ket][1], Orb_NN.Ind[2,P,J+1][Ket][2]
                    l_d, j_d = Orb[d].l, Orb[d].j

                    Bra_ab = O2b_temp_index(a,b,J,P,Orb_NN_t)
                    Orb_x = Orb_PreComp(a_max,j_d,l_d,Orb)

                    OppSum, OnnSum = 0.0, 0.0

                    @inbounds for i in Orb_x
                        Ket_ci = O2b_temp_index(c,i,J,P,Orb_NN_t)
                        OppME, OnnME = O2b_temp_transformation_ind4_T1_MEs(P,J,Bra_ab,Ket_ci,i,d,V,O_NN_t1)

                        OppSum += OppME
                        OnnSum += OnnME
                    end

                    Amp = -Float64((-1)^J) * sqrt(Float64((1 + kronecker_delta(a,b) * (-1)^J) * (1 + kronecker_delta(c,d) * (-1)^J)))

                    OppSum = Amp * OppSum
                    OnnSum = Amp * OnnSum

                    @views O_NN.qp40.pp[P,J+1][Ind] = OppSum
                    @views O_NN.qp40.nn[P,J+1][Ind] = OnnSum
                end

            end
        end

    end

    return O_NN
end

@inline function qpO2b_40_pp(a::Int64,b::Int64,c::Int64,d::Int64,J::Int64,P::Int64,O_NN::qpO2B,Orb::Vector{Orb1B},Orb_NN::Orb2B)
    if a < b
        a_s, b_s = b, a
        exp_ab = J + 1 + div(Orb[a].j + Orb[b].j, 2)
        Sign_ab = isodd(exp_ab) ? -1.0 : 1.0
    else
        a_s, b_s = a, b
        Sign_ab = 1
    end
    if c < d
        c_s, d_s = d, c
        exp_cd = J + 1 + div(Orb[c].j + Orb[d].j, 2)
        Sign_cd = isodd(exp_cd) ? -1.0 : 1.0
    else
        c_s, d_s = c, d
        Sign_cd = 1
    end
    Amp = Sign_ab * Sign_cd
    if (a < b) && (c < d)
        exp = div(Orb[a].j + Orb[b].j + Orb[c].j + Orb[d].j, 2)
        Amp = isodd(exp) ? -1.0 : 1.0
    end
    @inbounds Ind = O2b_index(a_s,b_s,c_s,d_s,J,P,1,Orb_NN)
    return Amp * O_NN.qp40.pp[P,J+1][Ind]
end

@inline function qpO2b_40_pn(a::Int64,b::Int64,c::Int64,d::Int64,J::Int64,P::Int64,O_NN::qpO2B,Orb_NN::Orb2B)
    @inbounds return O_NN.qp40.pn[P,J+1][O2b_index(a,b,c,d,J,P,0,Orb_NN)]
end

@inline function qpO2b_40_nn(a::Int64,b::Int64,c::Int64,d::Int64,J::Int64,P::Int64,O_NN::qpO2B,Orb::Vector{Orb1B},Orb_NN::Orb2B)
    if a < b
        a_s, b_s = b, a
        exp_ab = J + 1 + div(Orb[a].j + Orb[b].j, 2)
        Sign_ab = isodd(exp_ab) ? -1.0 : 1.0
    else
        a_s, b_s = a, b
        Sign_ab = 1
    end
    if c < d
        c_s, d_s = d, c
        exp_cd = J + 1 + div(Orb[c].j + Orb[d].j, 2)
        Sign_cd = isodd(exp_cd) ? -1.0 : 1.0
    else
        c_s, d_s = c, d
        Sign_cd = 1
    end
    Amp = Sign_ab * Sign_cd
    if (a < b) && (c < d)
        exp = div(Orb[a].j + Orb[b].j + Orb[c].j + Orb[d].j, 2)
        Amp = isodd(exp) ? -1.0 : 1.0
    end
    @inbounds Ind = O2b_index(a_s,b_s,c_s,d_s,J,P,1,Orb_NN)
    return Amp * O_NN.qp40.nn[P,J+1][Ind]
end
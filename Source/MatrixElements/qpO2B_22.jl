function qpO2b_22_allocate(Params::Parameters,Orb::Vector{Orb1B},Orb_NN::Orb2B,Q_NN::qpO2B,O_NN::O2B,U::O1B,V::O1B)
    # Preallocate arrays of allowed values of J & P ...
    JP = JP_initialize(Params.Calc.N2max + 1)

    println("\nAllocating (22) components of the given operator in the quasiparticle representation ...")

    # Initialize temporary V & W interaction arrays ...
    Q_NN_t1, Q_NN_t2, Orb_NN_t = qpO2b_22_initialize(Params,Orb)

    # Transformation of the 1st index ...
    @time Q_NN_t1 = qpO2b_22_allocate_ind1(Params,JP,Orb,Orb_NN,Orb_NN_t,O_NN,Q_NN_t1,U,V)

    # Transformation of the 2nd index ...
    @time Q_NN_t2 = qpO2b_22_allocate_ind2(Params,JP,Orb,Orb_NN_t,Q_NN_t1,Q_NN_t2,U,V)

    # Transformation of the 3rd index ...
    @time Q_NN_t1 = qpO2b_22_allocate_ind3(Params,JP,Orb,Orb_NN_t,Q_NN_t1,Q_NN_t2,U,V)

    # Transformation of the 4th index ...
    @time Q_NN = qpO2b_22_allocate_ind4(Params,JP,Orb,Orb_NN,Orb_NN_t,Q_NN,Q_NN_t1,Q_NN_t2,U,V)

    # Deallocate Q_NN_t ...
    Q_NN_t1 = nothing
    Q_NN_t2 = nothing

    # Perform the Garbage Collection ...
    GC.gc()

    println("\nAllocation of (22) components of the given 2-body NN quasiparticle operator done ...")

    return Q_NN
end

function qpO2b_22_initialize(Params::Parameters,Orb::Vector{Orb1B})
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

    # Initialize the 2-body NN operator matrices ...
        # Case of O_t1 ...
    O_pp_t1 = Vector{Matrix{Matrix{Float64}}}(undef,3)
    O_pn2002_t1 = Matrix{Matrix{Float64}}(undef,2,J_max+1)
    O_pn1111_t1 = Vector{Matrix{Matrix{Float64}}}(undef,4)
    O_pn0220_t1 = Matrix{Matrix{Float64}}(undef,2,J_max+1)
    O_nn_t1 = Vector{Matrix{Matrix{Float64}}}(undef,3)
        # Case of O_t2 ...
    O_pp_t2 = Vector{Matrix{Matrix{Float64}}}(undef,3)
    O_pn2002_t2 = Matrix{Matrix{Float64}}(undef,2,J_max+1)
    O_pn1111_t2 = Vector{Matrix{Matrix{Float64}}}(undef,4)
    O_pn0220_t2 = Matrix{Matrix{Float64}}(undef,2,J_max+1)
    O_nn_t2 = Vector{Matrix{Matrix{Float64}}}(undef,3)

    # Initialite sub-entries of O_t ...
    @inbounds for i in 1:3
        O_pp_t1[i] = Matrix{Matrix{Float64}}(undef,2,J_max+1)
        O_nn_t1[i] = Matrix{Matrix{Float64}}(undef,2,J_max+1)

        O_pp_t2[i] = Matrix{Matrix{Float64}}(undef,2,J_max+1)
        O_nn_t2[i] = Matrix{Matrix{Float64}}(undef,2,J_max+1)
    end
    @inbounds for i in 1:4
        O_pn1111_t1[i] = Matrix{Matrix{Float64}}(undef,2,J_max+1)

        O_pn1111_t2[i] = Matrix{Matrix{Float64}}(undef,2,J_max+1)
    end

 
    # Allocate entries of O ...
    @inbounds for P = 1:2
        @inbounds Threads.@threads for J=0:J_max
            Length = N_Orb_NN[P,J+1]

            @inbounds for i in 1:3
                O_pp_t1[i][P,J+1] = zeros(Float64,Length,Length)
                O_nn_t1[i][P,J+1] = zeros(Float64,Length,Length)
                O_pp_t2[i][P,J+1] = zeros(Float64,Length,Length)
                O_nn_t2[i][P,J+1] = zeros(Float64,Length,Length)
            end

            O_pn2002_t1[P,J+1] = zeros(Float64,Length,Length)
            O_pn0220_t1[P,J+1] = zeros(Float64,Length,Length)
            O_pn2002_t2[P,J+1] = zeros(Float64,Length,Length)
            O_pn0220_t2[P,J+1] = zeros(Float64,Length,Length)

            @inbounds for i in 1:4
                O_pn1111_t1[i][P,J+1] = zeros(Float64,Length,Length)
                O_pn1111_t2[i][P,J+1] = zeros(Float64,Length,Length)
            end


        end

    end

    # Allocate the NN operator in quasiparticle picture ...
    O_NN_t1 = O2B_22_Temp(O_pp_t1,O_pn2002_t1,O_pn1111_t1,O_pn0220_t1,O_nn_t1)
    O_NN_t2 = O2B_22_Temp(O_pp_t2,O_pn2002_t2,O_pn1111_t2,O_pn0220_t2,O_nn_t2)

    # Allocate the array for NN orbitals ...
    Orb_NN_t = Orb2B_Temp(Orb_NN_Dic,N_Orb_NN,Ind_Orb_NN)

    return O_NN_t1, O_NN_t2, Orb_NN_t
end

function qpO2b_22_allocate_ind1(Params::Parameters,JP::Vector{Vector{Int64}},Orb::Vector{Orb1B},Orb_NN::Orb2B,Orb_NN_t::Orb2B_Temp,O_NN::O2B,O_NN_t1::O2B_22_Temp,U::O1B,V::O1B)
    # Parameter initialization...
    N_max = Params.Calc.Nmax
    N_2max = 2*N_max
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Define function for the transformation MEs for the 1st index ...
    @inline function O2b_temp_transformation_ind1_MEs(J::Int64,P::Int64,i::Int64,a::Int64,b::Int64,c::Int64,d::Int64,U::O1B,V::O1B,O_NN::O2B,Orb::Vector{Orb1B},Orb_NN::Orb2B)
        return @views U.p[i,a] * O2b_pp(i,b,c,d,J,P,O_NN,Orb,Orb_NN), V.p[i,a] * O2b_pp(i,b,c,d,J,P,O_NN,Orb,Orb_NN),
                      U.p[i,a] * O2b_pn(i,b,c,d,J,P,O_NN,Orb_NN), V.p[i,a] * O2b_pn(i,b,c,d,J,P,O_NN,Orb_NN),
                      U.n[i,a] * O2b_nn(i,b,c,d,J,P,O_NN,Orb,Orb_NN), V.n[i,a] * O2b_nn(i,b,c,d,J,P,O_NN,Orb,Orb_NN)
    end

    # Transformation of the 1st index ...
    println("\nPerforming Quasiparticle Transformation of the residual 2-body NN interaction H22 in the 1st index...")
    @inbounds for i in JP
        J, P = i[1], i[2]
        if P == 1
            println("\tCalculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = +")
        else
            println("\tCalculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = -")
        end
        N = Orb_NN_t.N[P,J+1]
        @inbounds Threads.@threads for Bra in 1:N
            a = Orb_NN_t.Ind[P,J+1][Bra][1]
            b = Orb_NN_t.Ind[P,J+1][Bra][2]
            l_a, j_a = Orb[a].l, Orb[a].j

            Orb_x = Orb_PreComp(a_max,j_a,l_a,Orb)
            @inbounds for Ket in 1:N
                c = Orb_NN_t.Ind[P,J+1][Ket][1]
                d = Orb_NN_t.Ind[P,J+1][Ket][2]

                OppSum1, OppSum2, OppSum3 = 0.0, 0.0, 0.0
                Opn2002Sum, Opn1111Sum1, Opn1111Sum2, Opn1111Sum3, Opn1111Sum4, Opn0220Sum = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
                OnnSum1, OnnSum2, OnnSum3 = 0.0, 0.0, 0.0
                
                @inbounds for i in Orb_x
                    OppUME, OppVME, OpnUME, OpnVME, OnnUME, OnnVME = O2b_temp_transformation_ind1_MEs(J,P,i,a,b,c,d,U,V,O_NN,Orb,Orb_NN)

                    OppSum1 += OppUME
                    OppSum2 += OppUME
                    OppSum3 += OppVME

                    Opn2002Sum += OpnUME
                    Opn1111Sum1 += OpnUME
                    Opn1111Sum2 += OpnUME
                    Opn1111Sum3 += OpnVME
                    Opn1111Sum4 += OpnVME
                    Opn0220Sum += OpnVME

                    OnnSum1 += OnnUME
                    OnnSum2 += OnnUME
                    OnnSum3 += OnnVME
                end

                @views O_NN_t1.pp[1][P,J+1][Bra,Ket] = OppSum1
                @views O_NN_t1.pp[2][P,J+1][Bra,Ket] = OppSum2
                @views O_NN_t1.pp[3][P,J+1][Bra,Ket] = OppSum3

                @views O_NN_t1.pn2002[P,J+1][Bra,Ket] = Opn2002Sum
                @views O_NN_t1.pn1111[1][P,J+1][Bra,Ket] = Opn1111Sum1
                @views O_NN_t1.pn1111[2][P,J+1][Bra,Ket] = Opn1111Sum2
                @views O_NN_t1.pn1111[3][P,J+1][Bra,Ket] = Opn1111Sum3
                @views O_NN_t1.pn1111[4][P,J+1][Bra,Ket] = Opn1111Sum4
                @views O_NN_t1.pn0220[P,J+1][Bra,Ket] = Opn0220Sum

                @views O_NN_t1.nn[1][P,J+1][Bra,Ket] = OnnSum1
                @views O_NN_t1.nn[2][P,J+1][Bra,Ket] = OnnSum2
                @views O_NN_t1.nn[3][P,J+1][Bra,Ket] = OnnSum3
            end
        end
    end

    return O_NN_t1
end

function qpO2b_22_allocate_ind2(Params::Parameters,JP::Vector{Vector{Int64}},Orb::Vector{Orb1B},Orb_NN_t::Orb2B_Temp,O_NN_t1::O2B_22_Temp,O_NN_t2::O2B_22_Temp,U::O1B,V::O1B)
    # Parameter initialization...
    N_max = Params.Calc.Nmax
    N_2max = 2*N_max
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Define function for the transformation MEs for the 2nd index ...
    @inline function O2b_temp_transformation_ind2_MEs(P::Int64,J::Int64,Bra::Int64,Ket::Int64,i::Int64,b::Int64,T::Matrix{Float64},O_NN_t::Matrix{Matrix{Float64}})
        return @views T[i,b] * O_NN_t[P,J+1][Bra,Ket]
    end

    # Transformation of the 2nd index ...
    println("\nPerforming Quasiparticle Transformation of the residual 2-body NN interaction H22 in the 2nd index...")

    @inbounds for i in JP
        J, P = i[1], i[2]
        if P == 1
            println("\tCalculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = +")
        else
            println("\tCalculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = -")
        end
        N = Orb_NN_t.N[P,J+1]
        @inbounds Threads.@threads for Bra in 1:N
            a = Orb_NN_t.Ind[P,J+1][Bra][1]
            b = Orb_NN_t.Ind[P,J+1][Bra][2]
            l_b, j_b = Orb[b].l, Orb[b].j

            Orb_x = Orb_PreComp(a_max,j_b,l_b,Orb)
            @inbounds for Ket in 1:N
                c = Orb_NN_t.Ind[P,J+1][Ket][1]
                d = Orb_NN_t.Ind[P,J+1][Ket][2]

                OppSum1, OppSum2, OppSum3 = 0.0, 0.0, 0.0
                Opn2002Sum, Opn1111Sum1, Opn1111Sum2, Opn1111Sum3, Opn1111Sum4, Opn0220Sum = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
                OnnSum1, OnnSum2, OnnSum3 = 0.0, 0.0, 0.0
                
                @inbounds for i in Orb_x
                    Bra_ai = O2b_temp_index(a,i,J,P,Orb_NN_t)

                    OppME1 = O2b_temp_transformation_ind2_MEs(P,J,Bra_ai,Ket,i,b,U.p,O_NN_t1.pp[1])
                    OppME2 = O2b_temp_transformation_ind2_MEs(P,J,Bra_ai,Ket,i,b,V.p,O_NN_t1.pp[2])
                    OppME3 = O2b_temp_transformation_ind2_MEs(P,J,Bra_ai,Ket,i,b,V.p,O_NN_t1.pp[3])

                    Opn2002ME = O2b_temp_transformation_ind2_MEs(P,J,Bra_ai,Ket,i,b,V.n,O_NN_t1.pn2002)
                    Opn1111ME1 = O2b_temp_transformation_ind2_MEs(P,J,Bra_ai,Ket,i,b,U.n,O_NN_t1.pn1111[1])
                    Opn1111ME2 = O2b_temp_transformation_ind2_MEs(P,J,Bra_ai,Ket,i,b,V.n,O_NN_t1.pn1111[2])
                    Opn1111ME3 = O2b_temp_transformation_ind2_MEs(P,J,Bra_ai,Ket,i,b,U.n,O_NN_t1.pn1111[3])
                    Opn1111ME4 = O2b_temp_transformation_ind2_MEs(P,J,Bra_ai,Ket,i,b,V.n,O_NN_t1.pn1111[4])
                    Opn0220ME = O2b_temp_transformation_ind2_MEs(P,J,Bra_ai,Ket,i,b,U.n,O_NN_t1.pn0220)

                    OnnME1 = O2b_temp_transformation_ind2_MEs(P,J,Bra_ai,Ket,i,b,U.n,O_NN_t1.nn[1])
                    OnnME2 = O2b_temp_transformation_ind2_MEs(P,J,Bra_ai,Ket,i,b,V.n,O_NN_t1.nn[2])
                    OnnME3 = O2b_temp_transformation_ind2_MEs(P,J,Bra_ai,Ket,i,b,V.n,O_NN_t1.nn[3])

                    OppSum1 += OppME1
                    OppSum2 += OppME2
                    OppSum3 += OppME3

                    Opn2002Sum += Opn2002ME
                    Opn1111Sum1 += Opn1111ME1
                    Opn1111Sum2 += Opn1111ME2
                    Opn1111Sum3 += Opn1111ME3
                    Opn1111Sum4 += Opn1111ME4
                    Opn0220Sum += Opn0220ME

                    OnnSum1 += OnnME1
                    OnnSum2 += OnnME2
                    OnnSum3 += OnnME3

                end

                @views O_NN_t2.pp[1][P,J+1][Bra,Ket] = OppSum1
                @views O_NN_t2.pp[2][P,J+1][Bra,Ket] = OppSum2
                @views O_NN_t2.pp[3][P,J+1][Bra,Ket] = OppSum3

                @views O_NN_t2.pn2002[P,J+1][Bra,Ket] = Opn2002Sum
                @views O_NN_t2.pn1111[1][P,J+1][Bra,Ket] = Opn1111Sum1
                @views O_NN_t2.pn1111[2][P,J+1][Bra,Ket] = Opn1111Sum2
                @views O_NN_t2.pn1111[3][P,J+1][Bra,Ket] = Opn1111Sum3
                @views O_NN_t2.pn1111[4][P,J+1][Bra,Ket] = Opn1111Sum4
                @views O_NN_t2.pn0220[P,J+1][Bra,Ket] = Opn0220Sum

                @views O_NN_t2.nn[1][P,J+1][Bra,Ket] = OnnSum1
                @views O_NN_t2.nn[2][P,J+1][Bra,Ket] = OnnSum2
                @views O_NN_t2.nn[3][P,J+1][Bra,Ket] = OnnSum3

            end
        end
    end

    return O_NN_t2
end

function qpO2b_22_allocate_ind3(Params::Parameters,JP::Vector{Vector{Int64}},Orb::Vector{Orb1B},Orb_NN_t::Orb2B_Temp,O_NN_t1::O2B_22_Temp,O_NN_t2::O2B_22_Temp,U::O1B,V::O1B)
    # Parameter initialization...
    N_max = Params.Calc.Nmax
    N_2max = 2*N_max
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Define function for the transformation MEs for the 3rd index ...
    @inline function O2b_temp_transformation_ind3_MEs(P::Int64,J::Int64,Bra::Int64,Ket::Int64,i::Int64,c::Int64,T::Matrix{Float64},O_NN_t::Matrix{Matrix{Float64}})
        return @views T[i,c] * O_NN_t[P,J+1][Bra,Ket]
    end

    # Transformation of the 3rd index ...
    println("\nPerforming Quasiparticle (22) Transformation of the given 2-body NN operator in the 3rd index...")

    @inbounds for i in JP
        J, P = i[1], i[2]
        if P == 1
            println("\tCalculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = +")
        else
            println("\tCalculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = -")
        end
        N = Orb_NN_t.N[P,J+1]
        @inbounds Threads.@threads for Ket in 1:N
            c = Orb_NN_t.Ind[P,J+1][Ket][1]
            d = Orb_NN_t.Ind[P,J+1][Ket][2]
            l_c, j_c = Orb[c].l, Orb[c].j

            Orb_x = Orb_PreComp(a_max,j_c,l_c,Orb)

            @inbounds for Bra in 1:N

                OppSum1, OppSum2, OppSum3 = 0.0, 0.0, 0.0
                Opn2002Sum, Opn1111Sum1, Opn1111Sum2, Opn1111Sum3, Opn1111Sum4, Opn0220Sum = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
                OnnSum1, OnnSum2, OnnSum3 = 0.0, 0.0, 0.0

                @inbounds for i in Orb_x
                    Ket_id = O2b_temp_index(i,d,J,P,Orb_NN_t)

                    OppME1 = O2b_temp_transformation_ind3_MEs(P,J,Bra,Ket_id,i,c,U.p,O_NN_t2.pp[1])
                    OppME2 = O2b_temp_transformation_ind3_MEs(P,J,Bra,Ket_id,i,c,U.p,O_NN_t2.pp[2])
                    OppME3 = O2b_temp_transformation_ind3_MEs(P,J,Bra,Ket_id,i,c,V.p,O_NN_t2.pp[3])

                    Opn2002ME = O2b_temp_transformation_ind3_MEs(P,J,Bra,Ket_id,i,c,V.p,O_NN_t2.pn2002)
                    Opn1111ME1 = O2b_temp_transformation_ind3_MEs(P,J,Bra,Ket_id,i,c,U.p,O_NN_t2.pn1111[1])
                    Opn1111ME2 = O2b_temp_transformation_ind3_MEs(P,J,Bra,Ket_id,i,c,U.p,O_NN_t2.pn1111[2])
                    Opn1111ME3 = O2b_temp_transformation_ind3_MEs(P,J,Bra,Ket_id,i,c,V.p,O_NN_t2.pn1111[3])
                    Opn1111ME4 = O2b_temp_transformation_ind3_MEs(P,J,Bra,Ket_id,i,c,V.p,O_NN_t2.pn1111[4])
                    Opn0220ME = O2b_temp_transformation_ind3_MEs(P,J,Bra,Ket_id,i,c,U.p,O_NN_t2.pn0220)

                    OnnME1 = O2b_temp_transformation_ind3_MEs(P,J,Bra,Ket_id,i,c,U.n,O_NN_t2.nn[1])
                    OnnME2 = O2b_temp_transformation_ind3_MEs(P,J,Bra,Ket_id,i,c,U.n,O_NN_t2.nn[2])
                    OnnME3 = O2b_temp_transformation_ind3_MEs(P,J,Bra,Ket_id,i,c,V.n,O_NN_t2.nn[3])

                    OppSum1 += OppME1
                    OppSum2 += OppME2
                    OppSum3 += OppME3

                    Opn2002Sum += Opn2002ME
                    Opn1111Sum1 += Opn1111ME1
                    Opn1111Sum2 += Opn1111ME2
                    Opn1111Sum3 += Opn1111ME3
                    Opn1111Sum4 += Opn1111ME4
                    Opn0220Sum += Opn0220ME

                    OnnSum1 += OnnME1
                    OnnSum2 += OnnME2
                    OnnSum3 += OnnME3
                end

                @views O_NN_t1.pp[1][P,J+1][Bra,Ket] = OppSum1
                @views O_NN_t1.pp[2][P,J+1][Bra,Ket] = OppSum2
                @views O_NN_t1.pp[3][P,J+1][Bra,Ket] = OppSum3

                @views O_NN_t1.pn2002[P,J+1][Bra,Ket] = Opn2002Sum
                @views O_NN_t1.pn1111[1][P,J+1][Bra,Ket] = Opn1111Sum1
                @views O_NN_t1.pn1111[2][P,J+1][Bra,Ket] = Opn1111Sum2
                @views O_NN_t1.pn1111[3][P,J+1][Bra,Ket] = Opn1111Sum3
                @views O_NN_t1.pn1111[4][P,J+1][Bra,Ket] = Opn1111Sum4
                @views O_NN_t1.pn0220[P,J+1][Bra,Ket] = Opn0220Sum

                @views O_NN_t1.nn[1][P,J+1][Bra,Ket] = OnnSum1
                @views O_NN_t1.nn[2][P,J+1][Bra,Ket] = OnnSum2
                @views O_NN_t1.nn[3][P,J+1][Bra,Ket] = OnnSum3
            end
        end
    end

    return O_NN_t1
end

function qpO2b_22_allocate_ind4(Params::Parameters,JP::Vector{Vector{Int64}},Orb::Vector{Orb1B},Orb_NN::Orb2B,Orb_NN_t::Orb2B_Temp,O_NN::qpO2B,O_NN_t1::O2B_22_Temp,O_NN_t2::O2B_22_Temp,U::O1B,V::O1B)
    # Parameter initialization ...
    N_max = Params.Calc.Nmax
    N_2max = 2*N_max
    J_max = 2*N_max + 1
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Define function for the transformation MEs for the 4th index ...
    @inline function O2b_temp_transformation_ind4_MEs(P::Int64,J::Int64,Bra::Int64,Ket::Int64,i::Int64,d::Int64,T::Matrix{Float64},O_NN_t::Matrix{Matrix{Float64}})
        return @views T[i,d] * O_NN_t[P,J+1][Bra,Ket]
    end

    # Transformation of the 4th index ...
    println("\nPerforming Quasiparticle Transformation of the residual 2-body NN interaction H^(22) in the 4th index...")
        # Preallocate O_NN_t2 ... for the Pandya transformation ...
        println("\tPreallocating temporary array for the generalized Pandya transformation ...")
        @inbounds for i in JP
            J, P = i[1], i[2]
            if P == 1
                println("\t\tCalculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = +")
            else
                println("\t\tCalculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = -")
            end
            N = Orb_NN_t.N[P,J+1]
            @inbounds Threads.@threads for Ket in 1:N
                c = Orb_NN_t.Ind[P,J+1][Ket][1]
                d = Orb_NN_t.Ind[P,J+1][Ket][2]
                l_d, j_d = Orb[d].l, Orb[d].j
                
                Orb_x = Orb_PreComp(a_max,j_d,l_d,Orb)

                @inbounds for Bra in 1:N

                    OppSum = 0.0
                    Opn2002Sum, Opn1111Sum, Opn0220Sum = 0.0, 0.0, 0.0
                    OnnSum = 0.0

                    @inbounds for i in Orb_x
                        Ket_ci = O2b_temp_index(c,i,J,P,Orb_NN_t)

                        OppME = O2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,d,V.p,O_NN_t1.pp[2])

                        Opn2002ME = O2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,d,U.n,O_NN_t1.pn2002)
                        Opn1111ME2 = O2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,d,V.n,O_NN_t1.pn1111[2])
                        Opn1111ME3 = O2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,d,U.n,O_NN_t1.pn1111[3])
                        Opn0220ME = O2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,d,V.n,O_NN_t1.pn0220)

                        OnnME = O2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,d,V.n,O_NN_t1.nn[2])

                        OppSum += OppME

                        Opn2002Sum += Opn2002ME

                        Opn1111Sum += (Opn1111ME2 + Opn1111ME3)

                        Opn0220Sum += Opn0220ME

                        OnnSum += OnnME
                    end

                    @views O_NN_t2.pp[2][P,J+1][Bra,Ket] = OppSum

                    @views O_NN_t2.pn2002[P,J+1][Bra,Ket] = Opn2002Sum

                    @views O_NN_t2.pn1111[2][P,J+1][Bra,Ket] = Opn1111Sum

                    @views O_NN_t2.pn0220[P,J+1][Bra,Ket] = Opn0220Sum

                    @views O_NN_t2.nn[2][P,J+1][Bra,Ket] = OnnSum
                end
            end
        end

        # Allocate O^(22) ...
        println("\tAllocating the components of O^(22) ...")
        @inbounds for i in JP
            J, P = i[1], i[2]
            if P == 1
                println("\t\tCalculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = +")
            else
                println("\t\tCalculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = -")
            end

            N_T0, N_T1 = Orb_NN.N[1,P,J+1], Orb_NN.N[2,P,J+1]

            @inbounds Threads.@threads for Bra in 1:max(N_T0,N_T1)
                @inbounds for Ket in 1:Bra

                    # Case of pn entries ... T = 0
                    if Bra <= N_T0
                        Ind = Bra + (Ket - 1) * N_T0 - div(Ket * (Ket - 1),2)
                        a, b = Orb_NN.Ind[1,P,J+1][Bra][1], Orb_NN.Ind[1,P,J+1][Bra][2]
                        c, d = Orb_NN.Ind[1,P,J+1][Ket][1], Orb_NN.Ind[1,P,J+1][Ket][2]
                        j_a, j_b, j_c, j_d = Orb[a].j, Orb[b].j, Orb[c].j, Orb[d].j
                        l_a, l_b, l_c, l_d = Orb[a].l, Orb[b].l, Orb[c].l, Orb[d].l

                        Bra_ab = O2b_temp_index(a,b,J,P,Orb_NN_t)

                        Orb_x = Orb_PreComp(a_max,j_d,l_d,Orb)

                        Opn2002Sum, Opn1111Sum, Opn0220Sum = 0.0, 0.0, 0.0

                        @inbounds for i in Orb_x
                            Ket_ci = O2b_temp_index(c,i,J,P,Orb_NN_t)

                            Opn1111ME1 = O2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,d,U.n,O_NN_t1.pn1111[1])
                            Opn1111ME2 = O2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,d,V.n,O_NN_t1.pn1111[4])

                            Opn1111Sum += (Opn1111ME1 + Opn1111ME2)

                        end

                        Amp11 = Float64((-1)^(J))
                        Opn1111Sum = Amp11 * Opn1111Sum

                        @inbounds for I in 0:J_max

                            # pnO2002 && pnO0220 ...
                            if abs(j_a - j_c) <= 2*I && 2*I <= (j_a + j_c) && abs(j_b - j_d) <= 2*I && 2*I <= (j_b + j_d) &&
                                (rem(l_a + l_c,2) + 1) == P && (rem(l_b + l_d,2) + 1) == P
                            #if abs(j_a - j_b) <= 2*I && 2*I <= (j_a + j_b) && abs(j_c - j_d) <= 2*I && 2*I <= (j_c + j_d)

                                #=
                                Amp2002 = Float64((-1)^(J + I + div(j_b + j_c,2)) * (2*I + 1)) * f6j(j_a,j_c,2*J,j_d,j_b,2*I)
                                Amp0220 = Amp2002
            
                                Bra_ab, Ket_cd = O2b_temp_index(a,b,I,P,Orb_NN_t), O2b_temp_index(c,d,I,P,Orb_NN_t)

                                Opn2002Sum += Amp2002 * O_NN_t2.pn2002[P,I+1][Bra_ab,Ket_cd]
                                Opn0220Sum += Amp0220 * O_NN_t2.pn0220[P,I+1][Bra_ab,Ket_cd]
                                
                                # gives the same results as the one below ...
                                =#

                                
                                Amp2002 = Float64((-1)^(J + I + div(j_b + j_c,2)) * (2*I + 1)) * f6j(j_a,j_b,2*J,j_d,j_c,2*I) * sqrt(Float64((1 + kronecker_delta(a,c) * (-1)^I) * (1 + kronecker_delta(b,d) * (-1)^I)))
                                Amp0220 = Amp2002
            
                                Bra_ac, Ket_bd = O2b_temp_index(a,c,I,P,Orb_NN_t), O2b_temp_index(b,d,I,P,Orb_NN_t)

                                Opn2002Sum += Amp2002 * O_NN_t2.pn2002[P,I+1][Bra_ac,Ket_bd]
                                Opn0220Sum += Amp0220 * O_NN_t2.pn0220[P,I+1][Bra_ac,Ket_bd]
                                
                            end

                            # pnO1111 ...
                            if abs(j_a - j_d) <= 2*I && 2*I <= (j_a + j_d) && abs(j_c - j_b) <= 2*I &&
                            2*I <= (j_c + j_b) && (rem(l_a + l_d,2) + 1 == P) && (rem(l_b + l_c,2) + 1 == P)
                                Amp1111 = Float64((-1)^(J) * (2*I + 1)) * f6j(j_a,j_b,2*J,j_c,j_d,2*I)

                                Bra_ad, Ket_cb = O2b_temp_index(a,d,I,P,Orb_NN_t), O2b_temp_index(c,b,I,P,Orb_NN_t)

                                Opn1111Sum += Amp1111 * O_NN_t2.pn1111[2][P,I+1][Bra_ad,Ket_cb]
                            end

                        end

                        @views O_NN.qp22.pn2002[P,J+1][Ind] = Opn2002Sum
                        @views O_NN.qp22.pn1111[P,J+1][Ind] = Opn1111Sum
                        @views O_NN.qp22.pn0220[P,J+1][Ind] = Opn0220Sum
                    end

                    # Case of pp & nn entries ... T = 1
                    if Bra <= N_T1
                        Ind = Bra + (Ket - 1) * N_T1 - div(Ket * (Ket - 1),2)
                        a, b = Orb_NN.Ind[2,P,J+1][Bra][1], Orb_NN.Ind[2,P,J+1][Bra][2]
                        c, d = Orb_NN.Ind[2,P,J+1][Ket][1], Orb_NN.Ind[2,P,J+1][Ket][2]
                        j_a, j_b, j_c, j_d = Orb[a].j, Orb[b].j, Orb[c].j, Orb[d].j
                        l_a, l_b, l_c, l_d = Orb[a].l, Orb[b].l, Orb[c].l, Orb[d].l

                        Bra_ab = O2b_temp_index(a,b,J,P,Orb_NN_t)

                        Orb_x = Orb_PreComp(a_max,j_d,l_d,Orb)

                        OppSum, OnnSum = 0.0, 0.0

                        @inbounds for i in Orb_x
                            Ket_ci = O2b_temp_index(c,i,J,P,Orb_NN_t)

                            OppME1 = O2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,d,U.p,O_NN_t1.pp[1])
                            OppME3 = O2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,d,V.p,O_NN_t1.pp[3])

                            OnnME1 = O2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,d,U.n,O_NN_t1.nn[1])
                            OnnME3 = O2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,d,V.n,O_NN_t1.nn[3])

                            OppSum += (OppME1 + OppME3)
                            OnnSum += (OnnME1 + OnnME3)
                        end

                        Amp = Float64((-1)^(J + 1)) * sqrt(Float64((1 + kronecker_delta(a,b) * (-1)^J) * (1 + kronecker_delta(c,d) * (-1)^J))) # My og and Suhonen ...

                        OppSum = Amp * OppSum
                        OnnSum = Amp * OnnSum

                        # NNO22 ...
                        @inbounds for I in 0:J_max
                            if abs(j_a - j_d) <= 2*I && 2*I <= (j_a + j_d) && abs(j_b - j_c) <= 2*I &&
                            2*I <= (j_b + j_c) && (rem(l_a + l_d,2) + 1 == P) && (rem(l_b + l_c,2) + 1 == P)
                                AmpR = 4.0 * Float64((-1)^(J) * (2*I + 1)) * f6j(j_a,j_b,2*J,j_c,j_d,2*I) * sqrt(Float64((1 + kronecker_delta(a,d) * (-1)^I) * (1 + kronecker_delta(b,c) * (-1)^I))) # My og ... also in Suhonen ...

                                Bra_ad, Ket_cb = O2b_temp_index(a,d,I,P,Orb_NN_t), O2b_temp_index(c,b,I,P,Orb_NN_t)

                                OppSum += AmpR * O_NN_t2.pp[2][P,I+1][Bra_ad,Ket_cb]
                                OnnSum += AmpR * O_NN_t2.nn[2][P,I+1][Bra_ad,Ket_cb]
                            end
                        end

                        @views O_NN.qp22.pp[P,J+1][Ind] = OppSum
                        @views O_NN.qp22.nn[P,J+1][Ind] = OnnSum
                    end

                end
            end

        end

    return O_NN
end

@inline function qpO2b_22_pp(a::Int64,b::Int64,c::Int64,d::Int64,J::Int64,P::Int64,O_NN::qpO2B,Orb::Vector{Orb1B},Orb_NN::Orb2B)
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
    return Amp * O_NN.qp22.pp[P,J+1][Ind]
end

@inline function qpO2b_2002_pn(a::Int64,b::Int64,c::Int64,d::Int64,J::Int64,P::Int64,O_NN::qpO2B,Orb_NN::Orb2B)
    @inbounds return O_NN.qp22.pn2002[P,J+1][O2b_index(a,b,c,d,J,P,0,Orb_NN)]
end

@inline function qpO2b_1111_pn(a::Int64,b::Int64,c::Int64,d::Int64,J::Int64,P::Int64,O_NN::qpO2B,Orb_NN::Orb2B)
    @inbounds return O_NN.qp22.pn1111[P,J+1][O2b_index(a,b,c,d,J,P,0,Orb_NN)]
end

@inline function qpO2b_0220_pn(a::Int64,b::Int64,c::Int64,d::Int64,J::Int64,P::Int64,O_NN::qpO2B,Orb_NN::Orb2B)
    @inbounds return O_NN.qp22.pn0220[P,J+1][O2b_index(a,b,c,d,J,P,0,Orb_NN)]
end

@inline function qpO2b_22_nn(a::Int64,b::Int64,c::Int64,d::Int64,J::Int64,P::Int64,O_NN::qpO2B,Orb::Vector{Orb1B},Orb_NN::Orb2B)
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
    return Amp * O_NN.qp22.nn[P,J+1][Ind]
end
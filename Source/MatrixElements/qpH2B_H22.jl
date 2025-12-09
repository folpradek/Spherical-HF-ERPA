function qpH2b_allocate_H22(Params::Parameters,Orb::Vector{NOrb},Orb_NN::NNOrb,H_NN::qpH2B,V_NN::NNInt,U::pnMatrix,V::pnMatrix)
    println("\nAllocating H^(22) components of the quasiparticle Hamiltonian ...")

    # Preallocate arrays of allowed values of J & P ...
    JP = JP_initialize(Params.Calc.N2max + 1)

    # Initialize temporary V & W interaction arrays ...
    V_NN_t1, V_NN_t2, Orb_NN_t = qpH2b_initialize_H22(Params,Orb)

    # Transformation of the 1st index ...
    @time V_NN_t1 = qpH2b_allocate_H22_ind1(Params,JP,Orb,Orb_NN,Orb_NN_t,V_NN,V_NN_t1,U,V)

    # Transformation of the 2nd index ...
    @time V_NN_t2 = qpH2b_allocate_H22_ind2(Params,JP,Orb,Orb_NN_t,V_NN_t1,V_NN_t2,U,V)

    # Transformation of the 3rd index ...
    @time V_NN_t1 = qpH2b_allocate_H22_ind3(Params,JP,Orb,Orb_NN_t,V_NN_t1,V_NN_t2,U,V)

    # Transformation of the 4th index ...
    @time H_NN = qpH2b_allocate_H22_ind4(Params,JP,Orb,Orb_NN,Orb_NN_t,H_NN,V_NN_t1,V_NN_t2,U,V)

    # Deallocate V_NN_t & W_NN_t ...
    V_NN_t1, W_NN_t1 = nothing, nothing
    V_NN_t2, W_NN_t2 = nothing, nothing

    # Perform the Garbage Collection ...
    GC.gc()

    println("\nAllocation of H^(22) components of the quasiparticle Hamiltonian done ...")

    return H_NN
end

function qpH2b_initialize_H22(Params::Parameters,Orb::Vector{NOrb})
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
        # Case of V_t1 ...
    V_pp_t1 = Vector{Matrix{Matrix{Float64}}}(undef,3)
    V_pn2002_t1 = Matrix{Matrix{Float64}}(undef,2,J_max+1)
    V_pn1111_t1 = Vector{Matrix{Matrix{Float64}}}(undef,4)
    V_pn0220_t1 = Matrix{Matrix{Float64}}(undef,2,J_max+1)
    V_nn_t1 = Vector{Matrix{Matrix{Float64}}}(undef,3)
        # Case of V_t2 ...
    V_pp_t2 = Vector{Matrix{Matrix{Float64}}}(undef,3)
    V_pn2002_t2 = Matrix{Matrix{Float64}}(undef,2,J_max+1)
    V_pn1111_t2 = Vector{Matrix{Matrix{Float64}}}(undef,4)
    V_pn0220_t2 = Matrix{Matrix{Float64}}(undef,2,J_max+1)
    V_nn_t2 = Vector{Matrix{Matrix{Float64}}}(undef,3)

    # Initialite sub-entries of V_t ...
    @inbounds for i in 1:3
        V_pp_t1[i] = Matrix{Matrix{Float64}}(undef,2,J_max+1)
        V_nn_t1[i] = Matrix{Matrix{Float64}}(undef,2,J_max+1)

        V_pp_t2[i] = Matrix{Matrix{Float64}}(undef,2,J_max+1)
        V_nn_t2[i] = Matrix{Matrix{Float64}}(undef,2,J_max+1)
    end
    @inbounds for i in 1:4
        V_pn1111_t1[i] = Matrix{Matrix{Float64}}(undef,2,J_max+1)

        V_pn1111_t2[i] = Matrix{Matrix{Float64}}(undef,2,J_max+1)
    end

 
    # Allocate entries of V & W ...
    @inbounds for P = 1:2
        @inbounds Threads.@threads for J=0:J_max
            Length = N_Orb_NN[P,J+1]

            @inbounds for i in 1:3
                V_pp_t1[i][P,J+1] = zeros(Float64,Length,Length)
                V_nn_t1[i][P,J+1] = zeros(Float64,Length,Length)

                V_pp_t2[i][P,J+1] = zeros(Float64,Length,Length)
                V_nn_t2[i][P,J+1] = zeros(Float64,Length,Length)
            end

            V_pn2002_t1[P,J+1] = zeros(Float64,Length,Length)
            V_pn0220_t1[P,J+1] = zeros(Float64,Length,Length)

            V_pn2002_t2[P,J+1] = zeros(Float64,Length,Length)
            V_pn0220_t2[P,J+1] = zeros(Float64,Length,Length)

            @inbounds for i in 1:4
                V_pn1111_t1[i][P,J+1] = zeros(Float64,Length,Length)

                V_pn1111_t2[i][P,J+1] = zeros(Float64,Length,Length)
            end


        end

    end

    # Allocate the NN interaction in quasiparticle picture ...
    V_NN_t1 = V2B_H22_Temp(V_pp_t1,V_pn2002_t1,V_pn1111_t1,V_pn0220_t1,V_nn_t1)
    V_NN_t2 = V2B_H22_Temp(V_pp_t2,V_pn2002_t2,V_pn1111_t2,V_pn0220_t2,V_nn_t2)

    # Allocate the array for NN orbitals ...
    Orb_NN_t = NNOrb_Temp(Orb_NN_Dic,N_Orb_NN,Ind_Orb_NN)

    return V_NN_t1, V_NN_t2, Orb_NN_t
end

function qpH2b_allocate_H22_ind1(Params::Parameters,JP::Vector{Vector{Int64}},Orb::Vector{NOrb},Orb_NN::NNOrb,Orb_NN_t::NNOrb_Temp,V_NN::NNInt,V_NN_t1::V2B_H22_Temp,U::pnMatrix,V::pnMatrix)
    # Parameter initialization...
    N_max = Params.Calc.Nmax
    N_2max = 2*N_max
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Define function for the transformation MEs for the 1st index ...
    @inline function V2b_temp_transformation_ind1_MEs(J::Int64,i::Int64,a::Int64,b::Int64,c::Int64,d::Int64,U::pnMatrix,V::pnMatrix,V_NN::NNInt,Orb::Vector{NOrb},Orb_NN::NNOrb)
        return @views U.p[i,a] * V2B(i,b,c,d,J,1,V_NN.pp,Orb,Orb_NN), V.p[i,a] * V2B(i,b,c,d,J,1,V_NN.pp,Orb,Orb_NN),
                      U.p[i,a] * V2B(i,b,c,d,J,0,V_NN.pn,Orb,Orb_NN), V.p[i,a] * V2B(i,b,c,d,J,0,V_NN.pn,Orb,Orb_NN),
                      U.n[i,a] * V2B(i,b,c,d,J,1,V_NN.nn,Orb,Orb_NN), V.n[i,a] * V2B(i,b,c,d,J,1,V_NN.nn,Orb,Orb_NN)
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

                VppSum1, VppSum2, VppSum3 = 0.0, 0.0, 0.0
                Vpn2002Sum, Vpn1111Sum1, Vpn1111Sum2, Vpn1111Sum3, Vpn1111Sum4, Vpn0220Sum = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
                VnnSum1, VnnSum2, VnnSum3 = 0.0, 0.0, 0.0
                
                @inbounds for i in Orb_x
                    VppUME, VppVME, VpnUME, VpnVME, VnnUME, VnnVME = V2b_temp_transformation_ind1_MEs(J,i,a,b,c,d,U,V,V_NN,Orb,Orb_NN)

                    VppSum1 += VppUME
                    VppSum2 += VppUME
                    VppSum3 += VppVME

                    Vpn2002Sum += VpnUME
                    Vpn1111Sum1 += VpnUME
                    Vpn1111Sum2 += VpnUME
                    Vpn1111Sum3 += VpnVME
                    Vpn1111Sum4 += VpnVME
                    Vpn0220Sum += VpnVME

                    VnnSum1 += VnnUME
                    VnnSum2 += VnnUME
                    VnnSum3 += VnnVME

                end

                @views V_NN_t1.pp[1][P,J+1][Bra,Ket] = VppSum1
                @views V_NN_t1.pp[2][P,J+1][Bra,Ket] = VppSum2
                @views V_NN_t1.pp[3][P,J+1][Bra,Ket] = VppSum3

                @views V_NN_t1.pn2002[P,J+1][Bra,Ket] = Vpn2002Sum
                @views V_NN_t1.pn1111[1][P,J+1][Bra,Ket] = Vpn1111Sum1
                @views V_NN_t1.pn1111[2][P,J+1][Bra,Ket] = Vpn1111Sum2
                @views V_NN_t1.pn1111[3][P,J+1][Bra,Ket] = Vpn1111Sum3
                @views V_NN_t1.pn1111[4][P,J+1][Bra,Ket] = Vpn1111Sum4
                @views V_NN_t1.pn0220[P,J+1][Bra,Ket] = Vpn0220Sum

                @views V_NN_t1.nn[1][P,J+1][Bra,Ket] = VnnSum1
                @views V_NN_t1.nn[2][P,J+1][Bra,Ket] = VnnSum2
                @views V_NN_t1.nn[3][P,J+1][Bra,Ket] = VnnSum3

            end
        end
    end

    return V_NN_t1
end

function qpH2b_allocate_H22_ind2(Params::Parameters,JP::Vector{Vector{Int64}},Orb::Vector{NOrb},Orb_NN_t::NNOrb_Temp,V_NN_t1::V2B_H22_Temp,V_NN_t2::V2B_H22_Temp,U::pnMatrix,V::pnMatrix)
    # Parameter initialization...
    N_max = Params.Calc.Nmax
    N_2max = 2*N_max
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Define function for the transformation MEs for the 2nd index ...
        # Case of V ...
    @inline function V2b_temp_transformation_ind2_MEs(P::Int64,J::Int64,Bra::Int64,Ket::Int64,i::Int64,b::Int64,T::Matrix{Float64},V_NN_t::Matrix{Matrix{Float64}})
        return @views T[i,b] * V_NN_t[P,J+1][Bra,Ket]
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

                VppSum1, VppSum2, VppSum3 = 0.0, 0.0, 0.0
                Vpn2002Sum, Vpn1111Sum1, Vpn1111Sum2, Vpn1111Sum3, Vpn1111Sum4, Vpn0220Sum = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
                VnnSum1, VnnSum2, VnnSum3 = 0.0, 0.0, 0.0
                
                @inbounds for i in Orb_x
                    Bra_ai = V2b_temp_index(a,i,J,P,Orb_NN_t)

                    VppME1 = V2b_temp_transformation_ind2_MEs(P,J,Bra_ai,Ket,i,b,U.p,V_NN_t1.pp[1])
                    VppME2 = V2b_temp_transformation_ind2_MEs(P,J,Bra_ai,Ket,i,b,V.p,V_NN_t1.pp[2])
                    VppME3 = V2b_temp_transformation_ind2_MEs(P,J,Bra_ai,Ket,i,b,V.p,V_NN_t1.pp[3])

                    Vpn2002ME = V2b_temp_transformation_ind2_MEs(P,J,Bra_ai,Ket,i,b,V.n,V_NN_t1.pn2002)
                    Vpn1111ME1 = V2b_temp_transformation_ind2_MEs(P,J,Bra_ai,Ket,i,b,U.n,V_NN_t1.pn1111[1])
                    Vpn1111ME2 = V2b_temp_transformation_ind2_MEs(P,J,Bra_ai,Ket,i,b,U.n,V_NN_t1.pn1111[2])
                    Vpn1111ME3 = V2b_temp_transformation_ind2_MEs(P,J,Bra_ai,Ket,i,b,U.n,V_NN_t1.pn1111[3])
                    Vpn1111ME4 = V2b_temp_transformation_ind2_MEs(P,J,Bra_ai,Ket,i,b,V.n,V_NN_t1.pn1111[4])
                    Vpn0220ME = V2b_temp_transformation_ind2_MEs(P,J,Bra_ai,Ket,i,b,U.n,V_NN_t1.pn0220)

                    VnnME1 = V2b_temp_transformation_ind2_MEs(P,J,Bra_ai,Ket,i,b,U.n,V_NN_t1.nn[1])
                    VnnME2 = V2b_temp_transformation_ind2_MEs(P,J,Bra_ai,Ket,i,b,V.n,V_NN_t1.nn[2])
                    VnnME3 = V2b_temp_transformation_ind2_MEs(P,J,Bra_ai,Ket,i,b,V.n,V_NN_t1.nn[3])

                    VppSum1 += VppME1
                    VppSum2 += VppME2
                    VppSum3 += VppME3

                    Vpn2002Sum += Vpn2002ME
                    Vpn1111Sum1 += Vpn1111ME1
                    Vpn1111Sum2 += Vpn1111ME2
                    Vpn1111Sum3 += Vpn1111ME3
                    Vpn1111Sum4 += Vpn1111ME4
                    Vpn0220Sum += Vpn0220ME

                    VnnSum1 += VnnME1
                    VnnSum2 += VnnME2
                    VnnSum3 += VnnME3

                end

                @views V_NN_t2.pp[1][P,J+1][Bra,Ket] = VppSum1
                @views V_NN_t2.pp[2][P,J+1][Bra,Ket] = VppSum2
                @views V_NN_t2.pp[3][P,J+1][Bra,Ket] = VppSum3

                @views V_NN_t2.pn2002[P,J+1][Bra,Ket] = Vpn2002Sum
                @views V_NN_t2.pn1111[1][P,J+1][Bra,Ket] = Vpn1111Sum1
                @views V_NN_t2.pn1111[2][P,J+1][Bra,Ket] = Vpn1111Sum2
                @views V_NN_t2.pn1111[3][P,J+1][Bra,Ket] = Vpn1111Sum3
                @views V_NN_t2.pn1111[4][P,J+1][Bra,Ket] = Vpn1111Sum4
                @views V_NN_t2.pn0220[P,J+1][Bra,Ket] = Vpn0220Sum

                @views V_NN_t2.nn[1][P,J+1][Bra,Ket] = VnnSum1
                @views V_NN_t2.nn[2][P,J+1][Bra,Ket] = VnnSum2
                @views V_NN_t2.nn[3][P,J+1][Bra,Ket] = VnnSum3

            end
        end
    end

    return V_NN_t2
end

function qpH2b_allocate_H22_ind3(Params::Parameters,JP::Vector{Vector{Int64}},Orb::Vector{NOrb},Orb_NN_t::NNOrb_Temp,V_NN_t1::V2B_H22_Temp,V_NN_t2::V2B_H22_Temp,U::pnMatrix,V::pnMatrix)
    # Parameter initialization...
    N_max = Params.Calc.Nmax
    N_2max = 2*N_max
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Define function for the transformation MEs for the 3rd index ...
        # Case of V ...
    @inline function V2b_temp_transformation_ind3_MEs(P::Int64,J::Int64,Bra::Int64,Ket::Int64,i::Int64,c::Int64,T::Matrix{Float64},V_NN_t::Matrix{Matrix{Float64}})
        return @views T[i,c] * V_NN_t[P,J+1][Bra,Ket]
    end

    # Transformation of the 3rd index ...
    println("\nPerforming Quasiparticle Transformation of the residual 2-body NN interaction H22 in the 3rd index...")

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

                VppSum1, VppSum2, VppSum3 = 0.0, 0.0, 0.0
                Vpn2002Sum, Vpn1111Sum1, Vpn1111Sum2, Vpn1111Sum3, Vpn1111Sum4, Vpn0220Sum = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
                VnnSum1, VnnSum2, VnnSum3 = 0.0, 0.0, 0.0

                @inbounds for i in Orb_x
                    Ket_id = V2b_temp_index(i,d,J,P,Orb_NN_t)

                    VppME1 = V2b_temp_transformation_ind3_MEs(P,J,Bra,Ket_id,i,c,U.p,V_NN_t2.pp[1])
                    VppME2 = V2b_temp_transformation_ind3_MEs(P,J,Bra,Ket_id,i,c,U.p,V_NN_t2.pp[2])
                    VppME3 = V2b_temp_transformation_ind3_MEs(P,J,Bra,Ket_id,i,c,V.p,V_NN_t2.pp[3])

                    Vpn2002ME = V2b_temp_transformation_ind3_MEs(P,J,Bra,Ket_id,i,c,V.p,V_NN_t2.pn2002)
                    Vpn1111ME1 = V2b_temp_transformation_ind3_MEs(P,J,Bra,Ket_id,i,c,U.p,V_NN_t2.pn1111[1])
                    Vpn1111ME2 = V2b_temp_transformation_ind3_MEs(P,J,Bra,Ket_id,i,c,U.p,V_NN_t2.pn1111[2])
                    Vpn1111ME3 = V2b_temp_transformation_ind3_MEs(P,J,Bra,Ket_id,i,c,V.p,V_NN_t2.pn1111[3])
                    Vpn1111ME4 = V2b_temp_transformation_ind3_MEs(P,J,Bra,Ket_id,i,c,V.p,V_NN_t2.pn1111[4])
                    Vpn0220ME = V2b_temp_transformation_ind3_MEs(P,J,Bra,Ket_id,i,c,U.p,V_NN_t2.pn0220)

                    VnnME1 = V2b_temp_transformation_ind3_MEs(P,J,Bra,Ket_id,i,c,U.n,V_NN_t2.nn[1])
                    VnnME2 = V2b_temp_transformation_ind3_MEs(P,J,Bra,Ket_id,i,c,U.n,V_NN_t2.nn[2])
                    VnnME3 = V2b_temp_transformation_ind3_MEs(P,J,Bra,Ket_id,i,c,V.n,V_NN_t2.nn[3])

                    VppSum1 += VppME1
                    VppSum2 += VppME2
                    VppSum3 += VppME3

                    Vpn2002Sum += Vpn2002ME
                    Vpn1111Sum1 += Vpn1111ME1
                    Vpn1111Sum2 += Vpn1111ME2
                    Vpn1111Sum3 += Vpn1111ME3
                    Vpn1111Sum4 += Vpn1111ME4
                    Vpn0220Sum += Vpn0220ME

                    VnnSum1 += VnnME1
                    VnnSum2 += VnnME2
                    VnnSum3 += VnnME3

                end

                @views V_NN_t1.pp[1][P,J+1][Bra,Ket] = VppSum1
                @views V_NN_t1.pp[2][P,J+1][Bra,Ket] = VppSum2
                @views V_NN_t1.pp[3][P,J+1][Bra,Ket] = VppSum3

                @views V_NN_t1.pn2002[P,J+1][Bra,Ket] = Vpn2002Sum
                @views V_NN_t1.pn1111[1][P,J+1][Bra,Ket] = Vpn1111Sum1
                @views V_NN_t1.pn1111[2][P,J+1][Bra,Ket] = Vpn1111Sum2
                @views V_NN_t1.pn1111[3][P,J+1][Bra,Ket] = Vpn1111Sum3
                @views V_NN_t1.pn1111[4][P,J+1][Bra,Ket] = Vpn1111Sum4
                @views V_NN_t1.pn0220[P,J+1][Bra,Ket] = Vpn0220Sum

                @views V_NN_t1.nn[1][P,J+1][Bra,Ket] = VnnSum1
                @views V_NN_t1.nn[2][P,J+1][Bra,Ket] = VnnSum2
                @views V_NN_t1.nn[3][P,J+1][Bra,Ket] = VnnSum3

            end
        end
    end

    return V_NN_t1
end

function qpH2b_allocate_H22_ind4(Params::Parameters,JP::Vector{Vector{Int64}},Orb::Vector{NOrb},Orb_NN::NNOrb,Orb_NN_t::NNOrb_Temp,H_NN::qpH2B,V_NN_t1::V2B_H22_Temp,V_NN_t2::V2B_H22_Temp,U::pnMatrix,V::pnMatrix)
    # Parameter initialization...
    N_max = Params.Calc.Nmax
    N_2max = 2*N_max
    J_max = 2*N_max + 1
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Define function for the transformation MEs for the 4th index ...
    @inline function V2b_temp_transformation_ind4_MEs(P::Int64,J::Int64,Bra::Int64,Ket::Int64,i::Int64,d::Int64,T::Matrix{Float64},V_NN_t::Matrix{Matrix{Float64}})
        return @views T[i,d] * V_NN_t[P,J+1][Bra,Ket]
    end

    # Transformation of the 4th index ...
    println("\nPerforming Quasiparticle Transformation of the residual 2-body NN interaction H^(22) in the 4th index...")
        # Preallocate V_NN_t2 ... for the Pandya transformation ...
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

                    VppSum = 0.0
                    Vpn2002Sum, Vpn1111Sum, Vpn0220Sum = 0.0, 0.0, 0.0
                    VnnSum = 0.0

                    @inbounds for i in Orb_x
                        Ket_ci = V2b_temp_index(c,i,J,P,Orb_NN_t)

                        VppME = V2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,d,V.p,V_NN_t1.pp[2])

                        Vpn2002ME = V2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,d,U.n,V_NN_t1.pn2002)
                        Vpn1111ME2 = V2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,d,V.n,V_NN_t1.pn1111[2])
                        Vpn1111ME3 = V2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,d,U.n,V_NN_t1.pn1111[3])
                        Vpn0220ME = V2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,d,V.n,V_NN_t1.pn0220)

                        VnnME = V2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,d,V.n,V_NN_t1.nn[2])

                        VppSum += VppME

                        Vpn2002Sum += Vpn2002ME

                        Vpn1111Sum += (Vpn1111ME2 + Vpn1111ME3)

                        Vpn0220Sum += Vpn0220ME

                        VnnSum += VnnME

                    end

                    @views V_NN_t2.pp[2][P,J+1][Bra,Ket] = VppSum

                    @views V_NN_t2.pn2002[P,J+1][Bra,Ket] = Vpn2002Sum

                    @views V_NN_t2.pn1111[2][P,J+1][Bra,Ket] = Vpn1111Sum

                    @views V_NN_t2.pn0220[P,J+1][Bra,Ket] = Vpn0220Sum

                    @views V_NN_t2.nn[2][P,J+1][Bra,Ket] = VnnSum

                end
            end
        end

        # Allocate H^(22) ...
        println("\tAllocating the components of H^(22) ...")
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

                    # Case of pn interaction ... T = 0
                    if Bra <= N_T0
                        Ind = Bra + (Ket - 1) * N_T0 - div(Ket * (Ket - 1),2)
                        a, b = Orb_NN.Ind[1,P,J+1][Bra][1], Orb_NN.Ind[1,P,J+1][Bra][2]
                        c, d = Orb_NN.Ind[1,P,J+1][Ket][1], Orb_NN.Ind[1,P,J+1][Ket][2]
                        j_a, j_b, j_c, j_d = Orb[a].j, Orb[b].j, Orb[c].j, Orb[d].j
                        l_a, l_b, l_c, l_d = Orb[a].l, Orb[b].l, Orb[c].l, Orb[d].l

                        Bra_ab = V2b_temp_index(a,b,J,P,Orb_NN_t)

                        Orb_x = Orb_PreComp(a_max,j_d,l_d,Orb)

                        Hpn2002Sum, Hpn1111Sum, Hpn0220Sum = 0.0, 0.0, 0.0

                        @inbounds for i in Orb_x
                            Ket_ci = V2b_temp_index(c,i,J,P,Orb_NN_t)

                            Vpn1111ME1 = V2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,d,U.n,V_NN_t1.pn1111[1])
                            Vpn1111ME2 = V2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,d,V.n,V_NN_t1.pn1111[4])

                            Hpn1111Sum += (Vpn1111ME1 + Vpn1111ME2)

                        end

                        Amp11 = Float64((-1)^(J))
                        Hpn1111Sum = Amp11 * Hpn1111Sum

                        @inbounds for I in 0:J_max

                            # pnH2002 && pnH0220 ...
                            if abs(j_a - j_b) <= 2*I && 2*I <= (j_a + j_b) && abs(j_c - j_d) <= 2*I && 2*I <= (j_c + j_d)

                                Amp2002 = Float64((-1)^(J+div(j_b + j_c,2)) * (2*I + 1)) * f6j(j_a,j_c,2*J,j_d,j_b,2*I)
                                Amp0220 = Float64((-1)^(J+div(j_a + j_d,2)) * (2*I + 1)) * f6j(j_d,j_b,2*J,j_a,j_c,2*I)
            
                                Bra_ab, Ket_cd = V2b_temp_index(a,b,I,P,Orb_NN_t), V2b_temp_index(c,d,I,P,Orb_NN_t)

                                Hpn2002Sum += Amp2002 * V_NN_t2.pn2002[P,I+1][Bra_ab,Ket_cd]
                                Hpn0220Sum += Amp0220 * V_NN_t2.pn0220[P,I+1][Bra_ab,Ket_cd]
                            end

                            # pnH1111 ...
                            if abs(j_a - j_d) <= 2*I && 2*I <= (j_a + j_d) && abs(j_c - j_b) <= 2*I &&
                            2*I <= (j_c + j_b) && (rem(l_a + l_d,2) + 1 == P) && (rem(l_b + l_c,2) + 1 == P)

                                Amp1111 = Float64((-1)^(J) * (2*I + 1)) * f6j(j_a,j_b,2*J,j_c,j_d,2*I)

                                Bra_ad, Ket_cb = V2b_temp_index(a,d,I,P,Orb_NN_t), V2b_temp_index(c,b,I,P,Orb_NN_t)

                                Hpn1111Sum += Amp1111 * V_NN_t2.pn1111[2][P,I+1][Bra_ad,Ket_cb]
                            end

                        end

                        @views H_NN.H22.pn2002[P,J+1][Ind] = Hpn2002Sum
                        @views H_NN.H22.pn1111[P,J+1][Ind] = Hpn1111Sum
                        @views H_NN.H22.pn0220[P,J+1][Ind] = Hpn0220Sum

                    end

                    # Case of pp & nn interaction ... T = 1
                    if Bra <= N_T1
                        Ind = Bra + (Ket - 1) * N_T1 - div(Ket * (Ket - 1),2)
                        a, b = Orb_NN.Ind[2,P,J+1][Bra][1], Orb_NN.Ind[2,P,J+1][Bra][2]
                        c, d = Orb_NN.Ind[2,P,J+1][Ket][1], Orb_NN.Ind[2,P,J+1][Ket][2]
                        j_a, j_b, j_c, j_d = Orb[a].j, Orb[b].j, Orb[c].j, Orb[d].j
                        l_a, l_b, l_c, l_d = Orb[a].l, Orb[b].l, Orb[c].l, Orb[d].l

                        Bra_ab = V2b_temp_index(a,b,J,P,Orb_NN_t)

                        Orb_x = Orb_PreComp(a_max,j_d,l_d,Orb)

                        HppSum, HnnSum = 0.0, 0.0

                        @inbounds for i in Orb_x
                            Ket_ci = V2b_temp_index(c,i,J,P,Orb_NN_t)

                            VppME1 = V2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,d,U.p,V_NN_t1.pp[1])
                            VppME3 = V2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,d,V.p,V_NN_t1.pp[3])

                            VnnME1 = V2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,d,U.n,V_NN_t1.nn[1])
                            VnnME3 = V2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,d,V.n,V_NN_t1.nn[3])

                            HppSum += (VppME1 + VppME3)
                            HnnSum += (VnnME1 + VnnME3)

                        end

                        Amp = Float64((-1)^(J + 1))
                        HppSum = Amp * HppSum
                        HnnSum = Amp * HnnSum

                        # NNH22 ...
                        @inbounds for I in 0:J_max
                            if abs(j_a - j_d) <= 2*I && 2*I <= (j_a + j_d) && abs(j_c - j_b) <= 2*I &&
                            2*I <= (j_c + j_b) && (rem(l_a + l_d,2) + 1 == P) && (rem(l_b + l_c,2) + 1 == P)

                                AmpR = 4.0 * Float64((-1)^(J) * (2*I + 1)) * f6j(j_a,j_b,2*J,j_c,j_d,2*I)

                                Bra_ad, Ket_cb = V2b_temp_index(a,d,I,P,Orb_NN_t), V2b_temp_index(c,b,I,P,Orb_NN_t)

                                HppSum += AmpR * V_NN_t2.pp[2][P,I+1][Bra_ad,Ket_cb]
                                HnnSum += AmpR * V_NN_t2.nn[2][P,I+1][Bra_ad,Ket_cb]
                            end
                        end

                        @views H_NN.H22.pp[P,J+1][Ind] = HppSum
                        @views H_NN.H22.nn[P,J+1][Ind] = HnnSum

                    end

                end
            end

        end

    return H_NN
end

@inline function qpH2b22pp(a::Int64,b::Int64,c::Int64,d::Int64,J::Int64,t::Int64,H_NN::qpH2B,Orb::Vector{NOrb},Orb_NN::NNOrb)
    if a >= b && c >= d
        @inbounds Ind = V2B_Index(a,b,c,d,J,P,1,Orb_NN)
        @views h = H_NN.qpH22.pp[P,J+1][Ind]
        return h
    elseif a < b && c >= d
        Amp = Float64((-1)^(J + 1 + div(Orb[a].j + Orb[b].j,2)))
        @inbounds Ind = V2B_Index(b,a,c,d,J,P,1,Orb_NN)
        @views h = Amp * H_NN.qpH22.pp[P,J+1][Ind]
        return h
    elseif a >= b && c < d
        Amp = Float64((-1)^(J + 1 + div(Orb[c].j + Orb[d].j,2)))
        @inbounds Ind = V2B_Index(a,b,d,c,J,P,1,Orb_NN)
        @views h = Amp * H_NN.qpH22.pp[P,J+1][Ind]
        return h
    elseif a < b && c < d
        Amp = Float64((-1)^(div(Orb[a].j + Orb[b].j + Orb[c].j + Orb[d].j,2)))
        @inbounds Ind = V2B_Index(b,a,d,c,J,P,1,Orb_NN)
        @views h = Amp * H_NN.qpH22.pp[P,J+1][Ind]
        return h
    end
end

@inline function qpH2b2002pn(a::Int64,b::Int64,c::Int64,d::Int64,J::Int64,t::Int64,H_NN::qpH2B,Orb::Vector{NOrb},Orb_NN::NNOrb)
    @inbounds Ind = V2B_Index(a,b,c,d,J,P,0,Orb_NN)
    @views h = H_NN.qpH22.pn2002[P,J+1][Ind]
    return h
end

@inline function qpH2b1111pn(a::Int64,b::Int64,c::Int64,d::Int64,J::Int64,t::Int64,H_NN::qpH2B,Orb::Vector{NOrb},Orb_NN::NNOrb)
    @inbounds Ind = V2B_Index(a,b,c,d,J,P,0,Orb_NN)
    @views h = H_NN.qpH22.pn1111[P,J+1][Ind]
    return h
end

@inline function qpH2b0220pn(a::Int64,b::Int64,c::Int64,d::Int64,J::Int64,t::Int64,H_NN::qpH2B,Orb::Vector{NOrb},Orb_NN::NNOrb)
    @inbounds Ind = V2B_Index(a,b,c,d,J,P,0,Orb_NN)
    @views h = H_NN.qpH22.pn0220[P,J+1][Ind]
    return h
end

@inline function qpH2b22nn(a::Int64,b::Int64,c::Int64,d::Int64,J::Int64,t::Int64,H_NN::qpH2B,Orb::Vector{NOrb},Orb_NN::NNOrb)
    if a >= b && c >= d
        @inbounds Ind = V2B_Index(a,b,c,d,J,P,1,Orb_NN)
        @views h = H_NN.qpH22.nn[P,J+1][Ind]
        return h
    elseif a < b && c >= d
        Amp = Float64((-1)^(J + 1 + div(Orb[a].j + Orb[b].j,2)))
        @inbounds Ind = V2B_Index(b,a,c,d,J,P,1,Orb_NN)
        @views h = Amp * H_NN.qpH22.nn[P,J+1][Ind]
        return h
    elseif a >= b && c < d
        Amp = Float64((-1)^(J + 1 + div(Orb[c].j + Orb[d].j,2)))
        @inbounds Ind = V2B_Index(a,b,d,c,J,P,1,Orb_NN)
        @views h = Amp * H_NN.qpH22.nn[P,J+1][Ind]
        return h
    elseif a < b && c < d
        Amp = Float64((-1)^(div(Orb[a].j + Orb[b].j + Orb[c].j + Orb[d].j,2)))
        @inbounds Ind = V2B_Index(b,a,d,c,J,P,1,Orb_NN)
        @views h = Amp * H_NN.qpH22.nn[P,J+1][Ind]
        return h
    end
end


# TO BE REMOVED !!!!
function qpH2b_allocate_H22_ind4_BackUp(Params::Parameters,JP::Vector{Vector{Int64}},Orb::Vector{NOrb},Orb_NN::NNOrb,Orb_NN_t::NNOrb_Temp,H_NN::qpH2B,V_NN_t1::V2B_H22_Temp,U::pnMatrix,V::pnMatrix)
    # Parameter initialization...
    N_max = Params.Calc.Nmax
    N_2max = 2*N_max
    J_max = 2*N_max + 1
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Define function for the transformation MEs for the 4th index ...
    @inline function V2b_temp_transformation_ind4_MEs(P::Int64,J::Int64,Bra::Int64,Ket::Int64,i::Int64,d::Int64,T::Matrix{Float64},V_NN_t::Matrix{Matrix{Float64}})
        return @views T[i,d] * V_NN_t[P,J+1][Bra,Ket]
    end

    # Transformation of the 4th index ...
    println("\nPerforming Quasiparticle Transformation of the residual 2-body NN interaction H22 in the 4th index...")

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

                    Bra_ab = V2b_temp_index(a,b,J,P,Orb_NN_t)

                    Orb_x = Orb_PreComp(a_max,j_d,l_d,Orb)

                    Hpn2002Sum, Hpn1111Sum, Hpn0220Sum = 0.0, 0.0, 0.0

                    @inbounds for i in Orb_x
                        Ket_ci = V2b_temp_index(c,i,J,P,Orb_NN_t)

                        Vpn2002ME = V2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,c,U.n,V_NN_t1.pn2002)
                        Vpn1111ME1 = V2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,c,U.n,V_NN_t1.pn1111[1])
                        Vpn1111ME2 = V2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,c,V.n,V_NN_t1.pn1111[2])
                        Vpn1111ME3 = V2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,c,U.n,V_NN_t1.pn1111[3])
                        Vpn1111ME4 = V2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,c,V.n,V_NN_t1.pn1111[4])
                        Vpn0220ME = V2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,c,V.n,V_NN_t1.pn0220)

                        Hpn2002Sum += Vpn2002ME
                        Hpn1111Sum += (Vpn1111ME1 + Vpn1111ME2 + Vpn1111ME3 + Vpn1111ME4)
                        Hpn0220Sum += Vpn0220ME

                    end

                    Amp20 = Float64((-1)^J)

                    Hpn2002Sum = Amp20 * Hpn2002Sum
                    Hpn0220Sum = Amp20 * Hpn0220Sum

                    @views H_NN.H22.pn2002[P,J+1][Ind] = Hpn2002Sum
                    @views H_NN.H22.pn1111[P,J+1][Ind] = Hpn1111Sum
                    @views H_NN.H22.pn0220[P,J+1][Ind] = Hpn0220Sum

                end

                # Case of pp & nn interaction ... T = 1
                if Bra <= N_T1
                    Ind = Bra + (Ket - 1) * N_T1 - div(Ket * (Ket - 1),2)
                    a, b = Orb_NN.Ind[2,P,J+1][Bra][1], Orb_NN.Ind[2,P,J+1][Bra][2]
                    c, d = Orb_NN.Ind[2,P,J+1][Ket][1], Orb_NN.Ind[2,P,J+1][Ket][2]
                    l_d, j_d = Orb[d].l, Orb[d].j

                    Bra_ab = V2b_temp_index(a,b,J,P,Orb_NN_t)

                    Orb_x = Orb_PreComp(a_max,j_d,l_d,Orb)

                    HppSum20, HppSum11 = 0.0, 0.0
                    HnnSum20, HnnSum11 = 0.0, 0.0

                    @inbounds for i in Orb_x
                        Ket_ci = V2b_temp_index(c,i,J,P,Orb_NN_t)

                        VppME1 = V2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,c,U.p,V_NN_t1.pp[1])
                        VppME2 = V2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,c,U.p,V_NN_t1.pp[2])
                        VppME3 = V2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,c,V.p,V_NN_t1.pp[3])

                        VnnME1 = V2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,c,U.n,V_NN_t1.nn[1])
                        VnnME2 = V2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,c,U.n,V_NN_t1.nn[2])
                        VnnME3 = V2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,c,V.n,V_NN_t1.nn[3])

                        HppSum20 += (VppME1 + VppME3)
                        HnnSum20 += (VnnME1 + VnnME3)

                        HppSum11 += 4.0 * VppME2
                        HnnSum11 += 4.0 * VnnME2

                    end

                    Amp11 = 0.0
                    @inbounds for I in 0:J_max
                        Amp11 += 0.0 
                    end

                    @views H_NN.H22.pp[P,J+1][Ind] = HppSum
                    @views H_NN.H22.nn[P,J+1][Ind] = HnnSum

                end

            end
        end

    end

    @inbounds Threads.@threads for i in JP
        J, P = i[1], i[2]
        if P == 1
            println("\tCalculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = +")
        else
            println("\tCalculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = -")
        end
        N = Orb_NN.N[2,P,J+1]

        @inbounds for Bra in 1:N
            @inbounds for Ket in 1:Bra
                Ind = Bra + (Ket - 1) * N_T1 - div(Ket * (Ket - 1),2)
                a, b = Orb_NN.Ind[2,P,J+1][Bra][1], Orb_NN.Ind[2,P,J+1][Bra][2]
                c, d = Orb_NN.Ind[2,P,J+1][Ket][1], Orb_NN.Ind[2,P,J+1][Ket][2]
                l_d, j_d = Orb[d].l, Orb[d].j

                Bra_ab = V2b_temp_index(a,b,J,P,Orb_NN_t)

                Orb_x = Orb_PreComp(a_max,j_d,l_d,Orb)

                HppSum20, HppSum11 = 0.0, 0.0
                HnnSum20, HnnSum11 = 0.0, 0.0

                @inbounds for i in Orb_x
                    Ket_ci = V2b_temp_index(c,i,J,P,Orb_NN_t)

                    VppME1 = V2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,c,U.p,V_NN_t1.pp[1])
                    VppME2 = V2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,c,U.p,V_NN_t1.pp[2])
                    VppME3 = V2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,c,V.p,V_NN_t1.pp[3])

                    VnnME1 = V2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,c,U.n,V_NN_t1.nn[1])
                    VnnME2 = V2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,c,U.n,V_NN_t1.nn[2])
                    VnnME3 = V2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,c,V.n,V_NN_t1.nn[3])

                    HppSum20 += (VppME1 + VppME3)
                    HnnSum20 += (VnnME1 + VnnME3)

                    HppSum11 += 4.0 * VppME2
                    HnnSum11 += 4.0 * VnnME2

                end

                Amp11 = 0.0
                @inbounds for I in 0:J_max
                    Amp11 += 0.0 
                end

                @views H_NN.H22.pp[P,J+1][Ind] = HppSum
                @views H_NN.H22.nn[P,J+1][Ind] = HnnSum

            end
        end

    end

    return H_NN
end
function qpH2b_allocate_H31(Params::Parameters,Orb::Vector{NOrb},Orb_NN::NNOrb,H_NN::qpH2B,V_NN::NNInt,U::pnMatrix,V::pnMatrix)
    println("\nAllocating H^(31) components of the quasiparticle Hamiltonian ...")

    # Preallocate arrays of allowed values of J & P ...
    JP = JP_initialize(Params.Calc.N2max + 1)

    # Initialize temporary V & W interaction arrays ...
    V_NN_t1, V_NN_t2, Orb_NN_t = qpH2b_initialize_H31(Params,Orb)

    # Transformation of the 1st index ...
    @time V_NN_t1 = qpH2b_allocate_H31_ind1(Params,JP,Orb,Orb_NN,Orb_NN_t,V_NN,V_NN_t1,U,V)

    # Transformation of the 2nd index ...
    @time V_NN_t2 = qpH2b_allocate_H31_ind2(Params,JP,Orb,Orb_NN_t,V_NN_t1,V_NN_t2,U,V)

    # Transformation of the 3rd index ...
    @time V_NN_t1 = qpH2b_allocate_H31_ind3(Params,JP,Orb,Orb_NN_t,V_NN_t1,V_NN_t2,U,V)

    # Transformation of the 4th index ...
    @time H_NN = qpH2b_allocate_H31_ind4(Params,JP,Orb,Orb_NN,Orb_NN_t,H_NN,V_NN_t1,U,V)

    # Deallocate V_NN_t & W_NN_t ...
    V_NN_t1, W_NN_t1 = nothing, nothing
    V_NN_t2, W_NN_t2 = nothing, nothing

    # Perform the Garbage Collection ...
    GC.gc()

    println("\nAllocation of H^(31) components of the quasiparticle Hamiltonian done ...")

    return H_NN
end

function qpH2b_initialize_H31(Params::Parameters,Orb::Vector{NOrb})
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
    V_pp_t1 = Vector{Matrix{Matrix{Float64}}}(undef,2)
    V_pn2011_t1 = Vector{Matrix{Matrix{Float64}}}(undef,2)
    V_pn1120_t1 = Vector{Matrix{Matrix{Float64}}}(undef,2)
    V_nn_t1 = Vector{Matrix{Matrix{Float64}}}(undef,2)
        # Case of V_t2 ...
    V_pp_t2 = Vector{Matrix{Matrix{Float64}}}(undef,2)
    V_pn2011_t2 = Vector{Matrix{Matrix{Float64}}}(undef,2)
    V_pn1120_t2 = Vector{Matrix{Matrix{Float64}}}(undef,2)
    V_nn_t2 = Vector{Matrix{Matrix{Float64}}}(undef,2)

    # Initialite sub-entries of V_t ...
    @inbounds for i in 1:2
        V_pp_t1[i] = Matrix{Matrix{Float64}}(undef,2,J_max+1)
        V_pn2011_t1[i] = Matrix{Matrix{Float64}}(undef,2,J_max+1)
        V_pn1120_t1[i] = Matrix{Matrix{Float64}}(undef,2,J_max+1)
        V_nn_t1[i] = Matrix{Matrix{Float64}}(undef,2,J_max+1)

        V_pp_t2[i] = Matrix{Matrix{Float64}}(undef,2,J_max+1)
        V_pn2011_t2[i] = Matrix{Matrix{Float64}}(undef,2,J_max+1)
        V_pn1120_t2[i] = Matrix{Matrix{Float64}}(undef,2,J_max+1)
        V_nn_t2[i] = Matrix{Matrix{Float64}}(undef,2,J_max+1)
    end
 
    # Allocate entries of V & W ...
    @inbounds for P = 1:2
        @inbounds Threads.@threads for J=0:J_max
            Length = N_Orb_NN[P,J+1]

            @inbounds for i in 1:2
                V_pp_t1[i][P,J+1] = zeros(Float64,Length,Length)
                V_pn2011_t1[i][P,J+1] = zeros(Float64,Length,Length)
                V_pn1120_t1[i][P,J+1] = zeros(Float64,Length,Length)
                V_nn_t1[i][P,J+1] = zeros(Float64,Length,Length)

                V_pp_t2[i][P,J+1] = zeros(Float64,Length,Length)
                V_pn2011_t2[i][P,J+1] = zeros(Float64,Length,Length)
                V_pn1120_t2[i][P,J+1] = zeros(Float64,Length,Length)
                V_nn_t2[i][P,J+1] = zeros(Float64,Length,Length)
            end

        end

    end

    # Allocate the NN interaction in quasiparticle picture ...
    V_NN_t1 = V2B_H31_Temp(V_pp_t1,V_pn2011_t1,V_pn1120_t1,V_nn_t1)
    V_NN_t2 = V2B_H31_Temp(V_pp_t2,V_pn2011_t2,V_pn1120_t2,V_nn_t2)

    # Allocate the array for NN orbitals ...
    Orb_NN_t = NNOrb_Temp(Orb_NN_Dic,N_Orb_NN,Ind_Orb_NN)

    return V_NN_t1, V_NN_t2, Orb_NN_t
end

function qpH2b_allocate_H31_ind1(Params::Parameters,JP::Vector{Vector{Int64}},Orb::Vector{NOrb},Orb_NN::NNOrb,Orb_NN_t::NNOrb_Temp,V_NN::NNInt,V_NN_t1::V2B_H31_Temp,U::pnMatrix,V::pnMatrix)
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
    println("\nPerforming Quasiparticle Transformation of the residual 2-body NN interaction H31 in the 1st index...")
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

            Orb_x = Orb_PreComp(a_max,j_a,l_a,Orb)
            @inbounds for Ket in 1:N
                c = Orb_NN_t.Ind[P,J+1][Ket][1]
                d = Orb_NN_t.Ind[P,J+1][Ket][2]

                VppSum1, VppSum2 = 0.0, 0.0
                Vpn2011Sum1, Vpn2011Sum2, Vpn1120Sum1, Vpn1120Sum2 = 0.0, 0.0, 0.0, 0.0
                VnnSum1, VnnSum2 = 0.0, 0.0
                
                @inbounds for i in Orb_x
                    VppUME, VppVME, VpnUME, VpnVME, VnnUME, VnnVME = V2b_temp_transformation_ind1_MEs(J,i,a,b,c,d,U,V,V_NN,Orb,Orb_NN)

                    VppSum1 += VppUME
                    VppSum2 += VppVME

                    Vpn2011Sum1 += VpnUME
                    Vpn2011Sum2 += VpnUME
                    Vpn1120Sum1 += VpnUME
                    Vpn1120Sum2 += VpnVME

                    VnnSum1 += VnnUME
                    VnnSum2 += VnnVME

                end

                @views V_NN_t1.pp[1][P,J+1][Bra,Ket] = VppSum1
                @views V_NN_t1.pp[2][P,J+1][Bra,Ket] = VppSum2

                @views V_NN_t1.pn2011[1][P,J+1][Bra,Ket] = Vpn2011Sum1
                @views V_NN_t1.pn2011[2][P,J+1][Bra,Ket] = Vpn2011Sum2
                @views V_NN_t1.pn1120[1][P,J+1][Bra,Ket] = Vpn1120Sum1
                @views V_NN_t1.pn1120[2][P,J+1][Bra,Ket] = Vpn1120Sum2

                @views V_NN_t1.nn[1][P,J+1][Bra,Ket] = VnnSum1
                @views V_NN_t1.nn[2][P,J+1][Bra,Ket] = VnnSum2

            end
        end
    end

    return V_NN_t1
end

function qpH2b_allocate_H31_ind2(Params::Parameters,JP::Vector{Vector{Int64}},Orb::Vector{NOrb},Orb_NN_t::NNOrb_Temp,V_NN_t1::V2B_H31_Temp,V_NN_t2::V2B_H31_Temp,U::pnMatrix,V::pnMatrix)
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
    println("\nPerforming Quasiparticle Transformation of the residual 2-body NN interaction H31 in the 2nd index...")

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

                VppSum1, VppSum2 = 0.0, 0.0
                Vpn2011Sum1, Vpn2011Sum2, Vpn1120Sum1, Vpn1120Sum2 = 0.0, 0.0, 0.0, 0.0
                VnnSum1, VnnSum2 = 0.0, 0.0
                
                @inbounds for i in Orb_x
                    Bra_ai = V2b_temp_index(a,i,J,P,Orb_NN_t)

                    VppME1 = V2b_temp_transformation_ind2_MEs(P,J,Bra_ai,Ket,i,b,V.p,V_NN_t1.pp[1])
                    VppME2 = V2b_temp_transformation_ind2_MEs(P,J,Bra_ai,Ket,i,b,U.p,V_NN_t1.pp[2])

                    Vpn2011ME1 = V2b_temp_transformation_ind2_MEs(P,J,Bra_ai,Ket,i,b,U.n,V_NN_t1.pn2011[1])
                    Vpn2011ME2 = V2b_temp_transformation_ind2_MEs(P,J,Bra_ai,Ket,i,b,V.n,V_NN_t1.pn2011[2])
                    Vpn1120ME1 = V2b_temp_transformation_ind2_MEs(P,J,Bra_ai,Ket,i,b,U.n,V_NN_t1.pn1120[1])
                    Vpn1120ME2 = V2b_temp_transformation_ind2_MEs(P,J,Bra_ai,Ket,i,b,U.n,V_NN_t1.pn1120[2])

                    VnnME1 = V2b_temp_transformation_ind2_MEs(P,J,Bra_ai,Ket,i,b,V.n,V_NN_t1.nn[1])
                    VnnME2 = V2b_temp_transformation_ind2_MEs(P,J,Bra_ai,Ket,i,b,U.n,V_NN_t1.nn[2])

                    VppSum1 += VppME1
                    VppSum2 += VppME2

                    Vpn2011Sum1 += Vpn2011ME1
                    Vpn2011Sum2 += Vpn2011ME2
                    Vpn1120Sum1 += Vpn1120ME1
                    Vpn1120Sum2 += Vpn1120ME2

                    VnnSum1 += VnnME1
                    VnnSum2 += VnnME2

                end

                @views V_NN_t2.pp[1][P,J+1][Bra,Ket] = VppSum1
                @views V_NN_t2.pp[2][P,J+1][Bra,Ket] = VppSum2

                @views V_NN_t2.pn2011[1][P,J+1][Bra,Ket] = Vpn2011Sum1
                @views V_NN_t2.pn2011[2][P,J+1][Bra,Ket] = Vpn2011Sum2
                @views V_NN_t2.pn1120[1][P,J+1][Bra,Ket] = Vpn1120Sum1
                @views V_NN_t2.pn1120[2][P,J+1][Bra,Ket] = Vpn1120Sum2

                @views V_NN_t2.nn[1][P,J+1][Bra,Ket] = VnnSum1
                @views V_NN_t2.nn[2][P,J+1][Bra,Ket] = VnnSum2

            end
        end
    end

    return V_NN_t2
end

function qpH2b_allocate_H31_ind3(Params::Parameters,JP::Vector{Vector{Int64}},Orb::Vector{NOrb},Orb_NN_t::NNOrb_Temp,V_NN_t1::V2B_H31_Temp,V_NN_t2::V2B_H31_Temp,U::pnMatrix,V::pnMatrix)
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
    println("\nPerforming Quasiparticle Transformation of the residual 2-body NN interaction H31 in the 3rd index...")

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

                VppSum1, VppSum2 = 0.0, 0.0
                Vpn2011Sum1, Vpn2011Sum2, Vpn1120Sum1, Vpn1120Sum2 = 0.0, 0.0, 0.0, 0.0
                VnnSum1, VnnSum2 = 0.0, 0.0

                @inbounds for i in Orb_x
                    Ket_id = V2b_temp_index(i,d,J,P,Orb_NN_t)

                    VppME1 = V2b_temp_transformation_ind3_MEs(P,J,Bra,Ket_id,i,c,V.p,V_NN_t2.pp[1])
                    VppME2 = V2b_temp_transformation_ind3_MEs(P,J,Bra,Ket_id,i,c,U.p,V_NN_t2.pp[2])

                    Vpn2011ME1 = V2b_temp_transformation_ind3_MEs(P,J,Bra,Ket_id,i,c,V.p,V_NN_t2.pn2011[1])
                    Vpn2011ME2 = V2b_temp_transformation_ind3_MEs(P,J,Bra,Ket_id,i,c,V.p,V_NN_t2.pn2011[2])
                    Vpn1120ME1 = V2b_temp_transformation_ind3_MEs(P,J,Bra,Ket_id,i,c,U.p,V_NN_t2.pn1120[1])
                    Vpn1120ME2 = V2b_temp_transformation_ind3_MEs(P,J,Bra,Ket_id,i,c,V.p,V_NN_t2.pn1120[2])

                    VnnME1 = V2b_temp_transformation_ind3_MEs(P,J,Bra,Ket_id,i,c,V.n,V_NN_t2.nn[1])
                    VnnME2 = V2b_temp_transformation_ind3_MEs(P,J,Bra,Ket_id,i,c,U.n,V_NN_t2.nn[2])

                    VppSum1 += VppME1
                    VppSum2 += VppME2

                    Vpn2011Sum1 += Vpn2011ME1
                    Vpn2011Sum2 += Vpn2011ME2
                    Vpn1120Sum1 += Vpn1120ME1
                    Vpn1120Sum2 += Vpn1120ME2

                    VnnSum1 += VnnME1
                    VnnSum2 += VnnME2

                end

                @views V_NN_t1.pp[1][P,J+1][Bra,Ket] = VppSum1
                @views V_NN_t1.pp[2][P,J+1][Bra,Ket] = VppSum2

                @views V_NN_t1.pn2011[1][P,J+1][Bra,Ket] = Vpn2011Sum1
                @views V_NN_t1.pn2011[2][P,J+1][Bra,Ket] = Vpn2011Sum2
                @views V_NN_t1.pn1120[1][P,J+1][Bra,Ket] = Vpn1120Sum1
                @views V_NN_t1.pn1120[2][P,J+1][Bra,Ket] = Vpn1120Sum2

                @views V_NN_t1.nn[1][P,J+1][Bra,Ket] = VnnSum1
                @views V_NN_t1.nn[2][P,J+1][Bra,Ket] = VnnSum2

            end
        end
    end

    return V_NN_t1
end

function qpH2b_allocate_H31_ind4(Params::Parameters,JP::Vector{Vector{Int64}},Orb::Vector{NOrb},Orb_NN::NNOrb,Orb_NN_t::NNOrb_Temp,H_NN::qpH2B,V_NN_t1::V2B_H31_Temp,U::pnMatrix,V::pnMatrix)
    # Parameter initialization...
    N_max = Params.Calc.Nmax
    N_2max = 2*N_max
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Define function for the transformation MEs for the 4th index ...
    @inline function V2b_temp_transformation_ind4_MEs(P::Int64,J::Int64,Bra::Int64,Ket::Int64,i::Int64,d::Int64,T::Matrix{Float64},V_NN_t::Matrix{Matrix{Float64}})
        return @views T[i,d] * V_NN_t[P,J+1][Bra,Ket]
    end

    # Transformation of the 4th index ...
    println("\nPerforming Quasiparticle Transformation of the residual 2-body NN interaction H31 in the 4th index...")

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

                    Hpn2011Sum, Hpn1120Sum = 0.0, 0.0

                    @inbounds for i in Orb_x
                        Ket_ci = V2b_temp_index(c,i,J,P,Orb_NN_t)


                        Vpn2011ME1 = V2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,c,U.n,V_NN_t1.pn2011[1])
                        Vpn2011ME2 = V2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,c,V.n,V_NN_t1.pn2011[2])
                        Vpn1120ME1 = V2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,c,V.n,V_NN_t1.pn1120[1])
                        Vpn1120ME2 = V2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,c,V.n,V_NN_t1.pn1120[2])

                        Hpn2011Sum += (Vpn2011ME1 + Vpn2011ME2)
                        Hpn1120Sum += (Vpn1120ME1 + Vpn1120ME2)

                    end

                    Amp = Float64((-1)^J)

                    Hpn2011Sum = Amp * Hpn2011Sum
                    Hpn1120Sum = Amp * Hpn1120Sum

                    @views H_NN.H31.pn2011[P,J+1][Ind] = Hpn2011Sum
                    @views H_NN.H31.pn1120[P,J+1][Ind] = Hpn1120Sum

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

                        VppME1 = V2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,c,V.p,V_NN_t1.pp[1])
                        VppME2 = V2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,c,V.p,V_NN_t1.pp[2])

                        VnnME1 = V2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,c,V.n,V_NN_t1.nn[1])
                        VnnME2 = V2b_temp_transformation_ind4_MEs(P,J,Bra,Ket_ci,i,c,V.n,V_NN_t1.nn[2])

                        HppSum += -0.5 * (VppME1 + VppME2)
                        HnnSum += -0.5 * (VnnME1 + VnnME2)

                    end

                    Amp = 0.5 * Float64((-1)^J)

                    HppSum = Amp * HppSum
                    HnnSum = Amp * HnnSum

                    @views H_NN.H31.pp[P,J+1][Ind] = HppSum
                    @views H_NN.H31.nn[P,J+1][Ind] = HnnSum

                end

            end
        end

    end

    return H_NN
end

@inline function qpH2b31pp(a::Int64,b::Int64,c::Int64,d::Int64,J::Int64,t::Int64,H_NN::qpH2B,Orb::Vector{NOrb},Orb_NN::NNOrb)
    if a >= b && c >= d
        @inbounds Ind = V2B_Index(a,b,c,d,J,P,1,Orb_NN)
        @views h = H_NN.qpH31.pp[P,J+1][Ind]
        return h
    elseif a < b && c >= d
        Amp = Float64((-1)^(J + 1 + div(Orb[a].j + Orb[b].j,2)))
        @inbounds Ind = V2B_Index(b,a,c,d,J,P,1,Orb_NN)
        @views h = Amp * H_NN.qpH31.pp[P,J+1][Ind]
        return h
    elseif a >= b && c < d
        Amp = Float64((-1)^(J + 1 + div(Orb[c].j + Orb[d].j,2)))
        @inbounds Ind = V2B_Index(a,b,d,c,J,P,1,Orb_NN)
        @views h = Amp * H_NN.qpH31.pp[P,J+1][Ind]
        return h
    elseif a < b && c < d
        Amp = Float64((-1)^(div(Orb[a].j + Orb[b].j + Orb[c].j + Orb[d].j,2)))
        @inbounds Ind = V2B_Index(b,a,d,c,J,P,1,Orb_NN)
        @views h = Amp * H_NN.qpH31.pp[P,J+1][Ind]
        return h
    end
end

@inline function qpH2b2011pn(a::Int64,b::Int64,c::Int64,d::Int64,J::Int64,t::Int64,H_NN::qpH2B,Orb::Vector{NOrb},Orb_NN::NNOrb)
    @inbounds Ind = V2B_Index(a,b,c,d,J,P,0,Orb_NN)
    @views h = H_NN.qpH31.pn2011[P,J+1][Ind]
    return h
end

@inline function qpH2b1120pn(a::Int64,b::Int64,c::Int64,d::Int64,J::Int64,t::Int64,H_NN::qpH2B,Orb::Vector{NOrb},Orb_NN::NNOrb)
    @inbounds Ind = V2B_Index(a,b,c,d,J,P,0,Orb_NN)
    @views h = H_NN.qpH31.pn1120[P,J+1][Ind]
    return h
end

@inline function qpH2b31nn(a::Int64,b::Int64,c::Int64,d::Int64,J::Int64,t::Int64,H_NN::qpH2B,Orb::Vector{NOrb},Orb_NN::NNOrb)
    if a >= b && c >= d
        @inbounds Ind = V2B_Index(a,b,c,d,J,P,1,Orb_NN)
        @views h = H_NN.qpH31.nn[P,J+1][Ind]
        return h
    elseif a < b && c >= d
        Amp = Float64((-1)^(J + 1 + div(Orb[a].j + Orb[b].j,2)))
        @inbounds Ind = V2B_Index(b,a,c,d,J,P,1,Orb_NN)
        @views h = Amp * H_NN.qpH31.nn[P,J+1][Ind]
        return h
    elseif a >= b && c < d
        Amp = Float64((-1)^(J + 1 + div(Orb[c].j + Orb[d].j,2)))
        @inbounds Ind = V2B_Index(a,b,d,c,J,P,1,Orb_NN)
        @views h = Amp * H_NN.qpH31.nn[P,J+1][Ind]
        return h
    elseif a < b && c < d
        Amp = Float64((-1)^(div(Orb[a].j + Orb[b].j + Orb[c].j + Orb[d].j,2)))
        @inbounds Ind = V2B_Index(b,a,d,c,J,P,1,Orb_NN)
        @views h = Amp * H_NN.qpH31.nn[P,J+1][Ind]
        return h
    end
end
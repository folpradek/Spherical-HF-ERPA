function qpH2b_canonical_transformation(Params::Parameters,Orb::Vector{NOrb},Orb_NN::NNOrb,V_NN::NNInt,C::pnMatrix)
    # Transformation of V_NN & W_NN to the canonical HFB basis...
    println("\nPerforming the Canonical Transformation of residual NN V & W interactions into the target canonical basis ...")

    # Preallocate arrays of allowed values of J & P ...
    JP = JP_initialize(Params.Calc.N2max + 1)

    # Transformation of the 1st index ...
    V_NN_t1, Orb_NN_t = qpH2b_canonical_transformation_ind1(Params,JP,Orb,Orb_NN,V_NN,C)

    # Transformation of the 2nd index ...
    V_NN_t2 = qpH2b_canonical_transformation_ind2(Params,JP,Orb,Orb_NN_t,V_NN_t1,C)

    # Transformation of the 3rd index ...
    V_NN_t1 = qpH2b_canonical_transformation_ind3(Params,JP,Orb,Orb_NN_t,V_NN_t1,V_NN_t2,C)

    # Transformation of the 4th index ...
    V_NN = qpH2b_canonical_transformation_ind4(Params,JP,Orb,Orb_NN,Orb_NN_t,V_NN,V_NN_t1,C)

    # Deallocate temporary interaction arrays ...
    V_NN_t1 = nothing
    V_NN_t2 = nothing

    # Perform the Garbage collection ...
    GC.gc()
    
    println("\nResidual NN 2-body interaction V & W ready ...\n")

    return V_NN
end

function qpH2b_canonical_transformation_ind1(Params::Parameters,JP::Vector{Vector{Int64}},Orb::Vector{NOrb},Orb_NN::NNOrb,V_NN::NNInt,C::pnMatrix)
    # Parameter initialization...
    N_max = Params.Calc.Nmax
    N_2max = 2*N_max
    a_max = div((N_max + 1)*(N_max + 2),2)
    J_max = N_2max + 1

    # Make temporary arrays for NN interaction V ...
    V_NN_t1, Orb_NN_t = V2b_temp_initialize(Orb,N_max,Make_Orb_NN = true)

    # Define function for the transformation MEs for the 1st index ...
    @inline function V2b_temp_transformation_ind1_MEs(J::Int64,i::Int64,a::Int64,b::Int64,c::Int64,d::Int64,C::pnMatrix,V_NN::NNInt,Orb::Vector{NOrb},Orb_NN::NNOrb)
        return @views C.p[i,a] * V2B(i,b,c,d,J,1,V_NN.pp,Orb,Orb_NN), C.p[i,a] * V2B(i,b,c,d,J,0,V_NN.pn,Orb,Orb_NN),
                      C.n[i,a] * V2B(i,b,c,d,J,1,V_NN.nn,Orb,Orb_NN)
    end

    # Transformation of the 1st index ...
    println("\nPerforming Canonical Transformation of the residual 2-body NN interaction in the 1st index...")
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

                VppSum, VpnSum, VnnSum = 0.0, 0.0, 0.0
                
                @inbounds for i in Orb_x
                    VppME, VpnME, VnnME = V2b_temp_transformation_ind1_MEs(J,i,a,b,c,d,C,V_NN,Orb,Orb_NN)
                    VppSum += VppME
                    VpnSum += VpnME
                    VnnSum += VnnME

                end

                @views V_NN_t1.pp[P,J+1][Bra,Ket] = VppSum
                @views V_NN_t1.pn[P,J+1][Bra,Ket] = VpnSum
                @views V_NN_t1.nn[P,J+1][Bra,Ket] = VnnSum

            end
        end
    end

    return V_NN_t1, Orb_NN_t
end

function qpH2b_canonical_transformation_ind2(Params::Parameters,JP::Vector{Vector{Int64}},Orb::Vector{NOrb},Orb_NN_t::NNOrb_Temp,V_NN_t1::V2B_Temp,C::pnMatrix)
    # Parameter initialization...
    N_max = Params.Calc.Nmax
    N_2max = 2*N_max
    a_max = div((N_max + 1)*(N_max + 2),2)
    J_max = N_2max + 1

    # Make another temporary arrays for NN interaction V ...
    V_NN_t2, Orb_NN_t2 = V2b_temp_initialize(Orb,N_max,Make_Orb_NN = true)

    # Define function for the transformation MEs for the 2nd index ...
    @inline function V2b_temp_transformation_ind2_MEs(P::Int64,J::Int64,Bra::Int64,Ket::Int64,i::Int64,b::Int64,C::pnMatrix,V_NN_t::V2B_Temp)
        return @views C.p[i,b] * V_NN_t.pp[P,J+1][Bra,Ket], C.n[i,b] * V_NN_t.pn[P,J+1][Bra,Ket], C.n[i,b] * V_NN_t.nn[P,J+1][Bra,Ket]
    end

    # Transformation of the 2nd index ...
    println("\nPerforming Canonical Transformation of the residual 2-body NN interaction in the 2nd index...")

    @time @inbounds Threads.@threads for i in JP
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

                VppSum, VpnSum, VnnSum = 0.0, 0.0, 0.0
                
                @inbounds for i in Orb_x
                    Bra_ai = V2b_temp_index(a,i,J,P,Orb_NN_t)
                    VppME, VpnME, VnnME = V2b_temp_transformation_ind2_MEs(P,J,Bra_ai,Ket,i,b,C,V_NN_t1)
                    VppSum += VppME
                    VpnSum += VpnME
                    VnnSum += VnnME

                end

                @views V_NN_t2.pp[P,J+1][Bra,Ket] = VppSum
                @views V_NN_t2.pn[P,J+1][Bra,Ket] = VpnSum
                @views V_NN_t2.nn[P,J+1][Bra,Ket] = VnnSum

            end
        end
    end

    return V_NN_t2
end

function qpH2b_canonical_transformation_ind3(Params::Parameters,JP::Vector{Vector{Int64}},Orb::Vector{NOrb},Orb_NN_t::NNOrb_Temp,V_NN_t1::V2B_Temp,V_NN_t2::V2B_Temp,C::pnMatrix)
    # Parameter initialization...
    N_max = Params.Calc.Nmax
    N_2max = 2*N_max
    a_max = div((N_max + 1)*(N_max + 2),2)
    J_max = N_2max + 1

    # Define function for the transformation MEs for the 3rd index ...
    @inline function V2b_temp_transformation_ind3_MEs(P::Int64,J::Int64,Bra::Int64,Ket::Int64,i::Int64,c::Int64,C::pnMatrix,V_NN_t::V2B_Temp)
        return @views C.p[i,c] * V_NN_t.pp[P,J+1][Bra,Ket], C.p[i,c] * V_NN_t.pn[P,J+1][Bra,Ket], C.n[i,c] * V_NN_t.nn[P,J+1][Bra,Ket]
    end

    # Transformation of the 3rd index ...
    println("\nPerforming Canonical Transformation of the residual 2-body NN interaction in the 3rd index...")

    @time @inbounds Threads.@threads for i in JP
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
                VppSum, VpnSum, VnnSum = 0.0, 0.0, 0.0

                @inbounds for i in Orb_x
                    Ket_id = V2b_temp_index(i,d,J,P,Orb_NN_t)

                    VppME, VpnME, VnnME = V2b_temp_transformation_ind3_MEs(P,J,Bra,Ket_id,i,c,C,V_NN_t2)

                    VppSum += VppME
                    VpnSum += VpnME
                    VnnSum += VnnME

                end

                @views V_NN_t1.pp[P,J+1][Bra,Ket] = VppSum
                @views V_NN_t1.pn[P,J+1][Bra,Ket] = VpnSum
                @views V_NN_t1.nn[P,J+1][Bra,Ket] = VnnSum

            end
        end
    end

    return V_NN_t1
end

function qpH2b_canonical_transformation_ind4(Params::Parameters,JP::Vector{Vector{Int64}},Orb::Vector{NOrb},Orb_NN::NNOrb,Orb_NN_t::NNOrb_Temp,V_NN::NNInt,V_NN_t1::V2B_Temp,C::pnMatrix)
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

    # Index 4
    println("\nPerforming Canonical Transformation of the residual 2-body NN interaction in the 4th index...")

    @time @inbounds Threads.@threads for i in JP
        J, P = i[1], i[2]
        if P == 1
            println("\tCalculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = +")
        else
            println("\tCalculating   ...   J = " * string(J) * "/" * string(N_2max+1) * "\tP = -")
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

                    @inbounds for i in Orb_x
                        Ket_ci = V2b_temp_index(c,i,J,P,Orb_NN_t)

                        VpnME = V2b_temp_transformation_ind4_T0_MEs(P,J,Bra_ab,Ket_ci,i,d,C,V_NN_t1)

                        VpnSum += VpnME

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

                    @inbounds for i in Orb_x
                        Ket_ci = V2b_temp_index(c,i,J,P,Orb_NN_t)

                        VppME, VnnME = V2b_temp_transformation_ind4_T1_MEs(P,J,Bra_ab,Ket_ci,i,d,C,V_NN_t1)

                        VppSum += VppME
                        VnnSum += VnnME

                    end

                    @views V_NN.pp[P,J+1][Ind] = VppSum
                    @views V_NN.nn[P,J+1][Ind] = VnnSum

                end

            end
        end

    end

    return V_NN
end
function QTDA_diagonalize(Params::Parameters,Orb::Vector{Orb1B},Orb_2qp::qpOrb2B,A::Matrix{Matrix{Float64}},qpN::qpO1B,qpTrOp::qpTr1B)
    # Read calculation parameters ...
    N_max = Params.Calc.Nmax
    J_max = 2*N_max + 1
    Orthogon = Params.Calc.QTDA.Ortho

    # Make the list of values of J & P for iteration ...
    JP_list = JP_initialize(J_max)

    # Initialite the storage for QTDA solutions ...
    E_QTDA = Matrix{Vector{Float64}}(undef,2,J_max+1)
    X_QTDA = Matrix{Matrix{Float64}}(undef,2,J_max+1)

    # Perform the diagonalization of the QTDA matrix A ...
    println("\nDiagonalizing the QTDA matrix A ...")

    @inbounds for JP in JP_list
        J, P = JP[1], JP[2]
        A_JP = A[P,J+1]

        if Orthogon == false || !((J == 0 && P == 1) || (J == 1 && P == 2))
            E_QTDA_JP, X_QTDA_JP = eigen(A_JP, sortby=+)
            Sort = sortperm(E_QTDA_JP, by = x -> real(x))
            E_QTDA_JP, X_QTDA_JP = E_QTDA_JP[Sort], @views X_QTDA_JP[:,Sort]

        elseif Orthogon == true && J == 0 && P == 1

            # Read the dimension of 0+ 2-qp space ...
            N_qp = Orb_2qp.N[1,1]

            # Initialize the vector for orthogonalization of spurious 0+ level ...
            Spur_PN = QTDA_spurious_PN_initialize(Orb,Orb_2qp,qpN)

            V, Spur_ind = QTDA_spurious_orthogonalize(N_qp,Spur_PN)
            U = zeros(Float64,N_qp,N_qp)

            @inbounds for a in 1:(Spur_ind-1)
                @views U[:,a] = V[:,a]
            end

            @inbounds for a in (Spur_ind+1):N_qp
                @views U[:,a-1] = V[:,a]
            end

            A_JP_ort = zeros(Float64,N_qp - 1, N_qp - 1)

            Pm = zeros(Float64,N_qp,N_qp-1)
            @inbounds for a in 1:N_qp
                @inbounds for b in 1:(N_qp-1)
                    @inbounds for c in 1:N_qp
                        Pm[a,b] = Pm[a,b] + A_JP[a,c] * U[c,b]
                    end
                end
            end

            @inbounds for a in 1:(N_qp-1)
                @inbounds for b in 1:(N_qp-1)
                    @inbounds for c in 1:N_qp
                        A_JP_ort[a,b] = A_JP_ort[a,b] + U[c,a] * Pm[c,b]
                    end
                end
            end

            E_QTDA_JP, X_QTDA_JP = eigen(A_JP_ort, sortby=+)
            Sort_QTDA = sortperm(E_QTDA_JP, by = x -> real(x))
            E_QTDA_JP, X_QTDA_JP = @views E_QTDA_JP[Sort_QTDA], X_QTDA_JP[:,Sort_QTDA]

            Pm = zeros(Float64,N_qp,N_qp-1)
            @inbounds for a in 1:N_qp
                @inbounds for b in 1:(N_qp -1)
                    @inbounds for c in 1:(N_qp-1)
                        ME = U[a,c] * X_QTDA_JP[c,b]
                        Pm[a,b] += ME
                    end
                end
            end

            E_QTDA_JP_ort = zeros(Float64,N_qp)
            @inbounds for a in 2:N_qp
                E_QTDA_JP_ort[a] = E_QTDA_JP[a-1]
            end

            X_QTDA_JP_ort = zeros(Float64,N_qp,N_qp)
            @views X_QTDA_JP_ort[:,1] = Spur_PN
            @inbounds for a in 2:N_qp
                @views X_QTDA_JP_ort[:,a] = Pm[:,a-1]

            end

            X_QTDA_JP = X_QTDA_JP_ort
            E_QTDA_JP = E_QTDA_JP_ort
            
        elseif Orthogon == true && J == 1 && P == 2

            # Read the dimension of 1- 2-qp space ...
            N_qp = Orb_2qp.N[2,2]

            # Initialize the vector for orthogonalization of spurious 1- level ...
            Spur_CM = QTDA_spurious_CM_initialize(Orb_2qp,qpTrOp)

            V, Spur_ind = QTDA_spurious_orthogonalize(N_qp,Spur_CM)
            U = zeros(Float64,N_qp,N_qp)

            @inbounds for a in 1:(Spur_ind-1)
                @views U[:,a] = V[:,a]
            end

            @inbounds for a in (Spur_ind+1):N_qp
                @views U[:,a-1] = V[:,a]
            end

            A_JP_ort = zeros(Float64,N_qp - 1, N_qp - 1)

            Pm = zeros(Float64,N_qp,N_qp-1)
            @inbounds for a in 1:N_qp
                @inbounds for b in 1:(N_qp-1)
                    @inbounds for c in 1:N_qp
                        Pm[a,b] = Pm[a,b] + A_JP[a,c] * U[c,b]
                    end
                end
            end

            @inbounds for a in 1:(N_qp-1)
                @inbounds for b in 1:(N_qp-1)
                    @inbounds for c in 1:N_qp
                        A_JP_ort[a,b] = A_JP_ort[a,b] + U[c,a] * Pm[c,b]
                    end
                end
            end

            E_QTDA_JP, X_QTDA_JP = eigen(A_JP_ort, sortby=+)
            Sort_QTDA = sortperm(E_QTDA_JP, by = x -> real(x))
            E_QTDA_JP, X_QTDA_JP = @views E_QTDA_JP[Sort_QTDA], X_QTDA_JP[:,Sort_QTDA]

            Pm = zeros(Float64,N_qp,N_qp-1)
            @inbounds for a in 1:N_qp
                @inbounds for b in 1:(N_qp -1)
                    @inbounds for c in 1:(N_qp-1)
                        ME = U[a,c] * X_QTDA_JP[c,b]
                        Pm[a,b] += ME
                    end
                end
            end

            E_QTDA_JP_ort = zeros(Float64,N_qp)
            @inbounds for a in 2:N_qp
                E_QTDA_JP_ort[a] = E_QTDA_JP[a-1]
            end

            X_QTDA_JP_ort = zeros(Float64,N_qp,N_qp)
            @views X_QTDA_JP_ort[:,1] = Spur_CM
            @inbounds for a in 2:N_qp
                @views X_QTDA_JP_ort[:,a] = Pm[:,a-1]

            end

            X_QTDA_JP = X_QTDA_JP_ort
            E_QTDA_JP = E_QTDA_JP_ort

        end

        # Store the QTDA solutions ...
        E_QTDA[P,J+1] = E_QTDA_JP
        X_QTDA[P,J+1] = X_QTDA_JP
    end

    println("\nThe QTDA matrix A succesfully diagonalized ...")

    return E_QTDA, X_QTDA
end
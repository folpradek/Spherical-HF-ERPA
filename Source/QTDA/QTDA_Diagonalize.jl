function QTDA_diagonalize(Params::Parameters,Orb::Vector{Orb1B},Orb_2qp::qpOrb2B,A::Matrix{Matrix{Float64}})
    # Read calculation parameters ...
    N_max = Params.Calc.Nmax
    J_max = 2*N_max + 1
    #Orthogon = Params.Calc.QTDA.Ortho

    # Make the list of values of J & P for iteration ...
    JP_list = JP_initialize(J_max)

    # Initialite the storage for QTDA solutions ...
    E_QTDA = Matrix{Vector{Float64}}(undef,2,J_max+1)
    X_QTDA = Matrix{Matrix{Float64}}(undef,2,J_max+1)

    # Perform the diagonalization of the QTDA matrix A ...
    println("\nDiagonalizing the QTDA matrix A ...")

    @inbounds Threads.@threads for JP in JP_list
        J, P = JP[1], JP[2]
        A_JP = A[P,J+1]

        E_QTDA_JP, X_QTDA_JP = eigen(A_JP, sortby=+)
        Sort = sortperm(E_QTDA_JP, by = x -> real(x))
        E_QTDA_JP, X_QTDA_JP = E_QTDA_JP[Sort], @views X_QTDA_JP[:,Sort]

        # Store the QTDA solutions ...
        E_QTDA[P,J+1] = E_QTDA_JP
        X_QTDA[P,J+1] = X_QTDA_JP
    end

    println("\nThe QTDA matrix A succesfully diagonalized ...")

    return E_QTDA, X_QTDA
end
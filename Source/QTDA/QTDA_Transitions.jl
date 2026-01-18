function QTDA_rM(Params::Parameters,Orb_2qp::qpOrb2B,X_QTDA::Matrix{Matrix{Float64}},qpTrOp::qpTr1B)
    # Read calculation parameters ...
    Orthogon = Params.Calc.QTDA.Ortho

    # Evaluate the reduced transition metrix elements M^lambda ...
    println("\nCalculating the QTDA 1-phonon reduced transition matrix elements rM^lambda ...")

    # E0
    J, P = 0, 1
    N_qp = Orb_2qp.N[P,J+1]
    pM_E0 = Vector{ComplexF64}(undef,N_qp)
    nM_E0 = Vector{ComplexF64}(undef,N_qp)
    @inbounds for nu in 1:N_qp
        pM_E0Sum, nM_E0Sum = 0.0, 0.0
        @inbounds  for i_qp in 1:N_qp
            a, b, T_ab = Orb_2qp.i[P,J+1][i_qp].a, Orb_2qp.i[P,J+1][i_qp].b, Orb_2qp.i[P,J+1][i_qp].T
            if T_ab == -1
                ME = X_QTDA[P,J+1][i_qp,nu] * (qpTrOp.E0.qp20.p[a,b] + qpTrOp.E0.qp20.p[b,a]) / sqrt(1.0 + kronecker_delta(a,b))
                pM_E0Sum += ME
            end

            if T_ab == 1
                ME = X_QTDA[P,J+1][i_qp,nu] * (qpTrOp.E0.qp20.n[a,b] + qpTrOp.E0.qp20.n[b,a]) / sqrt(1.0 + kronecker_delta(a,b))
                nM_E0Sum += ME
            end
        end
        pM_E0[nu] = pM_E0Sum
        nM_E0[nu] = nM_E0Sum
    end

    # E1
    J, P = 1, 2
    N_qp = Orb_2qp.N[P,J+1]
    pM_E1 = Vector{ComplexF64}(undef,N_qp)
    nM_E1 = Vector{ComplexF64}(undef,N_qp)
    @inbounds for nu in 1:N_qp
        pM_E1Sum, nM_E1Sum = 0.0, 0.0
        @inbounds  for i_qp in 1:N_qp
            a, b, T_ab = Orb_2qp.i[P,J+1][i_qp].a, Orb_2qp.i[P,J+1][i_qp].b, Orb_2qp.i[P,J+1][i_qp].T
            if T_ab == -1
                ME = X_QTDA[P,J+1][i_qp,nu] * (qpTrOp.E1.qp20.p[a,b] + qpTrOp.E1.qp20.p[b,a]) / sqrt(1.0 + kronecker_delta(a,b))
                pM_E1Sum += ME
            end

            if T_ab == 1
                ME = X_QTDA[P,J+1][i_qp,nu] * (qpTrOp.E1.qp20.n[a,b] + qpTrOp.E1.qp20.n[b,a]) / sqrt(1.0 + kronecker_delta(a,b))
                nM_E1Sum += ME
            end
        end
        pM_E1[nu] = pM_E1Sum
        nM_E1[nu] = nM_E1Sum
    end

    # E2
    J, P = 2, 1
    N_qp = Orb_2qp.N[P,J+1]
    pM_E2 = Vector{ComplexF64}(undef,N_qp)
    nM_E2 = Vector{ComplexF64}(undef,N_qp)
    @inbounds for nu in 1:N_qp
        pM_E2Sum, nM_E2Sum = 0.0, 0.0
        @inbounds  for i_qp in 1:N_qp
            a, b, T_ab = Orb_2qp.i[P,J+1][i_qp].a, Orb_2qp.i[P,J+1][i_qp].b, Orb_2qp.i[P,J+1][i_qp].T
            if T_ab == -1
                ME = X_QTDA[P,J+1][i_qp,nu] * (qpTrOp.E2.qp20.p[a,b] + qpTrOp.E2.qp20.p[b,a]) / sqrt(1.0 + kronecker_delta(a,b))
                pM_E2Sum += ME
            end

            if T_ab == 1
                ME = X_QTDA[P,J+1][i_qp,nu] * (qpTrOp.E2.qp20.n[a,b] + qpTrOp.E2.qp20.n[b,a]) / sqrt(1.0 + kronecker_delta(a,b))
                nM_E2Sum += ME
            end
        end
        pM_E2[nu] = pM_E2Sum
        nM_E2[nu] = nM_E2Sum
    end

    # E3
    J = 3
    P = 2
    N_qp = Orb_2qp.N[P,J+1]
    pM_E3 = Vector{ComplexF64}(undef,N_qp)
    nM_E3 = Vector{ComplexF64}(undef,N_qp)
    @inbounds for nu in 1:N_qp
        pM_E3Sum, nM_E3Sum = 0.0, 0.0
        @inbounds  for i_qp in 1:N_qp
            a, b, T_ab = Orb_2qp.i[P,J+1][i_qp].a, Orb_2qp.i[P,J+1][i_qp].b, Orb_2qp.i[P,J+1][i_qp].T
            if T_ab == -1
                ME = X_QTDA[P,J+1][i_qp,nu] * (qpTrOp.E3.qp20.p[a,b] + qpTrOp.E3.qp20.p[b,a]) / sqrt(1.0 + kronecker_delta(a,b))
                pM_E3Sum += ME
            end

            if T_ab == 1
                ME =X_QTDA[P,J+1][i_qp,nu] * (qpTrOp.E3.qp20.n[a,b] + qpTrOp.E3.qp20.n[b,a]) / sqrt(1.0 + kronecker_delta(a,b))
                nM_E3Sum += ME
            end
        end
        pM_E3[nu] = pM_E3Sum
        nM_E3[nu] = nM_E3Sum
    end

    rM_E0 = pnCVector(Complex.(pM_E0),Complex.(nM_E0))
    rM_E1 = pnCVector(Complex.(pM_E1),Complex.(nM_E1))
    rM_E2 = pnCVector(Complex.(pM_E2),Complex.(nM_E2))
    rM_E3 = pnCVector(Complex.(pM_E3),Complex.(nM_E3))

    rM_QTDA = ReducedMultipole(rM_E0,rM_E1,rM_E2,rM_E3)

    println("\nQTDA 1-phonon reduced transition matrix elements rM^lambda succesfully calculated ...")

    return rM_QTDA
end

function QTDA_rB(Params::Parameters,Orb_2qp::qpOrb2B,rM::ReducedMultipole)
    # Read parameters ...
    Orthogon = Params.Calc.QTDA.Ortho

    println("\nCalculating the QTDA 1-phonon reduced transition intensities rB^lambda ...")

    # E0
    J, P = 0, 1
    N_qp = Orb_2qp.N[P,J+1]

    rB_phE0 = zeros(Float64,N_qp)
    rB_isE0 = zeros(Float64,N_qp)
    rB_ivE0 = zeros(Float64,N_qp)

    @inbounds for nu in 1:N_qp
        rB_phE0[nu] = abs(rM.E0.p[nu])^2
        rB_isE0[nu] = 0.25 * abs(rM.E0.p[nu] + rM.E0.n[nu])^2
        rB_ivE0[nu] = 0.25 * abs(rM.E0.p[nu] - rM.E0.n[nu])^2
    end

    # E1
    J, P = 1, 2
    N_qp = Orb_2qp.N[P,J+1]

    rB_phE1 = zeros(Float64,N_qp)
    rB_isE1 = zeros(Float64,N_qp)
    rB_ivE1 = zeros(Float64,N_qp)

    if Orthogon == true
        @inbounds for nu in 1:N_qp
            rB_phE1[nu] = abs(rM.E1.p[nu])^2
            rB_isE1[nu] = 0.25 * abs(rM.E1.p[nu] + rM.E1.n[nu])^2
            rB_ivE1[nu] = 0.25 * abs(rM.E1.p[nu] - rM.E1.n[nu])^2
        end
    else
        A, Z = Params.Calc.A, Params.Calc.Z
        e_p = Float64(A - Z) / Float64(A)
        e_n = Float64(Z) / Float64(A)
        @inbounds for nu in 1:N_qp
            rB_phE1[nu] = abs(rM.E1.p[nu])^2
            rB_isE1[nu] = 0.25 * abs(rM.E1.p[nu] + rM.E1.n[nu])^2
            rB_ivE1[nu] = abs(e_p * rM.E1.p[nu] - e_n * rM.E1.n[nu])^2
        end
    end

    # E2
    J, P = 2, 1
    N_qp = Orb_2qp.N[P,J+1]

    rB_phE2 = zeros(Float64,N_qp)
    rB_isE2 = zeros(Float64,N_qp)
    rB_ivE2 = zeros(Float64,N_qp)

    @inbounds for nu in 1:N_qp
        rB_phE2[nu] = abs(rM.E2.p[nu])^2
        rB_isE2[nu] = 0.25 * abs(rM.E2.p[nu] + rM.E2.n[nu])^2
        rB_ivE2[nu] = 0.25 * abs(rM.E2.p[nu] - rM.E2.n[nu])^2
    end

    # E3
    J, P = 3, 2
    N_qp = Orb_2qp.N[P,J+1]

    rB_phE3 = zeros(Float64,N_qp)
    rB_isE3 = zeros(Float64,N_qp)
    rB_ivE3 = zeros(Float64,N_qp)

    @inbounds for nu in 1:N_qp
        rB_phE3[nu] = abs(rM.E3.p[nu])^2
        rB_isE3[nu] = 0.25 * abs(rM.E3.p[nu] + rM.E3.n[nu])^2
        rB_ivE3[nu] = 0.25 * abs(rM.E3.p[nu] - rM.E3.n[nu])^2
    end

    rB_E0 = Transition(rB_phE0,rB_isE0,rB_ivE0)
    rB_E1 = Transition(rB_phE1,rB_isE1,rB_ivE1)
    rB_E2 = Transition(rB_phE2,rB_isE2,rB_ivE2)
    rB_E3 = Transition(rB_phE3,rB_isE3,rB_ivE3)

    rB = ReducedTransition(rB_E0,rB_E1,rB_E2,rB_E3)

    println("\nQTDA 1-phonon reduced transition intensities rB have been succesfully calculated ...")

    return rB
end
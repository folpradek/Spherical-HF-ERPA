function QTDA_rM(Params::Parameters,Orb::Vector{NOrb},Orb_2qp::qpOrb2B,X_QTDA::Matrix{Matrix{Float64}},qpTrOp::qpTr1B)
    # Read calculation parameters ...
    Orthogon = Params.Calc.QTDA.Ortho

    # Evaluate the reduced transition metrix elements M^lambda ...
    println("\nCalculating reduced transition matrix elements rM ...")

    # E0
    J, P = 0, 1
    N_qp = Orb_2qp.N[J+1,P]
    pM_E0 = Vector{ComplexF64}(undef,N_qp)
    nM_E0 = Vector{ComplexF64}(undef,N_qp)
    @inbounds for nu in 1:N_qp
        pM_E0Sum, nM_E0Sum = 0.0, 0.0
        @inbounds  for i_qp in 1:N_qp

            a, b, T_ab = Orb_2qp.qp[J+1,P][i_qp].a, Orb_2qp.qp[J+1,P][i_qp].b, Orb_2qp.qp[J+1,P][i_qp].T


            if t_ph == -1
                ME = qpTrOp.E0.qp22.p[a,b] * X_QTDA[J+1,P][i_qp,nu] / sqrt(1.0 + kron_delta(a,b))
                pM_E0Sum += ME
            end

            if t_ph == 1
                ME0_TDA = TrOp.E0.n[a_p,a_h] * X_TDA[J+1,P][Ind_ph,nu]
                nME0Sum_TDA += ME0_TDA

                ME0_RPA = TrOp.E0.n[a_p,a_h] * (X_RPA[J+1,P][Ind_ph,nu] + Y_RPA[J+1,P][Ind_ph,nu])
                nME0Sum_RPA += ME0_RPA
            end

        end
        prME0_TDA[nu] = pME0Sum_TDA
        nrME0_TDA[nu] = nME0Sum_TDA
        prME0_RPA[nu] = pME0Sum_RPA
        nrME0_RPA[nu] = nME0Sum_RPA
    end

    # E1
    J = 1
    P = 2
    N_ph = N_nu[J+1,P]
    prME1_TDA = Vector{ComplexF64}(undef,N_ph)
    prME1_RPA = Vector{ComplexF64}(undef,N_ph)
    nrME1_TDA = Vector{ComplexF64}(undef,N_ph)
    nrME1_RPA = Vector{ComplexF64}(undef,N_ph)
    @inbounds  for nu in 1:N_ph
        if (Orthogon == true && nu != 1) || (Orthogon == false)
            pME1Sum_TDA = ComplexF64(0.0)
            pME1Sum_RPA = ComplexF64(0.0)
            nME1Sum_TDA = ComplexF64(0.0)
            nME1Sum_RPA = ComplexF64(0.0)
            @inbounds  for Ind_ph in 1:N_ph
                ph = Orb_Phonon[J+1,P][Ind_ph]
                p, h = Phonon[ph].p, Phonon[ph].h
                t_ph = Phonon[ph].tz
                if t_ph == -1
                    a_p = Particle.p[p].a
                    a_h = Hole.p[h].a
                elseif t_ph == 1
                    a_p = Particle.n[p].a
                    a_h = Hole.n[h].a
                end

                if t_ph == -1
                    ME1_TDA = TrOp.E1.p[a_p,a_h] * X_TDA[J+1,P][Ind_ph,nu] * -1.0
                    pME1Sum_TDA += ME1_TDA

                    ME1_RPA = TrOp.E1.p[a_p,a_h] * (-1.0 * X_RPA[J+1,P][Ind_ph,nu] + Y_RPA[J+1,P][Ind_ph,nu])
                    pME1Sum_RPA += ME1_RPA
                end

                if t_ph == 1
                    ME1_TDA = TrOp.E1.n[a_p,a_h] * X_TDA[J+1,P][Ind_ph,nu] * -1.0
                    nME1Sum_TDA += ME1_TDA

                    ME1_RPA = TrOp.E1.n[a_p,a_h] * (-1.0 * X_RPA[J+1,P][Ind_ph,nu] + Y_RPA[J+1,P][Ind_ph,nu])
                    nME1Sum_RPA += ME1_RPA
                end

            end
            prME1_TDA[nu] = pME1Sum_TDA
            nrME1_TDA[nu] = nME1Sum_TDA
            prME1_RPA[nu] = pME1Sum_RPA
            nrME1_RPA[nu] = nME1Sum_RPA
        end
    end

    # E2
    J = 2
    P = 1
    N_ph = N_nu[J+1,P]
    prME2_TDA = Vector{ComplexF64}(undef,N_ph)
    prME2_RPA = Vector{ComplexF64}(undef,N_ph)
    nrME2_TDA = Vector{ComplexF64}(undef,N_ph)
    nrME2_RPA = Vector{ComplexF64}(undef,N_ph)
    @inbounds  for nu in 1:N_ph
        pME2Sum_TDA = ComplexF64(0.0)
        pME2Sum_RPA = ComplexF64(0.0)
        nME2Sum_TDA = ComplexF64(0.0)
        nME2Sum_RPA = ComplexF64(0.0)
        @inbounds  for Ind_ph in 1:N_ph
            ph = Orb_Phonon[J+1,P][Ind_ph]
            p, h = Phonon[ph].p, Phonon[ph].h
            t_ph = Phonon[ph].tz
            if t_ph == -1
                a_p = Particle.p[p].a
                a_h = Hole.p[h].a
            elseif t_ph == 1
                a_p = Particle.n[p].a
                a_h = Hole.n[h].a
            end

            if t_ph == -1
                ME2_TDA = TrOp.E2.p[a_p,a_h] * X_TDA[J+1,P][Ind_ph,nu]
                pME2Sum_TDA += ME2_TDA

                ME2_RPA = TrOp.E2.p[a_p,a_h] * (X_RPA[J+1,P][Ind_ph,nu] + Y_RPA[J+1,P][Ind_ph,nu])
                pME2Sum_RPA += ME2_RPA
            end

            if t_ph == 1
                ME2_TDA = TrOp.E2.n[a_p,a_h] * X_TDA[J+1,P][Ind_ph,nu]
                nME2Sum_TDA += ME2_TDA

                ME2_RPA = TrOp.E2.n[a_p,a_h] * (X_RPA[J+1,P][Ind_ph,nu] + Y_RPA[J+1,P][Ind_ph,nu])
                nME2Sum_RPA += ME2_RPA
            end

        end
        prME2_TDA[nu] = pME2Sum_TDA
        nrME2_TDA[nu] = nME2Sum_TDA
        prME2_RPA[nu] = pME2Sum_RPA
        nrME2_RPA[nu] = nME2Sum_RPA
    end

    # E3
    J = 3
    P = 2
    N_ph = N_nu[J+1,P]
    prME3_TDA = Vector{ComplexF64}(undef,N_ph)
    prME3_RPA = Vector{ComplexF64}(undef,N_ph)
    nrME3_TDA = Vector{ComplexF64}(undef,N_ph)
    nrME3_RPA = Vector{ComplexF64}(undef,N_ph)
    @inbounds  for nu in 1:N_ph
        pME3Sum_TDA = ComplexF64(0.0)
        pME3Sum_RPA = ComplexF64(0.0)
        nME3Sum_TDA = ComplexF64(0.0)
        nME3Sum_RPA = ComplexF64(0.0)
        @inbounds  for Ind_ph in 1:N_ph
            ph = Orb_Phonon[J+1,P][Ind_ph]
            p, h = Phonon[ph].p, Phonon[ph].h
            t_ph = Phonon[ph].tz
            if t_ph == -1
                a_p = Particle.p[p].a
                a_h = Hole.p[h].a
            elseif t_ph == 1
                a_p = Particle.n[p].a
                a_h = Hole.n[h].a
            end

            if t_ph == -1
                ME3_TDA = TrOp.E3.p[a_p,a_h] * X_TDA[J+1,P][Ind_ph,nu] * -1.0
                pME3Sum_TDA += ME3_TDA

                ME3_RPA = TrOp.E3.p[a_p,a_h] * (-1.0 * X_RPA[J+1,P][Ind_ph,nu] + Y_RPA[J+1,P][Ind_ph,nu])
                pME3Sum_RPA += ME3_RPA
            end

            if t_ph == 1
                ME3_TDA = TrOp.E3.n[a_p,a_h] * X_TDA[J+1,P][Ind_ph,nu] * -1.0
                nME3Sum_TDA += ME3_TDA

                ME3_RPA = TrOp.E3.n[a_p,a_h] * (-1.0 * X_RPA[J+1,P][Ind_ph,nu] + Y_RPA[J+1,P][Ind_ph,nu])
                nME3Sum_RPA += ME3_RPA
            end

        end
        prME3_TDA[nu] = pME3Sum_TDA
        nrME3_TDA[nu] = nME3Sum_TDA
        prME3_RPA[nu] = pME3Sum_RPA
        nrME3_RPA[nu] = nME3Sum_RPA
    end

    rM0_TDA = pnCVector(prME0_TDA,nrME0_TDA)
    rM1_TDA = pnCVector(prME1_TDA,nrME1_TDA)
    rM2_TDA = pnCVector(prME2_TDA,nrME2_TDA)
    rM3_TDA = pnCVector(prME3_TDA,nrME3_TDA)

    rM0_RPA = pnCVector(prME0_RPA,nrME0_RPA)
    rM1_RPA = pnCVector(prME1_RPA,nrME1_RPA)
    rM2_RPA = pnCVector(prME2_RPA,nrME2_RPA)
    rM3_RPA = pnCVector(prME3_RPA,nrME3_RPA)

    rM_TDA = ReducedMultipole(rM0_TDA,rM1_TDA,rM2_TDA,rM3_TDA)
    rM_RPA = ReducedMultipole(rM0_RPA,rM1_RPA,rM2_RPA,rM3_RPA)

    println("\nReduced transition matrix elements rM calculated ...")

    return rM_QTDA
end

function QTDA_rB(Params::Parameters,Orb_2qp::qpOrb2B,rM::ReducedMultipole)
    # Read parameters ...
    Orthogon = Params.Calc.QTDA.Ortho

    println("\nCalculating reduced transition intensities rB ...")

    # E0
    J = 0
    P = 1
    N_ph = N_nu[J+1,P]

    rB_phE0 = zeros(Float64,N_ph)
    rB_isE0 = zeros(Float64,N_ph)
    rB_ivE0 = zeros(Float64,N_ph)

    @inbounds for nu in 1:N_ph
        rB_phE0[nu] = abs(rM.E0.p[nu])^2
        rB_isE0[nu] = 0.25 * abs(rM.E0.p[nu] + rM.E0.n[nu])^2
        rB_ivE0[nu] = 0.25 * abs(rM.E0.p[nu] - rM.E0.n[nu])^2
    end

    # E1
    J = 1
    P = 2
    N_ph = N_nu[J+1,P]

    rB_phE1 = zeros(Float64,N_ph)
    rB_isE1 = zeros(Float64,N_ph)
    rB_ivE1 = zeros(Float64,N_ph)

    if Orthogon == true
        @inbounds for nu in 1:N_ph
            rB_phE1[nu] = abs(rM.E1.p[nu])^2
            rB_isE1[nu] = 0.25 * abs(rM.E1.p[nu] + rM.E1.n[nu])^2
            rB_ivE1[nu] = 0.25 * abs(rM.E1.p[nu] - rM.E1.n[nu])^2
        end
    else
        A = Params.Calc.A
        Z = Params.Calc.Z
        e_p = Float64(A - Z) / Float64(A)
        e_n = Float64(Z) / Float64(A)
        @inbounds for nu in 1:N_ph
            rB_phE1[nu] = abs(rM.E1.p[nu])^2
            rB_isE1[nu] = 0.25 * abs(rM.E1.p[nu] + rM.E1.n[nu])^2
            rB_ivE1[nu] = abs(e_p * rM.E1.p[nu] - e_n * rM.E1.n[nu])^2
        end
    end

    # E2
    J = 2
    P = 1
    N_ph = N_nu[J+1,P]

    rB_phE2 = zeros(Float64,N_ph)
    rB_isE2 = zeros(Float64,N_ph)
    rB_ivE2 = zeros(Float64,N_ph)

    @inbounds for nu in 1:N_ph
        rB_phE2[nu] = abs(rM.E2.p[nu])^2
        rB_isE2[nu] = 0.25 * abs(rM.E2.p[nu] + rM.E2.n[nu])^2
        rB_ivE2[nu] = 0.25 * abs(rM.E2.p[nu] - rM.E2.n[nu])^2
    end

    # E3
    J = 3
    P = 2
    N_ph = N_nu[J+1,P]

    rB_phE3 = zeros(Float64,N_ph)
    rB_isE3 = zeros(Float64,N_ph)
    rB_ivE3 = zeros(Float64,N_ph)

    @inbounds for nu in 1:N_ph
        rB_phE3[nu] = abs(rM.E3.p[nu])^2
        rB_isE3[nu] = 0.25 * abs(rM.E3.p[nu] + rM.E3.n[nu])^2
        rB_ivE3[nu] = 0.25 * abs(rM.E3.p[nu] - rM.E3.n[nu])^2
    end

    rB_E0 = Transition(rB_phE0,rB_isE0,rB_ivE0)
    rB_E1 = Transition(rB_phE1,rB_isE1,rB_ivE1)
    rB_E2 = Transition(rB_phE2,rB_isE2,rB_ivE2)
    rB_E3 = Transition(rB_phE3,rB_isE3,rB_ivE3)

    rB_ph = ReducedTransition(rB_E0,rB_E1,rB_E2,rB_E3)

    println("\nReduced transition intensities rB evaluated ...")

    return rB_ph
end
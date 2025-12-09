function HF_ERPA_OBDM_Iteration(Params::Parameters,N_nu::Matrix{Int64},N_Particle::pnInteger,Particle::pnSVector,N_Hole::pnInteger,Hole::pnSVector,X_RPA::Matrix{Matrix{ComplexF64}},Y_RPA::Matrix{Matrix{ComplexF64}},Rho_HF::pnMatrix,Rho_Ini::pnMatrix)
    # Read parameters ...
    N_max = Params.Calc.Nmax
    N_2max = 2*N_max
    J_max = N_2max + 1
    a_max = div((N_max + 1)*(N_max + 2),2)

    pU, nU = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)
    pU_ERPA, nU_ERPA = diagm(ones(Float64,a_max)), diagm(ones(Float64,a_max))

    # Initialize new density matrices ...
    pRho_Old, nRho_Old = deepcopy(Rho_Ini.p), deepcopy(Rho_Ini.n)
    pRho, nRho = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)

    # Iteration parameters ...
    Delta = 1.0
    Eps = 1e-6
    Cutting_Precision = 1e-8
    N_Iteration = 100
    Iteration = 0
    Delta_Odd = 1.0
    Delta_Even = 1.0

    # Make Particle-Hole orbitals list ...
    ParticleHole = Make_ParticleHole_List(J_max,N_Particle,Particle,N_Hole,Hole)

    # Make Particle-Particle & Hole-Hole lists for OBDM iteration ...
    Pairs = HF_ERPA_OBDM_Iteration_List(N_Particle,Particle,N_Hole,Hole)

    # Make J & P list for iteration ...
    JP_List = JP_Ini(J_max)

    # Initialize inner loop matrix ...
    pM = Matrix{Matrix{ComplexF64}}(undef,J_max+1,2)
    nM = Matrix{Matrix{ComplexF64}}(undef,J_max+1,2)
    @inbounds Threads.@threads for JP in JP_List
        J = JP[1]
        P = JP[2]
        J_ind = J + 1
        N_ph = N_nu[J_ind,P]
        pM[J_ind,P], nM[J_ind,P] = zeros(ComplexF64,N_ph,N_ph), zeros(ComplexF64,N_ph,N_ph)
    end

    println("\nStarting HF ERPA OBDM calculation ...")

    # OBDM iteration loop ...
    while (Delta > Eps) && (N_Iteration > Iteration)
        # Evaluate new OBDM matrix ...
        pRho, nRho = HF_ERPA_OBDM(N_nu,N_Particle,Particle,N_Hole,Hole,ParticleHole,Pairs,JP_List,pM,nM,X_RPA,Y_RPA,Rho_HF,pnMatrix(pRho_Old,nRho_Old))

        # Symmetrize OBDM ...
        pRho .= 0.5 .* (pRho .+ pRho')
        nRho .= 0.5 .* (nRho .+ nRho')

        # Remove off-diagonal elements from OBDM if diagonal approximation is set ...
        if Params.Calc.ERPA.ScOBH == false
            @inbounds for a in 1:a_max
                @inbounds for b in 1:a_max
                    if a != b
                        pRho[a,b] = 0.0
                        nRho[a,b] = 0.0
                    end
                end
            end
        end

        # Clean Rho ...
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                pr, nr = pRho[a,b], nRho[a,b]
                if abs(pr) < Cutting_Precision
                     pRho[a,b] = 0.0
                end
                if abs(nr) < Cutting_Precision
                    nRho[a,b] = 0.0
                end
            end
        end

        # Diagonalize OBDM ...
        pN, pU = eigen(pRho, sortby = nothing)
        nN, nU = eigen(nRho, sortby = nothing)

        # Clean U ...
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                pu, nu = pU[a,b], nU[a,b]
                if abs(pu) < Cutting_Precision
                     pU[a,b] = 0.0
                end
                if abs(nu) < Cutting_Precision
                    nU[a,b] = 0.0
                end
            end
        end

        # Reorder orbitals if swapped ...
        pN, pU = HF_ERPA_OBDM_Eigen_Reorder(a_max,pN,pU)
        nN, nU = HF_ERPA_OBDM_Eigen_Reorder(a_max,nN,nU)

        # Check positiveness of the s.p. orbits occupations probabilities differences ... the norm matrix ...
        Occupation = HF_ERPA_OBDM_Occupation_Check(pnVector(pN,nN),N_Particle,Particle,N_Hole,Hole)
        if Occupation == false || Iteration == 99
            println("\nOBDM occupation finished - Singular solution - ERPA iteration is about to finish with the last stable OBDM solution ...")
            return Rho_Ini, pnMatrix(diagm(ones(Float64,a_max)),diagm(ones(Float64,a_max)))
        end

        # Update the HF-ERPA NO transformation matrix ...
        pU_ERPA .= deepcopy(pU_ERPA * pU)
        nU_ERPA .= deepcopy(nU_ERPA * nU)

        # Clean U_ERPA ...
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                pu, nu = pU_ERPA[a,b], nU_ERPA[a,b]
                if abs(pu) < Cutting_Precision
                     pU_ERPA[a,b] = 0.0
                end
                if abs(nu) < Cutting_Precision
                    nU_ERPA[a,b] = 0.0
                end
            end
        end

        # Allocate new OBDM ...
        pRho, nRho = diagm(pN), diagm(nN)

        # Transform the old OBDM to the new HF-ERPA OBDM basis ...
        pRho_Old .= deepcopy(pU' * pRho_Old * pU)
        nRho_Old .= deepcopy(nU' * nRho_Old * nU)

        # Evaluate the difference of OBDMs ...
        dSum = 0.0
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                ME = (abs(pRho[a,b] - pRho_Old[a,b]) + abs(nRho[a,b] - nRho_Old[a,b])) / Float64(2*a_max^2)
                dSum += ME
            end
        end

        # Transform amplitudes X and Y into new HF-ERPA NO basis ...
        X_RPA, Y_RPA = HF_ERPA_OBDM_Transform_Amplitudes(N_nu,N_Particle,Particle,N_Hole,Hole,ParticleHole,JP_List,X_RPA,Y_RPA,pnMatrix(pU,nU))

        # Evaluate iteration ...
        Delta = dSum
        Iteration += 1

        # Check for degenerate solutions ... Tend to appear ...
        if rem(Iteration,2) == 0
            if abs(Delta - Delta_Even) < 1e-6
                println("\nOBDM iteration finished - cycled degenerate solution found ...")
                break
            end
            Delta_Even = Delta
        else
            if abs(Delta - Delta_Odd) < 1e-6
                println("\nOBDM iteration finished - cycled degenerate solution found ...")
                break
            end
            Delta_Odd = Delta
        end

        println("\nOBDM iteration:   Iteration Number = " * string(Iteration) * ",\tDelta = " * string(Delta))

        # Update old OBDM ...
        pRho_Old .= deepcopy(pRho)
        nRho_Old .= deepcopy(nRho)
    end

    println("\nERPA OBDM iteration has finished ...")

    return pnMatrix(pRho,nRho), pnMatrix(pU_ERPA,nU_ERPA)
end

function HF_ERPA_OBDM_Iteration_List(N_Particle::pnInteger,Particle::pnSVector,N_Hole::pnInteger,Hole::pnSVector)
    # Proton particle pairs ...
    Pair_N = 0
    @inbounds for p in 1:N_Particle.p
        a_p, l_p, j_p = OBDM_pParticle(Particle,p)
        @inbounds for q in 1:N_Particle.p
            a_q, l_q, j_q = OBDM_pParticle(Particle,q)
            if l_p == l_q && j_p == j_q && a_p >= a_q
                Pair_N += 1
            end
        end
    end

    # Proton hole pairs ...
    @inbounds for h in 1:N_Hole.p
        a_h, l_h, j_h = OBDM_pHole(Hole,h)
        @inbounds for g in 1:N_Hole.p
            a_g, l_g, j_g = OBDM_pHole(Hole,g)
            if l_h == l_g && j_h == j_g && a_h >= a_g
                Pair_N += 1
            end
        end
    end

    # Neutron particle pairs ...
    @inbounds for p in 1:N_Particle.n
        a_p, l_p, j_p = OBDM_nParticle(Particle,p)
        @inbounds for q in 1:N_Particle.n
            a_q, l_q, j_q = OBDM_nParticle(Particle,q)
            if l_p == l_q && j_p == j_q && a_p >= a_q
                Pair_N += 1
            end
        end
    end

    # Neutron hole pairs ...
    @inbounds for h in 1:N_Hole.n
        a_h, l_h, j_h = OBDM_nHole(Hole,h)
        @inbounds for g in 1:N_Hole.n
            a_g, l_g, j_g = OBDM_nHole(Hole,g)
            if l_h == l_g && j_h == j_g && a_h >= a_g
                Pair_N += 1
            end
        end
    end

    # Initialize vectors for Pair structure ...
    Pair_Type, Pair_a, Pair_b = zeros(Int64,Pair_N), zeros(Int64,Pair_N), zeros(Int64,Pair_N)

    # Initialize index counter ...
    Index = 0

    # Proton particle pairs ...
    @inbounds for p in 1:N_Particle.p
        a_p, l_p, j_p = OBDM_pParticle(Particle,p)
        @inbounds for q in 1:N_Particle.p
            a_q, l_q, j_q = OBDM_pParticle(Particle,q)
            if l_p == l_q && j_p == j_q && a_p >= a_q
                Index += 1
                Pair_Type[Index] = 1
                Pair_a[Index] = p
                Pair_b[Index] = q
            end
        end
    end

    # Proton hole pairs ...
    @inbounds for h in 1:N_Hole.p
        a_h, l_h, j_h = OBDM_pHole(Hole,h)
        @inbounds for g in 1:N_Hole.p
            a_g, l_g, j_g = OBDM_pHole(Hole,g)
            if l_h == l_g && j_h == j_g && a_h >= a_g
                Index += 1
                Pair_Type[Index] = 2
                Pair_a[Index] = h
                Pair_b[Index] = g
            end
        end
    end

    # Neutron particle pairs ...
    @inbounds for p in 1:N_Particle.n
        a_p, l_p, j_p = OBDM_nParticle(Particle,p)
        @inbounds for q in 1:N_Particle.n
            a_q, l_q, j_q = OBDM_nParticle(Particle,q)
            if l_p == l_q && j_p == j_q && a_p >= a_q
                Index += 1
                Pair_Type[Index] = 3
                Pair_a[Index] = p
                Pair_b[Index] = q
            end
        end
    end

    # Neutron hole pairs ...
    @inbounds for h in 1:N_Hole.n
        a_h, l_h, j_h = OBDM_nHole(Hole,h)
        @inbounds for g in 1:N_Hole.n
            a_g, l_g, j_g = OBDM_nHole(Hole,g)
            if l_h == l_g && j_h == j_g && a_h >= a_g
                Index += 1
                Pair_Type[Index] = 4
                Pair_a[Index] = h
                Pair_b[Index] = g
            end
        end
    end

    return OBDM_Iteration_Pairs(Pair_N,Pair_Type,Pair_a,Pair_b)
end

function HF_ERPA_OBDM_Transform_Amplitudes(N_nu::Matrix{Int64},N_Particle::pnInteger,Particle::pnSVector,N_Hole::pnInteger,Hole::pnSVector,ParticleHole::pnArray,JP_List::Vector{Vector{Int64}},X_RPA::Matrix{Matrix{ComplexF64}},Y_RPA::Matrix{Matrix{ComplexF64}},U::pnMatrix)
    # Read transformation matrices pU & nU ...
    pU, nU = U.p, U.n

    # Transform amplitudes X and Y into new HF-ERPA NO basis ...
    @inbounds Threads.@threads for JP in JP_List
        J = JP[1]
        P = JP[2]
        J_ind = J + 1
        N_ph = N_nu[J_ind,P]
        X_RPA_JP = zeros(ComplexF64,N_ph,N_ph)
        Y_RPA_JP = zeros(ComplexF64,N_ph,N_ph)
        @inbounds for nu in 1:N_ph
            # Protons
            @inbounds for h in 1:N_Hole.p
                a_h = Hole.p[h].a
                l_h = Hole.p[h].l
                j_h = Hole.p[h].j
                @inbounds for p in 1:N_Particle.p
                    a_p = Particle.p[p].a
                    l_p = Particle.p[p].l
                    j_p = Particle.p[p].j
                    XSum = ComplexF64(0.0)
                    YSum = ComplexF64(0.0)
                    if ( rem(l_p + l_h,2) + 1) == P && div(abs(j_p - j_h),2) <= J && J <= div(j_p + j_h,2)
                        ph = ParticleHole.p[J_ind,P,p,h]
                        XSum = ComplexF64(0.0)
                        YSum = ComplexF64(0.0)
                        @inbounds for g in 1:N_Hole.p
                            a_g = Hole.p[g].a
                            l_g = Hole.p[g].l
                            j_g = Hole.p[g].j
                            if l_g == l_h && j_g == j_h
                                @inbounds for q in 1:N_Particle.p
                                    a_q = Particle.p[q].a
                                    l_q = Particle.p[q].l
                                    j_q = Particle.p[q].j
                                    if l_q == l_p && j_q == j_p
                                        qg = ParticleHole.p[J_ind,P,q,g]
                                        YSum += ComplexF64(pU[a_q,a_p] * pU[a_g,a_h]) * Y_RPA[J_ind,P][qg,nu]
                                        XSum += ComplexF64(pU[a_q,a_p] * pU[a_g,a_h]) * X_RPA[J_ind,P][qg,nu]
                                    end
                                end
                            end
                        end
                        X_RPA_JP[ph,nu] += XSum
                        Y_RPA_JP[ph,nu] += YSum
                    end
                end
            end

            # Neutrons
            @inbounds for h in 1:N_Hole.n
                a_h = Hole.n[h].a
                l_h = Hole.n[h].l
                j_h = Hole.n[h].j
                @inbounds for p in 1:N_Particle.n
                    a_p = Particle.n[p].a
                    l_p = Particle.n[p].l
                    j_p = Particle.n[p].j
                    if ( rem(l_p + l_h,2) + 1) == P && div(abs(j_p - j_h),2) <= J && J <= div(j_p + j_h,2)
                        ph = ParticleHole.n[J_ind,P,p,h]
                        XSum = ComplexF64(0.0)
                        YSum = ComplexF64(0.0)
                        @inbounds for g in 1:N_Hole.n
                            a_g = Hole.n[g].a
                            l_g = Hole.n[g].l
                            j_g = Hole.n[g].j
                            if l_g == l_h && j_g == j_h
                                @inbounds for q in 1:N_Particle.n
                                    a_q = Particle.n[q].a
                                    l_q = Particle.n[q].l
                                    j_q = Particle.n[q].j
                                    if l_q == l_p && j_q == j_p
                                        qg = ParticleHole.n[J_ind,P,q,g]
                                        YSum += ComplexF64(nU[a_q,a_p] * nU[a_g,a_h]) * Y_RPA[J_ind,P][qg,nu]
                                        XSum += ComplexF64(nU[a_q,a_p] * nU[a_g,a_h]) * X_RPA[J_ind,P][qg,nu]
                                    end
                                end
                            end
                        end
                        X_RPA_JP[ph,nu] += XSum
                        Y_RPA_JP[ph,nu] += YSum
                    end
                end
            end
        end
        X_RPA[J_ind,P] .= X_RPA_JP
        Y_RPA[J_ind,P] .= Y_RPA_JP
    end

    return X_RPA, Y_RPA
end

function HF_ERPA_OBDM(N_nu::Matrix{Int64},N_Particle::pnInteger,Particle::pnSVector,N_Hole::pnInteger,Hole::pnSVector,ParticleHole::pnArray,Pairs::OBDM_Iteration_Pairs,JP_List::Vector{Vector{Int64}},pM::Matrix{Matrix{ComplexF64}},nM::Matrix{Matrix{ComplexF64}},X_RPA::Matrix{Matrix{ComplexF64}},Y_RPA::Matrix{Matrix{ComplexF64}},Rho_HF::pnMatrix,Rho_Old::pnMatrix)
    # Load the HF OBDM ...
    pRho, nRho = deepcopy(Rho_HF.p), deepcopy(Rho_HF.n)
    pRho_Old, nRho_Old = Rho_Old.p, Rho_Old.n
    
    # Precalculate inner loop ...
    @inbounds Threads.@threads for JP in JP_List
        J = JP[1]
        P = JP[2]
        J_ind = J + 1
        N_ph = N_nu[J_ind,P]
        pSum = ComplexF64(0.0)
        nSum = ComplexF64(0.0)
        @inbounds for nu in 1:N_ph
            @inbounds for mu in 1:N_ph
                pSum = ComplexF64(0.0)
                nSum = ComplexF64(0.0)
                @inbounds for q in 1:N_Particle.p
                    a_q, l_q, j_q = OBDM_pParticle(Particle,q)
                    @inbounds for g in 1:N_Hole.p
                        a_g, l_g, j_g = OBDM_pHole(Hole,g)
                        if (P == (rem(l_q + l_g,2) + 1)) && (abs(j_q - j_g) <= 2*J) && (2*J <= (j_q + j_g))
                            Amp_qg = OBDM_D(a_q,a_g,pRho_Old)
                            qg = ParticleHole.p[J_ind,P,q,g]
                            X_qg_nu = OBDM_X(X_RPA,J_ind,P,qg,nu)
                            X_qg_mu = OBDM_X(X_RPA,J_ind,P,qg,mu)
                            ME = OBDM_M(Amp_qg,X_qg_nu,X_qg_mu)
                            pSum += ME
                        end
                    end
                end
                @inbounds for q in 1:N_Particle.n
                    a_q, l_q, j_q = OBDM_nParticle(Particle,q)
                    @inbounds for g in 1:N_Hole.n
                        a_g, l_g, j_g = OBDM_nHole(Hole,g)
                        if (P == (rem(l_q + l_g,2) + 1)) && (abs(j_q - j_g) <= 2*J) && (2*J <= (j_q + j_g))
                            Amp_qg = OBDM_D(a_q,a_g,nRho_Old)
                            qg = ParticleHole.n[J_ind,P,q,g]
                            X_qg_nu = OBDM_X(X_RPA,J_ind,P,qg,nu)
                            X_qg_mu = OBDM_X(X_RPA,J_ind,P,qg,mu)
                            ME = OBDM_M(Amp_qg,X_qg_nu,X_qg_mu)
                            nSum += ME
                        end
                    end
                end
                pM[J_ind,P][nu,mu] = pSum
                nM[J_ind,P][nu,mu] = nSum
            end
        end
    end

    # OBDM iteration ...
    @inbounds Threads.@threads for Index in 1:Pairs.N
        Type = Pairs.T[Index]
        # pRho particles iteration ...
        if Type == 1
            p = Pairs.a[Index]
            q = Pairs.b[Index]
            a_p, l_p, j_p = OBDM_pParticle(Particle,p)
            a_q = Particle.p[q].a
            Sum = 0.0
            @inbounds for h in 1:N_Hole.p
                a_h, l_h, j_h = OBDM_pHole(Hole,h)
                Amp_ph, Amp_qh = OBDM_D(a_p,a_h,pRho_Old), OBDM_D(a_q,a_h,pRho_Old)
                P = rem(l_p + l_h,2) + 1
                @inbounds for J in div(abs(j_p - j_h),2):div(j_p + j_h,2)
                    J_ind = J + 1
                    N_ph = N_nu[J_ind,P]
                    ph, qh = ParticleHole.p[J_ind,P,p,h], ParticleHole.p[J_ind,P,q,h]
                    @inbounds for nu in 1:N_ph
                        Y_ph_nu = OBDM_Y(Y_RPA,J_ind,P,ph,nu)
                        Y_qh_nu = OBDM_Y(Y_RPA,J_ind,P,qh,nu)
                        @inbounds for mu in 1:N_ph
                            Y_qh_mu = OBDM_Y(Y_RPA,J_ind,P,qh,mu)
                            ME = Float64(2*J + 1) * OBDM_ME_Inner(pM,J_ind,P,nu,mu,Amp_ph,Amp_qh,Y_ph_nu,Y_qh_mu)
                            Sum += ME
                            if a_p > a_q
                                Y_ph_mu = OBDM_Y(Y_RPA,J_ind,P,ph,mu)
                                ME = Float64(2*J + 1) * OBDM_ME_Inner(pM,J_ind,P,mu,nu,Amp_ph,Amp_qh,Y_ph_mu,Y_qh_nu)
                                Sum += ME
                            end
                        end
                        ME = Float64(2*J + 1) * OBDM_ME_Outer(Amp_ph,Amp_qh,Y_ph_nu,Y_qh_nu)
                        Sum += ME
                    end
                end
            end
            Sum = Sum / Float64(j_p + 1)
            pRho[a_p,a_q] += Sum
            if a_p > a_q
                pRho[a_q,a_p] += Sum
            end
        # pRho holes iteration ...
        elseif Type == 2
            h = Pairs.a[Index]
            g = Pairs.b[Index]
            a_h, l_h, j_h = OBDM_pHole(Hole,h)
            a_g = Hole.p[g].a
            Sum = 0.0
            @inbounds for p in 1:N_Particle.p
                a_p, l_p, j_p = OBDM_pParticle(Particle,p)
                Amp_ph, Amp_pg = OBDM_D(a_p,a_h,pRho_Old), OBDM_D(a_p,a_g,pRho_Old)
                P = rem(l_p + l_h,2) + 1
                @inbounds for J in div(abs(j_p - j_h),2):div(j_p + j_h,2)
                    J_ind = J + 1
                    N_ph = N_nu[J_ind,P]
                    ph, pg = ParticleHole.p[J_ind,P,p,h], ParticleHole.p[J_ind,P,p,g]
                    @inbounds for nu in 1:N_ph
                        Y_ph_nu = OBDM_Y(Y_RPA,J_ind,P,ph,nu)
                        Y_pg_nu = OBDM_Y(Y_RPA,J_ind,P,pg,nu)
                        @inbounds for mu in 1:N_ph
                            Y_pg_mu = OBDM_Y(Y_RPA,J_ind,P,pg,mu)
                            ME = Float64(2*J + 1) * OBDM_ME_Inner(pM,J_ind,P,nu,mu,Amp_ph,Amp_pg,Y_ph_nu,Y_pg_mu)
                            Sum += ME
                            if a_h > a_g
                                Y_ph_mu = OBDM_Y(Y_RPA,J_ind,P,ph,mu)
                                ME = Float64(2*J + 1) * OBDM_ME_Inner(pM,J_ind,P,mu,nu,Amp_ph,Amp_pg,Y_ph_mu,Y_pg_nu)
                                Sum += ME
                            end
                        end
                        ME = Float64(2*J + 1) * OBDM_ME_Outer(Amp_ph,Amp_pg,Y_ph_nu,Y_pg_nu)
                        Sum += ME
                    end
                end
            end
            Sum = Sum / Float64(j_h + 1)
            pRho[a_h,a_g] -= Sum
            if a_h > a_g
                pRho[a_g,a_h] -= Sum
            end
        # nRho particles iteration ...
        elseif Type == 3
            p = Pairs.a[Index]
            q = Pairs.b[Index]
            a_p, l_p, j_p = OBDM_nParticle(Particle,p)
            a_q = Particle.n[q].a
            Sum = 0.0
            @inbounds for h in 1:N_Hole.n
                a_h, l_h, j_h = OBDM_nHole(Hole,h)
                Amp_ph, Amp_qh = OBDM_D(a_p,a_h,nRho_Old), OBDM_D(a_q,a_h,nRho_Old)
                P = rem(l_p + l_h,2) + 1
                @inbounds for J in div(abs(j_p - j_h),2):div(j_p + j_h,2)
                    J_ind = J + 1
                    N_ph = N_nu[J_ind,P]
                    ph, qh = ParticleHole.n[J_ind,P,p,h], ParticleHole.n[J_ind,P,q,h]
                    @inbounds for nu in 1:N_ph
                        Y_ph_nu = OBDM_Y(Y_RPA,J_ind,P,ph,nu)
                        Y_qh_nu = OBDM_Y(Y_RPA,J_ind,P,qh,nu)
                        @inbounds for mu in 1:N_ph
                            Y_qh_mu = OBDM_Y(Y_RPA,J_ind,P,qh,mu)
                            ME = Float64(2*J + 1) * OBDM_ME_Inner(nM,J_ind,P,nu,mu,Amp_ph,Amp_qh,Y_ph_nu,Y_qh_mu)
                            Sum += ME
                            if a_p > a_q
                                Y_ph_mu = OBDM_Y(Y_RPA,J_ind,P,ph,mu)
                                ME = Float64(2*J + 1) * OBDM_ME_Inner(nM,J_ind,P,mu,nu,Amp_ph,Amp_qh,Y_ph_mu,Y_qh_nu)
                                Sum += ME
                            end
                        end
                        ME = Float64(2*J + 1) * OBDM_ME_Outer(Amp_ph,Amp_qh,Y_ph_nu,Y_qh_nu)
                        Sum += ME
                    end
                end
            end
            Sum = Sum / Float64(j_p + 1)
            nRho[a_p,a_q] += Sum
            if a_p > a_q
                nRho[a_q,a_p] += Sum
            end
        # nRho holes iteration ...
        elseif Type == 4
            h = Pairs.a[Index]
            g = Pairs.b[Index]
            a_h, l_h, j_h = OBDM_nHole(Hole,h)
            a_g = Hole.n[g].a
            Sum = 0.0
            @inbounds for p in 1:N_Particle.n
                a_p, l_p, j_p = OBDM_nParticle(Particle,p)
                Amp_ph, Amp_pg = OBDM_D(a_p,a_h,nRho_Old), OBDM_D(a_p,a_g,nRho_Old)
                P = rem(l_p + l_h,2) + 1
                @inbounds for J in div(abs(j_p - j_h),2):div(j_p + j_h,2)
                    J_ind = J + 1
                    N_ph = N_nu[J_ind,P]
                    ph, pg = ParticleHole.n[J_ind,P,p,h], ParticleHole.n[J_ind,P,p,g]
                    @inbounds for nu in 1:N_ph
                        Y_ph_nu = OBDM_Y(Y_RPA,J_ind,P,ph,nu)
                        Y_pg_nu = OBDM_Y(Y_RPA,J_ind,P,pg,nu)
                        @inbounds for mu in 1:N_ph
                            Y_pg_mu = OBDM_Y(Y_RPA,J_ind,P,pg,mu)
                            ME = Float64(2*J + 1) * OBDM_ME_Inner(nM,J_ind,P,nu,mu,Amp_ph,Amp_pg,Y_ph_nu,Y_pg_mu)
                            Sum += ME
                            if a_h > a_g
                                Y_ph_mu = OBDM_Y(Y_RPA,J_ind,P,pg,mu)
                                ME = Float64(2*J + 1) * OBDM_ME_Inner(nM,J_ind,P,mu,nu,Amp_ph,Amp_pg,Y_ph_mu,Y_pg_nu)
                                Sum += ME
                            end
                        end
                        ME = Float64(2*J + 1) * OBDM_ME_Outer(Amp_ph,Amp_pg,Y_ph_nu,Y_pg_nu)
                        Sum += ME
                    end
                end
            end
            Sum = Sum / Float64(j_h + 1)
            nRho[a_h,a_g] -= Sum
            if a_h > a_g
                nRho[a_g,a_h] -= Sum
            end
        end
    end

    return pRho, nRho
end

function HF_ERPA_OBDM_Eigen_Reorder(a_max::Int64,N::Vector{Float64},U::Matrix{Float64})
    # Hungarian Algorithm based reordering of the HF-ERPA NO basis ...
    # Initialize the cost matrix for Hungarian algorithm ...
    C = ones(Float64,a_max,a_max) .- abs.(U)

    # Apply the Hungarian algorithm ....
    Reordering, Temp = hungarian(C)

    # Apply reordering ...
    U = U[:, Reordering]
    N = N[Reordering]

    @inbounds for a in 1:a_max
        if U[a,a] < 0.0
            U[:,a] .*= -1.0
        end
    end

    return N, U
end

function HF_ERPA_OBDM_Occupation_Check(N::pnVector,N_Particle::pnInteger,Particle::pnSVector,N_Hole::pnInteger,Hole::pnSVector)
    # Check proton s.p. orbits ...
    @inbounds for p in 1:N_Particle.p
        a_p = Particle.p[p].a
        @inbounds for h in 1:N_Hole.p
            a_h = Hole.p[h].a
            D_ph = N.p[a_h] - N.p[a_p]
            if D_ph < 0.0
                return false
            end
        end
    end
    # Check neutron s.p. orbits ...
    @inbounds for p in 1:N_Particle.n
        a_p = Particle.n[p].a
        @inbounds for h in 1:N_Hole.n
            a_h = Hole.n[h].a
            D_ph = N.n[a_h] - N.n[a_p]
            if D_ph < 0.0
                return false
            end
        end
    end
    return true
end


@inline function OBDM_D(a_p::Int64,a_h::Int64,Rho::Matrix{Float64})
    return Rho[a_h,a_h] - Rho[a_p,a_p]
end

@inline function OBDM_M(Amp_qg::Float64,X_qg_nu::ComplexF64,X_qg_mu::ComplexF64)
    return Float64(-0.5 * Amp_qg) * X_qg_mu * conj(X_qg_nu)
end

@inline function OBDM_ME_Inner(M::Matrix{Matrix{ComplexF64}},J_ind::Int64,P::Int64,nu::Int64,mu::Int64,Amp_ph::Float64,Amp_qg::Float64,Y_ph_nu::ComplexF64,Y_qg_mu::ComplexF64)
    return sqrt(Amp_ph * Amp_qg) * real(M[J_ind,P][nu,mu] * Y_ph_nu * conj(Y_qg_mu))
end

@inline function OBDM_ME_Outer(Amp_ph::Float64,Amp_qg::Float64,Y_ph_nu::ComplexF64,Y_qg_nu::ComplexF64)
    return sqrt(Amp_ph * Amp_qg) * real(Y_ph_nu * conj(Y_qg_nu))
end

@inline function OBDM_X(X_RPA::Matrix{Matrix{ComplexF64}},J_ind::Int64,P::Int64,ph::Int64,nu::Int64)
    return (@view X_RPA[J_ind,P][ph,nu])[1]
end

@inline function OBDM_Y(Y_RPA::Matrix{Matrix{ComplexF64}},J_ind::Int64,P::Int64,ph::Int64,nu::Int64)
    return (@view Y_RPA[J_ind,P][ph,nu])[1]
end

@inline function OBDM_pParticle(Particle::pnSVector,p::Int64)
    return Particle.p[p].a, Particle.p[p].l, Particle.p[p].j
end

@inline function OBDM_pHole(Hole::pnSVector,h::Int64)
    return Hole.p[h].a, Hole.p[h].l, Hole.p[h].j
end

@inline function OBDM_nParticle(Particle::pnSVector,p::Int64)
    return Particle.n[p].a, Particle.n[p].l, Particle.n[p].j
end

@inline function OBDM_nHole(Hole::pnSVector,h::Int64)
    return Hole.n[h].a, Hole.n[h].l, Hole.n[h].j
end
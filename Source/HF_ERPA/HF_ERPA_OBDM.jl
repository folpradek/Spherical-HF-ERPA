function HF_ERPA_OBDM(Params::Vector{Any},N_nu::Matrix{Int64},N_Particle::pnInteger,Particle::pnSVector,N_Hole::pnInteger,Hole::pnSVector,Orb::Vector{NOrb},X_RPA::Matrix{Matrix{ComplexF64}},Y_RPA::Matrix{Matrix{ComplexF64}})
    # Read parameters ...
    N_max = Params[4]
    N_2max = 2*N_max
    J_max = N_2max + 1
    a_max = div((N_max+1)*(N_max+2),2)

    # Make Particle-Hole orbitals list ...
    ParticleHole = Make_ParticleHole_List(J_max,N_Particle,Particle,N_Hole,Hole)
    
    # Make HF density matrices ...
    pRho_HF, nRho_HF = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)
    @inbounds for a in 1:a_max
        pRho_HF[a,a] = Orb[a].pO
        nRho_HF[a,a] = Orb[a].nO
    end

    # Initialize new density matrices ...
    pRho, nRho = deepcopy(pRho_HF), deepcopy(nRho_HF)
    pRho_old, nRho_old = deepcopy(pRho_HF), deepcopy(nRho_HF)

    Delta = 1.0
    Eps = 1e-6
    N_Iteration = 100
    Iteration = 0

    JP_List = JP_Ini(J_max)

    println("\nStarting HF ERPA OBDM calculation ...")

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

    @time while (Delta > Eps) && (N_Iteration > Iteration)

        pRho, nRho = deepcopy(pRho_HF), deepcopy(nRho_HF)

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
                                Amp_qg = OBDM_D(a_q,a_g,pRho_old)
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
                                Amp_qg = OBDM_D(a_q,a_g,nRho_old)
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

        Sum = zeros(Float64, Threads.nthreads())

        # pRho iteration ...
        @inbounds for p in 1:N_Particle.p
            a_p, l_p, j_p = OBDM_pParticle(Particle,p)
            @inbounds for h in 1:N_Hole.p
                a_h, l_h, j_h = OBDM_pHole(Hole,h)
                P = rem(l_p + l_h,2) + 1
                @inbounds for J in div(abs(j_p - j_h),2):div(j_p + j_h,2)
                    J_ind = J + 1
                    N_ph = N_nu[J_ind,P]
                    ph = ParticleHole.p[J_ind,P,p,h]
                    Amp_ph = Float64(2*J + 1) * OBDM_D(a_p,a_h,pRho_old)
                    fill!(Sum,0.0)
                    @inbounds Threads.@threads for nu in 1:N_ph
                        thread_id = Threads.threadid()
                        threadSum = 0.0
                        Y_ph_nu = OBDM_Y(Y_RPA,J_ind,P,ph,nu)
                        @inbounds for mu in 1:N_ph
                            Y_ph_mu = OBDM_Y(Y_RPA,J_ind,P,ph,mu)
                            ME = OBDM_ME_Inner(pM,J_ind,P,nu,mu,Amp_ph,Y_ph_nu,Y_ph_mu)
                            threadSum += ME
                        end
                        ME = OBDM_ME_Outer(Amp_ph,Y_ph_nu)
                        threadSum += ME
                        Sum[thread_id] += threadSum
                    end
                    
                    ME = sum(Sum)
                    pRho[a_h,a_h] -= ME / Float64(j_h + 1)
                    pRho[a_p,a_p] += ME / Float64(j_p + 1)
                end
            end
        end

        Sum = zeros(Float64, Threads.nthreads())

        # nRho iteration ...
        @inbounds for p in 1:N_Particle.n
            a_p, l_p, j_p = OBDM_nParticle(Particle,p)
            @inbounds for h in 1:N_Hole.n
                a_h, l_h, j_h = OBDM_nHole(Hole,h)
                P = rem(l_p + l_h,2) + 1
                @inbounds for J in div(abs(j_p - j_h),2):div(j_p + j_h,2)
                    J_ind = J + 1
                    N_ph = N_nu[J_ind,P]
                    ph = ParticleHole.n[J_ind,P,p,h]
                    Amp_ph = Float64(2*J + 1) * OBDM_D(a_p,a_h,nRho_old)
                    fill!(Sum,0.0)
                    @inbounds Threads.@threads for nu in 1:N_ph
                        thread_id = Threads.threadid()
                        threadSum = 0.0
                        Y_ph_nu = OBDM_Y(Y_RPA,J_ind,P,ph,nu)
                        @inbounds for mu in 1:N_ph
                            Y_ph_mu = OBDM_Y(Y_RPA,J_ind,P,ph,mu)
                            ME = OBDM_ME_Inner(nM,J_ind,P,nu,mu,Amp_ph,Y_ph_nu,Y_ph_mu)
                            threadSum += ME
                        end
                        ME = OBDM_ME_Outer(Amp_ph,Y_ph_nu)
                        threadSum += ME
                        Sum[thread_id] += threadSum
                    end
                    ME = sum(Sum)
                    nRho[a_h,a_h] -= ME / Float64(j_h + 1)
                    nRho[a_p,a_p] += ME / Float64(j_p + 1)
                end
            end
        end

        Sum = 0.0

        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                ME = (abs(pRho[a,b] - pRho_old[a,b]) + abs(nRho[a,b] - nRho_old[a,b])) / Float64(2*a_max^2)
                Sum += ME
            end
        end

        Delta = Sum
        Iteration += 1

        println("\nOBDM iteration:   Iteration Number = " * string(Iteration) * ",\tDelta = " * string(Delta))

        pRho_old, nRho_old = deepcopy(pRho), deepcopy(nRho)

    end

    println("\nERPA OBDM iteration has terminated ...")

    return pRho, nRho
end

@inline function OBDM_D(a_p::Int64,a_h::Int64,Rho::Matrix{Float64})
    return Rho[a_h,a_h] - Rho[a_p,a_p]
end

@inline function OBDM_M(Amp_qg::Float64,X_qg_nu::ComplexF64,X_qg_mu::ComplexF64)
    return Float64(-0.5 * Amp_qg) * X_qg_mu * conj(X_qg_nu)
end

@inline function OBDM_ME_Inner(M::Matrix{Matrix{ComplexF64}},J_ind::Int64,P::Int64,nu::Int64,mu::Int64,Amp_ph::Float64,Y_ph_nu::ComplexF64,Y_ph_mu::ComplexF64)
    return Amp_ph * real(M[J_ind,P][nu,mu] * Y_ph_nu * conj(Y_ph_mu))
end

@inline function OBDM_ME_Outer(Amp_ph::Float64,Y_ph_nu::ComplexF64)
    return Amp_ph * abs2(Y_ph_nu)
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
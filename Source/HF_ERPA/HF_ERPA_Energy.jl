function HF_ERPA_Energy(Params::Vector{Any},N_nu::Matrix{Int64},N_Particle::pnInteger,Particle::pnSVector,N_Hole::pnInteger,Hole::pnSVector,E_RPA::Matrix{Vector{ComplexF64}},Y_RPA::Matrix{Matrix{ComplexF64}})
    # Read calculation parameters ...
    N_max = Params[4]
    N_2max = 2*N_max
    J_max = N_2max + 1

    println("\nCalculating ERPA correlation energy ...")

    # Make Particle-Hole orbitals list ...
    ParticleHole = Make_ParticleHole_List(J_max,N_Particle,Particle,N_Hole,Hole)

    E_0 = 0.0

    # Proton ph iteration ...
    @inbounds for h in 1:N_Hole.p
        l_h = Hole.p[h].l
        j_h = Hole.p[h].j
        @inbounds for p in 1:N_Particle.p
            l_p = Particle.p[p].l
            j_p = Particle.p[p].j
            P = rem(l_p + l_h,2) + 1
            @inbounds for J in div(abs(j_p - j_h),2):div(j_p + j_h,2)
                N_ph = N_nu[J+1,P]
                ph = ParticleHole.p[J+1,P,p,h]
                Amp_J = Float64(2*J + 1)
                @inbounds for nu in 1:N_ph
                    Y_ph_nu = Y_RPA[J+1,P][ph,nu]
                    ME = -Amp_J * real(E_RPA[J+1,P][nu]) * abs(Y_ph_nu)^2
                    E_0 += ME
                end
            end
        end
    end

    # Neutron ph iteration ...
    @inbounds for h in 1:N_Hole.n
        l_h = Hole.n[h].l
        j_h = Hole.n[h].j
        @inbounds for p in 1:N_Particle.n
            l_p = Particle.n[p].l
            j_p = Particle.n[p].j
            P = rem(l_p + l_h,2) + 1
            @inbounds for J in div(abs(j_p - j_h),2):div(j_p + j_h,2)
                N_ph = N_nu[J+1,P]
                ph = ParticleHole.n[J+1,P,p,h]
                Amp_J = Float64(2*J + 1)
                @inbounds for nu in 1:N_ph
                    Y_ph_nu = Y_RPA[J+1,P][ph,nu]
                    ME = -Amp_J * real(E_RPA[J+1,P][nu]) * abs(Y_ph_nu)^2
                    E_0 += ME
                end
            end
        end
    end

    println("\nERPA ground-state correlation energy E_0^ERPA is ...     E_0 = " * string(round(E_0,digits=6)) * "\tMeV")
    println("\nTo get the total ground-state energy add the mean-field energy ...")

    return E_0
end
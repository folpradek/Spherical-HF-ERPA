function HF_ERPA_Energy(Params::Parameters,N_nu::Matrix{Int64},N_Particle::pnInteger,Particle::pnSVector,N_Hole::pnInteger,Hole::pnSVector,E_RPA::Matrix{Vector{ComplexF64}},Y_RPA::Matrix{Matrix{ComplexF64}})
    # Read calculation parameters ...
    N_max = Params.Calc.Nmax
    N_2max = 2*N_max
    J_max = N_2max + 1

    println("\nCalculating ERPA correlation energy ...")

    # Make Particle-Hole orbitals list ...
    ParticleHole = Make_ParticleHole_List(J_max,N_Particle,Particle,N_Hole,Hole)

    E_0 = 0.0

    # Proton ph iteration ...
    @inbounds for h in 1:N_Hole.p
        a_h = Hole.p[h].a
        l_h = Hole.p[h].l
        j_h = Hole.p[h].j
        @inbounds for p in 1:N_Particle.p
            a_p = Particle.p[p].a
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
        a_h = Hole.n[h].a
        l_h = Hole.n[h].l
        j_h = Hole.n[h].j
        @inbounds for p in 1:N_Particle.n
            a_p = Particle.n[p].a
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


# Evaluation of ERPA correlation energy ... for testing purposes!!!
function HF_ERPA_Energy_Alt(Params::Parameters,N_nu::Matrix{Int64},N_Particle::pnInteger,Particle::pnSVector,N_Hole::pnInteger,Hole::pnSVector,Orb::Vector{NOrb},H::pnMatrix,Orb_NN::NNOrb,VNN::NNInt,A::Matrix{Matrix{Float64}},E_RPA::Matrix{Vector{ComplexF64}},Y_RPA::Matrix{Matrix{ComplexF64}},Rho::pnMatrix)
    # Read calculation parameters ...
    N_max = Params.Calc.Nmax
    N_2max = 2*N_max
    J_max = N_2max + 1

    println("\nCalculating ERPA correlation energy ... Using the alternative 2-body density based formula ...")

    # Make Particle-Hole orbitals list ...
    ParticleHole = Make_ParticleHole_List(J_max,N_Particle,Particle,N_Hole,Hole)

    E_0 = 0.0
    E_0_Fermi = 0.0

    # Proton ph iteration ...
    @inbounds for h in 1:N_Hole.p
        a_h = Hole.p[h].a
        l_h = Hole.p[h].l
        j_h = Hole.p[h].j
        @inbounds for p in 1:N_Particle.p
            a_p = Particle.p[p].a
            l_p = Particle.p[p].l
            j_p = Particle.p[p].j
            P = rem(l_p + l_h,2) + 1
            @inbounds for J in div(abs(j_p - j_h),2):div(j_p + j_h,2)
                N_ph = N_nu[J+1,P]
                ph = ParticleHole.p[J+1,P,p,h]
                Amp_J = 0.5 * Float64(2*J + 1)

                @inbounds for nu in 1:N_ph
                    Y_ph_nu = Y_RPA[J+1,P][ph,nu]
                    ME = - Amp_J * real(E_RPA[J+1,P][nu]) * abs(Y_ph_nu)^2
                    E_0 += ME

                    @inbounds for g in 1:N_Hole.p
                        a_g = Hole.p[g].a
                        l_g = Hole.p[g].l
                        j_g = Hole.p[g].j
                        #if l_h == l_g && j_h == j_g
                            @inbounds for q in 1:N_Particle.p
                                a_q = Particle.p[q].a
                                l_q = Particle.p[q].l
                                j_q = Particle.p[q].j
                                if (rem(l_p + l_h,2) ==  rem(l_q + l_g,2)) && (abs(j_q - j_g) <= 2*J) && (2*J <= (j_q + j_g))
                                #if l_p == l_q && j_p == j_q
                                    qg = ParticleHole.p[J+1,P,q,g]
                                    Y_qg_nu = Y_RPA[J+1,P][qg,nu]
                                    #ME = - 0.5 * Amp_J * (sqrt((Rho.p[a_h,a_h] - Rho.p[a_p,a_p]) / (Rho.p[a_g,a_g] - Rho.p[a_q,a_q])) +
                                    #        sqrt((Rho.p[a_g,a_g] - Rho.p[a_q,a_q]) / (Rho.p[a_h,a_h] - Rho.p[a_p,a_p]))) * (H.p[a_p,a_q] *
                                    #        KroneckerDelta(a_h,a_g) - H.p[a_h,a_g] * KroneckerDelta(a_p,a_q)) * Y_ph_nu * Y_qg_nu
                                    ME = - Amp_J * A[J+1,P][ph,qg] * Y_ph_nu * Y_qg_nu
                                    E_0 += ME
                                end
                            end
                        #end
                    end

                end

                #ME = Amp_J * Rho.p[a_p,a_p] * V2B(a_p,a_h,a_p,a_h,J,1,VNN.pp,Orb,Orb_NN)
                #E_0_Fermi += ME
            end
        end
    end

    # Neutron ph iteration ...
    @inbounds for h in 1:N_Hole.n
        a_h = Hole.n[h].a
        l_h = Hole.n[h].l
        j_h = Hole.n[h].j
        @inbounds for p in 1:N_Particle.n
            a_p = Particle.n[p].a
            l_p = Particle.n[p].l
            j_p = Particle.n[p].j
            P = rem(l_p + l_h,2) + 1
            @inbounds for J in div(abs(j_p - j_h),2):div(j_p + j_h,2)
                N_ph = N_nu[J+1,P]
                ph = ParticleHole.n[J+1,P,p,h]
                Amp_J = 0.5 * Float64(2*J + 1)

                @inbounds for nu in 1:N_ph
                    Y_ph_nu = Y_RPA[J+1,P][ph,nu]
                    ME = - Amp_J * real(E_RPA[J+1,P][nu]) * abs(Y_ph_nu)^2
                    E_0 += ME

                    @inbounds for g in 1:N_Hole.n
                        a_g = Hole.n[g].a
                        l_g = Hole.n[g].l
                        j_g = Hole.n[g].j
                        #if l_h == l_g && j_h == j_g
                            @inbounds for q in 1:N_Particle.n
                                a_q = Particle.n[q].a
                                l_q = Particle.n[q].l
                                j_q = Particle.n[q].j
                                if (rem(l_p + l_h,2) ==  rem(l_q + l_g,2)) && (abs(j_q - j_g) <= 2*J) && (2*J <= (j_q + j_g))
                                #if l_p == l_q && j_p == j_q
                                    qg = ParticleHole.n[J+1,P,q,g]
                                    Y_qg_nu = Y_RPA[J+1,P][qg,nu]
                                    #ME = - 0.5 * Amp_J * (sqrt((Rho.n[a_h,a_h] - Rho.n[a_p,a_p]) / (Rho.n[a_g,a_g] - Rho.n[a_q,a_q])) +
                                    #        sqrt((Rho.n[a_g,a_g] - Rho.n[a_q,a_q]) / (Rho.n[a_h,a_h] - Rho.n[a_p,a_p]))) * (H.n[a_p,a_q] *
                                    #        KroneckerDelta(a_h,a_g) - H.n[a_h,a_g] * KroneckerDelta(a_p,a_q)) * Y_ph_nu * Y_qg_nu
                                    ME = - Amp_J * A[J+1,P][ph,qg] * Y_ph_nu * Y_qg_nu
                                    E_0 += ME
                                end
                            end
                        #end
                    end

                end

                #ME = Amp_J * Rho.n[a_p,a_p] * V2B(a_p,a_h,a_p,a_h,J,1,VNN.nn,Orb,Orb_NN)
                #E_0_Fermi += ME

            end
        end
    end

    println("ERPA ground-state correlation energy E_0^ERPA is ...     E_0_corr = " * string(round(real(E_0 + E_0_Fermi),digits=6)) * "\tMeV")
    #println("ERPA ground-state correlation energy due to Fermi sea depletion is ...     E_0_corr_dep = " * string(round(real(E_0_Fermi),digits=6)) * "\tMeV")
    println("To get the total ground-state energy add the mean-field energy ...")
    
    return
end

function HF_ERPA_Energy_Bosonic(Params::Parameters,N_nu::Matrix{Int64},N_Particle::pnInteger,Particle::pnSVector,N_Hole::pnInteger,Hole::pnSVector,E_RPA::Matrix{Vector{ComplexF64}},Y_RPA::Matrix{Matrix{ComplexF64}})
    # Read calculation parameters ...
    N_max = Params.Calc.Nmax
    N_2max = 2*N_max
    J_max = N_2max + 1

    println("\nCalculating ERPA correlation energy ... Using the bosonic algebra formula ...")

    # Make Particle-Hole orbitals list ...
    ParticleHole = Make_ParticleHole_List(J_max,N_Particle,Particle,N_Hole,Hole)

    E_0 = 0.0

    # Proton ph iteration ...
    @inbounds for h in 1:N_Hole.p
        a_h = Hole.p[h].a
        l_h = Hole.p[h].l
        j_h = Hole.p[h].j
        @inbounds for p in 1:N_Particle.p
            a_p = Particle.p[p].a
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
        a_h = Hole.n[h].a
        l_h = Hole.n[h].l
        j_h = Hole.n[h].j
        @inbounds for p in 1:N_Particle.n
            a_p = Particle.n[p].a
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

    println("ERPA ground-state correlation energy E_0^ERPA is ...     E_0 = " * string(round(E_0,digits=6)) * "\tMeV")
    println("To get the total ground-state energy add the mean-field energy ...")

    return
end
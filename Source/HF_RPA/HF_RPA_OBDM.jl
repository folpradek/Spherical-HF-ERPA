function HF_RPA_OBDM(Params::Vector{Any},Orb::Vector{NOrb},Orb_Phonon::Matrix{Vector{Int64}},Phonon::Vector{PhState},Particle::pnSVector,Hole::pnSVector,N_nu::Matrix{Int64},Y_RPA::Matrix{Matrix{ComplexF64}})
    # Read calculation parameters ..
    N_max = Params[4]
    N_2max = 2*N_max
    J_max = N_2max + 1
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Initialize density matrices ...
    pRho, nRho = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)

    @inbounds for a in 1:a_max
        pRho[a,a] = Orb[a].pO
        nRho[a,a] = Orb[a].nO
    end

    # Begin OBDM calculation ...
    @inbounds for J in 0:J_max
        @inbounds for P in 1:2
            N_ph = N_nu[J+1,P]
            @inbounds Threads.@threads for ph in 1:N_ph
                Ind_ph = Orb_Phonon[J+1,P][ph]
                p,h,t_ph,J_ph,P_ph,a_p,l_p,j_p,E_p,a_h,l_h,j_h,E_h = Phonon_Ind(Ind_ph,Phonon,Particle,Hole)
                Sum = 0.0
                @inbounds for nu in 1:N_ph
                    Y_norm = @views norm(Y_RPA[J+1,P][:,nu])^2
                    if abs(Y_norm) < 0.5
                        y = abs2(Y_RPA[J+1,P][ph,nu])
                        ME = y * Float64(2*J + 1)
                        Sum += ME
                    end
                end
                if t_ph == -1
                    pRho[a_h,a_h] -= Sum / Float64(j_h + 1)
                    pRho[a_p,a_p] += Sum / Float64(j_p + 1)
                elseif t_ph == 1
                    nRho[a_h,a_h] -= Sum / Float64(j_h + 1)
                    nRho[a_p,a_p] += Sum / Float64(j_p + 1)
                end
            end
        end
    end

    return pRho, nRho
end
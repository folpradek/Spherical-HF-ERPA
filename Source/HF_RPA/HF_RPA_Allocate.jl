function HF_RPA_Allocate(Params::Vector{Any},N_nu::Matrix{Int64},Orb_Phonon::Matrix{Vector{Int64}},Phonon::Vector{PhState},Particle::pnSVector,Hole::pnSVector,Orb::Vector{NOrb},Orb_NN::NNOrb,VNN::NNInt)
    # Read calculation parameters ...
    N_max = Params[4]
    N_2max = 2*N_max
    J_max = N_2max + 1

    JP_List = JP_Ini(J_max)
    
    println("\nAllocating TDA & RPA matrices A & B ...")

    #Allocate A & B matrices ...
    A = Matrix{Matrix{Float64}}(undef,J_max+1,2)
    B = Matrix{Matrix{Float64}}(undef,J_max+1,2)

    # Fill A & B matrices ...
    @inbounds Threads.@threads for JP in JP_List
        J = JP[1]
        P = JP[2]
        N_ph = N_nu[J+1,P]
        A_JP = zeros(Float64,N_ph,N_ph)
        B_JP = zeros(Float64,N_ph,N_ph)

        @inbounds for ind_ph in 1:N_ph
            ph = Orb_Phonon[J+1,P][ind_ph]
            p,h,t_ph,J,P,a_p,l_p,j_p,E_p,a_h,l_h,j_h,E_h = Phonon_Ind(ph,Phonon,Particle,Hole)
    
            @inbounds for ind_qg in 1:N_ph
                qg = Orb_Phonon[J+1,P][ind_qg]
                q,g,t_qg,J_qg,P_qg,a_q,l_q,j_q,E_q,a_g,l_g,j_g,E_g = Phonon_Ind(qg,Phonon,Particle,Hole)
    
                ASum = 0.0
                BSum = 0.0

                if t_ph == t_qg && ph == qg
                    ASum += (E_p - E_h)
                end

                if rem(l_p + l_g, 2) == rem(l_q + l_h, 2)
                    @inbounds for J_r in div(abs(j_p - j_g),2):div((j_p + j_g),2)
                        if (abs(j_q - j_h) <= 2*J_r) && (2*J_r <= (j_q + j_h))
                            if t_ph == -1 && t_qg == -1
                                phase_A = (-1)^(div(j_h + j_q,2) + J_r)
                                ME = phase_A * (2*J_r + 1) * f6j(j_p, j_h, 2*J, j_q, j_g, 2*J_r) * V2B(a_p,a_g,a_h,a_q,J_r,1,VNN.pp,Orb,Orb_NN)
                                ASum +=  ME
                            elseif t_ph == -1 && t_qg == 1
                                phase_A = (-1)^(div((j_h + j_q),2) + J_r)
                                ME = phase_A * (2*J_r + 1) * f6j(j_p, j_h, 2*J, j_q, j_g, 2*J_r) * V2B(a_p,a_g,a_h,a_q,J_r,0,VNN.pn,Orb,Orb_NN)
                                ASum += ME
                            elseif t_ph == 1 && t_qg == -1
                                phase_A = (-1)^(div((j_g + j_p),2) + J_r)
                                ME = phase_A * (2*J_r + 1) * f6j(j_q, j_g, 2*J, j_p, j_h, 2*J_r) * V2B(a_q,a_h,a_g,a_p,J_r,0,VNN.pn,Orb,Orb_NN)
                                ASum += ME
                            elseif t_ph == 1 && t_qg == 1
                                phase_A = (-1)^(div((j_h + j_q),2) + J_r)
                                ME = phase_A * (2*J_r + 1) * f6j(j_p, j_h, 2*J, j_q, j_g, 2*J_r) * V2B(a_p,a_g,a_h,a_q,J_r,1,VNN.nn,Orb,Orb_NN)
                                ASum +=  ME
                            end
                        end
                    end
                end

                if rem(l_p + l_q, 2) == rem(l_h + l_g, 2)
                    @inbounds for J_r in div(abs(j_p - j_q),2):div((j_p + j_q),2)
                        if (abs(j_h - j_g) <= 2*J_r) && (2*J_r <= (j_h + j_g))

                            if t_ph == -1 && t_qg == -1
                                Phase = Float64((-1)^(div(j_h + j_q,2) + J + J_r))
                                ME = Phase * Float64(2*J_r + 1) * f6j(j_p,j_h,2*J,j_g,j_q,2*J_r) * V2B(a_p,a_q,a_h,a_g,J_r,1,VNN.pp,Orb,Orb_NN)
                                BSum += ME
                            elseif t_ph == -1 && t_qg == 1
                                Phase = Float64((-1)^(div(j_h + j_q,2) + J + J_r))
                                ME = Phase * (2*J_r + 1) * f6j(j_p,j_h,2*J,j_g,j_q,2*J_r) * V2B(a_p,a_q,a_h,a_g,J_r,0,VNN.pn,Orb,Orb_NN)
                                BSum += ME
                            elseif t_ph == 1 && t_qg == -1
                                Phase = Float64((-1)^(div(j_g + j_p,2) + J + J_r))
                                ME = Phase * (2*J_r + 1) * f6j(j_q,j_g,2*J,j_h,j_p,2*J_r) * V2B(a_q,a_p,a_g,a_h,J_r,0,VNN.pn,Orb,Orb_NN)
                                BSum += ME
                            elseif t_ph == 1 && t_qg == 1
                                Phase = Float64((-1)^(div(j_h + j_q,2) + J + J_r))
                                ME = Phase * (2*J_r + 1) * f6j(j_p,j_h,2*J,j_g,j_q,2*J_r) * V2B(a_p,a_q,a_h,a_g,J_r,1,VNN.nn,Orb,Orb_NN)
                                BSum += ME
                            end

                        end
                    end
                end
                A_JP[ind_ph,ind_qg] = ASum
                B_JP[ind_ph,ind_qg] = BSum
            end
        end

        # For given J & P allocate A & B ...
        A[J+1,P] = A_JP
        B[J+1,P] = B_JP
    end

    println("\nMatrices A & B allocated ...")

    return A, B

end
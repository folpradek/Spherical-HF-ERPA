function HF_MBPT(Params::Parameters,Orb::Vector{Orb1B},Orb_NN::Orb2B,V_NN::O2B)
    # Make Particle-Hole orbitals ...
    N_Particle, Particle, N_Hole, Hole = orbitals_ph_make(Params,Orb)

    # Evaluate the HF-MBPT(2) ground-state energy correction ...
    HF_MBPT_energy(Params,Orb,Orb_NN,N_Particle,Particle,N_Hole,Hole,V_NN)

    # Evaluate the HF-MBPT(3) ground-state OBDM correction & radial densities ...
    HF_MBPT_density(Params,Orb,Orb_NN,N_Particle,Particle,N_Hole,Hole,V_NN)

    return
end

function HF_MBPT_energy(Params::Parameters,Orb::Vector{Orb1B},Orb_NN::Orb2B,N_Particle::pnInteger,Particle::pnSVector,N_Hole::pnInteger,Hole::pnSVector,V_NN::O2B)
    # Read parameters ...
    IO = Params.Calc.Path

    # Initialize the thread accumulators ...
    dE, dE_threads = 0.0, zeros(Float64,Threads.maxthreadid())

    # Calculate the MBPT(2) energy correction ...
    println("\nCalculating HF-MBPT(2) ground-state energy correction ... ")
    
    # proton-proton contribution ...
    @inbounds Threads.@threads :static for p in 1:N_Particle.p
        Sum, Tid = 0.0, Threads.threadid()
        a_p = Particle.p[p].a
        l_p = Particle.p[p].l
        j_p = Particle.p[p].j
        E_p = Particle.p[p].E
        @inbounds for q in 1:N_Particle.p
            a_q = Particle.p[q].a
            l_q = Particle.p[q].l
            j_q = Particle.p[q].j
            E_q = Particle.p[q].E
            P = rem(l_p + l_q, 2) + 1
            @inbounds for J in div(abs(j_p - j_q),2):div((j_p + j_q),2)
                @inbounds for h in 1:N_Hole.p
                    a_h = Hole.p[h].a
                    l_h = Hole.p[h].l
                    j_h = Hole.p[h].j
                    E_h = Hole.p[h].E
                    @inbounds for g in 1:N_Hole.p
                        a_g = Hole.p[g].a
                        l_g = Hole.p[g].l
                        j_g = Hole.p[g].j
                        E_g = Hole.p[g].E
                        if (P == (rem(l_h + l_g, 2) + 1)) && (abs(j_h - j_g) <= 2*J) && (2*J <= (j_h + j_g))
                            #Sum += Float64(2*J + 1) / 4.0 * abs(V2b(a_p,a_q,a_h,a_g,J,1,V_NN.pp,Orb,Orb_NN))^2 / (E_h + E_g - E_p - E_q)
                            Sum += Float64(2*J + 1) / 4.0 * abs(O2b_pp(a_p,a_q,a_h,a_g,J,P,V_NN,Orb,Orb_NN))^2 / (E_h + E_g - E_p - E_q)
                        end
                    end
                end
            end
        end
        dE_threads[Tid] += Sum 
    end
    dE += sum(dE_threads)

    # neutron-neutron contribution ...
    dE_threads = zeros(Float64,Threads.maxthreadid())
    @inbounds Threads.@threads :static for p in 1:N_Particle.n
        Sum, Tid = 0.0, Threads.threadid()
        a_p = Particle.n[p].a
        l_p = Particle.n[p].l
        j_p = Particle.n[p].j
        E_p = Particle.n[p].E
        @inbounds for q in 1:N_Particle.n
            a_q = Particle.n[q].a
            l_q = Particle.n[q].l
            j_q = Particle.n[q].j
            E_q = Particle.n[q].E
            P = rem(l_p + l_q, 2) + 1
            @inbounds for J in div(abs(j_p - j_q),2):div((j_p + j_q),2)
                @inbounds for h in 1:N_Hole.n
                    a_h = Hole.n[h].a
                    l_h = Hole.n[h].l
                    j_h = Hole.n[h].j
                    E_h = Hole.n[h].E
                    @inbounds for g in 1:N_Hole.n
                        a_g = Hole.n[g].a
                        l_g = Hole.n[g].l
                        j_g = Hole.n[g].j
                        E_g = Hole.n[g].E
                        if (P == (rem(l_h + l_g, 2) + 1)) && (abs(j_h - j_g) <= 2*J) && (2*J <= (j_h + j_g))
                            #Sum += Float64(2*J + 1) / 4.0 * abs(V2b(a_p,a_q,a_h,a_g,J,1,V_NN.nn,Orb,Orb_NN))^2 / (E_h + E_g - E_p - E_q)
                            Sum += Float64(2*J + 1) / 4.0 * abs(O2b_nn(a_p,a_q,a_h,a_g,J,P,V_NN,Orb,Orb_NN))^2 / (E_h + E_g - E_p - E_q)
                        end
                    end
                end
            end
        end
        dE_threads[Tid] += Sum
    end
    dE += sum(dE_threads)

    # proton-neutron contribution ...
    dE_threads = zeros(Float64,Threads.maxthreadid())
    @inbounds Threads.@threads :static for p in 1:N_Particle.p
        Sum, Tid = 0.0, Threads.threadid()
        a_p = Particle.p[p].a
        l_p = Particle.p[p].l
        j_p = Particle.p[p].j
        E_p = Particle.p[p].E
        @inbounds for q in 1:N_Particle.n
            a_q = Particle.n[q].a
            l_q = Particle.n[q].l
            j_q = Particle.n[q].j
            E_q = Particle.n[q].E
            P = rem(l_p + l_q, 2) + 1
            @inbounds for J in div(abs(j_p - j_q),2):div((j_p + j_q),2)
                @inbounds for h in 1:N_Hole.p
                    a_h = Hole.p[h].a
                    l_h = Hole.p[h].l
                    j_h = Hole.p[h].j
                    E_h = Hole.p[h].E
                    @inbounds for g in 1:N_Hole.n
                        a_g = Hole.n[g].a
                        l_g = Hole.n[g].l
                        j_g = Hole.n[g].j
                        E_g = Hole.n[g].E
                        if (P == (rem(l_h + l_g, 2) + 1)) && (abs(j_h - j_g) <= 2*J) && (2*J <= (j_h + j_g))
                            #Sum += Float64(2*J + 1) * abs(V2b(a_p,a_q,a_h,a_g,J,0,V_NN.pn,Orb,Orb_NN))^2 / (E_h + E_g - E_p - E_q)
                            Sum += Float64(2*J + 1) * abs(O2b_pn(a_p,a_q,a_h,a_g,J,P,V_NN,Orb_NN))^2 / (E_h + E_g - E_p - E_q)
                        end
                    end
                end
            end
        end
        dE_threads[Tid] += Sum
    end

    # Get the total HFB-BMBPT(2) energy from the thread accumulators ...
    dE += sum(dE_threads)

    # Print the total HF-MBPT(2) energy correction ...
    println("\tHF-MBPT(2) energy    ...   E^(2) = " * string(round(dE, sigdigits=9)) * " MeV")

    # Write the total HF-MBPT(2) energy correction to the summary file ...
    Summary =  open(string("IO/", IO, "/HF/HF_Summary.dat"), "a")
        println(Summary, "\nSpherical Hartree-Fock Leading Order Many-Body Perturbation Theory solution review:")
        println(Summary, "\nE_0^(2) = " * string(round(dE, sigdigits=9)) * "\t MeV \t\t ... \t LO HF-MBPT(2) correction to ground-state energy")
        println(Summary, "\nE_0^(2) / A = " * string(round(dE / Float64(Params.Calc.A), sigdigits=9)) * "\t MeV \t\t ... \t LO HF-MBPT(2) correction to ground-state energy per nucleon")
    close(Summary)

    return
end

function HF_MBPT_density(Params::Parameters,Orb::Vector{Orb1B},Orb_NN::Orb2B,N_Particle::pnInteger,Particle::pnSVector,N_Hole::pnInteger,Hole::pnSVector,V_NN::O2B)
    # Read calculation params ...
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Initialize density matrices ...
    dpRho = zeros(Float64,a_max,a_max)
    dnRho = zeros(Float64,a_max,a_max)

    # Start MBPT(3) LO density operator corrections ...
    println("\nEvaluating density LO corrections ...")

    # pRho 1p-1h pp
    @inbounds Threads.@threads for p in 1:N_Particle.p
        a_p = Particle.p[p].a
        l_p = Particle.p[p].l
        j_p = Particle.p[p].j
        E_p = Particle.p[p].E
        @inbounds for h in 1:N_Hole.p
            a_h = Hole.p[h].a
            l_h = Hole.p[h].l
            j_h = Hole.p[h].j
            E_h = Hole.p[h].E
            if j_p == j_h && l_p == l_h
                @inbounds for f in 1:N_Hole.p
                    a_f = Hole.p[f].a
                    l_f = Hole.p[f].l
                    j_f = Hole.p[f].j
                    E_f = Hole.p[f].E
                    P_pf, P_hf = rem(l_p + l_f,2) + 1, rem(l_h + l_f,2) + 1
                    @inbounds for J in div(abs(j_p - j_f),2):div((j_p + j_f),2)
                        @inbounds for r in 1:N_Particle.p
                            a_r = Particle.p[r].a
                            l_r = Particle.p[r].l
                            j_r = Particle.p[r].j
                            E_r = Particle.p[r].E
                            @inbounds for q in 1:N_Particle.p
                                a_q = Particle.p[q].a
                                l_q = Particle.p[q].l
                                j_q = Particle.p[q].j
                                E_q = Particle.p[q].E
                                if (abs(j_q - j_r) <= 2*J) && (2*J <= (j_q + j_r)) && (rem(l_q + l_r, 2) == rem(l_h + l_f, 2))
                                    ME = 0.5 * Float64(2*J + 1) / Float64(j_p + 1) * O2b_pp(a_p,a_f,a_q,a_r,J,P_pf,V_NN,Orb,Orb_NN) * 
                                        O2b_pp(a_q,a_r,a_h,a_f,J,P_hf,V_NN,Orb,Orb_NN) / (E_h - E_p) / (E_h + E_f - E_q - E_r)
                                    dpRho[a_p,a_h] += ME
                                    dpRho[a_h,a_p] += ME
                                end
                            end
                        end
                    end
                end

                @inbounds for r in 1:N_Particle.p
                    a_r = Particle.p[r].a
                    l_r = Particle.p[r].l
                    j_r = Particle.p[r].j
                    E_r = Particle.p[r].E
                    P_pr, P_hr = rem(l_p + l_r,2) + 1, rem(l_p + l_r,2) + 1
                    @inbounds for J in div(abs(j_h - j_r),2):div((j_h + j_r),2)
                        @inbounds for f in 1:N_Hole.p
                            a_f = Hole.p[f].a
                            l_f = Hole.p[f].l
                            j_f = Hole.p[f].j
                            E_f = Hole.p[f].E
                            @inbounds for g in 1:N_Hole.p
                                a_g = Hole.p[g].a
                                l_g = Hole.p[g].l
                                j_g = Hole.p[g].j
                                E_g = Hole.p[g].E
                                if (abs(j_g - j_f) <= 2*J) && (2*J <= (j_g + j_f)) && (rem(l_p + l_r, 2) == rem(l_g + l_f, 2))
                                    ME = - 0.5 * Float64(2*J + 1) / Float64(j_p + 1) * O2b_pp(a_p,a_r,a_g,a_f,J,P_pr,V_NN,Orb,Orb_NN) * 
                                            O2b_pp(a_g,a_f,a_h,a_r,J,P_hr,V_NN,Orb,Orb_NN) / (E_h - E_p) / (E_g + E_f - E_p - E_r)
                                    dpRho[a_p,a_h] += ME
                                    dpRho[a_h,a_p] += ME
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    # pRho 1p-1h pn
    @inbounds Threads.@threads for p in 1:N_Particle.p
        a_p = Particle.p[p].a
        l_p = Particle.p[p].l
        j_p = Particle.p[p].j
        E_p = Particle.p[p].E
        @inbounds for h in 1:N_Hole.p
            a_h = Hole.p[h].a
            l_h = Hole.p[h].l
            j_h = Hole.p[h].j
            E_h = Hole.p[h].E
            if j_p == j_h && l_p == l_h

                @inbounds for f in 1:N_Hole.n
                    a_f = Hole.n[f].a
                    l_f = Hole.n[f].l
                    j_f = Hole.n[f].j
                    E_f = Hole.n[f].E
                    P_pf, P_hf = rem(l_p + l_f,2) + 1, rem(l_h + l_f,2) + 1
                    @inbounds for J in div(abs(j_p - j_f),2):div((j_p + j_f),2)
                        @inbounds for r in 1:N_Particle.n
                            a_r = Particle.n[r].a
                            l_r = Particle.n[r].l
                            j_r = Particle.n[r].j
                            E_r = Particle.n[r].E
                            @inbounds for q in 1:N_Particle.p
                                a_q = Particle.p[q].a
                                l_q = Particle.p[q].l
                                j_q = Particle.p[q].j
                                E_q = Particle.p[q].E
                                if (abs(j_q - j_r) <= 2*J) && (2*J <= (j_q + j_r)) && (rem(l_q + l_r, 2) == rem(l_h + l_f, 2))
                                    ME = Float64(2*J + 1) / Float64(j_p + 1) * O2b_pn(a_p,a_f,a_q,a_r,J,P_pf,V_NN,Orb_NN) * 
                                            O2b_pn(a_q,a_r,a_h,a_f,J,P_hf,V_NN,Orb_NN) / (E_h - E_p) / (E_h + E_f - E_q - E_r)
                                    dpRho[a_p,a_h] += ME
                                    dpRho[a_h,a_p] += ME
                                end
                            end
                        end
                    end
                end

                @inbounds for r in 1:N_Particle.n
                    a_r = Particle.n[r].a
                    l_r = Particle.n[r].l
                    j_r = Particle.n[r].j
                    E_r = Particle.n[r].E
                    P_pr, P_hr = rem(l_p + l_r,2) + 1, rem(l_h + l_r,2) + 1
                    @inbounds for J in div(abs(j_h - j_r),2):div((j_h + j_r),2)
                        @inbounds for f in 1:N_Hole.n
                            a_f = Hole.n[f].a
                            l_f = Hole.n[f].l
                            j_f = Hole.n[f].j
                            E_f = Hole.n[f].E
                            @inbounds for g in 1:N_Hole.p
                                a_g = Hole.p[g].a
                                l_g = Hole.p[g].l
                                j_g = Hole.p[g].j
                                E_g = Hole.p[g].E
                                if (abs(j_g - j_f) <= 2*J) && (2*J <= (j_g + j_f)) && (rem(l_p + l_r, 2) == rem(l_g + l_f, 2))
                                    ME = - Float64(2*J + 1) / Float64(j_p + 1) * O2b_pn(a_p,a_r,a_g,a_f,J,P_pr,V_NN,Orb_NN) * 
                                            O2b_pn(a_g,a_f,a_h,a_r,J,P_hr,V_NN,Orb_NN) / (E_h - E_p) / (E_g + E_f - E_p - E_r)
                                    dpRho[a_p,a_h] += ME
                                    dpRho[a_h,a_p] += ME
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    # pRho 2p-2h pp
    @inbounds Threads.@threads for p in 1:N_Particle.p
        a_p = Particle.p[p].a
        l_p = Particle.p[p].l
        j_p = Particle.p[p].j
        E_p = Particle.p[p].E
        @inbounds for q in 1:N_Particle.p
            a_q = Particle.p[q].a
            l_q = Particle.p[q].l
            j_q = Particle.p[q].j
            E_q = Particle.p[q].E
            P_pq = rem(l_p + l_q,2) + 1
            @inbounds for J in div(abs(j_p - j_q),2):div((j_p + j_q),2)
                @inbounds for h in 1:N_Hole.p
                    a_h = Hole.p[h].a
                    l_h = Hole.p[h].l
                    j_h = Hole.p[h].j
                    E_h = Hole.p[h].E
                    @inbounds for g in 1:N_Hole.p
                        a_g = Hole.p[g].a
                        l_g = Hole.p[g].l
                        j_g = Hole.p[g].j
                        E_g = Hole.p[g].E
                        P_hg = rem(l_h + l_g,2) + 1
                        if (abs(j_h - j_g) <= 2*J) && (2*J <= (j_h + j_g)) && (rem(l_h + l_g, 2) == rem(l_p + l_q, 2))
                            @inbounds for a in 1:N_Particle.p
                                a_a = Particle.p[a].a
                                l_a = Particle.p[a].l
                                j_a = Particle.p[a].j
                                E_a = Particle.p[a].E
                                if (rem(l_h + l_g, 2) == rem(l_q + l_a, 2)) && (j_a == j_p) && (l_a == l_p)
                                    ME = 0.5 * Float64(2*J + 1) / Float64(j_p + 1) * O2b_pp(a_p,a_q,a_h,a_g,J,P_pq,V_NN,Orb,Orb_NN) *
                                            O2b_pp(a_h,a_g,a_a,a_q,J,P_hg,V_NN,Orb,Orb_NN) / (E_h + E_g - E_p - E_q) / (E_h + E_g - E_a - E_q)
                                    dpRho[a_p,a_a] += ME
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    @inbounds Threads.@threads for h in 1:N_Hole.p
        a_h = Hole.p[h].a
        l_h = Hole.p[h].l
        j_h = Hole.p[h].j
        E_h = Hole.p[h].E
        @inbounds for g in 1:N_Hole.p
            a_g = Hole.p[g].a
            l_g = Hole.p[g].l
            j_g = Hole.p[g].j
            E_g = Hole.p[g].E
            P_hg = rem(l_h + l_g,2) + 1
            @inbounds for J in div(abs(j_h - j_g),2):div((j_h + j_g),2)
                @inbounds for p in 1:N_Particle.p
                    a_p = Particle.p[p].a
                    l_p = Particle.p[p].l
                    j_p = Particle.p[p].j
                    E_p = Particle.p[p].E
                    @inbounds for q in 1:N_Particle.p
                        a_q = Particle.p[q].a
                        l_q = Particle.p[q].l
                        j_q = Particle.p[q].j
                        E_q = Particle.p[q].E
                        P_pq = rem(l_p + l_q,2) + 1
                        if (abs(j_p - j_q) <= 2*J) && (2*J <= (j_p + j_q)) && (rem(l_h + l_g, 2) == rem(l_p + l_q, 2))
                            @inbounds for c in 1:N_Hole.p
                                a_c = Hole.p[c].a
                                l_c = Hole.p[c].l
                                j_c = Hole.p[c].j
                                E_c = Hole.p[c].E
                                if (rem(l_g + l_c, 2) == rem(l_q + l_p, 2)) && (j_c == j_h) && (l_c == l_h)
                                    ME = - 0.5 * Float64(2*J + 1) / Float64(j_h + 1) * O2b_pp(a_p,a_q,a_h,a_g,J,P_hg,V_NN,Orb,Orb_NN) *
                                            O2b_pp(a_c,a_g,a_p,a_q,J,P_pq,V_NN,Orb,Orb_NN) / (E_h + E_g - E_p - E_q) / (E_c + E_g - E_p - E_q)
                                    dpRho[a_c,a_h] += ME
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    # pRho 2p-2h pn
    @inbounds Threads.@threads for p in 1:N_Particle.p
        a_p = Particle.p[p].a
        l_p = Particle.p[p].l
        j_p = Particle.p[p].j
        E_p = Particle.p[p].E
        @inbounds for q in 1:N_Particle.n
            a_q = Particle.n[q].a
            l_q = Particle.n[q].l
            j_q = Particle.n[q].j
            E_q = Particle.n[q].E
            P_pq = rem(l_p + l_q,2) + 1
            @inbounds for J in div(abs(j_p - j_q),2):div((j_p + j_q),2)
                @inbounds for h in 1:N_Hole.p
                    a_h = Hole.p[h].a
                    l_h = Hole.p[h].l
                    j_h = Hole.p[h].j
                    E_h = Hole.p[h].E
                    @inbounds for g in 1:N_Hole.n
                        a_g = Hole.n[g].a
                        l_g = Hole.n[g].l
                        j_g = Hole.n[g].j
                        E_g = Hole.n[g].E
                        P_hg = rem(l_h + l_g,2) + 1
                        if (abs(j_h - j_g) <= 2*J) && (2*J <= (j_h + j_g)) && (rem(l_h + l_g, 2) == rem(l_p + l_q, 2))
                            @inbounds for a in 1:N_Particle.p
                                a_a = Particle.p[a].a
                                l_a = Particle.p[a].l
                                j_a = Particle.p[a].j
                                E_a = Particle.p[a].E
                                if (rem(l_h + l_g, 2) == rem(l_q + l_a, 2)) && (j_a == j_p) && (l_a == l_p)
                                    ME = Float64(2*J + 1) / Float64(j_p + 1) * O2b_pn(a_p,a_q,a_h,a_g,J,P_pq,V_NN,Orb_NN) *
                                            O2b_pn(a_h,a_g,a_a,a_q,J,P_hg,V_NN,Orb_NN) / (E_h + E_g - E_p - E_q) / (E_h + E_g - E_a - E_q)
                                    dpRho[a_p,a_a] += ME
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    @inbounds Threads.@threads for h in 1:N_Hole.p
        a_h = Hole.p[h].a
        l_h = Hole.p[h].l
        j_h = Hole.p[h].j
        E_h = Hole.p[h].E
        @inbounds for g in 1:N_Hole.n
            a_g = Hole.n[g].a
            l_g = Hole.n[g].l
            j_g = Hole.n[g].j
            E_g = Hole.n[g].E
            @inbounds for J in div(abs(j_h - j_g),2):div((j_h + j_g),2)
                @inbounds for p in 1:N_Particle.p
                    a_p = Particle.p[p].a
                    l_p = Particle.p[p].l
                    j_p = Particle.p[p].j
                    E_p = Particle.p[p].E
                    @inbounds for q in 1:N_Particle.n
                        a_q = Particle.n[q].a
                        l_q = Particle.n[q].l
                        j_q = Particle.n[q].j
                        E_q = Particle.n[q].E
                        P_pq = rem(l_p + l_q,2) + 1
                        if (abs(j_p - j_q) <= 2*J) && (2*J <= (j_p + j_q)) && (rem(l_h + l_g, 2) == rem(l_p + l_q, 2))
                            @inbounds for c in 1:N_Hole.p
                                a_c = Hole.p[c].a
                                l_c = Hole.p[c].l
                                j_c = Hole.p[c].j
                                E_c = Hole.p[c].E
                                if (rem(l_g + l_c, 2) == rem(l_q + l_p, 2)) && (j_c == j_h) && (l_c == l_h)
                                    ME = - Float64(2*J + 1) / Float64(j_h + 1) * O2b_pn(a_p,a_q,a_h,a_g,J,P_pq,V_NN,Orb_NN) *
                                            O2b_pn(a_c,a_g,a_p,a_q,J,P_pq,V_NN,Orb_NN) / (E_h + E_g - E_p - E_q) / (E_c + E_g - E_p - E_q)
                                    dpRho[a_c,a_h] += ME
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    # nRho 1p-1h nn
    @inbounds Threads.@threads for p in 1:N_Particle.n
        a_p = Particle.n[p].a
        l_p = Particle.n[p].l
        j_p = Particle.n[p].j
        E_p = Particle.n[p].E
        @inbounds for h in 1:N_Hole.n
            a_h = Hole.n[h].a
            l_h = Hole.n[h].l
            j_h = Hole.n[h].j
            E_h = Hole.n[h].E
            if j_p == j_h && l_p == l_h
                @inbounds for f in 1:N_Hole.n
                    a_f = Hole.n[f].a
                    l_f = Hole.n[f].l
                    j_f = Hole.n[f].j
                    E_f = Hole.n[f].E
                    P_pf, P_hf = rem(l_p + l_f,2) + 1, rem(l_h + l_f,2) + 1
                    @inbounds for J in div(abs(j_p - j_f),2):div((j_p + j_f),2)
                        @inbounds for r in 1:N_Particle.n
                            a_r = Particle.n[r].a
                            l_r = Particle.n[r].l
                            j_r = Particle.n[r].j
                            E_r = Particle.n[r].E
                            @inbounds for q in 1:N_Particle.n
                                a_q = Particle.n[q].a
                                l_q = Particle.n[q].l
                                j_q = Particle.n[q].j
                                E_q = Particle.n[q].E
                                if (abs(j_q - j_r) <= 2*J) && (2*J <= (j_q + j_r)) && (rem(l_q + l_r, 2) == rem(l_h + l_f, 2))
                                    ME = 0.5 * Float64(2*J + 1) / Float64(j_p + 1) * O2b_nn(a_p,a_f,a_q,a_r,J,P_pf,V_NN,Orb,Orb_NN) * 
                                            O2b_nn(a_q,a_r,a_h,a_f,J,P_hf,V_NN,Orb,Orb_NN) / (E_h - E_p) / (E_h + E_f - E_q - E_r)
                                    dnRho[a_p,a_h] += ME
                                    dnRho[a_h,a_p] += ME
                                end
                            end
                        end
                    end
                end

                @inbounds for r in 1:N_Particle.n
                    a_r = Particle.n[r].a
                    l_r = Particle.n[r].l
                    j_r = Particle.n[r].j
                    E_r = Particle.n[r].E
                    P_pr, P_hr = rem(l_p + l_r,2) + 1, rem(l_h + l_r,2) + 1
                    @inbounds for J in div(abs(j_h - j_r),2):div((j_h + j_r),2)
                        @inbounds for f in 1:N_Hole.n
                            a_f = Hole.n[f].a
                            l_f = Hole.n[f].l
                            j_f = Hole.n[f].j
                            E_f = Hole.n[f].E
                            @inbounds for g in 1:N_Hole.n
                                a_g = Hole.n[g].a
                                l_g = Hole.n[g].l
                                j_g = Hole.n[g].j
                                E_g = Hole.n[g].E
                                if (abs(j_g - j_f) <= 2*J) && (2*J <= (j_g + j_f)) && (rem(l_p + l_r, 2) == rem(l_g + l_f, 2))
                                    ME = - 0.5 * Float64(2*J + 1) / Float64(j_p + 1) * O2b_nn(a_p,a_r,a_g,a_f,J,P_pr,V_NN,Orb,Orb_NN) * 
                                            O2b_nn(a_g,a_f,a_h,a_r,J,P_hr,V_NN,Orb,Orb_NN) / (E_h - E_p) / (E_g + E_f - E_p - E_r)
                                    dnRho[a_p,a_h] += ME
                                    dnRho[a_h,a_p] += ME
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    # nRho 1p-1h pn
    @inbounds Threads.@threads for p in 1:N_Particle.n
        a_p = Particle.n[p].a
        l_p = Particle.n[p].l
        j_p = Particle.n[p].j
        E_p = Particle.n[p].E
        @inbounds for h in 1:N_Hole.n
            a_h = Hole.n[h].a
            l_h = Hole.n[h].l
            j_h = Hole.n[h].j
            E_h = Hole.n[h].E
            if j_p == j_h && l_p == l_h

                @inbounds for f in 1:N_Hole.p
                    a_f = Hole.p[f].a
                    l_f = Hole.p[f].l
                    j_f = Hole.p[f].j
                    E_f = Hole.p[f].E
                    P_pf, P_hf = rem(l_p + l_f,2) + 1, rem(l_h + l_f,2) + 1
                    @inbounds for J in div(abs(j_p - j_f),2):div((j_p + j_f),2)
                        @inbounds for r in 1:N_Particle.p
                            a_r = Particle.p[r].a
                            l_r = Particle.p[r].l
                            j_r = Particle.p[r].j
                            E_r = Particle.p[r].E
                            @inbounds for q in 1:N_Particle.n
                                a_q = Particle.n[q].a
                                l_q = Particle.n[q].l
                                j_q = Particle.n[q].j
                                E_q = Particle.n[q].E
                                if (abs(j_q - j_r) <= 2*J) && (2*J <= (j_q + j_r)) && (rem(l_q + l_r, 2) == rem(l_h + l_f, 2))
                                    ME = Float64(2*J + 1) / Float64(j_p + 1) * O2b_pn(a_f,a_p,a_r,a_q,J,P_pf,V_NN,Orb_NN) * 
                                            O2b_pn(a_r,a_q,a_f,a_h,J,P_hf,V_NN,Orb_NN) / (E_h - E_p) / (E_h + E_f - E_q - E_r)
                                    dnRho[a_p,a_h] += ME
                                    dnRho[a_h,a_p] += ME
                                end
                            end
                        end
                    end
                end

                @inbounds for r in 1:N_Particle.p
                    a_r = Particle.p[r].a
                    l_r = Particle.p[r].l
                    j_r = Particle.p[r].j
                    E_r = Particle.p[r].E
                    P_pr, P_hr = rem(l_p + l_r,2) + 1, rem(l_h + l_r,2) + 1
                    @inbounds for J in div(abs(j_h - j_r),2):div((j_h + j_r),2)
                        @inbounds for f in 1:N_Hole.p
                            a_f = Hole.p[f].a
                            l_f = Hole.p[f].l
                            j_f = Hole.p[f].j
                            E_f = Hole.p[f].E
                            @inbounds for g in 1:N_Hole.n
                                a_g = Hole.n[g].a
                                l_g = Hole.n[g].l
                                j_g = Hole.n[g].j
                                E_g = Hole.n[g].E
                                if (abs(j_g - j_f) <= 2*J) && (2*J <= (j_g + j_f)) && (rem(l_p + l_r, 2) == rem(l_g + l_f, 2))
                                    ME = - Float64(2*J + 1) / Float64(j_p + 1) * O2b_pn(a_r,a_p,a_f,a_g,J,P_pr,V_NN,Orb_NN) * 
                                            O2b_pn(a_f,a_g,a_r,a_h,J,P_hr,V_NN,Orb_NN) / (E_h - E_p) / (E_g + E_f - E_p - E_r)
                                    dnRho[a_p,a_h] += ME
                                    dnRho[a_h,a_p] += ME
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    # nRho 2p-2h nn
    @inbounds Threads.@threads for p in 1:N_Particle.n
        a_p = Particle.n[p].a
        l_p = Particle.n[p].l
        j_p = Particle.n[p].j
        E_p = Particle.n[p].E
        @inbounds for q in 1:N_Particle.n
            a_q = Particle.n[q].a
            l_q = Particle.n[q].l
            j_q = Particle.n[q].j
            E_q = Particle.n[q].E
            P_pq = rem(l_p + l_q,2) + 1
            @inbounds for J in div(abs(j_p - j_q),2):div((j_p + j_q),2)
                @inbounds for h in 1:N_Hole.n
                    a_h = Hole.n[h].a
                    l_h = Hole.n[h].l
                    j_h = Hole.n[h].j
                    E_h = Hole.n[h].E
                    @inbounds for g in 1:N_Hole.n
                        a_g = Hole.n[g].a
                        l_g = Hole.n[g].l
                        j_g = Hole.n[g].j
                        E_g = Hole.n[g].E
                        P_hg = rem(l_h + l_g,2) + 1
                        if (abs(j_h - j_g) <= 2*J) && (2*J <= (j_h + j_g)) && (rem(l_h + l_g, 2) == rem(l_p + l_q, 2))
                            @inbounds for a in 1:N_Particle.n
                                a_a = Particle.n[a].a
                                l_a = Particle.n[a].l
                                j_a = Particle.n[a].j
                                E_a = Particle.n[a].E
                                if (rem(l_h + l_g, 2) == rem(l_q + l_a, 2)) && (j_a == j_p) && (l_a == l_p)
                                    ME = 0.5 * Float64(2*J + 1) / Float64(j_p + 1) * O2b_nn(a_p,a_q,a_h,a_g,J,P_pq,V_NN,Orb,Orb_NN) *
                                            O2b_nn(a_h,a_g,a_a,a_q,J,P_hg,V_NN,Orb,Orb_NN) / (E_h + E_g - E_p - E_q) / (E_h + E_g - E_a - E_q)
                                    dnRho[a_p,a_a] += ME
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    @inbounds Threads.@threads for h in 1:N_Hole.n
        a_h = Hole.n[h].a
        l_h = Hole.n[h].l
        j_h = Hole.n[h].j
        E_h = Hole.n[h].E
        @inbounds for g in 1:N_Hole.n
            a_g = Hole.n[g].a
            l_g = Hole.n[g].l
            j_g = Hole.n[g].j
            E_g = Hole.n[g].E
            @inbounds for J in div(abs(j_h - j_g),2):div((j_h + j_g),2)
                @inbounds for p in 1:N_Particle.n
                    a_p = Particle.n[p].a
                    l_p = Particle.n[p].l
                    j_p = Particle.n[p].j
                    E_p = Particle.n[p].E
                    @inbounds for q in 1:N_Particle.n
                        a_q = Particle.n[q].a
                        l_q = Particle.n[q].l
                        j_q = Particle.n[q].j
                        E_q = Particle.n[q].E
                        P_pq = rem(l_p + l_q,2) + 1
                        if (abs(j_p - j_q) <= 2*J) && (2*J <= (j_p + j_q)) && (rem(l_h + l_g, 2) == rem(l_p + l_q, 2))
                            @inbounds for c in 1:N_Hole.n
                                a_c = Hole.n[c].a
                                l_c = Hole.n[c].l
                                j_c = Hole.n[c].j
                                E_c = Hole.n[c].E
                                if (rem(l_g + l_c, 2) == rem(l_q + l_p, 2)) && (j_c == j_h) && (l_c == l_h)
                                    ME = - 0.5 * Float64(2*J + 1) / Float64(j_h + 1) * O2b_nn(a_p,a_q,a_h,a_g,J,P_pq,V_NN,Orb,Orb_NN) *
                                            O2b_nn(a_c,a_g,a_p,a_q,J,P_pq,V_NN,Orb,Orb_NN) / (E_h + E_g - E_p - E_q) / (E_c + E_g - E_p - E_q)
                                    dnRho[a_c,a_h] += ME
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    # nRho 2p-2h pn
    @inbounds Threads.@threads for p in 1:N_Particle.n
        a_p = Particle.n[p].a
        l_p = Particle.n[p].l
        j_p = Particle.n[p].j
        E_p = Particle.n[p].E
        @inbounds for q in 1:N_Particle.p
            a_q = Particle.p[q].a
            l_q = Particle.p[q].l
            j_q = Particle.p[q].j
            E_q = Particle.p[q].E
            @inbounds for J in div(abs(j_p - j_q),2):div((j_p + j_q),2)
                @inbounds for h in 1:N_Hole.n
                    a_h = Hole.n[h].a
                    l_h = Hole.n[h].l
                    j_h = Hole.n[h].j
                    E_h = Hole.n[h].E
                    @inbounds for g in 1:N_Hole.p
                        a_g = Hole.p[g].a
                        l_g = Hole.p[g].l
                        j_g = Hole.p[g].j
                        E_g = Hole.p[g].E
                        P_hg = rem(l_h + l_g,2) + 1
                        if (abs(j_h - j_g) <= 2*J) && (2*J <= (j_h + j_g)) && (rem(l_h + l_g, 2) == rem(l_p + l_q, 2))
                            @inbounds for a in 1:N_Particle.n
                                a_a = Particle.n[a].a
                                l_a = Particle.n[a].l
                                j_a = Particle.n[a].j
                                E_a = Particle.n[a].E
                                if (rem(l_h + l_g, 2) == rem(l_q + l_a, 2)) && (j_a == j_p) && (l_a == l_p)
                                    ME = Float64(2*J + 1) / Float64(j_p + 1) * O2b_pn(a_q,a_p,a_g,a_h,J,P_hg,V_NN,Orb_NN) *
                                            O2b_pn(a_g,a_h,a_q,a_a,J,P_hg,V_NN,Orb_NN) / (E_h + E_g - E_p - E_q) / (E_h + E_g - E_a - E_q)
                                    dnRho[a_p,a_a] += ME
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    @inbounds Threads.@threads for h in 1:N_Hole.n
        a_h = Hole.n[h].a
        l_h = Hole.n[h].l
        j_h = Hole.n[h].j
        E_h = Hole.n[h].E
        @inbounds for g in 1:N_Hole.p
            a_g = Hole.p[g].a
            l_g = Hole.p[g].l
            j_g = Hole.p[g].j
            E_g = Hole.p[g].E
            @inbounds for J in div(abs(j_h - j_g),2):div((j_h + j_g),2)
                @inbounds for p in 1:N_Particle.n
                    a_p = Particle.n[p].a
                    l_p = Particle.n[p].l
                    j_p = Particle.n[p].j
                    E_p = Particle.n[p].E
                    @inbounds for q in 1:N_Particle.p
                        a_q = Particle.p[q].a
                        l_q = Particle.p[q].l
                        j_q = Particle.p[q].j
                        E_q = Particle.p[q].E
                        P_pq = rem(l_p + l_q,2) + 1
                        if (abs(j_p - j_q) <= 2*J) && (2*J <= (j_p + j_q)) && (rem(l_h + l_g, 2) == rem(l_p + l_q, 2))
                            @inbounds for c in 1:N_Hole.n
                                a_c = Hole.n[c].a
                                l_c = Hole.n[c].l
                                j_c = Hole.n[c].j
                                E_c = Hole.n[c].E
                                if (rem(l_g + l_c, 2) == rem(l_q + l_p, 2)) && (j_c == j_h) && (l_c == l_h)
                                    ME = - Float64(2*J + 1) / Float64(j_h + 1) * O2b_pn(a_q,a_p,a_g,a_h,J,P_pq,V_NN,Orb_NN) *
                                            O2b_pn(a_g,a_c,a_q,a_p,J,P_pq,V_NN,Orb_NN) / (E_h + E_g - E_p - E_q) / (E_c + E_g - E_p - E_q)
                                    dnRho[a_c,a_h] += ME
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    # Read the transformation matrix C ... LHO -> HF ...
    C = O1b_import(Params,Orb,"IO/" * Params.Calc.Path * "/Bin/C_HF.bin")

    # Allocate the HF OBDM Rho ... in the LHO basis ...
    Rho_HF = HF_density_operator(a_max,C,Orb)

    # Transform Rho_HF to the HF basis ...
    Rho_HF = O1B(C.p' * Rho_HF.p * C.p, C.n' * Rho_HF.n * C.n)

    # Include HF-MBPT(3) corrections to Rho ...
    pRho, nRho = Rho_HF.p .+ dpRho, Rho_HF.n .+ dnRho

    # Symmetrize, clean & regularize Rho ...
        # Symmetrize density matrices ...
    pRho .= 0.5 .* (pRho .+ pRho')
    nRho .= 0.5 .* (nRho .+ nRho')
        # Eliminate any possible numerical noise spoiling block-diagonal structure ...
    @inbounds for a in 1:a_max
        @inbounds for b in 1:a_max
            if (Orb[a].j != Orb[b].j) || (Orb[a].l != Orb[b].l)
                pRho[a,b] = 0.0
                nRho[a,b] = 0.0
            end
        end
    end
        # Add a tiny deterministic diagonal splitting to lift degeneracies ...
    @inbounds for a in 1:a_max
        pRho[a,a] += 1e-10 * Float64(a)
        nRho[a,a] += 1e-10 * Float64(a)
    end

    # Diagonalize Rho ...
    pn, pD = eigen(Symmetric(pRho))
    nn, nD = eigen(Symmetric(nRho))

    # Reorder D & allocate new density matrix Rho ... in the HF-NAT basis ...
    Rho, D = HF_MBPT_reorder(Params,pnVector(pn,nn),O1B(pD,nD))

    # Calculate & export radial HF-MBPT(3) densities & radii ...
    Summary_File = "IO/" * Params.Calc.Path * "/HF/HF_Summary.dat"
    Densities_File = "IO/" * Params.Calc.Path * "/HF/Densities/HFMBPT_Radial_Densities.dat"
    OBDM_export(Params,Orb,Summary_File,Densities_File,O1B(C.p * D.p * Rho.p * D.p' * C.p', C.n * D.n * Rho.n * D.n' * C.n'),O1B(C.p * D.p, C.n * D.n))

    # Calculate & export HF-MBPT(3) occupation numbers ...
    HF_MBPT_occupation(Params,Rho_HF,Rho)

    return
end

function HF_MBPT_reorder(Params::Parameters,n::pnVector,C::O1B)
    # Read calculation parameters ...
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1) * (N_max + 2), 2)

    # Read needed arrays ...
    pC, pn = C.p, n.p
    nC, nn = C.n, n.n
 
    # Allocata temporary arrays for ordering ...
    pOrdering, nOrdering = zeros(Int64,a_max), zeros(Int64,a_max)
    V = diagm(ones(Float64,a_max))

    # Compute the largest overlaps ...
    @inbounds for a in 1:a_max
        pIndex, nIndex = 0, 0
        @views v = V[:,a]
        pMaximum, nMaximum = 0.0, 0.0
        @inbounds for b in 1:a_max
            @views pu = C.p[:,b]
            @views nu = C.n[:,b]
            pOverlap = abs(dot(v,pu))
            nOverlap = abs(dot(v,nu))
            if pOverlap > pMaximum
                pIndex = b
                pMaximum = pOverlap
            end
            if nOverlap > nMaximum
                nIndex = b
                nMaximum = nOverlap
            end
        end
        pOrdering[a] = pIndex
        nOrdering[a] = nIndex
    end

    # Reorder Occupations n & transformation matrix C ...
    @views pn .= pn[pOrdering]
    @views pC .= pC[:,pOrdering]

    @views nn .= nn[nOrdering]
    @views nC .= nC[:,nOrdering]

    # Define the diagonal density matrix Rho ...
    pRho, nRho = diagm(pn), diagm(nn)

    return O1B(pRho,nRho), O1B(pC,nC)
end

function HF_MBPT_occupation(Params::Parameters,Rho_HF::O1B,Rho_NAT::O1B)
    # Read parameters ...
    Output_File = Params.Calc.Path
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Set the export path ...
    Output_Path = "IO/" * Output_File * "/HF/Densities/HFMBPT_Occupation.dat"

    # Export the occupation probabilities to data file ...
    open(Output_Path, "w") do Write_File
        println(Write_File, "a\tpn_HF\tpn_MBPT\td_pn\tnn_HF\tnn_MBPT\td_nn")
        @inbounds for a in 1:a_max
            println(Write_File, string(a) * "\t" * string(Rho_HF.p[a,a]) * "\t" * string(Rho_NAT.p[a,a]) *
                    "\t" * string(round(Rho_HF.p[a,a] - Rho_NAT.p[a,a], digits = 5)) * "\t" * string(Rho_HF.n[a,a]) *
                    "\t" * string(Rho_NAT.n[a,a]) * "\t" * string(round(Rho_HF.n[a,a] - Rho_NAT.n[a,a], digits = 5)))
        end
    end

    return
end
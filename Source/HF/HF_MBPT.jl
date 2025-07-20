function HF_MBPT(Params::Parameters,Orb::Vector{NOrb},Orb_NN::NNOrb,VNN_res::NNInt)
    N_Particle, Particle, N_Hole, Hole = Make_ParticleHole_Orbitals(Params.Calc.Nmax,Params.Calc.Path,Orb)

    HFMBPT_Energy(Params,Orb,Orb_NN,N_Particle,Particle,N_Hole,Hole,VNN_res)

    HFMBPT_Density(Params,Orb,Orb_NN,N_Particle,Particle,N_Hole,Hole,VNN_res)

    return
end

function HFMBPT_Energy(Params::Parameters,Orb::Vector{NOrb},Orb_NN::NNOrb,N_Particle::pnInteger,Particle::pnSVector,N_Hole::pnInteger,Hole::pnSVector,VNN_res::NNInt)
    Output_File = Params.Calc.Path
    
    dE = 0.0

    dE_threads = Threads.Atomic{Float64}[Threads.Atomic{Float64}(0.0) for _ in 1:Threads.nthreads()]
    @inbounds Threads.@threads for p in 1:N_Particle.p
        thread_id = Threads.threadid()
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
                            @views dE_threads[thread_id][] += Float64(2*J + 1)/4.0 * abs(V2B(a_p,a_q,a_h,a_g,J,1,VNN_res.pp,Orb,Orb_NN))^2 / (E_h + E_g - E_p - E_q)
                        end
                    end
                end
            end
        end
    end
    dE += sum(x[] for x in dE_threads)

    dE_threads = Threads.Atomic{Float64}[Threads.Atomic{Float64}(0.0) for _ in 1:Threads.nthreads()]
    @inbounds Threads.@threads for p in 1:N_Particle.n
        thread_id = Threads.threadid()
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
                            @views dE_threads[thread_id][] += Float64(2*J + 1)/4.0 * abs(V2B(a_p,a_q,a_h,a_g,J,1,VNN_res.nn,Orb,Orb_NN))^2 / (E_h + E_g - E_p - E_q)
                        end
                    end
                end
            end
        end
    end
    dE += sum(x[] for x in dE_threads)

    dE_threads = Threads.Atomic{Float64}[Threads.Atomic{Float64}(0.0) for _ in 1:Threads.nthreads()]
    @inbounds Threads.@threads for p in 1:N_Particle.p
        thread_id = Threads.threadid()
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
                            @views dE_threads[thread_id][] += Float64(2*J + 1) * abs(V2B(a_p,a_q,a_h,a_g,J,0,VNN_res.pn,Orb,Orb_NN))^2 / (E_h + E_g - E_p - E_q)
                        end
                    end
                end
            end
        end
    end
    dE += sum(x[] for x in dE_threads)

    println("\nHF_MBPT(2) energy    ...   E^(2) = " * string(round(dE, sigdigits=9)) * " MeV")

    Summary_File =  open(string("IO/", Output_File, "/HF/HF_Summary.dat"), "a")
        println(Summary_File, "\nSpherical Hartree-Fock Leading Order Many-Body Perturbation Theory solution review:")
        println(Summary_File, "\nE_0^(2) = " * string(round(dE, sigdigits=9)) * "\t MeV \t\t ... \t LO HF-MBPT(2) correction to ground state energy")
    close(Summary_File)

    return
end

function HFMBPT_Density(Params::Parameters,Orb::Vector{NOrb},Orb_NN::NNOrb,N_Particle::pnInteger,Particle::pnSVector,N_Hole::pnInteger,Hole::pnSVector,VNN_res::NNInt)
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
                                    ME = 0.5 * Float64(2*J + 1) / Float64(j_p + 1) * V2B(a_p,a_f,a_q,a_r,J,1,VNN_res.pp,Orb,Orb_NN) * 
                                            V2B(a_q,a_r,a_h,a_f,J,1,VNN_res.pp,Orb,Orb_NN) / (E_h - E_p) / (E_h + E_f - E_q - E_r)
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
                                    ME = - 0.5 * Float64(2*J + 1) / Float64(j_p + 1) * V2B(a_p,a_r,a_g,a_f,J,1,VNN_res.pp,Orb,Orb_NN) * 
                                            V2B(a_g,a_f,a_h,a_r,J,1,VNN_res.pp,Orb,Orb_NN) / (E_h - E_p) / (E_g + E_f - E_p - E_r)
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
                                    ME = Float64(2*J + 1) / Float64(j_p + 1) * V2B(a_p,a_f,a_q,a_r,J,0,VNN_res.pn,Orb,Orb_NN) * 
                                            V2B(a_q,a_r,a_h,a_f,J,0,VNN_res.pn,Orb,Orb_NN) / (E_h - E_p) / (E_h + E_f - E_q - E_r)
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
                                    ME = - Float64(2*J + 1) / Float64(j_p + 1) * V2B(a_p,a_r,a_g,a_f,J,0,VNN_res.pn,Orb,Orb_NN) * 
                                            V2B(a_g,a_f,a_h,a_r,J,0,VNN_res.pn,Orb,Orb_NN) / (E_h - E_p) / (E_g + E_f - E_p - E_r)
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
                        if (abs(j_h - j_g) <= 2*J) && (2*J <= (j_h + j_g)) && (rem(l_h + l_g, 2) == rem(l_p + l_q, 2))
                            @inbounds for a in 1:N_Particle.p
                                a_a = Particle.p[a].a
                                l_a = Particle.p[a].l
                                j_a = Particle.p[a].j
                                E_a = Particle.p[a].E
                                if (rem(l_h + l_g, 2) == rem(l_q + l_a, 2)) && (j_a == j_p) && (l_a == l_p)
                                    ME = 0.5 * Float64(2*J + 1) / Float64(j_p + 1) * V2B(a_p,a_q,a_h,a_g,J,1,VNN_res.pp,Orb,Orb_NN) *
                                            V2B(a_h,a_g,a_a,a_q,J,1,VNN_res.pp,Orb,Orb_NN) / (E_h + E_g - E_p - E_q) / (E_h + E_g - E_a - E_q)
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
                        if (abs(j_p - j_q) <= 2*J) && (2*J <= (j_p + j_q)) && (rem(l_h + l_g, 2) == rem(l_p + l_q, 2))
                            @inbounds for c in 1:N_Hole.p
                                a_c = Hole.p[c].a
                                l_c = Hole.p[c].l
                                j_c = Hole.p[c].j
                                E_c = Hole.p[c].E
                                if (rem(l_g + l_c, 2) == rem(l_q + l_p, 2)) && (j_c == j_h) && (l_c == l_h)
                                    ME = - 0.5 * Float64(2*J + 1) / Float64(j_h + 1) * V2B(a_p,a_q,a_h,a_g,J,1,VNN_res.pp,Orb,Orb_NN) *
                                            V2B(a_c,a_g,a_p,a_q,J,1,VNN_res.pp,Orb,Orb_NN) / (E_h + E_g - E_p - E_q) / (E_c + E_g - E_p - E_q)
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
                        if (abs(j_h - j_g) <= 2*J) && (2*J <= (j_h + j_g)) && (rem(l_h + l_g, 2) == rem(l_p + l_q, 2))
                            @inbounds for a in 1:N_Particle.p
                                a_a = Particle.p[a].a
                                l_a = Particle.p[a].l
                                j_a = Particle.p[a].j
                                E_a = Particle.p[a].E
                                if (rem(l_h + l_g, 2) == rem(l_q + l_a, 2)) && (j_a == j_p) && (l_a == l_p)
                                    ME = Float64(2*J + 1) / Float64(j_p + 1) * V2B(a_p,a_q,a_h,a_g,J,0,VNN_res.pn,Orb,Orb_NN) *
                                            V2B(a_h,a_g,a_a,a_q,J,0,VNN_res.pn,Orb,Orb_NN) / (E_h + E_g - E_p - E_q) / (E_h + E_g - E_a - E_q)
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
                        if (abs(j_p - j_q) <= 2*J) && (2*J <= (j_p + j_q)) && (rem(l_h + l_g, 2) == rem(l_p + l_q, 2))
                            @inbounds for c in 1:N_Hole.p
                                a_c = Hole.p[c].a
                                l_c = Hole.p[c].l
                                j_c = Hole.p[c].j
                                E_c = Hole.p[c].E
                                if (rem(l_g + l_c, 2) == rem(l_q + l_p, 2)) && (j_c == j_h) && (l_c == l_h)
                                    ME = - Float64(2*J + 1) / Float64(j_h + 1) * V2B(a_p,a_q,a_h,a_g,J,0,VNN_res.pn,Orb,Orb_NN) *
                                            V2B(a_c,a_g,a_p,a_q,J,0,VNN_res.pn,Orb,Orb_NN) / (E_h + E_g - E_p - E_q) / (E_c + E_g - E_p - E_q)
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
                                    ME = 0.5 * Float64(2*J + 1) / Float64(j_p + 1) * V2B(a_p,a_f,a_q,a_r,J,1,VNN_res.nn,Orb,Orb_NN) * 
                                            V2B(a_q,a_r,a_h,a_f,J,1,VNN_res.nn,Orb,Orb_NN) / (E_h - E_p) / (E_h + E_f - E_q - E_r)
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
                                    ME = - 0.5 * Float64(2*J + 1) / Float64(j_p + 1) * V2B(a_p,a_r,a_g,a_f,J,1,VNN_res.nn,Orb,Orb_NN) * 
                                            V2B(a_g,a_f,a_h,a_r,J,1,VNN_res.nn,Orb,Orb_NN) / (E_h - E_p) / (E_g + E_f - E_p - E_r)
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
                                    ME = Float64(2*J + 1) / Float64(j_p + 1) * V2B(a_f,a_p,a_r,a_q,J,0,VNN_res.pn,Orb,Orb_NN) * 
                                            V2B(a_r,a_q,a_f,a_h,J,0,VNN_res.pn,Orb,Orb_NN) / (E_h - E_p) / (E_h + E_f - E_q - E_r)
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
                                    ME = - Float64(2*J + 1) / Float64(j_p + 1) * V2B(a_r,a_p,a_f,a_g,J,0,VNN_res.pn,Orb,Orb_NN) * 
                                            V2B(a_f,a_g,a_r,a_h,J,0,VNN_res.pn,Orb,Orb_NN) / (E_h - E_p) / (E_g + E_f - E_p - E_r)
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
                        if (abs(j_h - j_g) <= 2*J) && (2*J <= (j_h + j_g)) && (rem(l_h + l_g, 2) == rem(l_p + l_q, 2))
                            @inbounds for a in 1:N_Particle.n
                                a_a = Particle.n[a].a
                                l_a = Particle.n[a].l
                                j_a = Particle.n[a].j
                                E_a = Particle.n[a].E
                                if (rem(l_h + l_g, 2) == rem(l_q + l_a, 2)) && (j_a == j_p) && (l_a == l_p)
                                    ME = 0.5 * Float64(2*J + 1) / Float64(j_p + 1) * V2B(a_p,a_q,a_h,a_g,J,1,VNN_res.nn,Orb,Orb_NN) *
                                            V2B(a_h,a_g,a_a,a_q,J,1,VNN_res.nn,Orb,Orb_NN) / (E_h + E_g - E_p - E_q) / (E_h + E_g - E_a - E_q)
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
                        if (abs(j_p - j_q) <= 2*J) && (2*J <= (j_p + j_q)) && (rem(l_h + l_g, 2) == rem(l_p + l_q, 2))
                            @inbounds for c in 1:N_Hole.n
                                a_c = Hole.n[c].a
                                l_c = Hole.n[c].l
                                j_c = Hole.n[c].j
                                E_c = Hole.n[c].E
                                if (rem(l_g + l_c, 2) == rem(l_q + l_p, 2)) && (j_c == j_h) && (l_c == l_h)
                                    ME = - 0.5 * Float64(2*J + 1) / Float64(j_h + 1) * V2B(a_p,a_q,a_h,a_g,J,1,VNN_res.nn,Orb,Orb_NN) *
                                            V2B(a_c,a_g,a_p,a_q,J,1,VNN_res.nn,Orb,Orb_NN) / (E_h + E_g - E_p - E_q) / (E_c + E_g - E_p - E_q)
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
                        if (abs(j_h - j_g) <= 2*J) && (2*J <= (j_h + j_g)) && (rem(l_h + l_g, 2) == rem(l_p + l_q, 2))
                            @inbounds for a in 1:N_Particle.n
                                a_a = Particle.n[a].a
                                l_a = Particle.n[a].l
                                j_a = Particle.n[a].j
                                E_a = Particle.n[a].E
                                if (rem(l_h + l_g, 2) == rem(l_q + l_a, 2)) && (j_a == j_p) && (l_a == l_p)
                                    ME = Float64(2*J + 1) / Float64(j_p + 1) * V2B(a_q,a_p,a_g,a_h,J,0,VNN_res.pn,Orb,Orb_NN) *
                                            V2B(a_g,a_h,a_q,a_a,J,0,VNN_res.pn,Orb,Orb_NN) / (E_h + E_g - E_p - E_q) / (E_h + E_g - E_a - E_q)
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
                        if (abs(j_p - j_q) <= 2*J) && (2*J <= (j_p + j_q)) && (rem(l_h + l_g, 2) == rem(l_p + l_q, 2))
                            @inbounds for c in 1:N_Hole.n
                                a_c = Hole.n[c].a
                                l_c = Hole.n[c].l
                                j_c = Hole.n[c].j
                                E_c = Hole.n[c].E
                                if (rem(l_g + l_c, 2) == rem(l_q + l_p, 2)) && (j_c == j_h) && (l_c == l_h)
                                    ME = - Float64(2*J + 1) / Float64(j_h + 1) * V2B(a_q,a_p,a_g,a_h,J,0,VNN_res.pn,Orb,Orb_NN) *
                                            V2B(a_g,a_c,a_q,a_p,J,0,VNN_res.pn,Orb,Orb_NN) / (E_h + E_g - E_p - E_q) / (E_c + E_g - E_p - E_q)
                                    dnRho[a_c,a_h] += ME
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    HFMBPT_OBDM_Occupation(Params,Orb,pnMatrix(dpRho,dnRho))

    @time HFMBPT_Radial_Density(Params,Orb,pnMatrix(dpRho,dnRho))

    return
end

function HFMBPT_OBDM_Occupation(Params::Parameters,Orb::Vector{NOrb},dRho::pnMatrix)
    # Read parameters ...
    Orthogon = Params.Calc.Ortho
    Output_File = Params.Calc.Path
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1)*(N_max + 2),2)

    if Orthogon == true
        Output_Path = "IO/" * Output_File * "/HF/Densities/HFMBPT_Occupation_Ortho.dat"
    else
        Output_Path = "IO/" * Output_File * "/HF/Densities/HFMBPT_Occupation_Spur.dat"
    end

    # Initialize HF 1-body density matrix ...
    pRho_HF, nRho_HF = zeros(Float64,a_max,a_max),zeros(Float64,a_max,a_max)
    @inbounds for a in 1:a_max
        pRho_HF[a,a] = Orb[a].pO
        nRho_HF[a,a] = Orb[a].nO
    end

    # Construct the HF-MBPT OBDM ...
    pRho_MBPT, nRho_MBPT = pRho_HF .+ dRho.p, nRho_HF .+ dRho.n

    # Ensure Hermicity of OBDMs ...
    pRho_MBPT .= 0.5 .* (pRho_MBPT .+ pRho_MBPT')
    nRho_MBPT .= 0.5 .* (nRho_MBPT .+ nRho_MBPT')

    # Diagonalize OBDMs ...
    pN, pU = eigen(pRho_MBPT, sortby = nothing)
    nN, nU = eigen(pRho_MBPT, sortby = nothing)

    # Reorder according to HF orbital convention ...
    pN, pU = HFMBPT_OBDM_Reorder(a_max,pN,pU)
    nN, nU = HFMBPT_OBDM_Reorder(a_max,nN,nU)

    # Export the occupation probabilities to data file ...

    # Initialize occupation probabilities vectors ...
    pO_HF = Vector{Float64}(undef,a_max)
    nO_HF = Vector{Float64}(undef,a_max)
    pO_MBPT = Vector{Float64}(undef,a_max)
    nO_MBPT = Vector{Float64}(undef,a_max)
    pa = Vector{Int64}(undef,a_max)
    na = Vector{Int64}(undef,a_max)

    @inbounds for a in 1:a_max
        pO_HF[a] = pRho_HF[a,a]
        nO_HF[a] = nRho_HF[a,a]
        pO_MBPT[a] = round(pN[a], digits = 4)
        nO_MBPT[a] = round(nN[a], digits = 4)
        pa[a] = a
        na[a] = a
    end

    Sort_pInd = sortperm(pO_MBPT, rev = true)
    Sort_nInd = sortperm(nO_MBPT, rev = true)

    pO_HF = pO_HF[Sort_pInd]
    pO_MBPT = pO_MBPT[Sort_pInd]
    pa = pa[Sort_pInd]

    nO_HF = nO_HF[Sort_nInd]
    nO_MBPT = nO_MBPT[Sort_nInd]
    na = na[Sort_nInd]

    # HF-MBPT occupation export ...
    open(Output_Path, "w") do Write_File
        println(Write_File, "a\ta_p\tpN_HF\tpN_MBPT\tpN_dep\ta_n\tnN_HF\tnN_MBPT\tnN_dep")
        @inbounds for a in 1:a_max
            println(Write_File, string(a) * "\t" * string(pa[a]) * "\t" * string(pO_HF[a]) * "\t" * string(pO_MBPT[a]) * "\t" * string(round(pO_MBPT[a] - pO_HF[a], digits = 5)) * "\t" * string(na[a]) * "\t" * string(nO_HF[a]) * "\t" * string(nO_MBPT[a]) * "\t" * string(round(nO_MBPT[a] - nO_HF[a], digits = 5)))
        end
    end

    return
end

function HFMBPT_OBDM_Reorder(a_max::Int64,N::Vector{Float64},U::Matrix{Float64})
    Ordering = zeros(Int64,a_max)
    V = diagm(ones(Float64,a_max))
    @inbounds for a in 1:a_max
        Index = 0
        v = V[:,a]
        Maximum = 0.0
        @inbounds for b in 1:a_max
            u = U[:,b]
            Overlap = abs(dot(v,u))
            if Overlap > Maximum
                Index = b
                Maximum = Overlap
            end
        end
        Ordering[a] = Index
    end
    N = N[Ordering]
    U = U[:, Ordering]
    @inbounds for a in 1:a_max
        if U[a,a] < 0.0
            @inbounds for b in 1:a_max
                U[b,a] = -1.0 * U[b,a]
            end
        end
    end
    return N, U
end
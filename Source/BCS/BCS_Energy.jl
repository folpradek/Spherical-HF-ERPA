function BCS_Energy(Params::Parameters,Kappa::pnMatrix,Orb::Vector{NOrb},Orb_NN::NNOrb,VNN::NNInt)
   # Read parameters ...
    N_max = Params.Calc.Nmax
    N_2max = Params.Calc.N2max
    a_max = div((N_max + 1)*(N_max + 2),2)

    #Calculate the BCS ground-state energy ...
    println("\nCalculating total mean-field + BCS ground-state energy ...")

    E_BCS_partial = Threads.Atomic{Float64}[Threads.Atomic{Float64}(0.0) for _ in 1:Threads.nthreads()]

    @inbounds Threads.@threads for a = 1:a_max
        thread_id = Threads.threadid()
        n_a = Orb[a].n
        l_a = Orb[a].l
        j_a = Orb[a].j
        @inbounds for b = 1:a_max
            n_b = Orb[b].n
            l_b = Orb[b].l
            j_b = Orb[b].j
            if ((2*(n_a + n_b) + l_a + l_b) <= N_2max) && (l_a == l_b) && (j_a == j_b)
                @inbounds for d = 1:a_max
                    l_d = Orb[d].l
                    j_d = Orb[d].j
                    n_d = Orb[d].n
                    ja_jd_hat = sqrt(Float64((j_a + 1)*(j_d + 1)))
                    @inbounds for e = 1:a_max
                        n_e = Orb[e].n
                        l_e = Orb[e].l
                        j_e = Orb[e].j
                        if (l_d == l_e && j_d == j_e) && ((2*(n_d + n_e) + l_d + l_e) <= N_2max) && (rem(l_a + l_b,2) == rem(l_d + l_e,2))
                            @views E_BCS_partial[thread_id][] += 0.25 * ja_jd_hat * Kappa.p[a,b] * Kappa.p[d,e] * V2B(a,b,d,e,0,1,VNN.pp,Orb,Orb_NN)
                            @views E_BCS_partial[thread_id][] += 0.25 * ja_jd_hat * Kappa.n[a,b] * Kappa.n[d,e] * V2B(a,b,d,e,0,1,VNN.nn,Orb,Orb_NN)
                        end
                    end
                end
            end
        end
    end

    E_BCS = sum(x[] for x in E_BCS_partial)

    println("\nBCS ground-state pairing energy    ...   E_BCS = " * string(E_BCS) * " MeV")

    return E_BCS
end
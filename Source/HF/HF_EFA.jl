function HF_EFA(Params::Parameters,Orb::Vector{Orb1B},Orb_NNN::Orb3B,SPE::pnVector,C::O1B,V_NNN::Array{Vector{Vector{Float32}},4})
    println("\nStarting Spherical Equal-Filling Approximation Hartree-Fock calculation ... ")
    println("\nThe HF-EFA target nucleus is set to ...")
    println("A = " * string(Params.Calc.HF.EFA.A) * ", Z = " * string(Params.Calc.HF.EFA.Z) * "\n")
    println("\nThe HF-EFA reference nucleus is ...")
    println("A = " * string(Params.Calc.A) * ", Z = " * string(Params.Calc.Z))
    
    # Define new Parameters structure for the HF-EFA calculation ...
    Params_EFA = Parameters(Params.Int,
                Calculation_Parameters(
                A = Params.Calc.HF.EFA.A,
                Z = Params.Calc.HF.EFA.Z,
                hw = Params.Calc.hw,
                Nmax = Params.Calc.Nmax,
                N2max = Params.Calc.N2max,
                N3max = Params.Calc.N3max,                                                         
                CMS = Params.Calc.CMS,
                Path = Params.Calc.Path,
                HF = Params.Calc.HF,
                ))

    # 1-body kinetic operator ...
    T = T1b(Params_EFA,Orb)

    # 2-body bare NN interaction & Orbitals ...
    V_NN, Orb_NN = V2b_read(Params_EFA,Orb)

    # Calculate HF-EFA density operator Rho ...
    Rho = HF_EFA_density_operator(Params_EFA,SPE,C,Orb)

    # Evaluate the HF-EFA energy ...
    E_HF = HF_energy(Params_EFA,Rho,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN)

    # Evaluate the total HF-EFA kinetic energy ...
    T_HF = T1b_energy(Params_EFA,Rho,Orb,T)

    # Export the HF-EFA summary file ...
    HF_EFA_export_summary(Params_EFA,E_HF,T_HF)

    # Calculate & export radial HF-EFA densities & radii ...
    Summary_File = "IO/" * Params.Calc.Path * "/HF-EFA/HF-EFA_Summary.dat"
    Densities_File = "IO/" * Params.Calc.Path * "/HF-EFA/Densities/HF_Radial_Densities.dat"
    OBDM_export(Params_EFA,Orb,Summary_File,Densities_File,Rho,C)

    # Export information on the HF-EFA single-particle states ...
    HF_EFA_export_SPS_summary(Params_EFA,SPE,O1B(C.p' * Rho.p * C.p, C.n' * Rho.n * C.n),Orb)

    if Params.Calc.HF.EFA.BMF == true
        # Make density-dependent residual NN interaction ... NO2B approximation ...
        @time V_NN = V2b_residual_no2b(Params_EFA,Orb,Orb_NN,Orb_NNN,Rho,V_NN,V_NNN)

        # Transform the density-dependent NN interaction to the target HF basis ...
        @time W_NN = O2b_transformation(Params_EFA,Orb,Orb_NN,V_NN,C)

        # Evaluate the HF-MBPT(2) correlation energy ...
        HF_EFA_MBPT(Params_EFA,Orb,Orb_NN,SPE,O1B(C.p' * Rho.p * C.p, C.n' * Rho.n * C.n),V_NN)
    end

    # Deallocate the 2-body NN interaction ...
    V_NN = nothing
    W_NN = nothing

    # Perform the Garbage Collection ...
    GC.gc()

    println("\nThe Spherical Equal-Filling Approximation Hartree-Fock calculation finished ... ")

    return
end

function HF_EFA_density_operator(Params::Parameters,SPE::pnVector,C::O1B,Orb::Vector{Orb1B})
    # Read parameters ...
    Z = Params.Calc.HF.EFA.Z
    N = Params.Calc.HF.EFA.A - Z
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Initialize density Rho ...
    pRho, nRho = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)

    # Initialize the vector of occupation numbers ...
    pOcc, nOcc = zeros(Float64,a_max), zeros(Float64,a_max)

    # Initialize the counters for nucleon numbers ...
    Z_count, N_count = 0, 0

    # Determine the proton occupation numbers ...
    @inbounds for a in sortperm(SPE.p)
        j_a = Orb[a].j
        @inbounds for m_a in -j_a:2:j_a
            if (Z_count + 1) <= Z
                pOcc[a] += 1.0 / Float64(j_a + 1)
                Z_count += 1
            end
        end
        if Z_count == Z
            break
        end
    end

    # Determine the neutron occupation numbers ...
    @inbounds for a in sortperm(SPE.n)
        j_a = Orb[a].j
        @inbounds for m_a in -j_a:2:j_a
            if (N_count + 1) <= N
                nOcc[a] += 1.0 / Float64(j_a + 1)
                N_count += 1
            end
        end
        if N_count == N
            break
        end
    end

    # Evaluate the density matrix Rho ...
    @inbounds for k in 1:a_max
        @inbounds for l in 1:a_max
            if Orb[k].l == Orb[l].l && Orb[k].j == Orb[l].j
                pSum, nSum = 0.0, 0.0
                @inbounds for a in 1:a_max
                    if Orb[a].l == Orb[k].l && Orb[a].j == Orb[k].j
                        pSum += C.p[k,a] * C.p[l,a] * pOcc[a]
                    end
                    if Orb[a].l == Orb[k].l && Orb[a].j == Orb[k].j
                        nSum += C.n[k,a] * C.n[l,a] * nOcc[a]
                    end
                end
                pRho[k,l] = pSum
                nRho[k,l] = nSum
            end
        end
    end
    
    return O1B(pRho,nRho)
end

function HF_EFA_MBPT(Params::Parameters,Orb::Vector{Orb1B},Orb_NN::Orb2B,SPE::pnVector,Rho::O1B,V_NN::O2B)
    # Read parameters ...
    IO = Params.Calc.Path

    # Local function to generate lists of particle & hole indices for MBPT iteration ...
    function HF_EFA_MBPT_indices(Params::Parameters,Rho::O1B)
        # Read parameters ...
        N_max = Params.Calc.Nmax
        a_max = div((N_max + 1)*(N_max + 2),2)

        # Initialize the counters for particle & hole states ...
        pN_p, pN_h = 0, 0
        nN_p, nN_h = 0, 0

        # Count the number of particle & hole states ...
        @inbounds for a in 1:a_max
            if abs(Rho.p[a,a] - 1.0) > 1e-3
                pN_p += 1
            end
            if abs(Rho.p[a,a]) > 1e-3
                pN_h += 1
            end
            if abs(Rho.n[a,a] - 1.0) > 1e-3
                nN_p += 1
            end
            if abs(Rho.n[a,a]) > 1e-3
                nN_h += 1
            end
        end

        # Initialize the lists of particle & hole states ...
        pParticle_list = Vector{Int64}(undef,pN_p)
        pHole_list = Vector{Int64}(undef,pN_h)
        nParticle_list = Vector{Int64}(undef,nN_p)
        nHole_list = Vector{Int64}(undef,nN_h)

        # Reset the counters ...
        pN_p, pN_h = 0, 0
        nN_p, nN_h = 0, 0

        @inbounds for a in 1:a_max
            if abs(Rho.p[a,a] - 1.0) > 1e-3
                pN_p += 1
                pParticle_list[pN_p] = a
            end
            if abs(Rho.p[a,a]) > 1e-3
                pN_h += 1
                pHole_list[pN_h] = a
            end
            if abs(Rho.n[a,a] - 1.0) > 1e-3
                nN_p += 1
                nParticle_list[nN_p] = a
            end
            if abs(Rho.n[a,a]) > 1e-3
                nN_h += 1
                nHole_list[nN_h] = a
            end
        end
        return (pParticle_list,pHole_list,nParticle_list,nHole_list)
    end

    # Get the lists of particle & hole states ...
    Orb_2p2h = HF_EFA_MBPT_indices(Params,Rho)

    # Initialize the thread accumulators ...
    dE, dE_threads = 0.0, zeros(Float64,Threads.maxthreadid())

    # Calculate the MBPT(2) energy correction ...
    println("\nCalculating HF-EFA-MBPT(2) ground-state correlation energy ... ")
    
    # proton-proton contribution ...
    @inbounds Threads.@threads :static for p in Orb_2p2h[1]
        Sum, Tid = 0.0, Threads.threadid()
        l_p, j_p, E_p = Orb[p].l, Orb[p].j, SPE.p[p]
        @inbounds for q in Orb_2p2h[1]
            l_q, j_q, E_q = Orb[q].l, Orb[q].j, SPE.p[q]
            P = rem(l_p + l_q, 2) + 1
            @inbounds for J in div(abs(j_p - j_q),2):div((j_p + j_q),2)
                @inbounds for h in Orb_2p2h[2]
                    l_h, j_h, E_h = Orb[h].l, Orb[h].j, SPE.p[h]
                    @inbounds for g in Orb_2p2h[2]
                        l_g, j_g, E_g = Orb[g].l, Orb[g].j, SPE.p[g]
                        Denom = E_h + E_g - E_p - E_q
                        if (P == (rem(l_h + l_g, 2) + 1)) && (abs(j_h - j_g) <= 2*J) && (2*J <= (j_h + j_g)) && abs(Denom) > 1e-3
                            Sum += Float64(2*J + 1) / 4.0 * abs(O2b_pp(p,q,h,g,J,P,V_NN,Orb,Orb_NN))^2 / Denom * (1.0 - Rho.p[p,p]) * (1.0 - Rho.p[q,q]) * Rho.p[h,h] * Rho.p[g,g]
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
    @inbounds Threads.@threads :static for p in Orb_2p2h[3]
        Sum, Tid = 0.0, Threads.threadid()
        l_p, j_p, E_p = Orb[p].l, Orb[p].j, SPE.n[p]
        @inbounds for q in Orb_2p2h[3]
            l_q, j_q, E_q = Orb[q].l, Orb[q].j, SPE.n[q]
            P = rem(l_p + l_q, 2) + 1
            @inbounds for J in div(abs(j_p - j_q),2):div((j_p + j_q),2)
                @inbounds for h in Orb_2p2h[4]
                    l_h, j_h, E_h = Orb[h].l, Orb[h].j, SPE.n[h]
                    @inbounds for g in Orb_2p2h[4]
                        l_g, j_g, E_g = Orb[g].l, Orb[g].j, SPE.n[g]
                        Denom = E_h + E_g - E_p - E_q
                        if (P == (rem(l_h + l_g, 2) + 1)) && (abs(j_h - j_g) <= 2*J) && (2*J <= (j_h + j_g)) && abs(Denom) > 1e-3
                            Sum += Float64(2*J + 1) / 4.0 * abs(O2b_nn(p,q,h,g,J,P,V_NN,Orb,Orb_NN))^2 / Denom * (1.0 - Rho.n[p,p]) * (1.0 - Rho.n[q,q]) * Rho.n[h,h] * Rho.n[g,g]
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
    @inbounds Threads.@threads :static for p in Orb_2p2h[1]
        Sum, Tid = 0.0, Threads.threadid()
        l_p, j_p, E_p = Orb[p].l, Orb[p].j, SPE.p[p]
        @inbounds for q in Orb_2p2h[3]
            l_q, j_q, E_q = Orb[q].l, Orb[q].j, SPE.n[q]
            P = rem(l_p + l_q, 2) + 1
            @inbounds for J in div(abs(j_p - j_q),2):div((j_p + j_q),2)
                @inbounds for h in Orb_2p2h[2]
                    l_h, j_h, E_h = Orb[h].l, Orb[h].j, SPE.p[h]
                    @inbounds for g in Orb_2p2h[4]
                        l_g, j_g, E_g = Orb[g].l, Orb[g].j, SPE.n[g]
                        Denom = E_h + E_g - E_p - E_q
                        if (P == (rem(l_h + l_g, 2) + 1)) && (abs(j_h - j_g) <= 2*J) && (2*J <= (j_h + j_g)) && abs(Denom) > 1e-3
                            Sum += Float64(2*J + 1) * abs(O2b_pn(p,q,h,g,J,P,V_NN,Orb_NN))^2 / Denom * (1.0 - Rho.p[p,p]) * (1.0 - Rho.n[q,q]) * Rho.p[h,h] * Rho.n[g,g]
                        end
                    end
                end
            end
        end
        dE_threads[Tid] += Sum 
    end
    dE += sum(dE_threads)

    # Print the total HF-EFA-MBPT(2) energy correction ...
    println("\tHF-EFA-MBPT(2) energy    ...   E^(2) = " * string(round(dE, sigdigits=9)) * " MeV")

    # Write the total HF-EFA-MBPT(2) energy correction to the summary file ...
    Summary =  open(string("IO/", IO, "/HF-EFA/HF-EFA_Summary.dat"), "a")
        println(Summary, "\nSpherical Equal-Filling Approximation Hartree-Fock Leading Order Many-Body Perturbation Theory solution review:")
        println(Summary, "\nE_0^(2) = " * string(round(dE, sigdigits=9)) * "\t MeV \t\t ... \t LO HF-EFA-MBPT(2) correction to ground-state energy")
        println(Summary, "\nE_0^(2) / A = " * string(round(dE / Float64(Params.Calc.HF.EFA.A), sigdigits=9)) * "\t MeV \t\t ... \t LO HF-EFA-MBPT(2) correction to ground-state energy per nucleon")
    close(Summary)

    return
end

function HF_EFA_export_summary(Params::Parameters,E_HF::Float64,T_HF::Float64)
    # Read parameters ...
    hw = Params.Calc.hw
    A = Params.Calc.HF.EFA.A
    Z = Params.Calc.HF.EFA.Z
    N_max = Params.Calc.Nmax
    N_2max = Params.Calc.N2max
    N_3max = Params.Calc.N3max
    CMS = Params.Calc.CMS
    Output_File = Params.Calc.Path

    # Export the HF calculation summary ...
    println("\nExporting HF-EFA calculation summary ...")
    Summary =  open(string("IO/", Output_File, "/HF-EFA/HF-EFA_Summary.dat"), "a")
        println(Summary, "Nuclid data:    A = " * string(A) * ", Z = " * string(Z))
        println(Summary, "\nCalculation data:    hw = " * string(hw) * " , N_max = " * string(N_max) *
                " , N_2max = " * string(N_2max) * " , N_3max = " * string(N_3max) * ", J-basis size = " *
                string(div((N_max+1)*(N_max+2),2)) * ", M-basis size = " * string(div((N_max+1)*(N_max+2)*(N_max+3),6)))
        println(Summary, "\t\t\tCenter-of-Mass System correction is set to:   CMS = ''" * string(CMS) * "''")
        if CMS == "CMS1+2B"
            println(Summary, "\n\t\t\tCombined 1-body + 2-body Center-of-Mass System (CMS) motion correction is set ...")
        elseif CMS == "CMS2B"
            println(Summary, "\n\t\t\tOnly pure 2-body Center-of-Mass System (CMS) motion correction is set ...")
        else
            println(Summary, "\n\t\t\tNo Center-of-Mass System (CMS) motion correction is set ...")
        end
        println(Summary, "\nSpherical Equal-Filling Approximation Hartree-Fock solution review:")
        println(Summary, "\nE_HF = " * string(round(E_HF, sigdigits=9)) * "\t MeV \t\t ... \t Mean-field ground state energy")
        println(Summary, "T_HF = " * string(round(T_HF, sigdigits=9)) * "\t MeV \t\t ... \t Mean-field total kinetic energy")
        println(Summary, "\nE_HF / A = " * string(round(E_HF / Float64(A), sigdigits=9)) * "\t MeV / Nucleon \t\t ... \t Mean-field ground state energy")
        println(Summary, "T_HF / A = " * string(round(T_HF / Float64(A), sigdigits=9)) * "\t MeV / Nucleon \t\t ... \t Mean-field total kinetic energy")
    close(Summary)

    return
end

function HF_EFA_export_SPS_summary(Params::Parameters,SPE::pnVector,Rho::O1B,Orb::Vector{Orb1B})
    # Read parameters
    N_max = Params.Calc.Nmax
    N_2max = Params.Calc.N2max
    Output_File = Params.Calc.Path
    a_max = div((N_max + 1)*(N_max + 2),2)
    J_max = N_2max + 1

    # Initialize array for radial number n ...
    pn, nn = zeros(Int64,N_max+1,J_max+1), zeros(Int64,N_max+1,J_max+1)

    # Precompute energy ordered s.p. states for print ...
    Sort_pOrbs = sortperm(SPE.p)
    Sort_nOrbs = sortperm(SPE.n)

    # Export orbitals ...
    println("\nExporting HF-EFA calculation summary ...")
    Summary =  open(string("IO/", Output_File, "/HF-EFA/HF-EFA_Summary.dat"), "a")
        println(Summary,"\nProton single-particle states")
        println(Summary,"____________________________________________________________________")
        @printf(Summary, "%4s %4s %4s %8s %14s", "n", "l", "2j", "O_a", "E\n")
        @inbounds for Ind_a in 1:a_max
            a = Sort_pOrbs[Ind_a]
            l_a = Orb[a].l
            j_a = Orb[a].j
            n_a = pn[l_a+1,j_a+1]
            O_a = abs(Rho.p[a,a])
            E_a = SPE.p[a]
            @printf(Summary, "%4d %4d %4d %8.3f %14.5f\tMeV\n",n_a, l_a, j_a, O_a, E_a)
            pn[l_a+1,j_a+1] += 1
        end
        println(Summary,"____________________________________________________________________")

        println(Summary,"\nNeutron single-particle states")
        println(Summary,"____________________________________________________________________")
        @printf(Summary, "%4s %4s %4s %8s %14s", "n", "l", "2j", "O_a", "E\n")
        for Ind_a in 1:a_max
            a = Sort_nOrbs[Ind_a]
            l_a = Orb[a].l
            j_a = Orb[a].j
            n_a = nn[l_a+1,j_a+1]
            O_a = abs(Rho.n[a,a])
            E_a = SPE.n[a]
        @printf(Summary, "%4d %4d %4d %8.3f %14.5f\tMeV\n",n_a, l_a, j_a, O_a, E_a)
            nn[l_a+1,j_a+1] += 1
        end
        println(Summary,"____________________________________________________________________")

    close(Summary)

    return
end
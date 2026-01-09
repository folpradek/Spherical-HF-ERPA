function HFB_summary(Params::Parameters,E_HFB::Vector{Float64},T_HFB::Float64,Lambda::pnFloat,dA::pnFloat,Convergence::Bool,Iteration::Int64)
    # Read parameters ...
    hw = Params.Calc.hw
    A, Z = Params.Calc.A, Params.Calc.Z
    dZ, dN = dA.p, dA.n
    N_max = Params.Calc.Nmax
    N_2max = Params.Calc.N2max
    N_3max = Params.Calc.N3max
    CMS = Params.Calc.CMS
    Output_File = Params.Calc.Path
    epsilon = Params.Calc.HFB.Tol

    # Export the HF calculation summary ...
    println("\nExporting HFB calculation summary ...")
    Summary =  open(string("IO/", Output_File, "/HFB/HFB_Summary.dat"), "a")
        println(Summary, "Spherical Hartree-Fock-Bogoliubov mean-field calculation summary")
        println(Summary, "\nTarget nuclid:       A = " * string(A) * ", Z = " * string(Z))
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
        println(Summary, "\nPairing interaction parameters:")
        println(Summary, "cP2N     = " * string(Params.Int.cP2N))
        println(Summary, "cP3N     = " * string(Params.Int.cP3N))
        # To be finished ...
        #println(Summary, "\n\tHFB calculation parameters:")
        #println(Summary, "\t\tp2     = " * string(Params.Calc.Pairing.p2))
        #println(Summary, "\t\tp3     = " * string(Params.Calc.Pairing.p3))
        if Convergence == true
            println(Summary, "\nHFB equations converged in " * string(Iteration) * " iterations with precision epsilon = " * string(epsilon) * " ...")
        else
            println(Summary, "\n\nHFB equations did NOT converge in " * string(Iteration) * " iterations with precision epsilon = " * string(epsilon) * " ...")
        end
        println(Summary, "\nSpherical HFB solution review:")
        println(Summary, "\nNumber of HFB Iterations = " * string(Iteration) * ", Convergence precision = " * string(round(epsilon, digits=8)) * " (dE <-> MeV, dN <-> Nucleons)")
        println(Summary, "\nE_HFB      = " * string(round(sum(E_HFB), sigdigits=9)) * "\t MeV \t\t ... \t HFB mean-field ground-state energy")
        println(Summary, "E_MF       = " * string(round(E_HFB[1], sigdigits=9)) * "\t MeV \t\t ... \t mean-field ground-state energy")
        println(Summary, "E_Par      = " * string(round(E_HFB[2] + E_HFB[3], sigdigits=9)) * "\t MeV \t\t ... \t pairing ground-state  energy")
        println(Summary, "E_Par (2N) = " * string(round(E_HFB[2], sigdigits=9)) * "\t MeV \t\t ... \t pairing ground-state  energy")
        println(Summary, "E_Par (3N) = " * string(round(E_HFB[3], sigdigits=9)) * "\t MeV \t\t ... \t pairing ground-state  energy")
        println(Summary, "T_HFB      = " * string(round(T_HFB, sigdigits=9)) * "\t MeV \t\t ... \t HFB mean-field ground-state kinetic energy")
        println(Summary, "\nE_HFB / A      = " * string(round(sum(E_HFB) / A, digits=9)) * "\t MeV / Nucleon \t\t ... \t HFB mean-field ground-state energy per nucleon")
        println(Summary, "E_MF  / A      = " * string(round(E_HFB[1] / A, sigdigits=9)) * "\t MeV / Nucleon \t\t ... \t mean-field ground-state energy per nucleon")
        println(Summary, "E_Par / A      = " * string(round((E_HFB[2] + E_HFB[3]) / A, sigdigits=9)) * "\t MeV / Nucleon \t\t ... \t pairing ground-state  energy per nucleon")
        println(Summary, "E_Par / A (2N) = " * string(round(E_HFB[2] / A, sigdigits=9)) * "\t MeV / Nucleon \t\t ... \t pairing ground-state  energy")
        println(Summary, "E_Par / A (3N) = " * string(round(E_HFB[3] / A, sigdigits=9)) * "\t MeV / Nucleon \t\t ... \t pairing ground-state  energy")
        println(Summary, "T_HFB / A      = " * string(round(T_HFB / A, digits=9)) * "\t MeV / Nucleon\t\t ... \t HFB mean-field ground-state kinetic energy per nucleon")
        println(Summary, "\n\nValues of HFB chemical potentials:")
        println(Summary, "\npLambda =  " * string(Lambda.p) * " MeV")
        println(Summary, "nLambda =  " * string(Lambda.n) * " MeV")
        println(Summary, "\nDispersion of HFB particle numbers dZ & dN:")
        println(Summary, "\ndZ = " * string(round(dZ, sigdigits=6)))
        println(Summary, "dN = " * string(round(dN, sigdigits=6)))
        println(Summary, "\ndZ / Z = " * string(round(dZ / Float64(Z), sigdigits=6)))
        println(Summary, "dN / N = " * string(round(dN / Float64(A - Z), sigdigits=6)))
    close(Summary)

    return
end

function HFB_summary_SQS(Params::Parameters,SQE::pnVector,SQE_C::pnVector,Rho_C::O1B,Orb::Vector{Orb1B})
    # Read parameters
    N_max = Params.Calc.Nmax
    N_2max = Params.Calc.N2max
    Output_File = Params.Calc.Path
    a_max = div((N_max + 1)*(N_max + 2),2)
    J_max = N_2max + 1

    # Initialize temporary counters for radial numbers ...
    pn, nn = zeros(Int64,N_max+1,J_max+1), zeros(Int64,N_max+1,J_max+1)

    # Precompute energy ordered s.p. states for print ...
    Sort_pOrbs = sortperm(SQE.p, by = x -> real(x))
    Sort_nOrbs = sortperm(SQE.n, by = x -> real(x))

    # Export orbitals ...
    println("\nExporting information on single-quasiparticle orbits ...")
    Summary =  open(string("IO/", Output_File, "/HFB/HFB_Summary.dat"), "a")
        println(Summary,"\nProton single-quasiparticle states")
        println(Summary,"_______________________________________________________________________________________________________")
        @printf(Summary, "%4s %4s %4s %10s %14s %14s", "n", "l", "2j", "O_a", "E_C", "E_Q\n")
        @inbounds for Ind_a in 1:a_max
            a = Sort_pOrbs[Ind_a]
            l_a = Orb[a].l
            j_a = Orb[a].j
            n_a = pn[l_a+1,j_a+1]
            @printf(Summary, "%4d %4d %4d %10.3f %14.5f %14.5f\tMeV\n",
                    n_a, l_a, j_a,Rho_C.p[a,a],SQE_C.p[a],SQE.p[a])
            pn[l_a+1,j_a+1] += 1
        end
        println(Summary,"_______________________________________________________________________________________________________")

        println(Summary,"\nNeutron single-quasiparticle states")
        println(Summary,"_______________________________________________________________________________________________________")
        @printf(Summary, "%4s %4s %4s %10s %14s %14s", "n", "l", "2j", "O_a", "E_C", "E_Q\n")
        @inbounds for Ind_a in 1:a_max
            a = Sort_nOrbs[Ind_a]
            l_a = Orb[a].l
            j_a = Orb[a].j
            n_a = nn[l_a+1,j_a+1]
            @printf(Summary, "%4d %4d %4d %10.3f %14.5f %14.5f\tMeV\n",
                    n_a, l_a, j_a,Rho_C.n[a,a],SQE_C.n[a],SQE.n[a])
            nn[l_a+1,j_a+1] += 1
        end
        println(Summary,"_______________________________________________________________________________________________________")

    close(Summary)

    return
end

function HFB_export(Params::Parameters,Orb::Vector{Orb1B},Orb_NN::Orb2B,C::O1B,U::O1B,V::O1B,H_N::qpO1B,H_NN::qpO2B)
    # Read parameters
    Output_File = Params.Calc.Path

    # Export HFB canonical basis matrix C ...
    C_Export_Path = "IO/" * Output_File * "/Bin/C.bin"
    O1b_export(Params,Orb,C,C_Export_Path)

    # Export HFB amplitudes U & V ... in the canonical basis ...
        # Case of U ...
    U_Export_Path = "IO/" * Output_File * "/Bin/U.bin"
    O1b_export(Params,Orb,U,U_Export_Path)

        # Case of V ...
    V_Export_Path = "IO/" * Output_File * "/Bin/V.bin"
    O1b_export(Params,Orb,V,V_Export_Path)

    # Export HFB 1-body Hamiltonian H_N ...
    H1B_Export_Path = "IO/" * Output_File * "/Bin/qpH1B.bin"
    qpO1B_export(Params,Orb,H_N,H1B_Export_Path)

    # Export HFB 2-body Hamiltonian H_NN ..
    H2B_Export_Path = "IO/" * Output_File * "/Bin/qpH2B.bin"
    qpO2b_export(Params,Orb_NN,H_NN,H2B_Export_Path)

    return
end
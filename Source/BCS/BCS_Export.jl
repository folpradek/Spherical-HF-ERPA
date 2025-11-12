function BCS_Summary(Params::Parameters,Params_ref::Parameters,E_MF::Float64,E_BCS::Float64,T_BCS::Float64,Lambda::pnFloat,dA::pnFloat,epsilon::Float64,Iteration::Int64)
    # Read parameters ...
    hw = Params.Calc.hw
    A, Z = Params.Calc.A, Params.Calc.Z
    A_ref, Z_ref = Params_ref.Calc.A, Params_ref.Calc.Z
    dZ, dN = dA.p, dA.n
    N_max = Params.Calc.Nmax
    N_2max = Params.Calc.N2max
    N_3max = Params.Calc.N3max
    CMS = Params.Calc.CMS
    Output_File = Params.Calc.Path

    # Export the HF calculation summary ...
    println("\nExporting BCS calculation summary ...")
    Summary =  open(string("IO/", Output_File, "/BCS/BCS_Summary.dat"), "a")
        println(Summary, "Target nuclid:       A = " * string(A) * ", Z = " * string(Z))
        println(Summary, "Reference nuclid:    A = " * string(A_ref) * ", Z = " * string(Z_ref))
        println(Summary, "\nCalculation data:    hw = " * string(hw) * " , N_max = " * string(N_max) *
                " , N_2max = " * string(N_2max) * " , N_3max = " * string(N_3max) * ", J-basis size = " *
                string(div((N_max+1)*(N_max+2),2)) * ", M-basis size = " * string(div((N_max+1)*(N_max+2)*(N_max+3),6)))
        println(Summary, "                     Center of mass correction option is set to:    " * string(CMS))
        if CMS == "CMS1+2B"
            println(Summary, "\nCombined 1-body + 2-body center of mass motion correction is included ...")
        elseif CMS == "CMS2B"
            println(Summary, "\nOnly pure 2-body center of mass motion correction is included ...")
        else
            println(Summary, "\nNo center of mass motion correction is included ...")
        end
        println(Summary, "\nSpherical BCS solution review:")
        println(Summary, "\nNumber of BCS Iterations = " * string(Iteration) * ", Convergence precision = " * string(round(epsilon, digits=8)) * " Nucleons")
        println(Summary, "\nE_MF  = " * string(round(E_MF, sigdigits=9)) * "\t MeV \t\t ... \t Mean-field ground state energy")
        println(Summary, "E_BCS = " * string(round(E_BCS, sigdigits=9)) * "\t MeV \t\t ... \t BCS pairing ground state energy")
        println(Summary, "T_BCS = " * string(round(T_BCS, sigdigits=9)) * "\t MeV \t\t ... \t BCS mean-field ground-state kinetic energy")
        println(Summary, "\n\nE_MF / A  = " * string(round(E_MF / A, digits=9)) * "\t MeV/Nucleon \t\t ... \t Mean-field ground state energy per nucleon")
        println(Summary, "E_BCS / A = " * string(round(E_BCS / A, digits=9)) * "\t MeV/Nucleon\t\t ... \t BCS pairing ground state energy per nucleon")
        println(Summary, "T_BCS / A = " * string(round(T_BCS / A, digits=9)) * "\t MeV / Nucleon\t\t ... \t BCS mean-field ground-state kinetic energy per nucleon")
        println(Summary, "\n\nValues of BCS chemical potentials ...")
        println(Summary, "pLambda =  " * string(Lambda.p) * " MeV")
        println(Summary, "nLambda =  " * string(Lambda.n) * " MeV")
    close(Summary)

    return
end

function BCS_SQS_Summary(Params::Parameters,SPE::pnVector,SQE::pnVector,U::pnVector,V::pnVector,Orb::Vector{NOrb})
    # Read parameters
    N_max = Params.Calc.Nmax
    N_2max = Params.Calc.N2max
    Output_File = Params.Calc.Path

    a_max = div((N_max + 1)*(N_max + 2),2)
    J_max = N_2max + 1

    pn, nn = zeros(Int64,N_max+1,J_max+1), zeros(Int64,N_max+1,J_max+1)

    # Precompute energy ordered s.p. states for print ...
    Sort_pOrbs = sortperm(SPE.p, by = x -> real(x))
    Sort_nOrbs = sortperm(SPE.n, by = x -> real(x))

    # Export orbitals ...
    println("\nExporting information on single-(quasi)particle orbits ...")
    Summary_File =  open(string("IO/", Output_File, "/BCS/BCS_Summary.dat"), "a")
        println(Summary_File,"\nProton single-(quasi)particle states")
        println(Summary_File,"____________________________________________________________________")
        println(Summary_File,"\tn\t\tl\t\t2j\t\tOcc\t\t\te\t\t\t\t\tU^2\t\t\t\tV^2\t\tE")
        for Ind_a in 1:a_max
            a = Sort_pOrbs[Ind_a]
            l_a = Orb[a].l
            j_a = Orb[a].j
            n_a = pn[l_a+1,j_a+1]
            if Orb[a].pO == 1
                O_a = "h"
            else
                O_a = "p"
            end
            Row = "\t" * string(n_a) * "\t\t" * string(l_a) * "\t\t" * string(j_a) * "\t\t" *  O_a * "\t\t" * string(round(SPE.p[a], sigdigits = 9)) * "\t\t\tMeV\t\t" *
                    string(round(U.p[a]^2, sigdigits = 5)) * "\t\t" * string(round(V.p[a]^2, sigdigits = 5)) * "\t\t" * string(round(SQE.p[a], sigdigits = 9)) * "\t\t\tMeV"
            println(Summary_File, Row)
            pn[l_a+1,j_a+1] += 1
        end
        println(Summary_File,"____________________________________________________________________")

        println(Summary_File,"\nNeutron single-(quasi)particle states")
        println(Summary_File,"____________________________________________________________________")
        println(Summary_File,"\tn\t\tl\t\t2j\t\tOcc\t\t\te\t\t\t\t\tU^2\t\t\t\tV^2\t\tE")
        for Ind_a in 1:a_max
            a = Sort_nOrbs[Ind_a]
            l_a = Orb[a].l
            j_a = Orb[a].j
            n_a = nn[l_a+1,j_a+1]
            if Orb[a].nO == 1
                O_a = "h"
            else
                O_a = "p"
            end
            Row = "\t" * string(n_a) * "\t\t" * string(l_a) * "\t\t" * string(j_a) * "\t\t" *  O_a * "\t\t" * string(round(SPE.n[a], sigdigits = 9)) * "\t\t\tMeV\t\t" *
                    string(round(U.n[a]^2, sigdigits = 5)) * "\t\t" * string(round(V.n[a]^2, sigdigits = 5)) * "\t\t" * string(round(SQE.n[a], sigdigits = 9)) * "\t\t\tMeV"
            println(Summary_File, Row)
            nn[l_a+1,j_a+1] += 1
        end
        println(Summary_File,"____________________________________________________________________")

    close(Summary_File)

    return
end
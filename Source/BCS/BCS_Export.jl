function BCS_Summary(Params::Parameters,Params_ref::Parameters,E_MF::Float64,E_BCS::Float64,lambda::pnFloat,dA::pnFloat,epsilon::Float64,Iteration::Int64)
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
    Summary_File =  open(string("IO/", Output_File, "/BCS/BCS_Summary.dat"), "a")
        println(Summary_File, "Target nuclid:       A = " * string(A) * ", Z = " * string(Z))
        println(Summary_File, "Reference nuclid:    A = " * string(A_ref) * ", Z = " * string(Z_ref))
        println(Summary_File, "\nCalculation data:    hw = " * string(hw) * " , N_max = " * string(N_max) *
                " , N_2max = " * string(N_2max) * " , N_3max = " * string(N_3max) * ", J-basis size = " *
                string(div((N_max+1)*(N_max+2),2)) * ", M-basis size = " * string(div((N_max+1)*(N_max+2)*(N_max+3),6)))
        println(Summary_File, "                     Center of mass correction option is set to:    " * string(CMS))
        if CMS == "CMS1+2B"
            println(Summary_File, "\nCombined 1-body + 2-body center of mass motion correction is included ...")
        elseif CMS == "CMS2B"
            println(Summary_File, "\nOnly pure 2-body center of mass motion correction is included ...")
        else
            println(Summary_File, "\nNo center of mass motion correction is included ...")
        end
        println(Summary_File, "\nSpherical BCS solution review:")
        println(Summary_File, "\nNumber of BCS Iterations = " * string(Iteration) * ", Convergence precision = " * string(round(epsilon, digits=8)) * " Nucleons")
        println(Summary_File, "\nE_MF  = " * string(round(E_MF, sigdigits=9)) * "\t MeV \t\t ... \t Mean-field ground state energy")
        println(Summary_File, "E_BCS = " * string(round(E_BCS, sigdigits=9)) * "\t MeV \t\t ... \t BCS pairing ground state energy")
        println(Summary_File, "\n\nE_MF / A  = " * string(round(E_MF / A, digits=9)) * "\t MeV/Nucleon \t\t ... \t Mean-field ground state energy per nucleon")
        println(Summary_File, "E_BCS / A = " * string(round(E_BCS / A, digits=9)) * "\t MeV/Nucleon\t\t ... \t BCS pairing ground state energy per nucleon")
        println(Summary_File, "\n\nValues of BCS chemical potentials ...")
        println(Summary_File, "pLambda =  " * string(lambda.p) * " MeV")
        println(Summary_File, "nLambda =  " * string(lambda.n) * " MeV")
    close(Summary_File)

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
                    string(round(U.p[a], sigdigits = 5)) * "\t\t" * string(round(V.p[a], sigdigits = 5)) * "\t\t" * string(round(SQE.p[a], sigdigits = 9)) * "\t\t\tMeV"
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
                    string(round(U.n[a], sigdigits = 5)) * "\t\t" * string(round(V.n[a], sigdigits = 5)) * "\t\t" * string(round(SQE.n[a], sigdigits = 9)) * "\t\t\tMeV"
            println(Summary_File, Row)
            nn[l_a+1,j_a+1] += 1
        end
        println(Summary_File,"____________________________________________________________________")

    close(Summary_File)

    return
end
# To be refined ...
function BCS_Export(Params::Parameters,U::pnMatrix,SPEnergies::pnVector)
    # Read parameters
    N_max = Params.Calc.Nmax
    Output_File = Params.Calc.Path

    a_max = div((N_max + 1)*(N_max + 2),2)

    # Export densities ...
    pU_Export = "IO/" * Output_File * "/Bin/pU_HF.bin"
    open(pU_Export, "w") do Export_File
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                ME = @views U.p[a,b]
                write(Export_File, Float64(ME))
            end
        end
    end

    nU_Export = "IO/" * Output_File * "/Bin/nU_HF.bin"
    open(nU_Export, "w") do Export_File
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                ME = @views U.n[a,b]
                write(Export_File, Float64(ME))
            end
        end
    end

    pE_Export = "IO/" * Output_File * "/Bin/pE_HF.bin"
    open(pE_Export, "w") do Export_File
        @inbounds for a in 1:a_max
            ME = @views SPEnergies.p[a]
            write(Export_File, Float64(ME))
        end
    end

    nE_Export = "IO/" * Output_File * "/Bin/nE_HF.bin"
    open(nE_Export, "w") do Export_File
        @inbounds for a in 1:a_max
            ME = @views SPEnergies.n[a]
            write(Export_File, Float64(ME))
        end
    end

    return
end
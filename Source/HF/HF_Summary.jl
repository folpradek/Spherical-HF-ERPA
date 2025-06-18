function HF_Summary(Params::Vector{Any},E_HF::Float64,T_HF::Float64,Epsilon::Float64,Iteration::Int64)
    HbarOmega = Params[1]
    A = Params[5]
    Z = Params[6]
    N_max = Params[7]
    N_2max = Params[8]
    N_3max = Params[9]
    CMS = Params[10]
    Output_File = Params[13]

    println("\nExporting HF calculation summary ...")
    Summary_File =  open(string("IO/", Output_File, "/HF/HF_Summary.dat"), "a")
        println(Summary_File, "Nuclid data:    A = " * string(A) * ", Z = " * string(Z))
        println(Summary_File, "\nCalculation data:    HbarOmega = " * string(HbarOmega) * " , N_max = " * string(N_max) *
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
        println(Summary_File, "\nSpherical Hartree-Fock solution review:")
        println(Summary_File, "\nNumber of Iterations = " * string(Iteration) * ", Precision = " * string(round(Epsilon, digits=8)) * " MeV")
        println(Summary_File, "\nE_HF = " * string(round(E_HF, sigdigits=9)) * "\t MeV \t\t ... \t Mean-field ground state energy")
        println(Summary_File, "T_HF = " * string(round(T_HF, sigdigits=9)) * "\t MeV \t\t ... \t Mean-field total kinetic energy")
    close(Summary_File)

    return
end

function HF_SPS_Summary(Params::Vector{Any},pSPEnergies::Vector{Float64},nSPEnergies::Vector{Float64},Orb::Vector{NOrb})
    N_max = Params[7]
    N_2max = Params[8]
    Output_File = Params[13]

    a_max = div((N_max + 1)*(N_max + 2),2)
    J_max = N_2max + 1

    pn = zeros(Int64,N_max+1,J_max+1)
    nn = zeros(Int64,N_max+1,J_max+1)

    # Precompute energy ordered s.p. states for print ...
    Sort_pOrbs = sortperm(pSPEnergies, by = x -> real(x))
    Sort_nOrbs = sortperm(nSPEnergies, by = x -> real(x))


    println("\nExporting HF calculation summary ...")
    Summary_File =  open(string("IO/", Output_File, "/HF/HF_Summary.dat"), "a")
        println(Summary_File,"\nProton single-particle states")
        println(Summary_File,"____________________________________________________________________")
        println(Summary_File,"\tn\t\tl\t\t2j\t\tOcc\t\t\tE")
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
            E_a = pSPEnergies[a]
            Row = "\t" * string(n_a) * "\t\t" * string(l_a) * "\t\t" * string(j_a) * "\t\t" *  O_a * "\t\t" * string(round(E_a, sigdigits = 9)) * "\t\t\tMeV"
            println(Summary_File, Row)
            pn[l_a+1,j_a+1] += 1
        end
        println(Summary_File,"____________________________________________________________________")

        println(Summary_File,"\nNeutron single-particle states")
        println(Summary_File,"____________________________________________________________________")
        println(Summary_File,"\tn\t\tl\t\t2j\t\tOcc\t\t\tE")
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
            E_a = nSPEnergies[a]
            Row = "\t" * string(n_a) * "\t\t" * string(l_a) * "\t\t" * string(j_a) * "\t\t" *  O_a * "\t\t" * string(round(E_a, sigdigits = 9)) * "\t\t\tMeV"
            println(Summary_File, Row)
            nn[l_a+1,j_a+1] += 1
        end
        println(Summary_File,"____________________________________________________________________")

    close(Summary_File)

    return
end
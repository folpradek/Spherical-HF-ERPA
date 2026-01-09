function HF_summary(Params::Parameters,E_HF::Float64,T_HF::Float64,Iteration::Int64)
    # Read parameters ...
    hw = Params.Calc.hw
    A = Params.Calc.A
    Z = Params.Calc.Z
    N_max = Params.Calc.Nmax
    N_2max = Params.Calc.N2max
    N_3max = Params.Calc.N3max
    CMS = Params.Calc.CMS
    Output_File = Params.Calc.Path
    epsilon = Params.Calc.HF.Tol

    # Export the HF calculation summary ...
    println("\nExporting HF calculation summary ...")
    Summary =  open(string("IO/", Output_File, "/HF/HF_Summary.dat"), "a")
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
        println(Summary, "\nSpherical Hartree-Fock solution review:")
        println(Summary, "\nNumber of Iterations = " * string(Iteration) * ", Precision = " * string(round(epsilon, digits=8)) * " MeV")
        println(Summary, "\nE_HF = " * string(round(E_HF, sigdigits=9)) * "\t MeV \t\t ... \t Mean-field ground state energy")
        println(Summary, "T_HF = " * string(round(T_HF, sigdigits=9)) * "\t MeV \t\t ... \t Mean-field total kinetic energy")
        println(Summary, "\nE_HF / A = " * string(round(E_HF / Float64(A), sigdigits=9)) * "\t MeV / Nucleon \t\t ... \t Mean-field ground state energy")
        println(Summary, "T_HF / A = " * string(round(T_HF / Float64(A), sigdigits=9)) * "\t MeV / Nucleon \t\t ... \t Mean-field total kinetic energy")
    close(Summary)

    return
end

function HF_SPS_summary(Params::Parameters,h::O1B,Orb::Vector{Orb1B})
    # Read parameters
    N_max = Params.Calc.Nmax
    N_2max = Params.Calc.N2max
    Output_File = Params.Calc.Path
    a_max = div((N_max + 1)*(N_max + 2),2)
    J_max = N_2max + 1

    # Read the diagonal part of h ...
    pSPE, nSPE = diag(h.p), diag(h.n)

    # Initialize array for radial number n ...
    pn, nn = zeros(Int64,N_max+1,J_max+1), zeros(Int64,N_max+1,J_max+1)

    # Precompute energy ordered s.p. states for print ...
    Sort_pOrbs = sortperm(pSPE, by = x -> real(x))
    Sort_nOrbs = sortperm(nSPE, by = x -> real(x))

    # Export orbitals ...
    println("\nExporting HF calculation summary ...")
    Summary =  open(string("IO/", Output_File, "/HF/HF_Summary.dat"), "a")
        println(Summary,"\nProton single-particle states")
        println(Summary,"____________________________________________________________________")
        println(Summary,"\tn\t\tl\t\t2j\t\tOcc\t\t\tE")
        @inbounds for Ind_a in 1:a_max
            a = Sort_pOrbs[Ind_a]
            l_a, j_a = Orb[a].l, Orb[a].j
            n_a = pn[l_a+1,j_a+1]
            if Orb[a].pO == 1
                O_a = "h"
            else
                O_a = "p"
            end
            E_a = pSPE[a]
            Row = "\t" * string(n_a) * "\t\t" * string(l_a) * "\t\t" * string(j_a) * "\t\t" *  O_a * "\t\t" * string(round(E_a, sigdigits = 9)) * "\t\t\tMeV"
            println(Summary, Row)
            pn[l_a+1,j_a+1] += 1
        end
        println(Summary,"____________________________________________________________________")

        println(Summary,"\nNeutron single-particle states")
        println(Summary,"____________________________________________________________________")
        println(Summary,"\tn\t\tl\t\t2j\t\tOcc\t\t\tE")
        for Ind_a in 1:a_max
            a = Sort_nOrbs[Ind_a]
            l_a, j_a = Orb[a].l, Orb[a].j
            n_a = nn[l_a+1,j_a+1]
            if Orb[a].nO == 1
                O_a = "h"
            else
                O_a = "p"
            end
            E_a = nSPE[a]
            Row = "\t" * string(n_a) * "\t\t" * string(l_a) * "\t\t" * string(j_a) * "\t\t" *  O_a * "\t\t" * string(round(E_a, sigdigits = 9)) * "\t\t\tMeV"
            println(Summary, Row)
            nn[l_a+1,j_a+1] += 1
        end
        println(Summary,"____________________________________________________________________")

    close(Summary)

    return
end

function HF_export(Params::Parameters,C::O1B,h::O1B,Orb::Vector{Orb1B})
    # Read parameters
    N_max = Params.Calc.Nmax
    Output_File = Params.Calc.Path
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Export HF orbitals matrix C ...
    C_Export = "IO/" * Output_File * "/Bin/C_HF.bin"
    open(C_Export, "w") do Export_File
        # Proton orbitals ...
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                if Orb[a].l == Orb[b].l &&Orb[a].j == Orb[b].j
                    ME = @views C.p[a,b]
                    write(Export_File, Float64(ME))
                end
            end
        end

        # Neutron orbitals ...
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                if Orb[a].l == Orb[b].l &&Orb[a].j == Orb[b].j
                    ME = @views C.n[a,b]
                    write(Export_File, Float64(ME))
                end
            end
        end
    end

    # Export HF 1-body Hamiltonian h ... 
    h_Export = "IO/" * Output_File * "/Bin/h_HF.bin"
    open(h_Export, "w") do Export_File
        # Proton Hamiltonian ...
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                if Orb[a].l == Orb[b].l &&Orb[a].j == Orb[b].j
                    ME = @views h.p[a,b]
                    write(Export_File, Float64(ME))
                end
            end
        end

        # Neutron Hamiltonian ...
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                if Orb[a].l == Orb[b].l &&Orb[a].j == Orb[b].j
                    ME = @views h.n[a,b]
                    write(Export_File, Float64(ME))
                end
            end
        end

    end

    return
end
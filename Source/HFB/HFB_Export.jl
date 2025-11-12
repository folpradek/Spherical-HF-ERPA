function HFB_Summary(Params::Parameters,E_HFB::Float64,T_HFB::Float64,Lambda::pnFloat,dA::pnFloat,epsilon::Float64,Iteration::Int64)
    # Read parameters ...
    hw = Params.Calc.hw
    A, Z = Params.Calc.A, Params.Calc.Z
    dZ, dN = dA.p, dA.n
    N_max = Params.Calc.Nmax
    N_2max = Params.Calc.N2max
    N_3max = Params.Calc.N3max
    CMS = Params.Calc.CMS
    Output_File = Params.Calc.Path

    # Export the HF calculation summary ...
    println("\nExporting HFB calculation summary ...")
    Summary =  open(string("IO/", Output_File, "/HFB/HFB_Summary.dat"), "a")
        println(Summary, "Target nuclid:       A = " * string(A) * ", Z = " * string(Z))
        println(Summary, "\nCalculation data:    hw = " * string(hw) * " , N_max = " * string(N_max) *
                " , N_2max = " * string(N_2max) * " , N_3max = " * string(N_3max) * ", J-basis size = " *
                string(div((N_max+1)*(N_max+2),2)) * ", M-basis size = " * string(div((N_max+1)*(N_max+2)*(N_max+3),6)))
        println(Summary, "                     Center of mass correction option is set to:    " * string(CMS))
        if CMS == "CMS1+2B"
            println(Summary, "\nCombined 1-body + 2-body Center-of-Mass (CM) motion correction is included ...")
        elseif CMS == "CMS2B"
            println(Summary, "\nOnly pure 2-body Center-of-Mass (CM) motion correction is included ...")
        else
            println(Summary, "\nNo Center-of-Mass (CM) motion correction is included ...")
        end
        println(Summary, "\nSpherical HFB solution review:")
        println(Summary, "\nNumber of HFB Iterations = " * string(Iteration) * ", Convergence precision = " * string(round(epsilon, digits=8)) * " (dE <-> MeV, dN <-> Nucleons)")
        println(Summary, "\nE_HFB  = " * string(round(E_HFB, sigdigits=9)) * "\t MeV \t\t ... \t HFB mean-field ground-state energy")
        println(Summary, "T_HFB = " * string(round(T_HFB, sigdigits=9)) * "\t MeV \t\t ... \t HFB mean-field ground-state kinetic energy")
        println(Summary, "\n\nE_HFB / A  = " * string(round(E_HFB / A, digits=9)) * "\t MeV / Nucleon \t\t ... \t HFB mean-field ground-state energy per nucleon")
        println(Summary, "T_HFB / A = " * string(round(T_HFB / A, digits=9)) * "\t MeV / Nucleon\t\t ... \t HFB mean-field ground-state kinetic energy per nucleon")
        println(Summary, "\n\nValues of HFB chemical potentials:")
        println(Summary, "\npLambda =  " * string(Lambda.p) * " MeV")
        println(Summary, "nLambda =  " * string(Lambda.n) * " MeV")
        println(Summary, "\nDispersion of HFB particle numbers dZ & dN:")
        println(Summary, "\ndZ = " * string(round(dA.p, sigdigits=6)))
        println(Summary, "dN = " * string(round(dA.n, sigdigits=6)))
        println(Summary, "\ndZ / Z = " * string(round(dA.p / Float64(Z), sigdigits=6)))
        println(Summary, "dN / N = " * string(round(dA.n / Float64(A - Z), sigdigits=6)))
    close(Summary)

    return
end

function HFB_SQS_Summary(Params::Parameters,SQE::pnVector,SQE_C::pnVector,Rho_C::pnMatrix,Orb::Vector{NOrb})
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
        println(Summary,"\tn\t\tl\t\t2j\t\tO_a\t\t\tE_C & E_Q")
        @inbounds for Ind_a in 1:a_max
            a = Sort_pOrbs[Ind_a]
            l_a = Orb[a].l
            j_a = Orb[a].j
            n_a = pn[l_a+1,j_a+1]
            Row = "\t" * string(n_a) * "\t\t" * string(l_a) * "\t\t" * string(j_a) * "\t\t" * string(round(Rho_C.p[a,a], digits = 3)) *
                  "\t\t" * string(round(SQE_C.p[a], sigdigits = 5)) * "\t\t" * string(round(SQE.p[a], sigdigits = 5)) * "\t\t\tMeV"
            println(Summary, Row)
            pn[l_a+1,j_a+1] += 1
        end
        println(Summary,"_______________________________________________________________________________________________________")

        println(Summary,"\nNeutron single-quasiparticle states")
        println(Summary,"_______________________________________________________________________________________________________")
        println(Summary,"\tn\t\tl\t\t2j\t\tO_a\t\t\tE_C & E_Q")
        @inbounds for Ind_a in 1:a_max
            a = Sort_nOrbs[Ind_a]
            l_a = Orb[a].l
            j_a = Orb[a].j
            n_a = nn[l_a+1,j_a+1]
            Row = "\t" * string(n_a) * "\t\t" * string(l_a) * "\t\t" * string(j_a) * "\t\t" * string(round(Rho_C.n[a,a], digits = 3)) *
                  "\t\t" * string(round(SQE_C.n[a], sigdigits = 5)) * "\t\t" * string(round(SQE.n[a], sigdigits = 5)) * "\t\t\tMeV"
            println(Summary, Row)
            nn[l_a+1,j_a+1] += 1
        end
        println(Summary,"_______________________________________________________________________________________________________")

    close(Summary)

    return
end

# Yet to be defined ...
function HFB_Export(Params::Parameters,C::pnMatrix,U::pnMatrix,V::pnMatrix,SQE::pnVector)
    # Read parameters
    N_max = Params.Calc.Nmax
    Output_File = Params.Calc.Path

    a_max = div((N_max + 1)*(N_max + 2),2)

    # Export densities ...
    pU_Export = "IO/" * Output_File * "/Bin/pC.bin"
    open(pU_Export, "w") do Export_File
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                ME = @views U.p[a,b]
                write(Export_File, Float64(ME))
            end
        end
    end

    nU_Export = "IO/" * Output_File * "/Bin/nC.bin"
    open(nU_Export, "w") do Export_File
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                ME = @views U.n[a,b]
                write(Export_File, Float64(ME))
            end
        end
    end

    pE_Export = "IO/" * Output_File * "/Bin/pE.bin"
    open(pE_Export, "w") do Export_File
        @inbounds for a in 1:a_max
            ME = @views SPEnergies.p[a]
            write(Export_File, Float64(ME))
        end
    end

    nE_Export = "IO/" * Output_File * "/Bin/nE.bin"
    open(nE_Export, "w") do Export_File
        @inbounds for a in 1:a_max
            ME = @views SPEnergies.n[a]
            write(Export_File, Float64(ME))
        end
    end

    return
end
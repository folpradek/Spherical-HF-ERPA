function HFB_Summary(Params::Parameters,E_HFB::Vector{Float64},T_HFB::Float64,Lambda::pnFloat,dA::pnFloat,epsilon::Float64,Iteration::Int64)
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
        println(Summary, "\nE_HFB  = " * string(round(sum(E_HFB), sigdigits=9)) * "\t MeV \t\t ... \t HFB mean-field ground-state energy")
        println(Summary, "\tE_MF  = " * string(round(E_HFB[1], sigdigits=9)) * "\t MeV \t\t ... \t mean-field ground-state energy")
        println(Summary, "\tE_Par = " * string(round(E_HFB[2], sigdigits=9)) * "\t MeV \t\t ... \t pairing ground-state  energy")
        println(Summary, "\tT_HFB = " * string(round(T_HFB, sigdigits=9)) * "\t MeV \t\t ... \t HFB mean-field ground-state kinetic energy")
        println(Summary, "\n\nE_HFB / A  = " * string(round(sum(E_HFB) / A, digits=9)) * "\t MeV / Nucleon \t\t ... \t HFB mean-field ground-state energy per nucleon")
        println(Summary, "\tE_MF  / A = " * string(round(E_HFB[1] / A, sigdigits=9)) * "\t MeV \t\t ... \t mean-field ground-state energy per nucleon")
        println(Summary, "\tE_Par / A = " * string(round(E_HFB[2] / A, sigdigits=9)) * "\t MeV \t\t ... \t pairing ground-state  energy per nucleon")
        println(Summary, "\tT_HFB / A = " * string(round(T_HFB / A, digits=9)) * "\t MeV / Nucleon\t\t ... \t HFB mean-field ground-state kinetic energy per nucleon")
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
        #println(Summary,"\tn\t\tl\t\t2j\t\tO_a\t\t\tE_C & E_Q")
        @printf(Summary, "%4s %4s %4s %10s %14s %14s", "n", "l", "2j", "O_a", "E_C", "E_Q\n")
        @inbounds for Ind_a in 1:a_max
            a = Sort_pOrbs[Ind_a]
            l_a = Orb[a].l
            j_a = Orb[a].j
            n_a = pn[l_a+1,j_a+1]
            #Row = "\t" * string(n_a) * "\t\t" * string(l_a) * "\t\t" * string(j_a) * "\t\t" * string(round(Rho_C.p[a,a], digits = 3)) *
            #      "\t\t" * string(round(SQE_C.p[a], sigdigits = 5)) * "\t\t" * string(round(SQE.p[a], sigdigits = 5)) * "\t\t\tMeV"
            #println(Summary, Row)
            @printf(Summary, "%4d %4d %4d %10.3f %14.5f %14.5f\tMeV\n",
                    n_a, l_a, j_a,Rho_C.p[a,a],SQE_C.p[a],SQE.p[a])
            pn[l_a+1,j_a+1] += 1
        end
        println(Summary,"_______________________________________________________________________________________________________")

        println(Summary,"\nNeutron single-quasiparticle states")
        println(Summary,"_______________________________________________________________________________________________________")
        #println(Summary,"\tn\tl\t2j\t\tO_a\t\t\tE_C\t\t\tE_Q")
        @printf(Summary, "%4s %4s %4s %10s %14s %14s", "n", "l", "2j", "O_a", "E_C", "E_Q\n")
        @inbounds for Ind_a in 1:a_max
            a = Sort_nOrbs[Ind_a]
            l_a = Orb[a].l
            j_a = Orb[a].j
            n_a = nn[l_a+1,j_a+1]
            #Row = "\t" * string(n_a) * "\t\t" * string(l_a) * "\t\t" * string(j_a) * "\t\t" * string(round(Rho_C.n[a,a], digits = 3)) *
            #      "\t\t" * string(round(SQE_C.n[a], sigdigits = 5)) * "\t\t" * string(round(SQE.n[a], sigdigits = 5)) * "\t\t\tMeV"
            #println(Summary, Row)
            @printf(Summary, "%4d %4d %4d %10.3f %14.5f %14.5f\tMeV\n",
                    n_a, l_a, j_a,Rho_C.n[a,a],SQE_C.n[a],SQE.n[a])
            nn[l_a+1,j_a+1] += 1
        end
        println(Summary,"_______________________________________________________________________________________________________")

    close(Summary)

    return
end

# To be finished on ...
function HFB_Export(Params::Parameters,Orb_NN::NNOrb,C::pnMatrix,U::pnMatrix,V::pnMatrix,H_N::qpH1B,H_NN::qpH2B)
    # Read parameters
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1)*(N_max + 2),2)
    Output_File = Params.Calc.Path

    # Export HFB canonical basis matrix C ...
    C_Export = "IO/" * Output_File * "/Bin/C.bin"
    open(C_Export, "w") do Export_File
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a
                ME = @views C.p[a,b]
                write(Export_File, Float64(ME))
            end
        end
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a
                ME = @views C.n[a,b]
                write(Export_File, Float64(ME))
            end
        end
    end

    # Export HFB amplitudes U & V ... in the canonical basis ...
        # Case of U ...
    U_Export = "IO/" * Output_File * "/Bin/U.bin"
    open(U_Export, "w") do Export_File
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                ME = @views U.p[a,b]
                write(Export_File, Float64(ME))
            end
        end
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                ME = @views V.p[a,b]
                write(Export_File, Float64(ME))
            end
        end
    end
        # Case of V ...
    V_Export = "IO/" * Output_File * "/Bin/V.bin"
    open(V_Export, "w") do Export_File
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                ME = @views U.n[a,b]
                write(Export_File, Float64(ME))
            end
        end
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                ME = @views V.n[a,b]
                write(Export_File, Float64(ME))
            end
        end
    end

    # Export HFB 1-body Hamiltonian H1B ...
    H1B_Export = "IO/" * Output_File * "/Bin/H_N.bin"
    open(H1B_Export, "w") do Export_File
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a
                ME = @views H_N.H11.p[a,b]
                write(Export_File, Float64(ME))
            end
        end
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a
                ME = @views H_N.H11.n[a,b]
                write(Export_File, Float64(ME))
            end
        end
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a
                ME = @views H_N.H20.p[a,b]
                write(Export_File, Float64(ME))
            end
        end
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a
                ME = @views H_N.H20.n[a,b]
                write(Export_File, Float64(ME))
            end
        end
    end

    #qpH2b_export(Params,Orb_NN,H_NN)

    return
end

function qpH2b_export(Params::Parameters,Orb_NN::NNOrb,H_NN::qpH2B)
    # Read parameters ...
    N_max = Params.Calc.Nmax
    N_2max = 2*N_max
    J_max = N_2max + 1
    Output_File = Params.Calc.Path

    H2B_Export = "IO/" * Output_File * "/Bin/H2B.bin"

    # Export residual 2-body interaction ...

    println("\nExporting residual 2-body quasiparticle Hamiltonian H_NN into output binary file ... ''" * string(H2B_Export) * "''")
    open(H2B_Export, "w") do Export_File
        # ppH40
        @inbounds for J in 0:J_max
            @inbounds for P in 1:2
                N_T1 = Orb_NN.N[2,P,J+1]
                @inbounds for Bra in 1:N_T1
                    @inbounds for Ket in 1:Bra
                        Ind = Bra + (Ket - 1) * N_T1 - div(Ket * (Ket - 1),2)
                        ME = @views H_NN.H40.pp[P,J+1][Ind]
                        write(Export_File, Float64(ME))
                    end
                end
            end
        end

        # pnH40
        @inbounds for J in 0:J_max
            @inbounds for P in 1:2
                N_T0 = Orb_NN.N[1,P,J+1]
                @inbounds for Bra in 1:N_T0
                    @inbounds for Ket in 1:N_T0
                        Ind = Bra + (Ket - 1) * N_T0 - div(Ket * (Ket - 1),2)
                        ME = @views H_NN.H40.pn[P,J+1][Ind]
                        write(Export_File, Float64(ME))
                    end
                end
            end
        end

        # nnH40
        @inbounds for J in 0:J_max
            @inbounds for P in 1:2
                N_T1 = Orb_NN.N[2,P,J+1]
                @inbounds for Bra in 1:N_T1
                    @inbounds for Ket in 1:Bra
                        Ind = Bra + (Ket - 1) * N_T1 - div(Ket * (Ket - 1),2)
                        ME = @views H_NN.H40.nn[P,J+1][Ind]
                        write(Export_File, Float64(ME))
                    end
                end
            end
        end

        # ppH31
        @inbounds for J in 0:J_max
            @inbounds for P in 1:2
                N_T1 = Orb_NN.N[2,P,J+1]
                @inbounds for Bra in 1:N_T1
                    @inbounds for Ket in 1:Bra
                        Ind = Bra + (Ket - 1) * N_T1 - div(Ket * (Ket - 1),2)
                        ME = @views H_NN.H40.pp[P,J+1][Ind]
                        write(Export_File, Float64(ME))
                    end
                end
            end
        end

    end

    println("\nResidual 2-body interaction succesfully exported into file V2B_res_HF.bin ...")

    return
end
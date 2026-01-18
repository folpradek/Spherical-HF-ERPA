function QTDA_export(Params::Parameters,Orb::Vector{Orb1B},Orb_2qp::qpOrb2B,E_QTDA::Matrix{Vector{Float64}},X_QTDA::Matrix{Matrix{Float64}},rB_QTDA::ReducedTransition)

    # Export of QTDA summary ...
    @time QTDA_summary(Params,Orb_2qp)

    # Export of QTDA spectra ...
    @time QTDA_spectrum_export(Params,Orb_2qp,E_QTDA)

    # Export of QTDA plot-ready spectra ...
    @time QTDA_spectrum_export_plot(Params,Orb_2qp,E_QTDA)

    # Export of QTDA |X|^2 amplitudes, not plot ready ...
    #@time QTDA_amplitudes_export(Params,N_nu,E_RPA,X_RPA,Y_RPA)

    # Export oQTDA electric transitions ...
    @time QTDA_transitions_export(Params,Orb_2qp,E_QTDA,rB_QTDA)

    # Export dimensions of phonon subspaces ...
    @time QTDA_phonon_space_export(Params,Orb_2qp)
    
    return
end

function QTDA_summary(Params::Parameters,Orb_2qp::qpOrb2B)
    # Read parameters ...
    A = Params.Calc.A
    Z = Params.Calc.Z
    hw = Params.Calc.hw
    N_max = Params.Calc.Nmax
    N_2max = 2*N_max
    J_max = N_2max + 1
    CMS = Params.Calc.CMS
    #Orthogon = Params.Calc.QTDA.Ortho
    Output_File = Params.Calc.Path

    # Set the path for "HF_RPA_Summary.dat" output file ...
    Output_Path = "IO/" * Output_File * "/QTDA/QTDA_Summary.dat"

    println("\nPreparing QTDA summary ...")
    Summary =  open(Output_Path, "w")
        println(Summary, "Spherical quasiparticle mean-field Quasiparticle Tamm-Dancoff Approximation review:")
        println(Summary, "\nNuclid data:    A = " * string(A) * ", Z = " * string(Z))
        println(Summary, "\nCalculation data:    hw = " * string(hw) * " , N_max = " * string(N_max) *
                ", s.p. J-basis size = " * string(div((N_max+1)*(N_max+2),2)) * ", s.p. M-basis size = " * string(div((N_max+1)*(N_max+2)*(N_max+3),6)))
        println(Summary, "\t\t\tCenter-of-Mass System correction is set to:   CMS = ''" * string(CMS) * "''")
        if CMS == "CMS1+2B"
            println(Summary, "\n\t\t\tCombined 1-body + 2-body Center-of-Mass System (CMS) motion correction is set ...")
        elseif CMS == "CMS2B"
            println(Summary, "\n\t\t\tOnly pure 2-body Center-of-Mass System (CMS) motion correction is set ...")
        else
            println(Summary, "\n\t\t\tNo Center-of-Mass System (CMS) motion correction is set ...")
        end
        #if Orthogon == true
        #    println(Summary, "\nOrthogonalization of 1- spurious state is included ...")
        #else
        #    println(Summary, "\nOrthogonalization of 1- spurious state is NOT included ...")
        #end
        N_qp = 0
        @inbounds for J in 0:J_max
            @inbounds for P in 1:2
                n_qp = Orb_2qp.N[P,J+1] * (2*J + 1)
                N_qp += n_qp
            end
        end
        println(Summary, "\nTotal number of 1-phonon quasiparticle states in J-scheme = " * string(sum(Orb_2qp.N)))
        println(Summary, "Total number of 1-phonon quasiparticle states in M-scheme = " * string(N_qp))
    close(Summary)

    println("\nQTDA summary exported ...")

    return
end

function QTDA_spectrum_export(Params::Parameters,Orb_2qp::qpOrb2B,E_QTDA::Matrix{Vector{Float64}})
    # Read parameters ...
    N_max = Params.Calc.Nmax
    J_max = 2*N_max + 1
    Output_File = Params.Calc.Path

    # Set the out-put file paths ...
    Output_Path = "IO/" * Output_File * "/QTDA/Spectra/QTDA.dat"

    # QTDA spectrum export
    println("\nPreparing export of QTDA spectra ...")
    open(Output_Path, "w") do Write_File
        println(Write_File, "J\tP\tE")
        @inbounds for J in 0:J_max
            @inbounds for P in 1:2
                N_qp = Orb_2qp.N[P,J+1]
                @inbounds for nu in 1:N_qp
                    if P == 1
                        println(Write_File, string(J) * "\t" * "+" * "\t" * string(E_QTDA[P,J+1][nu]))
                    else
                        println(Write_File, string(J) * "\t" * "-" * "\t" * string(E_QTDA[P,J+1][nu]))
                    end
                end
                println(Write_File,"\n")
            end
        end
    end

    println("\nExport of QTDA spectra succesfully finished ...")
    return
end

function QTDA_spectrum_export_plot(Params::Parameters,Orb_2qp::qpOrb2B,E_QTDA::Matrix{Vector{Float64}})
    # Read parameters ...
    N_max = Params.Calc.Nmax
    J_max = 2*N_max + 1
    IO = Params.Calc.Path

    # Set the out-put file paths ...
    Output_Path = "IO/" * IO * "/QTDA/Spectra/QTDA_Plot.dat"

    # Spectra label gap parameter ...
    Delta = 1.0

    # Prepare the QTDA plot-ready spectra for export ...
    println("\nPreparing plot-ready export of QTDA spectra ...")

    N_qp = Int64(sum(Orb_2qp.N))

    QTDA_export = Matrix{Float64}(undef,4,N_qp+1)

    # 0+ ground state addition ...
    QTDA_export[1,1] = 0.0
    QTDA_export[2,1] = 1.0
    QTDA_export[3,1] = 0.0
    QTDA_export[4,1] = 0.0

    qp_count = 1
    @inbounds for J in 0:J_max
        @inbounds for P in 1:2
            N_ab = Orb_2qp.N[P,J+1]
            @inbounds for nu in 1:N_ab
                qp_count += 1
                QTDA_export[1,qp_count] = Float64(J)
                QTDA_export[2,qp_count] = Float64(P)
                QTDA_export[3,qp_count] = E_QTDA[P,J+1][nu]
                QTDA_export[4,qp_count] = E_QTDA[P,J+1][nu]
            end
        end
    end

    E_QTDA_export = @views QTDA_export[3,:]

    Sort = sortperm(E_QTDA_export, by = x -> real(x))

    QTDA_export .= @views QTDA_export[:,Sort]

    # Adjust the positions of J & P labels ...
    @inbounds for nu in 2:(N_qp+1)
        E_1 = QTDA_export[4,nu-1]
        E_2 = QTDA_export[4,nu]
        if E_2 < E_1
            QTDA_export[4,nu] = E_1 + 1.0
        elseif abs(E_2 - E_1) < Delta
            QTDA_export[4,nu] = E_1 + 1.0
        end
    end

    # QTDA spectrum plot data export ...
    println("\nPerforming the export of plot-ready QTDA spectra ...")
    open(Output_Path, "w") do Write_File
        println(Write_File, "J\tP\tE\tm")
        @inbounds for nu in 1:N_qp
            J = Int64(round(QTDA_export[1,nu]))
            P = "P"
            if abs(QTDA_export[2,nu] - 1.0) < 1e-3
                P = "+"
            else
                P = "-"
            end
            E = QTDA_export[3,nu]
            E_m = QTDA_export[4,nu]
            Row = string(string(J) * "\t" * P * "\t" * string(round(E,sigdigits = 5)) * "\t" * string(round(E_m,sigdigits = 5)))
            println(Write_File, Row)
        end
    end

    println("\nQTDA plot-ready spectra successfully exported ...")

    return
end

function QTDA_phonon_space_export(Params::Parameters,Orb_2qp::qpOrb2B)
    # Read parameters ...
    N_max = Params.Calc.Nmax
    N_2max = 2*N_max
    J_max = N_2max + 1
    Output_File = Params.Calc.Path

    # Output file name & path ...
    Output_Path = "IO/" * Output_File * "/QTDA/QTDA_Phonon_Subspace_Size.dat"

    # Export dimensiones of 1-phonon quasiparticle subspaces ...
    println("\nExporting dimensions of QTDA 1-phonon quasiparticle subspaces ...")
    open(Output_Path, "w") do Write_File
        println(Write_File, "J\tN_phonon")
        @inbounds for J in 0:J_max
            N_qp_Jp, N_qp_Jm = Orb_2qp.N[1,J+1], Orb_2qp.N[2,J+1]
            Row = string(string(J) * "\t" * string(N_qp_Jp + N_qp_Jm))
            println(Write_File, Row)
        end
        println(Write_File, "\nSubspace\tP = +")
        @inbounds for J in 0:J_max
            N_qp = Orb_2qp.N[1,J+1]
            Row = string(string(J) * "\t" *string(N_qp))
            println(Write_File, Row)
        end
        println(Write_File, "\nSubspace\tP = -")
        @inbounds for J in 0:J_max
            N_qp = Orb_2qp.N[2,J+1]
            Row = string(string(J) * "\t" *string(N_qp))
            println(Write_File, Row)
        end
    end

    println("\nDimensions of QTDA 1-phonon subspaces successfully exported ...")

    return
end




# WIP

function QTDA_transitions_export(Params::Parameters,Orb_2qp::qpOrb2B,E_QTDA::Matrix{Vector{Float64}},rB_QTDA::ReducedTransition)
    # Read parameters ...
    Orthogon = Params.Calc.QTDA.Ortho

    # Set the export path ...
    Output_File = "IO/" * Params.Calc.Path

    println("\nPreparing export of QTDA 1-phonon electromagnetic transition intensities ...")

    # QTDA E0 export ...
    open(Output_File * "/QTDA/Transitions/E0/QTDA_E0.dat", "w") do Write_File
        J, P = 0, 1
        println(Write_File, "E\tB_ph\tB_is\tB_iv")
        @inbounds for nu = 1:Orb_2qp.N[P,J+1]
            println(Write_File, string(round(E_QTDA[P,J+1][nu], sigdigits=4)) * "\t" * string(round(rB_QTDA.E0.ph[nu],digits = 4)) * "\t" * string(round(rB_QTDA.E0.is[nu],digits = 4)) * "\t" * string(round(rB_QTDA.E0.iv[nu],digits = 4)))
        end
    end

    if Orthogon == true && 1 == 2
        # QTDA E1 export ...
        open(Output_File * "/QTDA/Transitions/E1/QTDA_E1.dat", "w") do Write_File
            J, P = 1, 2
            println(Write_File, "E\tB_ph\tB_is\tB_iv")
            @inbounds for nu in 1:Orb_2qp.N[P,J+1]
                println(Write_File, string(round(E_QTDA[P,J+1][nu], sigdigits=4)) * "\t" * string(round(rB_QTDA.E1.ph[nu],digits = 4)) * "\t" * string(round(rB_QTDA.E1.is[nu],digits = 4)) * "\t" * string(round(rB_QTDA.E1.iv[nu],digits = 4)))
            end
        end
    else
        # QTDA E1 export ...
        open(Output_File * "/QTDA/Transitions/E1/QTDA_E1.dat", "w") do Write_File
            J, P = 1, 2
            println(Write_File, "E\tB_E1_ph\tB_E1_is\tB_E1_iv")
            @inbounds for nu in 1:Orb_2qp.N[P,J+1]
                println(Write_File, string(round(E_QTDA[P,J+1][nu], sigdigits=4)) * "\t" * string(round(rB_QTDA.E1.ph[nu],digits = 4)) * "\t" * string(round(rB_QTDA.E1.is[nu],digits = 4)) * "\t" * string(round(rB_QTDA.E1.iv[nu],digits = 4)))
            end
        end
    end

    # QTDA E2 export ...
    open(Output_File * "/QTDA/Transitions/E2/QTDA_E2.dat", "w") do Write_File
        J, P = 2, 1
        println(Write_File, "E\tB_ph\tB_is\tB_iv")
        @inbounds for nu = 1:Orb_2qp.N[P,J+1]
            println(Write_File, string(round(E_QTDA[P,J+1][nu], sigdigits=4)) * "\t" * string(round(rB_QTDA.E2.ph[nu],digits = 4)) * "\t" * string(round(rB_QTDA.E2.is[nu],digits = 4)) * "\t" * string(round(rB_QTDA.E2.iv[nu],digits = 4)))
        end
    end

    # QTDA E3 export ...
    open(Output_File * "/QTDA/Transitions/E3/QTDA_E3.dat", "w") do Write_File
        J, P = 3, 2
        println(Write_File, "E\tB_ph\tB_is\tB_iv")
        @inbounds for nu = 1:1:Orb_2qp.N[P,J+1]
            println(Write_File, string(round(E_QTDA[P,J+1][nu], sigdigits=4)) * "\t" * string(round(rB_QTDA.E3.ph[nu],digits = 4)) * "\t" * string(round(rB_QTDA.E3.is[nu],digits = 4)) * "\t" * string(round(rB_QTDA.E3.iv[nu],digits = 4)))
        end
    end

    println("\nQTDA 1-phonon electromagnetic transition intensities succesfully exported ...")

    return
end
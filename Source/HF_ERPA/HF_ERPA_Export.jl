function HF_ERPA_Export(Params::Parameters,Orb::Vector{NOrb},N_nu::Matrix{Int64},E_0_RPA::Float64,E_RPA::Matrix{Vector{ComplexF64}},X_RPA::Matrix{Matrix{ComplexF64}},Y_RPA::Matrix{Matrix{ComplexF64}},CI_RPA::Matrix{Vector{Float64}},rB_RPA::ReducedTransition,Rho::pnMatrix,U::pnMatrix)
    
    # Export ERPA summary file ...
    HF_ERPA_Summary(Params,N_nu,E_0_RPA)

    # Export of ERPA spectra ...
    HF_ERPA_Spectrum_Export(Params,N_nu,E_RPA,CI_RPA)

    # Export of ERPA plot-ready spectra ...
    HF_ERPA_Plot_Spectrum_Export(Params,N_nu,E_RPA)

    # Export single-particle level occupations ...
    HF_ERPA_Occupation_Export(Params,Orb,Rho)

    # Export ERPA radial densities & radii ...
    HF_ERPA_Radial_Density(Params,Orb,Rho,U)

    # Export of ERPA |X|^2 & |Y|^2 amplitudes, not plot ready ...
    HF_ERPA_Amplitudes_Export(Params,N_nu,E_RPA,X_RPA,Y_RPA)

    # Export of ERPA electric transitions ...
    HF_ERPA_Transition_Export(Params,N_nu,E_RPA,rB_RPA)

    # Export dimensions of phonon subspaces ...
    HF_ERPA_Phonon_Space_Export(Params,N_nu)

    return
end

function HF_ERPA_Summary(Params::Parameters,N_nu::Matrix{Int64},E_0::Float64)
    # Read calculation parameters ...
    A = Params.Calc.A
    Z = Params.Calc.Z
    HbarOmega = Params.Calc.hw
    N_max = Params.Calc.Nmax
    N_2max = 2*N_max
    J_max = N_2max + 1
    Orthogon = Params.Calc.Ortho
    Output_File = Params.Calc.Path

    # Export name & path ...
    if Orthogon == true
        Output_Path = "IO/" * Output_File * "/ERPA/HF_ERPA_Summary_Ortho.dat"
    else
        Output_Path = "IO/" * Output_File * "/ERPA/HF_ERPA_Summary_Spur.dat"
    end
    
    println("\nPreparing ERPA summary ...")
    Summary_File =  open(Output_Path, "w")
        println(Summary_File, "Spherical Hartree-Fock Extended-Random-Phase Approximation review:")
        println(Summary_File, "\nNuclid data:    A = " * string(A) * ", Z = " * string(Z))
        println(Summary_File, "\nCalculation data:    HbarOmega = " * string(HbarOmega) * " , N_max = " * string(N_max) *
                ", s.p. J-basis size = " * string(div((N_max+1)*(N_max+2),2)) * ", s.p. M-basis size = " * string(div((N_max+1)*(N_max+2)*(N_max+3),6)))
        if Orthogon == true
            println(Summary_File, "\nOrthogonalization of 1- spurious state included ...")
        else
            println(Summary_File, "\nOrthogonalization of 1- spurious state not included ...")
        end
        N_M_tot = 0
        @inbounds for J in 0:J_max
            @inbounds for P in 1:2
                n = N_nu[J+1,P] * (2*J + 1)
                N_M_tot += n
            end
        end
        println(Summary_File, "\nTotal number of 1 phonon states in J-scheme = " * string(sum(N_nu)))
        println(Summary_File, "Total number of 1 phonon states in M-scheme = " * string(N_M_tot))
        println(Summary_File, "\nE_Corr_ERPA = " * string(round(E_0, digits = 8)) * "\t MeV \t...\tERPA ground-state correlation energy")
        println(Summary_File, "\n\tTo get the total correct ground-state energy add the HF-mean field energy ...")
    close(Summary_File)

    println("\nERPA summary exported ...")

    return

end

function HF_ERPA_Spectrum_Export(Params::Parameters,N_nu::Matrix{Int64},E_RPA::Matrix{Vector{ComplexF64}},CI_RPA::Matrix{Vector{Float64}})
    # Read calculation parameters ...
    N_max = Params.Calc.Nmax
    N_2max = 2*N_max
    J_max = N_2max + 1
    Orthogon = Params.Calc.Ortho
    Output_File = Params.Calc.Path

    # Export name & path ...
    if Orthogon == true
        Output_Path_RPA = "IO/" * Output_File * "/ERPA/Spectra/ERPA_Ortho.dat"
    else
        Output_Path_RPA = "IO/" * Output_File * "/ERPA/Spectra/ERPA_Spur.dat"
    end

    println("\nPreparing ERPA spectra export ...")

    # ERPA spectrum export
    open(Output_Path_RPA, "w") do Write_File
        println(Write_File, "CI\tJ\tP\tE")
        @inbounds for J in 0:J_max
            @inbounds for P in 1:2
                N_ph = N_nu[J+1,P]
                @inbounds for nu in 1:N_ph
                    if P == 1
                        println(Write_File, string(round(CI_RPA[J+1,P][nu],digits=3)) * "\t" * string(J) * "\t" * "+" * "\t" * string(E_RPA[J+1,P][nu]))
                    else
                        println(Write_File, string(round(CI_RPA[J+1,P][nu],digits=3)) * "\t" * string(J) * "\t" * "-" * "\t" * string(E_RPA[J+1,P][nu]))
                    end
                end
                println(Write_File,"\n")
            end
        end
    end

    println("\nERPA spectra export finished ...")

    return
end

function HF_ERPA_Plot_Spectrum_Export(Params::Parameters,N_nu::Matrix{Int64},E_RPA::Matrix{Vector{ComplexF64}})
    # Read calculation parameters ...
    N_max = Params.Calc.Nmax
    N_2max = 2*N_max
    J_max = N_2max + 1
    Orthogon = Params.Calc.Ortho
    Output_File = Params.Calc.Path

    # Export name & path ...
    if Orthogon == true
        Output_Path_RPA = "IO/" * Output_File * "/ERPA/Spectra/ERPA_Plot_Ortho.dat"
    else
        Output_Path_RPA = "IO/" * Output_File * "/ERPA/Spectra/ERPA_Plot_Spur.dat"
    end

    # Spectra label gap parameter ...
    Delta = 1.0

    println("\nPreparing plot-ready export of RPA solutions ...")

    N_ph = Int64(sum(N_nu))

    RPA_Solution = Matrix{Float64}(undef,4,N_ph+1)

    # 0+ ground state addition ...

    RPA_Solution[1,1] = 0.0
    RPA_Solution[2,1] = 1.0
    RPA_Solution[3,1] = 0.0
    RPA_Solution[4,1] = 0.0

    nu_count = 1
    @inbounds for J in 0:J_max
        @inbounds for P in 1:2
            N_qg = N_nu[J+1,P]
            @inbounds for nu in 1:N_qg
                nu_count += 1
                RPA_Solution[1,nu_count] = Float64(J)
                RPA_Solution[2,nu_count] = Float64(P)
                RPA_Solution[3,nu_count] = real(E_RPA[J+1,P][nu])
                RPA_Solution[4,nu_count] = real(E_RPA[J+1,P][nu])
            end
        end
    end

    E_RPA_full = @views RPA_Solution[3,:]

    Sort_RPA = sortperm(E_RPA_full, by = x -> real(x))

    RPA_Solution .= @views RPA_Solution[:,Sort_RPA]

    # Adjust the JP label positions ...
    @inbounds for nu in 2:(N_ph+1)

        e_RPA_1 = RPA_Solution[4,nu-1]
        e_RPA_2 = RPA_Solution[4,nu]
        if e_RPA_2 < e_RPA_1
            RPA_Solution[4,nu] = e_RPA_1 + 1.0
        elseif abs(e_RPA_2 - e_RPA_1) < Delta
            RPA_Solution[4,nu] = e_RPA_1 + 1.0
        end
    end

    # RPA spectrum plot data export ...
    open(Output_Path_RPA, "w") do Write_File
        println(Write_File, "J\tP\tE\tm")
        @inbounds for nu in 1:N_ph
            J = Int64(round(RPA_Solution[1,nu]))
            P = "P"
            if abs(RPA_Solution[2,nu] - 1.0) < 1e-3
                P = "+"
            else
                P = "-"
            end
            E = RPA_Solution[3,nu]
            E_m = RPA_Solution[4,nu]
            Row = string(string(J) * "\t" * P * "\t" * string(round(E,sigdigits = 5)) * "\t" * string(round(E_m,sigdigits = 5)))
            println(Write_File, Row)
        end
    end

    println("\nRPA & TDA plot-ready solutions exported ...")

    return
end

function HF_ERPA_Occupation_Export(Params::Parameters,Orb::Vector{NOrb},Rho::pnMatrix)
    # Read parameters ...
    N_max = Params.Calc.Nmax
    Orthogon = Params.Calc.Ortho
    Output_File = Params.Calc.Path

    a_max = div((N_max+1)*(N_max+2),2)

    pRho, nRho = Rho.p, Rho.n

    if Orthogon == true
        Output_Path = "IO/" * Output_File * "/ERPA/Densities/ERPA_Occupation_Ortho.dat"
    else
        Output_Path = "IO/" * Output_File * "/ERPA/Densities/ERPA_Occupation_Spur.dat"
    end
    
    # Make HF density matrices ...
    pRho_HF, nRho_HF = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)
    @inbounds for a in 1:a_max
        pRho_HF[a,a] = Orb[a].pO
        nRho_HF[a,a] = Orb[a].nO
    end

    pO_HF = Vector{Float64}(undef,a_max)
    nO_HF = Vector{Float64}(undef,a_max)
    pO_ERPA = Vector{Float64}(undef,a_max)
    nO_ERPA = Vector{Float64}(undef,a_max)
    pa = Vector{Int64}(undef,a_max)
    na = Vector{Int64}(undef,a_max)

    @inbounds for a in 1:a_max
        pO_HF[a] = pRho_HF[a,a]
        nO_HF[a] = nRho_HF[a,a]
        pO_ERPA[a] = round(pRho[a,a], digits = 4)
        nO_ERPA[a] = round(nRho[a,a], digits = 4)
        pa[a] = a
        na[a] = a
    end

    Sort_pInd = sortperm(pO_ERPA, rev = true)
    Sort_nInd = sortperm(nO_ERPA, rev = true)

    pO_HF = pO_HF[Sort_pInd]
    pO_ERPA = pO_ERPA[Sort_pInd]
    pa = pa[Sort_pInd]

    nO_HF = nO_HF[Sort_nInd]
    nO_ERPA = nO_ERPA[Sort_nInd]
    na = na[Sort_nInd]

    # ERPA spectrum export
    open(Output_Path, "w") do Write_File
        println(Write_File, "a\ta_p\tpN_HF\tpN_ERPA\tpN_dep\ta_n\tnN_HF\tnN_ERPA\tnN_dep")
        @inbounds for a in 1:a_max
            println(Write_File, string(a) * "\t" * string(pa[a]) * "\t" * string(pO_HF[a]) * "\t" * string(pO_ERPA[a]) * "\t" * string(round(pO_ERPA[a] - pO_HF[a], digits = 5)) * "\t" * string(na[a]) * "\t" * string(nO_HF[a]) * "\t" * string(nO_ERPA[a]) * "\t" * string(round(nO_ERPA[a] - nO_HF[a], digits = 5)))
        end
    end

end

function HF_ERPA_Amplitudes_Export(Params::Parameters,N_nu::Matrix{Int64},E_RPA::Matrix{Vector{ComplexF64}},X_RPA::Matrix{Matrix{ComplexF64}},Y_RPA::Matrix{Matrix{ComplexF64}})
    # Read parameters ...
    N_max = Params.Calc.Nmax
    N_2max = 2*N_max
    J_max = N_2max + 1
    Orthogon = Params.Calc.Ortho
    Output_File = Params.Calc.Path

    # First a systematic export according to JP sectors ...
    #______________________________________________________

    # Export name & path ...
    if Orthogon == true
        Output_Path = "IO/" * Output_File * "/ERPA/Amplitudes/HF_ERPA_Amplitudes_Ordered_Ortho.dat"
    else
        Output_Path = "IO/" * Output_File * "/ERPA/Amplitudes/HF_ERPA_Amplitudes_Ordered_Spur.dat"
    end

    println("\nPreparing export of ERPA amplitudes X & Y ...")

    open(Output_Path, "w") do Write_File
        println(Write_File, "E\t|X|^2\t\t|Y|^2")
        @inbounds for J in 0:J_max
            @inbounds for P in 1:2
                if P == 1
                    println(Write_File, "\nJ = " * string(J) * ",\tP = +")
                else
                    println(Write_File, "\nJ = " * string(J) * ",\tP = -")
                end
                N_ph = N_nu[J+1,P]
                if Orthogon == true
                    if (J == 1 && P == 2)
                        @inbounds for nu in 2:N_ph
                            x = @views norm(X_RPA[J+1,P][:,nu])^2
                            y = @views norm(Y_RPA[J+1,P][:,nu])^2
                            println(Write_File, string(round(real(E_RPA[J+1,P][nu]),digits=4)) * "\t" * string(round(x,digits=7)) * "\t" * string(round(y,digits=7)))
                        end
                    else
                        @inbounds for nu in 1:N_ph
                            x = @views norm(X_RPA[J+1,P][:,nu])^2
                            y = @views norm(Y_RPA[J+1,P][:,nu])^2
                            println(Write_File, string(round(real(E_RPA[J+1,P][nu]),digits=4)) * "\t" * string(round(x,digits=7)) * "\t" * string(round(y,digits=7)))
                        end
                    end
                else
                    @inbounds for nu in 1:N_ph
                        x = @views norm(X_RPA[J+1,P][:,nu])^2
                        y = @views norm(Y_RPA[J+1,P][:,nu])^2
                        println(Write_File, string(round(real(E_RPA[J+1,P][nu]),digits=4)) * "\t" * string(round(x,digits=7)) * "\t" * string(round(y,digits=7)))
                    end
                end
            end
        end
    end

    # Second an energy ordered export suitable for analysis ...
    #___________________________________________________________

    println("\nExported norms of ERPA amplitudes |X|^2 & |Y|^2 ...")

    # Allocate new export path ...
    Output_File = Params.Calc.Path

    if Orthogon == true
        Output_Path = "IO/" * Output_File * "/ERPA/Amplitudes/HF_ERPA_Amplitudes_Ortho.dat"
    else
        Output_Path = "IO/" * Output_File * "/ERPA/Amplitudes/HF_ERPA_Amplitudes_Spur.dat"
    end

    println("\nPreparing export of energy-ordered RPA |X|^2 & |Y|^2 amplitudes ...")

    N_ph = sum(N_nu)
    E_RPA_ord = Vector{Float64}(undef,N_ph)
    y_RPA = Vector{Float64}(undef,N_ph)
    x_RPA = Vector{Float64}(undef,N_ph)
    J_RPA = Vector{Int64}(undef,N_ph)
    P_RPA = Vector{Int64}(undef,N_ph)

    ph_count = 0
    @inbounds for J in 0:J_max
        @inbounds for P in 1:2
            N_ph = N_nu[J+1,P]
            @inbounds for nu in 1:N_ph
                ph_count += 1
                E_RPA_ord[ph_count] = real(E_RPA[J+1,P][nu])
                xSum = 0.0
                ySum = 0.0
                @inbounds for ph in 1:N_ph
                    x = abs2(X_RPA[J+1,P][ph,nu])
                    y = abs2(Y_RPA[J+1,P][ph,nu])
                    ySum += y
                    xSum += x
                end
                x_RPA[ph_count] = xSum
                y_RPA[ph_count] = ySum
                J_RPA[ph_count] = J + 1 
                P_RPA[ph_count] = P
            end
        end
    end

    Sort_Ind = sortperm(E_RPA_ord)
    E_RPA_ord = E_RPA_ord[Sort_Ind]
    x_RPA = x_RPA[Sort_Ind]
    y_RPA = y_RPA[Sort_Ind]
    J_RPA = J_RPA[Sort_Ind]
    P_RPA = P_RPA[Sort_Ind]

    open(Output_Path, "w") do Write_File
        N_ph = sum(N_nu)
        println(Write_File, "Ind\tE\tJ\tP\t|Y|^2\t|X|^2")
        @inbounds for ph in 1:N_ph
            if P_RPA[ph] == 1
                println(Write_File, string(ph) * "\t" * string(round(E_RPA_ord[ph],sigdigits=5)) * "\t" * string(J_RPA[ph] - 1) * "\t+" * "\t" * string(round(y_RPA[ph],digits=5)) * "\t" * string(round(x_RPA[ph],digits=5)))
            elseif P_RPA[ph] == 2
                println(Write_File, string(ph) * "\t" * string(round(E_RPA_ord[ph],sigdigits=5)) * "\t" * string(J_RPA[ph] - 1) * "\t-" * "\t" * string(round(y_RPA[ph],digits=5)) * "\t" * string(round(x_RPA[ph],digits=5)))
            end
        end
    end

    return
end

function HF_ERPA_Transition_Export(Params::Parameters,N_nu::Matrix{Int64},E_RPA::Matrix{Vector{ComplexF64}},rB_RPA::ReducedTransition)
    # Read parameters ...
    Orthogon = Params.Calc.Ortho
    Output_File = "IO/" * Params.Calc.Path

    println("\nPreparing export of ERPA transition intensities ...")


    # ERPA E0 export
    open(Output_File * "/ERPA/Transitions/E0/ERPA_E0.dat", "w") do Write_File
        J, P = 0, 1
        println(Write_File, "E\tB_ph\tB_is\tB_iv")
        @inbounds for nu = 1:N_nu[J+1,P]
            println(Write_File, string(round(real(E_RPA[J+1,P][nu]), sigdigits=4)) * "\t" * string(round(rB_RPA.E0.ph[nu],sigdigits = 4)) * "\t" * string(round(rB_RPA.E0.is[nu],sigdigits = 4)) * "\t" * string(round(rB_RPA.E0.iv[nu],sigdigits = 4)))
        end
    end


    if Orthogon == true
        # ERPA E1 export
        open(Output_File * "/ERPA/Transitions/E1/ERPA_E1_Ortho.dat", "w") do Write_File
            J, P = 1, 2
            println(Write_File, "E\tB_ph\tB_is\tB_iv")
            @inbounds for nu in 1:N_nu[J+1,P]
                println(Write_File, string(round(real(E_RPA[J+1,P][nu]), sigdigits=4)) * "\t" * string(round(rB_RPA.E1.ph[nu],digits = 4)) * "\t" * string(round(rB_RPA.E1.is[nu],digits = 4)) * "\t" * string(round(rB_RPA.E1.iv[nu],digits = 4)))
            end
        end
    else
        # ERPA E1 export
        open(Output_File * "/ERPA/Transitions/E1/ERPA_E1_Spur.dat", "w") do Write_File
            J, P = 1, 2
            println(Write_File, "E\tB_E1_ph\tB_E1_is\tB_E1_iv")
            @inbounds for nu in 1:N_nu[J+1,P]
                println(Write_File, string(round(real(E_RPA[J+1,P][nu]), sigdigits=4)) * "\t" * string(round(rB_RPA.E1.ph[nu],digits = 4)) * "\t" * string(round(rB_RPA.E1.is[nu],digits = 4)) * "\t" * string(round(rB_RPA.E1.iv[nu],digits = 4)))
            end
        end
    end

    # ERPA E2 export
    open(Output_File * "/ERPA/Transitions/E2/ERPA_E2.dat", "w") do Write_File
        J, P = 2, 1
        println(Write_File, "E\tB_ph\tB_is\tB_iv")
        @inbounds for nu = 1:N_nu[J+1,P]
            println(Write_File, string(round(real(E_RPA[J+1,P][nu]), sigdigits=4)) * "\t" * string(round(rB_RPA.E2.ph[nu],digits = 4)) * "\t" * string(round(rB_RPA.E2.is[nu],digits = 4)) * "\t" * string(round(rB_RPA.E2.iv[nu],digits = 4)))
        end
    end

    # ERPA E3 export
    open(Output_File * "/ERPA/Transitions/E3/ERPA_E3.dat", "w") do Write_File
        J, P = 3, 2
        println(Write_File, "E\tB_ph\tB_is\tB_iv")
        @inbounds for nu = 1:N_nu[J+1,P]
            println(Write_File, string(round(real(E_RPA[J+1,P][nu]), sigdigits=4)) * "\t" * string(round(rB_RPA.E3.ph[nu],digits = 4)) * "\t" * string(round(rB_RPA.E3.is[nu],digits = 4)) * "\t" * string(round(rB_RPA.E3.iv[nu],digits = 4)))
        end
    end

    println("\nERPA transition intensities succesfully exported ...")

    return
end

function HF_ERPA_Phonon_Space_Export(Params::Parameters,N_nu::Matrix{Int64})
    # Read parameters ...
    N_max = Params.Calc.Nmax
    N_2max = 2*N_max
    J_max = N_2max + 1
    Output_File = Params.Calc.Path

    println("\nPreparing export of dimensions of ERPA phonon subspaces ...")

    # Output file name & path ...
    Output_Path = "IO/" * Output_File * "/ERPA/HF_ERPA_Phonon_Subspace_Size.dat"

    # Export of 1-phonon subspace dimensions ...
    open(Output_Path, "w") do Write_File
        println(Write_File, "J\tN_phonon")
        @inbounds for J in 0:J_max
            Row = string(string(J) * "\t" * string(N_nu[J+1,1] + N_nu[J+1,2]))
            println(Write_File, Row)
        end
        println(Write_File, "\nSubspace\tP = +")
        @inbounds for J in 0:J_max
            Row = string(string(J) * "\t" *string(N_nu[J+1,1]))
            println(Write_File, Row)
        end
        println(Write_File, "\nSubspace\tP = -")
        @inbounds for J in 0:J_max
            Row = string(string(J) * "\t" *string(N_nu[J+1,2]))
            println(Write_File, Row)
        end
    end

    println("\nDimensions of ERPA phonon subspaces exported ...")

    return
end
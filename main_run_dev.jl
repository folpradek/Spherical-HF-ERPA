using DelimitedFiles, LinearAlgebra, CGcoefficient, BenchmarkTools, Hungarian
include("Source/Import_Structures.jl")
include("Source/Functions/Functions.jl")
include("Source/MatrixElements/Orb.jl")
include("Source/MatrixElements/TrOp.jl")
include("Source/MatrixElements/TrOp.jl")
include("Source/MatrixElements/T1B.jl")
include("Source/MatrixElements/T2B.jl")
include("Source/MatrixElements/V2B.jl")
include("Source/MatrixElements/V3B_NO2B.jl")
include("Source/MatrixElements/V2B_Res.jl")
include("Source/HF/HF.jl")
include("Source/HF/HF_Solver.jl")
include("Source/HF/HF_Density_Operator.jl")
include("Source/HF/HF_Orbital_Ordering.jl")
include("Source/HF/HF_Radial_Density.jl")
include("Source/HF/HF_Radial_Potential.jl")
include("Source/HF/HF_Energy.jl")
include("Source/HF/HF_Summary.jl")
include("Source/HF/HF_Export.jl")
include("Source/HF/HF_MBPT.jl")
include("Source/HF/HF_Radial_MBPT.jl")
include("Source/HF_RPA/HF_RPA_Solver.jl")
include("Source/HF_RPA/HF_RPA.jl")
include("Source/HF_RPA/HF_RPA_Phonon_Count.jl") 
include("Source/HF_ERPA/HF_ERPA_Hamiltonian.jl")
include("Source/HF_RPA/HF_RPA_Allocate.jl")
include("Source/HF_RPA/HF_RPA_Spurious.jl")
include("Source/HF_RPA/HF_RPA_Diagonalize.jl")
include("Source/HF_RPA/HF_RPA_Corr_Energy.jl")
include("Source/HF_RPA/HF_RPA_OBDM.jl")
include("Source/HF_RPA/HF_RPA_Collectivity.jl")
include("Source/HF_RPA/HF_RPA_Radial_Density.jl")
include("Source/HF_RPA/HF_RPA_Transitions.jl")
include("Source/HF_RPA/HF_RPA_Export.jl")
include("Source/HF_ERPA/HF_ERPA_Solver.jl")
include("Source/HF_ERPA/HF_ERPA.jl")
include("Source/HF_ERPA/HF_ERPA_Phonon_Count.jl")
include("Source/HF_ERPA/HF_ERPA_Iteration.jl")
include("Source/HF_ERPA/HF_ERPA_OBDM.jl")
include("Source/HF_ERPA/HF_ERPA_Allocate.jl")
include("Source/HF_ERPA/HF_ERPA_Diagonalize.jl")
include("Source/HF_ERPA/HF_ERPA_Spurious.jl")
include("Source/HF_ERPA/HF_ERPA_Transitions.jl")
include("Source/HF_ERPA/HF_ERPA_Energy.jl")
include("Source/HF_ERPA/HF_ERPA_Radial_Density.jl")
include("Source/HF_ERPA/HF_ERPA_Collectivity.jl")
include("Source/HF_ERPA/HF_ERPA_Export.jl")


# pHF RPA ...
include("Source/pHF_RPA/pHF_RPA_Solver.jl")
include("Source/pHF_RPA/pHF_RPA.jl")
include("Source/pHF_RPA/pHF_RPA_Allocate.jl")
include("Source/pHF_RPA/pHF_RPA_Diagonalize.jl")

# BCS
include("Source/BCS/BCS.jl")
include("Source/BCS/BCS_Solver.jl")
include("Source/BCS/BCS_Export.jl")

# HFB
include("Source/HFB/HFB.jl")
include("Source/HFB/HFB_Solver.jl")
include("Source/HFB/HFB_Export.jl")

function run_main()

    # HF, HF-RPA, HF_ERPA runs ...
    IntParams = Interaction_Parameters(
                NN_File = "IO/NN.bin",
                NNN_File = "IO/NNN.bin",
                hw = 16.0,
                Nmax = 3,
                N2max = 6,
                N3max = 9
                )

    CalcParams = Calculation_Parameters(
                A = 16,
                Z = 8,
                hw = 16.0,
                Nmax = 3,
                N2max = 6,
                N3max = 9,
                Ortho = true,
                CMS = "CMS1+2B",
                BMF = true,
                Path = "pHF_RPA_A16_Z8_hw16.0_Nmax3_N2max6_N3max9_CMS1+2B",
                ERPA = Parameters_ERPA(OBDM = "Full", ScOBH = true, Sc3N = true),
                Format = "Bin",
                #cV_res = 0.01
                )

    #HF_Solver(Parameters(IntParams,CalcParams))

    #HF_RPA_Solver(Parameters(IntParams,CalcParams))

    #HF_ERPA_Solver(Parameters(IntParams,CalcParams))

    # BCS testing ...
    IntParams = Interaction_Parameters(
                NN_File = "IO/NN.bin",
                NNN_File = "IO/NNN.bin",
                hw = 16.0,
                Nmax = 3,
                N2max = 6,
                N3max = 9
                )


    CalcParams = Calculation_Parameters(
                A = 18,
                Z = 8,
                hw = 16.0,
                Nmax = 3,
                N2max = 6,
                N3max = 9,
                Ortho = true
                )

    BCS(Parameters(IntParams,CalcParams))
    #HFB(Parameters(IntParams,CalcParams))


end

@time run_main()
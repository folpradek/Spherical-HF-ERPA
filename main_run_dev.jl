using DelimitedFiles, LinearAlgebra, CGcoefficient, BenchmarkTools, Hungarian, Mmap
include("Source/Import_Structures.jl")
include("Source/Functions/Functions.jl")
include("Source/MatrixElements/Orb.jl")
include("Source/MatrixElements/OBDM.jl")
include("Source/MatrixElements/TrOp.jl")
include("Source/MatrixElements/T1B.jl")
include("Source/MatrixElements/T2B.jl")
include("Source/MatrixElements/V2B.jl")
include("Source/MatrixElements/V3B_NO2B.jl")
include("Source/MatrixElements/V2B_Res.jl")
include("Source/HF/HF.jl")
include("Source/HF/HF_Solver.jl")
include("Source/HF/HF_Allocate.jl")
include("Source/HF/HF_Density_Operator.jl")
include("Source/HF/HF_Ordering.jl")
include("Source/HF/HF_Radial_Potential.jl")
include("Source/HF/HF_Energy.jl")
include("Source/HF/HF_Export.jl")
include("Source/HF/HF_MBPT.jl")
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
include("Source/BCS/BCS_Residual_Interaction.jl")
include("Source/BCS/BCS_Density_Operator.jl")
include("Source/BCS/BCS_Chemical_Potential.jl")
include("Source/BCS/BCS_Particle_Number.jl")
include("Source/BCS/BCS_Allocate.jl")
include("Source/BCS/BCS_Energy.jl")
include("Source/BCS/BCS_Export.jl")

# HFB
include("Source/HFB/HFB.jl")
include("Source/HFB/HFB_Solver.jl")
include("Source/HFB/HFB_Density_Operator.jl")
include("Source/HFB/HFB_Chemical_Potential.jl")
include("Source/HFB/HFB_Particle_Number.jl")
include("Source/HFB/HFB_Broyden.jl")
include("Source/HFB/HFB_Ordering.jl")
include("Source/HFB/HFB_Allocate.jl")
include("Source/HFB/HFB_Canonical_Basis.jl")
include("Source/HFB/HFB_Energy.jl")
include("Source/HFB/HFB_Export.jl")

function run_main()

    # Interaction parameters structure ...
    IntParams = Interaction_Parameters(
                NN_File = "IO/NN.bin",          #   Path to the NN interaction binary file
                NNN_File = "IO/NNN.bin",        #   Path to the NNN interaction binary file
                hw = 16.0,
                Nmax = 3,
                N2max = 6,
                N3max = 9
                )

    # Calculation parameters structure ...
    CalcParams = Calculation_Parameters(
                A = 18,
                Z = 8,
                hw = 16.0,
                Nmax = 3,
                N2max = 6,
                N3max = 9,
                Ortho = true,
                CMS = "CMS1+2B",
                BMF = true,
                Path = "A18_Z8_hw16.0_Nmax3_N2max6_N3max9",
                #Format = "Bin",
                #ERPA = Parameters_ERPA(OBDM = "Full", ScOBH = true, Sc3N = true),
                #cV_res = 0.01
                #Pairing = Parameters_Pairing(ScBCS = false, pL0 = -5.0, nL0 = -5.0, pK0 = 0.5, nK0 = 0.5, pdN0 = 2.0, ndN0 = 2.0),
                )

    # HF calculation call ...
    #HF(Parameters(IntParams,CalcParams))

    # HF-RPA calculation call ...
    #HF_RPA_Solver(Parameters(IntParams,CalcParams))

    # HF-ERPA calculation call ...
    #HF_ERPA_Solver(Parameters(IntParams,CalcParams))

    # HF-BCS calculation call ...
    #BCS(Parameters(IntParams,CalcParams))

    # HFB calculation call ...
    HFB(Parameters(IntParams,CalcParams))

end

@time run_main()

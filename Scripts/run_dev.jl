using DelimitedFiles, LinearAlgebra, CGcoefficient, BenchmarkTools, Hungarian, Mmap, Printf
include("../Source/Import_Structures.jl")
include("../Source/Functions/Functions.jl")
include("../Source/MatrixElements/Orb.jl")
include("../Source/MatrixElements/O1B.jl")
include("../Source/MatrixElements/T1B.jl")
include("../Source/MatrixElements/Tr1B.jl")
include("../Source/MatrixElements/OBDM.jl")
include("../Source/MatrixElements/T2B.jl")
include("../Source/MatrixElements/O2B.jl")
include("../Source/MatrixElements/V2B.jl")
include("../Source/MatrixElements/V3B_NO2B.jl")
include("../Source/MatrixElements/qpO2B.jl")
include("../Source/MatrixElements/qpO2B_40.jl")
include("../Source/MatrixElements/qpO2B_31.jl")
include("../Source/MatrixElements/qpO2B_22.jl")

 # HF + MBPT ...
include("../Source/HF/HF.jl")
include("../Source/HF/HF_Solver.jl")
include("../Source/HF/HF_Allocate.jl")
include("../Source/HF/HF_Density_Operator.jl")
include("../Source/HF/HF_Ordering.jl")
include("../Source/HF/HF_Radial_Potential.jl")
include("../Source/HF/HF_Energy.jl")
include("../Source/HF/HF_Export.jl")
include("../Source/HF/HF_MBPT.jl")
# HF-RPA ...
include("../Source/HF_RPA/HF_RPA.jl")
include("../Source/HF_RPA/HF_RPA_Solver.jl")
include("../Source/HF_RPA/HF_RPA_Phonon_Count.jl") 
include("../Source/HF_RPA/HF_RPA_Allocate.jl")
include("../Source/HF_RPA/HF_RPA_Spurious.jl")
include("../Source/HF_RPA/HF_RPA_Diagonalize.jl")
include("../Source/HF_RPA/HF_RPA_Energy.jl")
include("../Source/HF_RPA/HF_RPA_OBDM.jl")
include("../Source/HF_RPA/HF_RPA_Collectivity.jl")
include("../Source/HF_RPA/HF_RPA_Transitions.jl")
include("../Source/HF_RPA/HF_RPA_Export.jl")
# HF-RRPA
include("../Source/HF_RRPA/HF_RRPA.jl")
include("../Source/HF_RRPA/HF_RRPA_Solver.jl")
include("../Source/HF_RRPA/HF_RRPA_Phonon_Count.jl")
include("../Source/HF_RRPA/HF_RRPA_Hamiltonian.jl")
include("../Source/HF_RRPA/HF_RRPA_OBDM.jl")
include("../Source/HF_RRPA/HF_RRPA_Allocate.jl")
include("../Source/HF_RRPA/HF_RRPA_Diagonalize.jl")
include("../Source/HF_RRPA/HF_RRPA_Spurious.jl")
include("../Source/HF_RRPA/HF_RRPA_Transitions.jl")
include("../Source/HF_RRPA/HF_RRPA_Energy.jl")
include("../Source/HF_RRPA/HF_RRPA_Collectivity.jl")
include("../Source/HF_RRPA/HF_RRPA_Export.jl")

# HF-ScRPA ...
include("../Source/Temp/HF_ScRPA/HF_ScRPA_TBDM.jl")
include("../Source/Temp/HF_ScRPA/HF_ScRPA_Solver.jl")
include("../Source/Temp/HF_ScRPA/HF_ScRPA.jl")
include("../Source/Temp/HF_ScRPA/HF_ScRPA_Allocate.jl")
include("../Source/Temp/HF_ScRPA/HF_ScRPA_Diagonalize.jl")

# BCS
include("../Source/BCS/BCS.jl")
include("../Source/BCS/BCS_Solver.jl")
include("../Source/BCS/BCS_Residual_Interaction.jl")
include("../Source/BCS/BCS_Density_Operator.jl")
include("../Source/BCS/BCS_Chemical_Potential.jl")
include("../Source/BCS/BCS_Particle_Number.jl")
include("../Source/BCS/BCS_Allocate.jl")
include("../Source/BCS/BCS_Energy.jl")
include("../Source/BCS/BCS_Export.jl")

# HFB
include("../Source/HFB/HFB.jl")
include("../Source/HFB/HFB_Solver.jl")
include("../Source/HFB/HFB_Density_Operator.jl")
include("../Source/HFB/HFB_Chemical_Potential.jl")
include("../Source/HFB/HFB_Particle_Number.jl")
include("../Source/HFB/HFB_Broyden.jl")
include("../Source/HFB/HFB_Ordering.jl")
include("../Source/HFB/HFB_Allocate.jl")
include("../Source/HFB/HFB_Canonical_Basis.jl")
include("../Source/HFB/HFB_Energy.jl")
include("../Source/HFB/HFB_Export.jl")

# QTDA solver ...
include("../Source/QTDA/QTDA.jl")
include("../Source/QTDA/QTDA_Solver.jl")
include("../Source/QTDA/QTDA_Allocate.jl")
include("../Source/QTDA/QTDA_Diagonalize.jl")
include("../Source/QTDA/QTDA_Export.jl")

function main_run_dev()

    # Interaction parameters structure ...
    IntParams = Interaction_Parameters(
                NN_File = "IO/NN.bin",
                NNN_File = "IO/NNN.bin",
                hw = 16.0,
                Nmax = 3,
                N2max = 6,
                N3max = 9,
                #cRes = 1.0,
                #cP2N = 1.0,
                #cP3N = 1.0
                )

    # Calculation parameters structure ...
    CalcParams = Calculation_Parameters(
                A = 16,
                Z = 8,
                hw = 16.0,
                Nmax = 3,
                N2max = 6,
                N3max = 9,
                CMS = "CMS1+2B",
                Path = "A16_Z8_hw16.0_Nmax3_N2max6_N3max9",
                HF = HF_Parameters(Tol = 1e-7, IMax = 100, BMF = true, HRF = false),
                RPA = RPA_Parameters(Ortho = true),
                RRPA = RRPA_Parameters(Ortho = true, Tol = 1e-7, IMax = 100, ScV3N = true, ScOBH = true, dOBDM = false)
                )

    # HF calculation call ...
    #HF(Parameters(IntParams,CalcParams))

    # HF-RPA calculation call ...
    #HF_RPA(Parameters(IntParams,CalcParams))

    # HF-RRPA calculation call ...
    #HF_RRPA(Parameters(IntParams,CalcParams))

    CalcParams = Calculation_Parameters(
                A = 18,
                Z = 8,
                hw = 16.0,
                Nmax = 3,
                N2max = 6,
                N3max = 9,
                CMS = "CMS1+2B",
                Path = "BCS_A18_Z8_hw16.0_Nmax3_N2max6_N3max9",
                BCS = BCS_Parameters(ScBCS = false, Tol = 1e-6, IMax = 500, q = 0.05, pD0 = 0.5, nD0 = 0.5, BMF = true, HRF = false),
                HFB = HFB_Parameters(Tol = 1e-6, IMax = 150, dLmax = 0.5, BMF = true, HRF = true)
                )

    # HF-BCS calculation call ...
    BCS(Parameters(IntParams,CalcParams))

    CalcParams = Calculation_Parameters(
                A = 18,
                Z = 8,
                hw = 16.0,
                Nmax = 3,
                N2max = 6,
                N3max = 9,
                CMS = "CMS1+2B",
                Path = "HFB_A18_Z8_hw16.0_Nmax3_N2max6_N3max9",
                BCS = BCS_Parameters(ScBCS = false, Tol = 1e-6, IMax = 500, q = 0.05, pD0 = 0.5, nD0 = 0.5, BMF = true, HRF = false),
                HFB = HFB_Parameters(Tol = 1e-6, IMax = 150, dLmax = 0.5, BMF = true, HRF = true),
                QTDA = QTDA_Parameters(Ortho = true)
                )

    # HFB calculation call ...
    HFB(Parameters(IntParams,CalcParams))

    # QTDA calculation call ...
    QTDA(Parameters(IntParams,CalcParams))

end

@time main_run_dev()

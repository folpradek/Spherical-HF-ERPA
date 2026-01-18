module NuJuliet
    # Import all needed packages ...
    using BenchmarkTools, DelimitedFiles, LinearAlgebra, CGcoefficient, Hungarian, Mmap, Printf

    # Matrix elements & general functions ...
    include("Import_Structures.jl")
    include("Functions/Functions.jl")
    include("MatrixElements/Orb.jl")
    include("MatrixElements/O1B.jl")
    include("MatrixElements/T1B.jl")
    include("MatrixElements/Tr1B.jl")
    include("MatrixElements/OBDM.jl")
    include("MatrixElements/T2B.jl")
    include("MatrixElements/O2B.jl")
    include("MatrixElements/V2B.jl")
    include("MatrixElements/V3B_NO2B.jl")
    include("MatrixElements/qpO2B.jl")
    include("MatrixElements/qpO2B_40.jl")
    include("MatrixElements/qpO2B_31.jl")
    include("MatrixElements/qpO2B_22.jl")

    # HF + MBPT ...
    include("HF/HF.jl")
    include("HF/HF_Solver.jl")
    include("HF/HF_Allocate.jl")
    include("HF/HF_Density_Operator.jl")
    include("HF/HF_Ordering.jl")
    include("HF/HF_Radial_Potential.jl")
    include("HF/HF_Energy.jl")
    include("HF/HF_Export.jl")
    include("HF/HF_MBPT.jl")

    # HF-RPA ...
    include("HF_RPA/HF_RPA_Solver.jl")
    include("HF_RPA/HF_RPA.jl")
    include("HF_RPA/HF_RPA_Phonon_Count.jl") 
    include("HF_RPA/HF_RPA_Allocate.jl")
    include("HF_RPA/HF_RPA_Spurious.jl")
    include("HF_RPA/HF_RPA_Diagonalize.jl")
    include("HF_RPA/HF_RPA_Energy.jl")
    include("HF_RPA/HF_RPA_OBDM.jl")
    include("HF_RPA/HF_RPA_Collectivity.jl")
    include("HF_RPA/HF_RPA_Transitions.jl")
    include("HF_RPA/HF_RPA_Export.jl")

    # HF-Renormalized-RPA ...
    include("HF_RRPA/HF_RRPA.jl")
    include("HF_RRPA/HF_RRPA_Solver.jl")
    include("HF_RRPA/HF_RRPA_Phonon_Count.jl")
    include("HF_RRPA/HF_RRPA_Hamiltonian.jl")
    include("HF_RRPA/HF_RRPA_OBDM.jl")
    include("HF_RRPA/HF_RRPA_Allocate.jl")
    include("HF_RRPA/HF_RRPA_Diagonalize.jl")
    include("HF_RRPA/HF_RRPA_Spurious.jl")
    include("HF_RRPA/HF_RRPA_Transitions.jl")
    include("HF_RRPA/HF_RRPA_Energy.jl")
    include("HF_RRPA/HF_RRPA_Collectivity.jl")
    include("HF_RRPA/HF_RRPA_Export.jl")

    # BCS ...
    include("BCS/BCS.jl")
    include("BCS/BCS_Solver.jl")
    include("BCS/BCS_Residual_Interaction.jl")
    include("BCS/BCS_Density_Operator.jl")
    include("BCS/BCS_Chemical_Potential.jl")
    include("BCS/BCS_Particle_Number.jl")
    include("BCS/BCS_Allocate.jl")
    include("BCS/BCS_Energy.jl")
    include("BCS/BCS_BMBPT.jl")
    include("BCS/BCS_Export.jl")


    # HFB solver ...
    include("HFB/HFB.jl")
    include("HFB/HFB_Solver.jl")
    include("HFB/HFB_Density_Operator.jl")
    include("HFB/HFB_Chemical_Potential.jl")
    include("HFB/HFB_Particle_Number.jl")
    include("HFB/HFB_Broyden.jl")
    include("HFB/HFB_Ordering.jl")
    include("HFB/HFB_Allocate.jl")
    include("HFB/HFB_Canonical_Basis.jl")
    include("HFB/HFB_Energy.jl")
    include("HFB/HFB_BMBPT.jl")
    include("HFB/HFB_Export.jl")

    # QTDA solver ...
    include("QTDA/QTDA.jl")
    include("QTDA/QTDA_Solver.jl")
    include("QTDA/QTDA_Allocate.jl")
    include("QTDA/QTDA_Diagonalize.jl")
    include("QTDA/QTDA_Transitions.jl")
    include("QTDA/QTDA_Export.jl")


    # pHF-RPA ... for testing purposes !!!
        # placement is really temporary ...
    #include("Temp/pHF_RPA/HF_ScRPA_TBDM.jl")
    #include("Temp/pHF_RPA/pHF_RPA_Solver.jl")
    #include("Temp/pHF_RPA/pHF_RPA.jl")
    #include("Temp/pHF_RPA/pHF_RPA_Allocate.jl")
    #include("Temp/pHF_RPA/pHF_RPA_Diagonalize.jl")


    for n in names(@__MODULE__; all=true)
        if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
            @eval export $n
        end
    end
    
end
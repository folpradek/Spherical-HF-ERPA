module HF_RPA_Module
    # Import all needed packages ...
    using BenchmarkTools, DelimitedFiles, LinearAlgebra, CGcoefficient

    # Matrix elements & general functions ...
    include("Import_Structures.jl")
    include("Functions/Functions.jl")
    include("MatrixElements/Orb.jl")
    include("MatrixElements/TrOp.jl")
    include("MatrixElements/T1B.jl")
    include("MatrixElements/T2B.jl")
    include("MatrixElements/V2B.jl")
    include("MatrixElements/V3B_NO2B.jl")
    include("MatrixElements/V2B_Res.jl")

    # HF + LO HF-MBPT solver ...
    include("HF/HF_Solver.jl")
    include("HF/HF_NO2B.jl")
    include("HF/HF_Density_Operator.jl")
    include("HF/HF_Orbital_Ordering.jl")
    include("HF/HF_Radial_Density.jl")
    include("HF/HF_Radial_Potential.jl")
    include("HF/HF_Energy.jl")
    include("HF/HF_Summary.jl")
    include("HF/HF_Export.jl")
    include("HF/HF_MBPT.jl")
    include("HF/HF_Radial_MBPT.jl")

    # HF-RPA solver ...
    include("HF_RPA/HF_RPA_Solver.jl")
    include("HF_RPA/HF_RPA.jl")
    include("HF_RPA/HF_RPA_Phonon_Count.jl") 
    include("HF_RPA/HF_RPA_Allocate.jl")
    include("HF_RPA/HF_RPA_Spurious.jl")
    include("HF_RPA/HF_RPA_Diagonalize.jl")
    include("HF_RPA/HF_RPA_Corr_Energy.jl")
    include("HF_RPA/HF_RPA_OBDM.jl")
    include("HF_RPA/HF_RPA_Collectivity.jl")
    include("HF_RPA/HF_RPA_Radial_Density.jl")
    include("HF_RPA/HF_RPA_Transitions.jl")
    include("HF_RPA/HF_RPA_Export.jl")

    # HF-Extended-RPA solver ... A.K.A. renormalized-RPA ...
    include("HF_ERPA/HF_ERPA_Solver.jl")
    include("HF_ERPA/HF_ERPA.jl")
    include("HF_ERPA/HF_ERPA_Phonon_Count.jl")
    include("HF_ERPA/HF_ERPA_Iteration.jl")
    include("HF_ERPA/HF_ERPA_OBDM.jl")
    include("HF_ERPA/HF_ERPA_Allocate.jl")
    include("HF_ERPA/HF_ERPA_Diagonalize.jl")
    include("HF_ERPA/HF_ERPA_Spurious.jl")
    include("HF_ERPA/HF_ERPA_Transitions.jl")
    include("HF_ERPA/HF_ERPA_Energy.jl")
    include("HF_ERPA/HF_ERPA_Radial_Density.jl")
    include("HF_ERPA/HF_ERPA_Collectivity.jl")
    include("HF_ERPA/HF_ERPA_Export.jl")

    for n in names(@__MODULE__; all=true)
        if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
            @eval export $n
        end
    end
    
end
# Parameter structures ...
    Base.@kwdef struct Interaction_Parameters
        NN_File::String = "none"
        NNN_File::String = "none"
        hw::Float64 = 0.0
        Nmax::Int64 = 0
        N2max::Int64 = 0
        N3max::Int64 = 0
    end

    Base.@kwdef struct Parameters_ERPA
        Sc3N::Bool = true
        ScOBH::Bool = true
        OBDM::String = "Full"
    end

    Base.@kwdef struct Parameters_Pairing
        ScBCS::Bool = false
        s2::Float64 = 1.0
        s3::Float64 = 1.0
        pL0::Float64 = -5.0
        nL0::Float64 = -5.0
        pK0::Float64 = 5.0
        nK0::Float64 = 5.0
        pa0::Float64 = sqrt(2.0)
        na0::Float64 = sqrt(2.0)
    end

    Base.@kwdef struct Calculation_Parameters
        A::Int64 = 0
        Z::Int64 = 0
        hw::Float64 = 0.0
        Nmax::Int64 = 0
        N2max::Int64 = 0
        N3max::Int64 = 0
        CMS::String = "CMS1+2B"
        Ortho::Bool = true
        BMF::Bool = true
        cV_res::Float64 = 1.0
        Path::String = "A" * string(A) * "_Z" * string(Z) *
                       "_hw" * string(hw) * "_Nmax" * string(Nmax) *
                       "_N2max" * string(N2max) * "_N3max" * string(N3max) * "_" * CMS
        Format::String = "Bin"
        ERPA::Parameters_ERPA = Parameters_ERPA(OBDM = "Full", ScOBH = true, Sc3N = true)
        Pairing::Parameters_Pairing = Parameters_Pairing(ScBCS = false, s2 = 1.0, s3 = 1.0)
    end

    struct Parameters
        Int::Interaction_Parameters
        Calc::Calculation_Parameters
    end

# Orbitals, Single-Particle, Particle-Hole & Phonon State structures ...
    struct NOrb
        a::Int64
        n::Int64
        l::Int64
        j::Int64
        pO::Int64
        nO::Int64
    end

    struct SPState
        a::Int64
        l::Int64
        j::Int64
        E::Float64
    end

    struct pnSVector
        p::Vector{SPState}
        n::Vector{SPState}
    end

    struct pnInteger
        p::Int64
        n::Int64
    end

    struct pnFloat
        p::Float64
        n::Float64
    end

    struct pnVector
        p::Vector{Float64}
        n::Vector{Float64}
    end

    struct pnMatrix
        p::Matrix{Float64}
        n::Matrix{Float64}
    end

    struct pnArray
        p::Array
        n::Array
    end

    struct pnCVector
        p::Vector{ComplexF64}
        n::Vector{ComplexF64}
    end

    struct PhState
        p::Int64
        h::Int64
        tz::Int64
        J::Int64
        P::Int64
    end

# Interaction array structures ...
    # Standard bare NN & NNN interaction structures ...
    struct NNInt
        pp::Matrix{Vector{Float64}}
        pn::Matrix{Vector{Float64}}
        nn::Matrix{Vector{Float64}}
    end

    struct NNOrb
        Dic::Dict{Tuple{Int8,Int8,Int8,Int16,Int16},Int32}
        N::Array{Int64}
        Ind::Array{Vector{Vector{Int64}},3}
    end

    struct NNNOrb
        Dic::Dict{Tuple{Int8,Int8,Int8,Int8,Int16,Int16,Int8,Int8,Int8},Int32}
        N::Array{Vector{Int64}}
    end

    # Residual NN interaction structures ...
    struct NNInt_Res
        pp::Matrix{Matrix{Float64}}
        pn::Matrix{Matrix{Float64}}
        nn::Matrix{Matrix{Float64}}
    end

    struct NNOrb_Res
        Dic::Dict{Tuple{Int8,Int8,Int16,Int16},Int32}
        N::Matrix{Int64}
        Ind::Matrix{Vector{Vector{Int64}}}
    end

    struct V2B_Temp
        pp::Matrix{Matrix{Float64}}
        pn::Matrix{Matrix{Float64}}
        nn::Matrix{Matrix{Float64}}
    end

    struct NNOrb_Temp
        Dic::Dict{Tuple{Int8,Int8,Int16,Int16},Int32}
        N::Matrix{Int64}
        Ind::Matrix{Vector{Vector{Int64}}}
    end

    # Residual quasiparticle NN interaction structures ...
    struct qpH40
        pp::Matrix{Vector{Float64}}
        pn::Matrix{Vector{Float64}}
        nn::Matrix{Vector{Float64}}
    end

    struct qpH31
        pp::Matrix{Vector{Float64}}
        pn2011::Matrix{Vector{Float64}}
        pn1120::Matrix{Vector{Float64}}
        nn::Matrix{Vector{Float64}}
    end

    struct qpH22
        pp::Matrix{Vector{Float64}}
        pn2002::Matrix{Vector{Float64}}
        pn1111::Matrix{Vector{Float64}}
        pn0220::Matrix{Vector{Float64}}
        nn::Matrix{Vector{Float64}}
    end

    struct qpH2B
        H40::qpH40
        H31::qpH31
        H22::qpH22
    end

    struct qpH1B
        H11::pnMatrix
        H20::pnMatrix
    end

    struct V2B_H40_Temp
        pp::Matrix{Matrix{Float64}}
        pn::Matrix{Matrix{Float64}}
        nn::Matrix{Matrix{Float64}}
    end

    struct V2B_H31_Temp
        pp::Vector{Matrix{Matrix{Float64}}}
        pn2011::Vector{Matrix{Matrix{Float64}}}
        pn1120::Vector{Matrix{Matrix{Float64}}}
        nn::Vector{Matrix{Matrix{Float64}}}
    end

    struct V2B_H22_Temp
        pp::Vector{Matrix{Matrix{Float64}}}
        pn2002::Matrix{Matrix{Float64}}
        pn1111::Vector{Matrix{Matrix{Float64}}}
        pn0220::Matrix{Matrix{Float64}}
        nn::Vector{Matrix{Matrix{Float64}}}
    end

# Electromagnetic transition operator structures ...
    struct TranOper
        E0::pnMatrix
        E1::pnMatrix
        E2::pnMatrix
        E3::pnMatrix
        M1::pnMatrix
        M2::pnMatrix
        M3::pnMatrix
    end

    struct ReducedMultipole
        E0::pnCVector
        E1::pnCVector
        E2::pnCVector
        E3::pnCVector
    end

    struct Transition
        ph::Vector{Float64}
        is::Vector{Float64}
        iv::Vector{Float64}
    end

    struct ReducedTransition
        E0::Transition
        E1::Transition
        E2::Transition
        E3::Transition
    end

# ERPA auxilliary structures ...
    struct OBDM_Iteration_Pairs
        N::Int64
        T::Vector{Int64}
        a::Vector{Int64}
        b::Vector{Int64}
    end

# HFB structures ...
    # Strcutures for the modified Broyden's algorithm ...
    struct HFB_Broyden_Key
        a::Int64
        b::Int64
    end

    struct HFB_Broyden
        Map::Matrix{Int64}
        Key::Vector{HFB_Broyden_Key}
        M::Int64
        m::Int64
    end

    mutable struct HFB_Broyden_Vector
        x_Rho::pnVector
        y_Rho::pnVector
        z_Rho::pnVector
        X_Rho::pnMatrix
        r_Rho::pnVector
        s_Rho::pnVector
        R_Rho::pnMatrix
        x_Kappa::pnVector
        y_Kappa::pnVector
        z_Kappa::pnVector
        X_Kappa::pnMatrix
        r_Kappa::pnVector
        s_Kappa::pnVector
        R_Kappa::pnMatrix
    end
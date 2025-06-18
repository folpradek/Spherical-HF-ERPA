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
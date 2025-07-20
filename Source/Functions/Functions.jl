function KroneckerDelta(a::Int,b::Int)
    if a == b
        return 1
    else
        return 0
    end
end

pretty_summarysize(x) = Base.format_bytes(Base.summarysize(x))

function lorentzian(x::Float64,w::Float64)
    return w / (x^2 + w^2 / 4.0) / (2.0 * pi)
end

function GeneralizedLaguerre(n::Int,k::Number,x::AbstractArray)
    if n == 0
        return ones(eltype(x), size(x))
    elseif n == 1
        return 1 .+ k .- x
    else
        GL1 = ones(eltype(x), size(x))
        GL2 = 1 .+ k .- x
        GL = GL2
        @inbounds for j in 2:n
            GL = ((k + 2*j - 1 .- x) .* GL2 .+ (-k - j + 1) .* GL1) ./ j
            GL1 = GL2
            GL2 = GL
        end
        return GL
    end
end

function GeneralizedLaguerre(n::Int,k::Number,x::Number)
    if n == 0
        return 1
    elseif n == 1
        return 1 + k - x
    else
        GL1 = 1
        GL2 = 1 + k - x
        GL = GL2
        @inbounds for j in 2:n
            GL = ((k + 2*j - 1 - x) * GL2 + (-k - j + 1) * GL1) / j
            GL1 = GL2
            GL2 = GL
        end
        return GL
    end
end

function doublefactorial(number::Integer)
    fact = one(number)
    @inbounds for m in iseven(number)+1:2:number
        fact *= m
    end
    return fact
end

function Psi_rad_LHO(r::Float64,n::Int64,l::Int64,nu::Float64)
    N = sqrt(sqrt(2 * nu^3 / π) * 2^(n + 2*l + 3) * factorial(n) * nu^l / doublefactorial(2*n + 2*l + 1))
    Psi =  N * exp(-nu * r^2) * r^l * GeneralizedLaguerre(n, l + 0.5, 2 * nu * r^2)
    return Psi
end

function Integrate_Trap(x::Vector{Float64},y::Vector{Float64})
    n = length(x) - 1
    I = 0.0
    @inbounds for i in 1:n
        dx = x[i+1] - x[i]
        I += 0.5 * (y[i] + y[i+1]) * dx
    end
    return I
end

function Radial_Moment_Matrix_LHO(lambda::Int64,HbarOmega::Float64,Orb::Vector{NOrb})
    HbarC = 197.326980
    m_n = 939.565346
    m_p = 938.272013
    b_osc = 1.0/ (HbarC / sqrt(0.5 * (m_p + m_n) * HbarOmega))
    
    a_max = length(Orb)
    Radial_lambda = zeros(Float64,a_max,a_max)

    @inbounds Threads.@threads for a in 1:a_max
        n_a = Orb[a].n
        l_a = Orb[a].l
        @inbounds for b in 1:a_max
            n_b = Orb[b].n
            l_b = Orb[b].l
            ME = Radial_Moment_LHO(lambda,n_a,l_a,n_b,l_b,b_osc)
            Radial_lambda[a,b] = ME
        end
    end

    return Radial_lambda
end

function Radial_Moment_LHO(k::Int64,n_a::Int64,l_a::Int64,n_b::Int64,l_b::Int64,b_osc::Float64)
    gamma = Vector{Float64}(undef,171)
    gamma[1] = 1.0
    gamma[2] = sqrt(pi)
    gamma[3] = 1.0
    @inbounds for n in 3:170
        gamma[n+1] = (Float64(n)/2.0 - 1.0)*gamma[n-1]
    end

    dfactorial = Vector{Float64}(undef,171)
    dfactorial[1] = 1.0
    dfactorial[2] = 1.0
    @inbounds for n in 2:170
        dfactorial[n+1] = Float64(n)*dfactorial[n-1]
    end

    gi = Vector{Float64}(undef,171)
    gi[1] = sqrt(pi) / 2.0
    gi[2] = 0.5
    @inbounds for n in 2:170
        if rem(n,2) == 0
            gi[n+1] = sqrt(pi) * dfactorial[n] / Float64(2^(1+div(n,2)))
        elseif rem(n,2) == 1
            gi[n+1] = dfactorial[n] / Float64(2^(1+div(n-1,2)))
        end
    end

    Amp = 4.0 / (b_osc^(k+0) * sqrt(pi)) * sqrt(Float64(factorial(n_a) * factorial(n_b) * 2^(n_a + n_b + l_a + l_b)) /
                dfactorial[2*n_a + 2*l_a + 2] / dfactorial[2*n_b + 2*l_b + 2])

    Sum = 0.0

    @inbounds for m_a in 0:n_a
        @inbounds for m_b in 0:n_b
            Sum += Float64((-1)^(m_a + m_b) * gamma[2*n_a + 2*l_a + 3 + 1] * gamma[2*n_b + 2*l_b + 3 + 1] *
                    gi[k + 0 + 2 + l_a + l_b + 2*m_a + 2*m_b + 1] / factorial(m_a) / factorial(m_b) /
                    factorial(n_a - m_a) / factorial(n_b - m_b) / gamma[2*l_a + 2*m_a + 3 + 1] / gamma[2*l_b + 2*m_b + 3 + 1])
        end
    end

    I = Sum * Amp

    return I
end

function logm(A::Matrix{Float64})
    D, U = eigen(A)
    l = length(D)
    @inbounds for i in 1:l
        if D[i] < 0.0
            D[i] = -1.0 * D[i]
        end
    end
    D = diagm(log.(D))
    A_log = U * D * U'
    return A_log
end

function JP_Ini(J_max::Int64)
    JP = Vector{Vector{Int64}}(undef,2*(J_max+1))
    JP_count = 0
    @inbounds for J in 0:J_max
        @inbounds for P in 1:2
            JP_count += 1
            JP[JP_count] = [J, P]
        end
    end
    return JP
end
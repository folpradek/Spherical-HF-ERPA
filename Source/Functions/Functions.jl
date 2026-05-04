function KroneckerDelta(a::Int,b::Int)
    if a == b
        return 1
    else
        return 0
    end
end

@inline function kronecker_delta(a::Int,b::Int)
    if a == b
        return 1
    else
        return 0
    end
end

@inline function phase(Arg::Int)
    Phase = isodd(Arg) ? -1.0 : 1.0
    return Phase
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

function integrate_trap(x::Vector{Float64},y::Vector{Float64})
    n = length(x) - 1
    I = 0.0
    @inbounds for i in 1:n
        dx = x[i+1] - x[i]
        I += 0.5 * (y[i] + y[i+1]) * dx
    end
    return I
end

function radial_moment_matrix_LHO(lambda::Int64,hw::Float64,Orb::Vector{Orb1B})
    HbarC = 197.326980
    m_n = 939.565346
    m_p = 938.272013
    b_osc = 1.0/ (HbarC / sqrt(0.5 * (m_p + m_n) * hw))
    
    a_max = length(Orb)
    Radial_lambda = zeros(Float64,a_max,a_max)

    @inbounds Threads.@threads for a in 1:a_max
        n_a = Orb[a].n
        l_a = Orb[a].l
        @inbounds for b in 1:a_max
            n_b = Orb[b].n
            l_b = Orb[b].l
            ME = radial_moment_LHO(lambda,n_a,l_a,n_b,l_b,b_osc)
            Radial_lambda[a,b] = ME
        end
    end

    return Radial_lambda
end

function radial_moment_LHO(k::Int64,n_a::Int64,l_a::Int64,n_b::Int64,l_b::Int64,b_osc::Float64)
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

function JP_initialize(J_max::Int64)
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

function integrate_quadrature(x::Vector{Float64},f::Vector{Float64};M::Int64=0)
    # Read size of field x ...
    N = length(x)

    # Check if the Romberg quadrature is possible ...
    if N < 2
        error("Romberg integration requires at least two points.")
    end

    if N != length(f)
        error("Romberg integration requires arrays of equal sizes!")
    end

    if mod((N-1),2) != 0
        error("Romberg integration requires array size to be of form (2^k + 1)!")
    end

    if M > 0 && M < Int(floor(log2(N-1)))
        n = M
    else
        n = Int(floor(log2(N-1)))
    end

    # Check if at least one Richardson's extrapolation for the Romberg iteration is possible
    # if not, performs trapezoidal rule ...
    if n < 1
        h = x[end] - x[1]
        return h * (f[1] + f[end]) / 2.0
    end

    # Initialize the Romberg table collumn ...
    R = zeros(Float64,n+1)

    # Calculate the Romberg table ... R_j,0 ...
    @inbounds for j in 0:n
        t_j = 2^j
        h_j = (x[end] - x[1]) / t_j

        # Trapezoidal rule ...
        I = (f[1] + f[end]) / 2.0

        s = div((N - 1), t_j)
        @inbounds for k in 1:(t_j - 1)
            I += f[1+k*s]
        end
        R[j+1] = I * h_j

    end

    # Extrapolation for the rest of the Romberg's table ...
    @inbounds for k in 1:n
        @inbounds for j in n:-1:k
            f_j = 4.0^j
            R[j+1] = (f_j * R[j + 1] - R[j]) / (f_j - 1.0)
        end
    end

    return R[n+1]
end

function radial_phase_test(Params::Parameters)
    # Read parameters ...
    A = Params.Calc.A
    hw = Params.Calc.hw
    N_max = Params.Calc.Nmax
    a_max = div((N_max+1)*(N_max+2),2)

    # Initialite the angular momentum algebra ...
    wigner_init_float(75, "Jmax", 9)

    # Define basic constants ...
    hc = 197.326980
    m_n = 939.565346
    m_p = 938.272013
    nu_proton = 0.5 * m_p * hw / hc^2
    nu_neutron = 0.5 * m_n * hw / hc^2
    b_osc = sqrt(0.5 * (939.565346 + 938.272013) * hw) / 197.326980

    # Define the properties of radial grid ...
    r_min, r_max = 0.0 + 1e-8, 2.5 * 1.2 * Float64(A)^(1/3)
    N_grid = 2^14 + 1 # 8192 + 1 grid points ...

    # Make s.p. orbitals ...
    Orb = orbitals_make(Params)

    # Define some local auxiliary functions ...
    @inline function rME_nabla_Omega(l_a::Int64,l_b::Int64)
        rME = 0.0
        if l_a == (l_b + 1)
            rME += - l_b * sqrt(Float64(l_b + 1))
        elseif l_a == (l_b - 1)
            rME += - (l_b + 1) * sqrt(Float64(l_b))
        end
        return rME
    end

    @inline function rME_n(l_a::Int64,l_b::Int64)
        if abs(l_b - 1) <= l_a && l_a <= (l_b + 1)
            rME = sqrt(Float64(2*l_b + 1)) * fCG(2*l_b,2,2*l_a,0,0,0)
            return rME
        else
            return 0.0
        end
    end

    @inline function rME_rN(n_a::Int64,l_a::Int64,n_b::Int64,l_b::Int64,N::Int64,hw::Float64)
        b_osc = sqrt(0.5 * (939.565346 + 938.272013) * hw) / 197.326980
        rME = radial_moment_LHO(N,n_a,l_a,n_b,l_b,b_osc)
        return rME
    end

    @inline function rME_rN_dr(n_a::Int64,l_a::Int64,n_b::Int64,l_b::Int64,N::Int64,hw::Float64)
        b_osc = sqrt(0.5 * (939.565346 + 938.272013) * hw) / 197.326980
        rME = Float64(l_b) * radial_moment_LHO(N-1,n_a,l_a,n_b,l_b,b_osc) - b_osc^2 * radial_moment_LHO(N+1,n_a,l_a,n_b,l_b,b_osc)
        if 1 <= n_b
            rME -= 2.0 * b_osc * sqrt(Float64(n_b)) * radial_moment_LHO(N,n_a,l_a,n_b-1,l_b+1,b_osc)
        end
        return rME
    end

    @inline function dr_Psi_rad_LHO(r::Float64,n::Int64,l::Int64,nu::Float64,h::Float64)
        Psi = Psi_rad_LHO(r,n,l,nu)
        Psi_ph = Psi_rad_LHO(r+h,n,l,nu)
        Psi_p2h = Psi_rad_LHO(r+2.0*h,n,l,nu)
        dPsi = (-3.0 * Psi + 4.0 * Psi_ph - Psi_p2h) / (2.0 * h)
        return dPsi
    end

    @inline function d2r_Psi_rad_LHO(r::Float64,n::Int64,l::Int64,nu::Float64,h::Float64)
        Psi = Psi_rad_LHO(r,n,l,nu)
        Psi_ph = Psi_rad_LHO(r+h,n,l,nu)
        Psi_p2h = Psi_rad_LHO(r+2.0*h,n,l,nu)
        dPsi = (Psi - 2.0 * Psi_ph + Psi_p2h) / (h^2)
        return dPsi
    end

    # Initialize the radial grid ...
    r_grid = range(r_min, stop = r_max, length = N_grid)
    r_grid = collect(r_grid)

    # Initialize the radial grid for single-particle orbitals ...
    pPsi_grid = Matrix{Float64}(undef,a_max,N_grid)
    nPsi_grid = Matrix{Float64}(undef,a_max,N_grid)
    pdPsi_grid = Matrix{Float64}(undef,a_max,N_grid)
    ndPsi_grid = Matrix{Float64}(undef,a_max,N_grid)
    pd2Psi_grid = Matrix{Float64}(undef,a_max,N_grid)
    nd2Psi_grid = Matrix{Float64}(undef,a_max,N_grid)

    # Set the numerical epsilon for differentiation ...
    h = 1e-4
    M = 6

    # Precalculate the radial wave functions and their derivatives on the grid ...
    @inbounds Threads.@threads for i in 1:N_grid
        r = r_grid[i]
        @inbounds for a in 1:a_max
            l_a, j_a, n_a = Orb[a].l, Orb[a].j, Orb[a].n
            pPsi = Psi_rad_LHO(r,n_a,l_a,nu_proton)
            nPsi = Psi_rad_LHO(r,n_a,l_a,nu_neutron)
            pdPsi = dr_Psi_rad_LHO(r,n_a,l_a,nu_proton,h)
            ndPsi = dr_Psi_rad_LHO(r,n_a,l_a,nu_neutron,h)
            pd2Psi = d2r_Psi_rad_LHO(r,n_a,l_a,nu_proton,h)
            nd2Psi = d2r_Psi_rad_LHO(r,n_a,l_a,nu_neutron,h)
            pPsi_grid[a,i] = pPsi
            nPsi_grid[a,i] = nPsi
            pdPsi_grid[a,i] = pdPsi
            ndPsi_grid[a,i] = ndPsi
            pd2Psi_grid[a,i] = pd2Psi
            nd2Psi_grid[a,i] = nd2Psi
        end
    end

    # Initialize the matrix for matrix elements ...

    T_Miy = zeros(Float64,a_max,a_max)
    T_Num = zeros(Float64,a_max,a_max)
    T_Num2 = zeros(Float64,a_max,a_max)
    Grad_Miy = zeros(Float64,a_max,a_max)
    Grad_Num = zeros(Float64,a_max,a_max)
    Grad_Num2 = zeros(Float64,a_max,a_max)
    r2p_Miy = zeros(Float64,a_max,a_max)
    r2p_Num = zeros(Float64,a_max,a_max)
    r2p_Num2 = zeros(Float64,a_max,a_max)

    # Calculate the 1-body kinetic operator ... Miyagi convention ...
    T = T1b(Params,Orb)

    # Allocate the 1-body kinetic operator matrix elements in Miyagi convention ...
    T_Miy .= T.p

    # Calculate all the other operators ...
    @inbounds for a in 1:a_max
        n_a, l_a, j_a = Orb[a].n, Orb[a].l, Orb[a].j
        @inbounds for b in 1:a_max
            n_b, l_b, j_b = Orb[b].n, Orb[b].l, Orb[b].j
            if rem(l_a + l_b + 1,2) == 0 && abs(j_a - j_b) <= 2 && (j_a + j_b) >= 2
                Grad_Miy[a,b] = rGrad(a,b,Orb)
                Grad_Num2[a,b] = (integrate_quadrature(r_grid, r_grid.^2 .* pPsi_grid[a,:] .* pdPsi_grid[b,:];M=M) * rME_n(l_a,l_b) + 
                                  integrate_quadrature(r_grid, r_grid .* pPsi_grid[a,:] .* pPsi_grid[b,:];M=M) * rME_nabla_Omega(l_a,l_b)) *
                                 f6j(j_a,j_b,2,2*l_b,2*l_a,1) * sqrt(Float64((j_a + 1) * (j_b + 1))) * Float64((-1)^(l_a + 1 + div(j_b + 1,2))) / b_osc * phase(n_a + n_b)
                Grad_Num[a,b] = (integrate_quadrature(r_grid, r_grid.^2 .* pPsi_grid[a,:] .* pdPsi_grid[b,:];M=M) * rME_n(l_a,l_b) + 
                                  integrate_quadrature(r_grid, r_grid .* pPsi_grid[a,:] .* pPsi_grid[b,:];M=M) * rME_nabla_Omega(l_a,l_b)) *
                                 f6j(j_a,j_b,2,2*l_b,2*l_a,1) * sqrt(Float64((j_a + 1) * (j_b + 1))) * Float64((-1)^(l_a + 1 + div(j_b + 1,2))) / b_osc
            end

            if l_a == l_b && j_a == j_b
                T_Num2[a,b] = -0.5 * hc^2 / m_p * (integrate_quadrature(r_grid, r_grid.^2 .* pPsi_grid[a,:] .* pd2Psi_grid[b,:];M=M) +
                                            2.0 * integrate_quadrature(r_grid, r_grid .* pPsi_grid[a,:] .* pdPsi_grid[b,:];M=M) -
                                l_a * (l_a + 1) * integrate_quadrature(r_grid, pPsi_grid[a,:] .* pPsi_grid[b,:];M=M)) * phase(n_a + n_b)
                T_Num[a,b] = -0.5 * hc^2 / m_p * (integrate_quadrature(r_grid, r_grid.^2 .* pPsi_grid[a,:] .* pd2Psi_grid[b,:];M=M) +
                                            2.0 * integrate_quadrature(r_grid, r_grid .* pPsi_grid[a,:] .* pdPsi_grid[b,:];M=M) -
                                l_a * (l_a + 1) * integrate_quadrature(r_grid, pPsi_grid[a,:] .* pPsi_grid[b,:];M=M))
            end

            if j_a == j_b && l_a == l_b
                r2p_Num2[a,b] = integrate_quadrature(r_grid, r_grid.^4 .* pPsi_grid[a,:] .* pPsi_grid[b,:];M=M) * b_osc^2 * phase(n_a + n_b)
                r2p_Num[a,b] = integrate_quadrature(r_grid, r_grid.^4 .* pPsi_grid[a,:] .* pPsi_grid[b,:];M=M) * b_osc^2
            end

            if a == b
                r2p_Miy[a,b] += Float64(2*n_a + l_a) + 1.5
            end

            if l_a == l_b && j_a == j_b && n_a == (n_b + 1)
                r2p_Miy[a,b] -= sqrt(n_a * (n_a + l_a + 0.5))
            end

            if l_a == l_b && j_a == j_b && n_b == (n_a + 1)
                r2p_Miy[a,b] -= sqrt(n_b * (n_b + l_b + 0.5))
            end

        end
    end

    # Display the results ...
    println("\n\n\n___________________________________________________________________")
    println("1-body kinetic operator matrix elements in Miyagi convention:")
    display(T_Miy)
    println("___________________________________________________________________")
    println("1-body kinetic operator matrix elements numerical calculation (original phase):")
    display(T_Num)
    println("___________________________________________________________________")
    println("1-body kinetic operator matrix elements numerical calculation (updated phase):")
    display(T_Num2)
    println("\n\n\n___________________________________________________________________")
    println("1-body gradient operator matrix elements in Miyagi convention:")
    display(Grad_Miy)
    println("___________________________________________________________________")
    println("1-body gradient operator matrix elements numerical calculation (original phase):")
    display(Grad_Num)
    println("___________________________________________________________________")
    println("1-body gradient operator matrix elements numerical calculation (updated phase):")
    display(Grad_Num2)
    println("\n\n\n___________________________________________________________________")
    println("1-body r^2 p operator matrix elements in Miyagi convention:")
    display(r2p_Miy)
    println("___________________________________________________________________")
    println("1-body r^2 p operator matrix elements numerical calculation (original phase):")
    display(r2p_Num)
    println("___________________________________________________________________")
    println("1-body r^2 p operator matrix elements numerical calculation (updated phase):")
    display(r2p_Num2)



    return
end
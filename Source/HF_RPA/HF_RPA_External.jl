function HF_RPA_import_binary(Params::Parameters,N_nu::Matrix{Int64})
    # Parameters / dimensions
    N_max = Params.Calc.Nmax
    J_max = 2*N_max + 1
    OutputFile = Params.Calc.Path

    # Initialite the list of values of J & P for iteration ...
    JP_list = JP_initialize(J_max)

    println("\nImporting TDA & RPA solutions from binary format ...")

    # Initialize amplitudes & energies ...
        # Case of TDA ...
    X_TDA = Matrix{Matrix{Float64}}(undef,J_max+1,2)
    E_TDA = Matrix{Vector{Float64}}(undef,J_max+1,2)
        # Case of RPA ...
    X_RPA = Matrix{Matrix{ComplexF64}}(undef,J_max+1,2)
    Y_RPA = Matrix{Matrix{ComplexF64}}(undef,J_max+1,2)
    E_RPA = Matrix{Vector{ComplexF64}}(undef,J_max+1,2)

    # Initialite each J & P block ...
    @inbounds for JP in JP_list
        J, P = JP[1], JP[2]
        N_ph = N_nu[J+1,P]
        X_TDA[J+1,P] = Matrix{Float64}(undef,N_ph,N_ph)
        E_TDA[J+1,P] = Vector{Float64}(undef,N_ph)
        X_RPA[J+1,P] = Matrix{ComplexF64}(undef,N_ph,N_ph)
        Y_RPA[J+1,P] = Matrix{ComplexF64}(undef,N_ph,N_ph)
        E_RPA[J+1,P] = Vector{ComplexF64}(undef,N_ph)
    end

    # Read TDA amplitudes X ...
    Import_Path_TDA_X = "IO/$OutputFile/Bin/TDA_X.bin"
    open(Import_Path_TDA_X, "r") do Import_File
        @inbounds for JP in JP_list
            J, P = JP[1], JP[2]
            N_ph = N_nu[J+1,P]
            N_ph == 0 && continue
            Buffer = Vector{Float64}(undef,N_ph*N_ph)
            read!(Import_File,Buffer)
            X_TDA[J+1,P] .= reshape(Buffer,N_ph,N_ph)
        end
    end

    # Read TDA energies E ...
    Import_Path_TDA_X = "IO/$OutputFile/Bin/TDA_E.bin"
    open(Import_Path_TDA_X, "r") do Import_File
        @inbounds for JP in JP_list
            J, P = JP[1], JP[2]
            N_ph = N_nu[J+1,P]
            N_ph == 0 && continue
            Buffer = Vector{Float64}(undef,N_ph)
            read!(Import_File,Buffer)
            E_TDA[J+1,P] .= Buffer
        end
    end

    # Read RPA amplitudes X ...
    Import_Path_RPA_X = "IO/$OutputFile/Bin/RPA_X.bin"
    open(Import_Path_RPA_X, "r") do Import_File
        @inbounds for JP in JP_list
            J, P = JP[1], JP[2]
            N_ph = N_nu[J+1,P]
            N_ph == 0 && continue
            Buffer = Vector{ComplexF64}(undef,N_ph*N_ph)
            read!(Import_File,Buffer)
            X_RPA[J+1,P] .= reshape(Buffer,N_ph,N_ph)
        end
    end

    # Read RPA amplitudes Y ...
    Import_Path_RPA_Y = "IO/$OutputFile/Bin/RPA_Y.bin"
    open(Import_Path_RPA_Y, "r") do Import_File
        @inbounds for JP in JP_list
            J, P = JP[1], JP[2]
            N_ph = N_nu[J+1,P]
            N_ph == 0 && continue
            Buffer = Vector{ComplexF64}(undef,N_ph*N_ph)
            read!(Import_File,Buffer)
            Y_RPA[J+1,P] .= reshape(Buffer,N_ph,N_ph)
        end
    end

    # Read RPA energies E ...
    Import_Path_RPA_X = "IO/$OutputFile/Bin/RPA_E.bin"
    open(Import_Path_RPA_X, "r") do Import_File
        @inbounds for JP in JP_list
            J, P = JP[1], JP[2]
            N_ph = N_nu[J+1,P]
            N_ph == 0 && continue
            Buffer = Vector{ComplexF64}(undef,N_ph)
            read!(Import_File,Buffer)
            E_RPA[J+1,P] .= Buffer
        end
    end

    println("\tTDA & RPA solutions have been succesufully imported ...")

    return X_TDA, E_TDA, X_RPA, Y_RPA, E_RPA
end

function HF_RPA_transition_density(Params::Parameters;JP::String,nu_list::Vector{Int64},File_Name::String = "Phonon_Densities",Weighted::Bool = false)
    # Read parameters ...
    A = Params.Calc.A
    hw = Params.Calc.hw
    N_max = Params.Calc.Nmax
    a_max = div((N_max+1)*(N_max+2),2)

    # Define basic constants ...
    hbarc = 197.326980
    m_n = 939.565346
    m_p = 938.272013
    nu_proton = 0.5 * m_p * hw / hbarc^2
    nu_neutron = 0.5 * m_n * hw / hbarc^2

    # Define the properties of radial grid ...
    r_min, r_max = 0.0 + 1e-8, 2.5 * 1.2 * Float64(A)^(1/3)
    N_grid = 2^13 + 1 # 8192 + 1 grid points ...

    # Read the values of J & P ...
    P_char = string(JP[end])
    J = parse(Int64,JP[1:end-1])

    # Check input formating ...
    if (P_char != "+") && (P_char != "-")
        error("Invalid JP = \"$JP\"; use J+ or J-.")
    end

    # Set internal value of parity P ...
    P = 1
    if P_char == "-"
        P = 2
    end

    # Make s.p. orbitals ...
    Orb = orbitals_make(Params)

    # Prepare Particle & Hole orbitals ...
    N_Particle, Particle, N_Hole, Hole = orbitals_ph_make(Params,Orb)

    # Prepare 1p-1h phonon states ...
    N_Phonon, Phonon = orbitals_one_phonon_make(N_Particle,Particle,N_Hole,Hole)

    # Count & pre-index all phonon states in JP subspaces ...
    N_nu, Orb_Phonon = HF_RPA_phonon_count(Params,N_Phonon,Phonon)

    # Import transformation matrix mapping LHO & reference basis ...
    @time C = O1b_import(Params,Orb,"IO/" * Params.Calc.Path * "/Bin/C_HF.bin")

    # Import TDA & RPA solutions ...
    @time X_TDA, E_TDA, X_RPA, Y_RPA, E_RPA = HF_RPA_import_binary(Params,N_nu)

    # Determine the number of phonon levels & the total number of phonons ...
    N_phonon = length(nu_list)
    N_ph = N_nu[J+1,P]

    # Initialize the radial grid ...
    r_grid = range(r_min, stop = r_max, length = N_grid)
    r_grid = collect(r_grid)

    # Initialize the radial grid for single-particle orbitals ...
    pOrb_grid = Matrix{Float64}(undef,a_max,N_grid)
    nOrb_grid = Matrix{Float64}(undef,a_max,N_grid)

    # Initialize the radial grid for 1-body transition operators ...
    pM_grid = ones(Float64,a_max,N_grid)
    nM_grid = ones(Float64,a_max,N_grid)

    # Calculate the radial representation of the reference basis orbitals ...
    @inbounds Threads.@threads for i in 1:N_grid
        r = r_grid[i]
        @inbounds for a in 1:a_max
            l_a, j_a = Orb[a].l, Orb[a].j
            pSum, nSum = 0.0, 0.0
            @inbounds for k in 1:a_max
                n_k, l_k, j_k = Orb[k].n, Orb[k].l, Orb[k].j
                if l_a != l_k || j_a != j_k
                    continue
                end
                pPsi = C.p[k,a] * Psi_rad_LHO(r,n_k,l_k,nu_proton)
                nPsi = C.p[k,a] * Psi_rad_LHO(r,n_k,l_k,nu_neutron)
                pSum += pPsi
                nSum += nPsi
            end
            pOrb_grid[a,i] = pSum
            nOrb_grid[a,i] = nSum
        end
    end

    # If set to true, calculate the radial representation of 1-body transition operator ...
    if Weighted == true
        println("to be implemented ...")
    end

    # Initialize radial grids for given phonon levels ...
        # Case of TDA ...
    pRho_nu_TDA = Matrix{Float64}(undef,N_phonon,N_grid)
    nRho_nu_TDA = Matrix{Float64}(undef,N_phonon,N_grid)
        # Case of RPA ...
    pRho_nu_RPA = Matrix{Float64}(undef,N_phonon,N_grid)
    nRho_nu_RPA = Matrix{Float64}(undef,N_phonon,N_grid)

    # For given phonon levels evaluate their radial representation ...
    @inbounds Threads.@threads for i in 1:N_grid
        r = r_grid[i]
        @inbounds for (nu_ind, nu) in enumerate(nu_list)
            pTDASum, nTDASum = 0.0, 0.0
            pRPASum, nRPASum = 0.0, 0.0
            @inbounds for ph in 1:N_ph
                i_ph = Orb_Phonon[J+1,P][ph]
                p, h, t_ph = Phonon[i_ph].p, Phonon[i_ph].h, Phonon[i_ph].tz
                x_RPA, y_RPA = real(X_RPA[J+1,P][ph,nu]), real(Y_RPA[J+1,P][ph,nu])
                x_TDA = X_TDA[J+1,P][ph,nu]
                if t_ph == -1
                    pRho_TDA = Float64((-1)^J) * x_TDA * pOrb_grid[p,i] * pOrb_grid[h,i] * pM_grid[p,i] * pM_grid[h,i]
                    pRho_RPA = (Float64((-1)^J) * x_RPA + y_RPA) * pOrb_grid[p,i] * pOrb_grid[h,i] * pM_grid[p,i] * pM_grid[h,i]
                    pTDASum += pRho_TDA
                    pRPASum += pRho_RPA
                elseif t_ph == 1
                    nRho_TDA = Float64((-1)^J) * x_TDA * nOrb_grid[p,i] * nOrb_grid[h,i] * nM_grid[p,i] * nM_grid[h,i]
                    nRho_RPA = (Float64((-1)^J) * x_RPA + y_RPA) * nOrb_grid[p,i] * nOrb_grid[h,i] * nM_grid[p,i] * nM_grid[h,i]
                    nTDASum += nRho_TDA
                    nRPASum += nRho_RPA
                end
            end
            pRho_nu_TDA[nu_ind,i] = pTDASum * r^2 / Float64(2*J + 1)
            nRho_nu_TDA[nu_ind,i] = nTDASum * r^2 / Float64(2*J + 1)
            pRho_nu_RPA[nu_ind,i] = pRPASum * r^2 / Float64(2*J + 1)
            nRho_nu_RPA[nu_ind,i] = nRPASum * r^2 / Float64(2*J + 1)
        end
    end

    # Perform the export of radial phonon densities ...
    @inbounds for (nu_ind, nu) in enumerate(nu_list)
        e_TDA, e_RPA = E_TDA[J+1,P][nu], real(E_RPA[J+1,P][nu])
        sE_TDA, sE_RPA = string(round(e_TDA,digits=3)), string(round(e_RPA,digits=3))

        # Export TDA transition densities ...
        Output_File_TDA = "IO/" * Params.Calc.Path * "/RPA/Densities/TDA_" * File_Name * "_JP$JP _nu$nu _E$sE_TDA .dat"
        open(Output_File_TDA, "w") do Write_File
            @printf(Write_File, "%-20s %-20s %-20s %-20s %-20s\n", "r", "rho_p", "rho_n", "rho_is", "rho_iv")
            @inbounds for i in 1:N_grid
                r, pRho, nRho = r_grid[i], pRho_nu_TDA[nu_ind,i], nRho_nu_TDA[nu_ind,i]
                isRho = 0.5 * (pRho + nRho)
                ivRho = 0.5 * (pRho - nRho)
                @printf(Write_File, "%-20.8f %-20.8f %-20.8f %-20.8f %-20.8f\n", r, pRho, nRho, isRho, ivRho)
            end
        end

        # Export TDA transition densities ...
        Output_File_RPA = "IO/" * Params.Calc.Path * "/RPA/Densities/RPA_" * File_Name * "_JP$JP _nu$nu _E$sE_RPA .dat"
        open(Output_File_RPA, "w") do Write_File
            @printf(Write_File, "%-20s %-20s %-20s %-20s %-20s\n", "r", "rho_p", "rho_n", "rho_is", "rho_iv")
            @inbounds for i in 1:N_grid
                r, pRho, nRho = r_grid[i], pRho_nu_RPA[nu_ind,i], nRho_nu_RPA[nu_ind,i]
                isRho = 0.5 * (pRho + nRho)
                ivRho = 0.5 * (pRho - nRho)
                @printf(Write_File, "%-20.8f %-20.8f %-20.8f %-20.8f %-20.8f\n", r, pRho, nRho, isRho, ivRho)
            end
        end

    end







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

    @inline function rME_YL(l_a::Int64,j_a::Int64,l_b::Int64,j_b::Int64,L::Int64)
        if rem(l_a + l_b + L,2) != 0
            return 0.0
        end
        rME = phase(L + div(j_a - 1,2)) * sqrt(Float64((j_a + 1) * (j_b + 1)) / (4.0 * pi)) * fCG(j_a,j_b,2*L,1,-1,0)
        return rME
    end

    @inline function rME_YL_cross_n(l_a::Int64,l_b::Int64,L::Int64,J::Int64)
        rME = 0.0
        if abs(J - l_b) > l_a || l_a > (J + l_b)
            return rME
        end
        Amp = phase(l_a + l_b + J) * sqrt(Float64((2*J + 1))) 
        @inbounds for l in max(abs(L - l_a),abs(1 - l_b)):min(L + l_a, 1 + l_b)
            if rem(l_a + l + L,2) == 0 && rem(l_b + l + 1,2) == 0
                ME = f6j(2*L,2,2*J,2*l_b,2*l_a,2*l) * rME_YL(l_a,l,L) * rME_n(l,l_b)
                rME += ME
            end
        end
        rME = Amp * rME
        return rME
    end

    @inline function rME_n_cross_YL(l_a::Int64,l_b::Int64,L::Int64,J::Int64)
        rME = 0.0
        if abs(J - l_b) > l_a || l_a > (J + l_b)
            return rME
        end
        Amp = phase(l_a + l_b + J) * sqrt(Float64((2*J + 1))) 
        @inbounds for l in max(abs(1 - l_a),abs(L - l_b)):min(1 + l_a, L + l_b)
            if rem(l_a + l + 1,2) == 0 && rem(l_b + l + L,2) == 0
                ME = f6j(2,2*L,2*J,2*l_b,2*l_a,2*l) * rME_YL(l,l_b,L) * rME_n(l_a,l)
                rME += ME
            end
        end
        rME = Amp * rME
        return rME
    end

    @inline function rME_YL_cross_nabla_Omega(l_a::Int64,l_b::Int64,L::Int64,J::Int64)
        rME = 0.0
        if abs(J - l_b) > l_a || l_a > (J + l_b)
            return rME
        end
        Amp = phase(l_a + l_b + J) * sqrt(Float64((2*J + 1) * (2*L + 1)) / (4.0 * pi))
        @inbounds for l in max(abs(L - l_a), abs(1 - l_b)):min(L + l_a, 1 + l_b)
            if rem(l_a + l + L, 2) == 0 && rem(l + l_b + 1, 2) == 0
                ME = sqrt(Float64(2*l + 1)) * f6j(2*L,2,2*J,2*l_b,2*l_a,2*l) * fCG(2*l,2*L,2*l_a,0,0,0) * rME_nabla_Omega(l,l_b)
                rME += ME
            end
        end
        rME = Amp * rME
        return rME
    end

    @inline function rME_nabla_Omega_cross_YL(l_a::Int64,l_b::Int64,L::Int64,J::Int64)
        rME = 0.0
        if abs(J - l_b) > l_a || l_a > (J + l_b)
            return rME
        end
        Amp = phase(l_a + l_b + J) * sqrt(Float64((2*J + 1) * (2*L + 1)) / (4.0 * pi))
        @inbounds for l in max(abs(1 - l_a), abs(L - l_b)):min(1 + l_a, L + l_b)
            if rem(l_a + l + 1, 2) == 0 && rem(l + l_b + L, 2) == 0
                ME = sqrt(Float64(2*l_b + 1)) * f6j(2,2*L,2*J,2*l_b,2*l_a,2*l) * fCG(2*l_b,2*L,2*l,0,0,0) * rME_nabla_Omega(l_a,l)
                rME += ME
            end
        end
        rME = Amp * rME
        return rME
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

    @inline function rME_YL(l_a::Int64,l_b::Int64,L::Int64)
        rME = sqrt(Float64((2*l_b + 1) * (2*L + 1)) / (4.0 * pi)) * fCG(2*l_b,2*L,2*l_a,0,0,0)
        return rME
    end

    @inline function rME_YJL(l_a::Int64,l_b::Int64,L::Int64,J::Int64)
        rME = 0.0
        Amp = sqrt(Float64(2*J + 1)) * phase(l_a + l_b + J)
        @inbounds for l in max(abs(L - l_a),abs(1 - l_b)):min(L + l_a, 1 + l_b)
            if rem(l_a + l + L,2) == 0 && rem(l_b + l + 1,2) == 0
                ME = f6j(2*L,2,2*J,2*l_b,2*l_a,2*l) * rME_YL(l_a,l,L) * rME_n(l,l_b)
                rME += ME
            end
        end
        rME *= Amp
        return rME
    end

    @inline function rME_nabla_dot_rN_YJL(a::Int64,b::Int64,N::Int64,L::Int64,J::Int64,hw::Float64,Orb::Vector{Orb1B})
        n_a, l_a, j_a = Orb[a].n, Orb[a].l, Orb[a].j
        n_b, l_b, j_b = Orb[b].n, Orb[b].l, Orb[b].j

        rME = phase(l_a + J + div(j_b + 1,2)) *  sqrt(Float64((j_a + 1) * (j_b + 1))) * f6j(2*l_b,j_b,1,j_a,2*l_a,2*J) *
                ((rME_rN_dr(n_a,l_a,n_b,l_b,N,hw) - rME_rN_dr(n_b,l_b,n_a,l_a,N,hw)) * rME_YL(l_a,l_b,J) *
                 (kronecker_delta(L,J-1) * sqrt(J / (2*J + 1)) - kronecker_delta(L,J+1) * sqrt((J + 1) / (2*J + 1))) +
                  rME_rN(n_a,l_a,n_b,l_b,N-1,hw) * (rME_YL_cross_nabla_Omega(l_a,l_b,L,J) + rME_YL_cross_nabla_Omega(l_b,l_a,L,J)))

        return rME
    end

    @inline function rME_nabla_cross_S_dot_rN_YJL(a::Int64,b::Int64,N::Int64,L::Int64,J::Int64,hw::Float64,Orb::Vector{Orb1B})
        n_a, l_a, j_a = Orb[a].n, Orb[a].l, Orb[a].j
        n_b, l_b, j_b = Orb[b].n, Orb[b].l, Orb[b].j

        rME = 0.0

        if L != 0 && N != 0
            rME = 0.5 * sqrt(6.0 * Float64((j_a + 1) * (2*J + 1) * (j_b + 1))) * f9j(2*l_a,1,j_a,2*l_b,1,j_b,2*J,2,2*J) *
                    rME_rN(n_a,l_a,n_b,l_b,N-1,hw) * rME_YL(l_a,l_b,J) * (kronecker_delta(L,J-1) * sqrt(Float64(J + 1) / Float64(2*J + 1)) * Float64(N - J + 1) -
                    kronecker_delta(L,J+1) * sqrt(Float64(J) / Float64(2*J + 1)) * Float64(N + J + 2))
        end

        return rME
    end

    @inline function rME_rN_YL(a::Int64,b::Int64,N::Int64,L::Int64,hw::Float64,Orb::Vector{Orb1B})
        l_a, l_b = Orb[a].l, Orb[b].l
        if rem(l_a + l_b + L,2) == 0
            n_a, j_a = Orb[a].n, Orb[a].j
            n_b, j_b = Orb[b].n, Orb[b].j
            rY = phase(L + div(j_a - 1,2)) * sqrt(Float64((j_a + 1) * (j_b + 1)) / (4.0 * pi)) * fCG(j_a,j_b,2*L,1,-1,0)

            b_osc = sqrt(0.5 * (939.565346 + 938.272013) * hw) / 197.326980
            rN = radial_moment_LHO(N,n_a,l_a,n_b,l_b,b_osc)

            rME = rN * rY
            return rME
        else
            return 0.0
        end
    end

    wigner_init_float(75, "Jmax", 9)


    @inbounds Threads.@threads for i in 1:N_grid
        r = r_grid[i]
        @inbounds for a in 1:a_max
            l_a, j_a, n_a = Orb[a].l, Orb[a].j, Orb[a].n
            pPsi = Psi_rad_LHO(r,n_a,l_a,nu_proton)
            nPsi = Psi_rad_LHO(r,n_a,l_a,nu_neutron)
            pOrb_grid[a,i] = pPsi
            nOrb_grid[a,i] = nPsi
        end
    end

    rm1 = zeros(Float64,a_max,a_max)

    @inbounds for a in 1:a_max
        n_a, l_a, j_a = Orb[a].n, Orb[a].l, Orb[a].j
        @inbounds for b in 1:a_max
            n_b, l_b, j_b = Orb[b].n, Orb[b].l, Orb[b].j
            @views r = r_grid[:]
            @views pO_a = pOrb_grid[a,:]
            @views pO_b = pOrb_grid[b,:]

            prm1 = integrate_trap(r_grid, pO_a .* pO_b .* r)
            rm1[a,b] = prm1
        end
    end


    Grad1 = zeros(Float64,a_max,a_max)
    Grad2 = zeros(Float64,a_max,a_max)

    b_osc = sqrt(0.5 * (939.565346 + 938.272013) * hw) / 197.326980

    @inbounds for a in 1:a_max
        n_a, l_a, j_a = Orb[a].n, Orb[a].l, Orb[a].j
        @inbounds for b in 1:a_max
            n_b, l_b, j_b = Orb[b].n, Orb[b].l, Orb[b].j
            if rem(l_a + l_b + 1,2) == 0 && abs(j_a - j_b) <= 2 && (j_a + j_b) >= 2
                Grad1[a,b] = rGrad(a,b,Orb) * b_osc
                Grad2[a,b] = phase(l_a + 1 + div(j_b + 1,2)) * f6j(j_a,2*l_a,1,2*l_b,j_b,2) * sqrt(Float64((j_a + 1) * (j_b + 1))) * 
                            (rME_rN_dr(n_a,l_a,n_b,l_b,0,hw) * rME_n(l_a,l_b) + rm1[a,b] * rME_nabla_Omega(l_a,l_b))
                            #(rME_rN_dr(n_a,l_a,n_b,l_b,0,hw) * rME_n(l_a,l_b) + 0.0*radial_moment_LHO(-1,n_a,l_a,n_b,l_b,hw) * rME_nabla_Omega(l_a,l_b))
            end
        end
    end

    display(Grad1)
    display(Grad2)

    return
end
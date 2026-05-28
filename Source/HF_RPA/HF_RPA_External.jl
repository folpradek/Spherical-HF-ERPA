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

function HF_RPA_transition_densities(Params::Parameters;JP::String,nu_list::Vector{Int64},File_Name::String = "Phonon_Densities",Weighted::Bool = false)
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

        # Export RPA transition densities ...
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

    return
end

function HF_RPA_transition_currents(Params::Parameters;JP::String,nu_list::Vector{Int64},File_Name::String = "Phonon_Currents",M::Int64=-1)
    # Initialize the angular momentum algebra ...
    wigner_init_float(75, "Jmax", 9)
    
    # Read parameters ...
    A = Params.Calc.A
    Z = Params.Calc.Z
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
    r_min, r_max = 0.01, 2.5 * 1.2 * Float64(A)^(1/3)
    rN_grid = 2^13 + 1 # 8192 + 1 grid points ...
    tN_grid = 256

    # Read the values of J & P ...
    P_char = string(JP[end])
    J = parse(Int64,JP[1:end-1])

    # Check input formating ...
    if (P_char != "+") && (P_char != "-")
        error("Invalid JP = \"$JP\"; use J+ or J- ...")
    end

    # Set internal value of parity P ...
    P = 1
    if P_char == "-"
        P = 2
    end

    # Set the value of projection M ...
    if M == -1
        M = J
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

    # Precalculate the reduced matrix elements of Electrogmanetic 1-body transition operators ...
        # Iniztialize the 1-body Electromagnetic transition operators ...
    TrOp = Tr1b_initialize(Params,Orb,1.0)
         # Transform the 1-body transition operators to the reference basis ...
    TrOp = Tr1b_transformation(Params,Orb,C,TrOp)

    # Precalculate the matrix elements of 1-body transition operators ...
        # Initialize the matrix elements for given phonon levels ...
    pM, nM = Vector{Float64}(undef,N_phonon), Vector{Float64}(undef,N_phonon)
    isM, ivM = Vector{Float64}(undef,N_phonon), Vector{Float64}(undef,N_phonon)
        # Calculate the matrix elements for given phonon levels ...
    if (J,P) in [(0,1),(2,1),(3,2)]
        @inbounds for (nu_ind, nu) in enumerate(nu_list)
            pSum, nSum, isSum, ivSum = 0.0, 0.0, 0.0, 0.0

            @inbounds for ph in 1:N_ph
                i_ph = Orb_Phonon[J+1,P][ph]
                p, h, t_ph = Phonon[i_ph].p, Phonon[i_ph].h, Phonon[i_ph].tz
                x_TDA = X_TDA[J+1,P][ph,nu]

                MME = phase(J) / sqrt(Float64(2*J + 1)) * x_TDA

                if t_ph == -1
                    a_p, a_h = Particle.p[p].a, Hole.p[h].a

                    if J == 0 && P == 1
                        MME = TrOp.E0.p[a_p,a_h] * MME

                    elseif J == 2 && P == 1
                        MME = TrOp.E2.p[a_p,a_h] * MME

                    elseif J == 3 && P == 2
                        MME = TrOp.E3.p[a_p,a_h] * MME

                    end

                    pSum += MME
                    isSum += 0.5 * MME
                    ivSum += 0.5 * MME

                elseif t_ph == 1
                    a_p, a_h = Particle.n[p].a, Hole.n[h].a

                    if J == 0 && P == 1
                        MME = TrOp.E0.n[a_p,a_h] * MME

                    elseif J == 2 && P == 1
                        MME = TrOp.E2.n[a_p,a_h] * MME

                    elseif J == 3 && P == 2
                        MME = TrOp.E3.n[a_p,a_h] * MME

                    end

                    nSum += MME
                    isSum += 0.5 * MME
                    ivSum -= 0.5 * MME

                end

            end

            pM[nu_ind] = pSum
            nM[nu_ind] = nSum
            isM[nu_ind] = isSum
            ivM[nu_ind] = ivSum
        end

    elseif (J,P) == (1,2)
        @inbounds for (nu_ind, nu) in enumerate(nu_list)
            pSum, nSum, isSum, ivSum = 0.0, 0.0, 0.0, 0.0

            @inbounds for ph in 1:N_ph
                i_ph = Orb_Phonon[J+1,P][ph]
                p, h, t_ph = Phonon[i_ph].p, Phonon[i_ph].h, Phonon[i_ph].tz
                x_TDA = X_TDA[J+1,P][ph,nu]

                Amp = phase(J) / sqrt(Float64(2*J + 1)) * x_TDA

                if t_ph == -1
                    a_p, a_h = Particle.p[p].a, Hole.p[h].a

                    pMME = TrOp.E1.p[a_p,a_h] * Amp
                    isMME = 0.5 * TrOp.E1_C.p[a_p,a_h] * Amp
                    ivMME = TrOp.E1.p[a_p,a_h] * Amp * Float64(A - Z) / Float64(A)

                    pSum += pMME
                    isSum += isMME
                    ivSum += ivMME

                elseif t_ph == 1
                    a_p, a_h = Particle.n[p].a, Hole.n[h].a

                    nMME = TrOp.E1.n[a_p,a_h] * Amp
                    isMME = 0.5 * TrOp.E1_C.n[a_p,a_h] * Amp
                    ivMME = - TrOp.E1.n[a_p,a_h] * Amp * Float64(Z) / Float64(A)

                    nSum += nMME
                    isSum += isMME
                    ivSum += ivMME

                end

            end

            pM[nu_ind] = pSum
            nM[nu_ind] = nSum
            isM[nu_ind] = isSum
            ivM[nu_ind] = ivSum
        end

    end

    # Initialize the radial & angular grid ...
        # Case of the radial grid ...
    r_grid = range(r_min, stop = r_max, length = rN_grid)
    r_grid = collect(r_grid)
        # Case of the angular theta grid ...
    theta_grid = range(0.0, stop = 2.0*pi, length = tN_grid)
    theta_grid = collect(theta_grid)

    # Initialize the radial grid for single-particle orbitals ...
        # Case of wavefunctions ...
    pR_grid = Matrix{Float64}(undef,a_max,rN_grid)
    nR_grid = Matrix{Float64}(undef,a_max,rN_grid)
        # Case of wavefunction derivatives ...
    pdR_grid = Matrix{Float64}(undef,a_max,rN_grid)
    ndR_grid = Matrix{Float64}(undef,a_max,rN_grid)

    # Calculate the radial representation of the reference basis orbitals ...
    @inbounds Threads.@threads for i in 1:rN_grid
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
                nPsi = C.n[k,a] * Psi_rad_LHO(r,n_k,l_k,nu_neutron)
                pSum += pPsi
                nSum += nPsi
            end
            pR_grid[a,i] = pSum
            nR_grid[a,i] = nSum
        end
    end

    # Calculate the radial representation of the radial derivative of reference basis orbitals ...
    @inbounds Threads.@threads for a in 1:a_max
        pR_a = collect(@view pR_grid[a,:])
        nR_a = collect(@view nR_grid[a,:])
        pdR_a = differentiate_Oh6(r_grid,pR_a)
        ndR_a = differentiate_Oh6(r_grid,nR_a)
        @views pdR_grid[a,:] .= pdR_a
        @views ndR_grid[a,:] .= ndR_a
    end

    # Prepare the current map grid ...
        # Grid limits & step size ...
    x_min, x_max, xN_grid = -0.8 * r_max, 0.8 * r_max, 40
    z_min, z_max, zN_grid = -0.8 * r_max, 0.8 * r_max, 40
    dx, dz = (x_max - x_min) / Float64(xN_grid), (z_max - z_min) / Float64(zN_grid)
        # Grid size ...
    N_map = xN_grid * zN_grid
        # Map coordinate & index grid ...
    i_map = Matrix{Int64}(undef,2,N_map)
    x_map = Matrix{Float64}(undef,2,N_map)
        # Single state current maps ...
    pj_map = Array{Float64}(undef,3,N_ph,N_map)
    nj_map = Array{Float64}(undef,3,N_ph,N_map)
    pJ_map = Array{Float64}(undef,3,N_phonon,N_map)
    nJ_map = Array{Float64}(undef,3,N_phonon,N_map)

    # Prepare the current map grid ...
        # Precalculate radial & angular grid step size ...
    dr = (r_max - r_min) / (rN_grid - 1)
    dtheta = 2.0 * pi / (tN_grid - 1)  # if theta_grid spans 0..2*pi inclusive
        # Reset the grid size counter ...
    N_map = 0
        # Prepare the cartesian-like grid from initial radial x angular grid ...
    @inbounds for ind_x in 1:xN_grid
        x = x_min + Float64(ind_x) * dx
        @inbounds for ind_z in 1:zN_grid
            z = z_min + Float64(ind_z) * dz

            r = sqrt(x^2 + z^2)
            theta = atan(x,z)
            if theta < 0.0
                theta += 2.0*pi
            end

            n = Int64(round((r - r_min) / dr)) + 1
            m = Int64(round((theta - 0.0) / dtheta)) + 1

            n = clamp(n, 1, rN_grid)
            m = clamp(m, 1, tN_grid)

            N_map += 1
            i_map[1,N_map] = n
            i_map[2,N_map] = m
            x_map[1,N_map] = r_grid[n]
            x_map[2,N_map] = theta_grid[m]
        end
    end

    # Evaluate the grid of the 1-body transition current matrix elements ...
    @time @inbounds Threads.@threads for N in 1:N_map
        n, m = i_map[1,N], i_map[2,N]
        r, theta = r_grid[n], theta_grid[m]

        @inbounds for ph in 1:N_ph
            i_ph = Orb_Phonon[J+1,P][ph]
            p, h, t_ph = Phonon[i_ph].p, Phonon[i_ph].h, Phonon[i_ph].tz

            if t_ph == -1
                a_p, a_h = Particle.p[p].a, Hole.p[h].a
            elseif t_ph == 1
                a_p, a_h = Particle.n[p].a, Hole.n[h].a
            end

            l_p, j_p = Orb[a_p].l, Orb[a_p].j
            l_h, j_h = Orb[a_h].l, Orb[a_h].j
            j_x, j_y, j_z = 0.0, 0.0, 0.0

            if t_ph == -1
                R_p, R_h = pR_grid[a_p,n], pR_grid[a_h,n]
                dR_p, dR_h = pdR_grid[a_p,n], pdR_grid[a_h,n]
            elseif t_ph == 1
                R_p, R_h = nR_grid[a_p,n], nR_grid[a_h,n]
                dR_p, dR_h = ndR_grid[a_p,n], ndR_grid[a_h,n]
            end

            @inbounds for m_p in -j_p:2:j_p
                @inbounds for m_h in -j_h:2:j_h

                    if (m_p - m_h) != 2*M
                        continue
                    end

                    Amp = fCG(j_p,j_h,2*J,m_p,-m_h,2*M) * phase(j_h + m_h)

                    (jx,jy,jz) = J1b(r,theta,0.0,l_p,j_p,m_p,R_p,dR_p,l_h,j_h,-m_h,R_h,dR_h)

                    j_x += Amp * real(jx)
                    j_y += Amp * real(jy)
                    j_z += Amp * real(jz)

                end
            end

            if t_ph == -1
                pj_map[1,ph,N] = j_x
                pj_map[2,ph,N] = j_y
                pj_map[3,ph,N] = j_z
            elseif t_ph == 1
                nj_map[1,ph,N] = j_x
                nj_map[2,ph,N] = j_y
                nj_map[3,ph,N] = j_z
            end

        end

    end

    # Evaluate the grid of transition current matrix elements for given phonon levels ...
    @inbounds Threads.@threads for N in 1:N_map
        @inbounds for (nu_ind, nu) in enumerate(nu_list)
            pJxSum, pJySum, pJzSum = 0.0, 0.0, 0.0
            nJxSum, nJySum, nJzSum = 0.0, 0.0, 0.0

            @inbounds for ph in 1:N_ph
                i_ph = Orb_Phonon[J+1,P][ph]
                p, h, t_ph = Phonon[i_ph].p, Phonon[i_ph].h, Phonon[i_ph].tz

                if t_ph == -1
                    a_p, a_h = Particle.p[p].a, Hole.p[h].a
                elseif t_ph == 1
                    a_p, a_h = Particle.n[p].a, Hole.n[h].a
                end

                x_TDA = X_TDA[J+1,P][ph,nu]

                if t_ph == -1
                    pJxSum += x_TDA * pj_map[1,ph,N]
                    pJySum += x_TDA * pj_map[2,ph,N]
                    pJzSum += x_TDA * pj_map[3,ph,N]

                elseif t_ph == 1
                    nJxSum += x_TDA * nj_map[1,ph,N]
                    nJySum += x_TDA * nj_map[2,ph,N]
                    nJzSum += x_TDA * nj_map[3,ph,N]
                end

            end

            pJ_map[1,nu_ind,N] = pJxSum
            pJ_map[2,nu_ind,N] = pJySum
            pJ_map[3,nu_ind,N] = pJzSum
            nJ_map[1,nu_ind,N] = nJxSum
            nJ_map[2,nu_ind,N] = nJySum
            nJ_map[3,nu_ind,N] = nJzSum
        end
    
    end

    # Perform the export of of individual phonon level transition current maps ...
    @inbounds for (nu_ind, nu) in enumerate(nu_list)
        e_TDA, e_RPA = E_TDA[J+1,P][nu], real(E_RPA[J+1,P][nu])
        sE_TDA, sE_RPA = string(round(e_TDA,digits=3)), string(round(e_RPA,digits=3))

        # Export TDA transition densities ...
        Output_File_TDA = "IO/" * Params.Calc.Path * "/RPA/Densities/TDA_" * File_Name * "_JP$(JP)_nu$(nu)_E$(sE_TDA).dat"
        open(Output_File_TDA, "w") do Write_File
            @printf(Write_File, "%-20s %-20s %-20s %-20s %-20s %-20s %-20s %-20s %-20s %-20s\n", "x", "z", "pJ_x", "pJ_z", "nJ_x", "nJ_z", "isJ_x", "isJ_z", "ivJ_x", "ivJ_z")
            @inbounds for N in 1:N_map
                r, theta = x_map[1,N], x_map[2,N]
                x = r * sin(theta)
                z = r * cos(theta)
                pJ_x = pJ_map[1,nu_ind,N]
                pJ_z = pJ_map[3,nu_ind,N]
                nJ_x = nJ_map[1,nu_ind,N]
                nJ_z = nJ_map[3,nu_ind,N]
                isJ_x = 0.5 * (pJ_x + nJ_x)
                isJ_z = 0.5 * (pJ_z + nJ_z)
                ivJ_x = 0.5 * (pJ_x - nJ_x)
                ivJ_z = 0.5 * (pJ_z - nJ_z)
                @printf(Write_File, "%-20.8f %-20.8f %-20.8f %-20.8f %-20.8f %-20.8f %-20.8f %-20.8f %-20.8f %-20.8f\n", x, z, pJ_x, pJ_z, nJ_x, nJ_z, isJ_x, isJ_z, ivJ_x, ivJ_z)
            end
        end

        #=
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
        =#
    end

    display("Export of individual phonon level transition current maps successfuly completed ...")

    if ((J,P) in [(0,1),(1,2),(2,1),(3,2)]) == false
        return
    end

    # Perform the export of of averaged transition current maps ...
        # Export TDA averged transition densities ...
    Output_File_TDA_M = "IO/" * Params.Calc.Path * "/RPA/Densities/TDA_" * File_Name * "_JP$(JP)_Averaged.dat"
    open(Output_File_TDA_M, "w") do Write_File
        @printf(Write_File, "%-20s %-20s %-20s %-20s %-20s %-20s %-20s %-20s %-20s %-20s %-20s %-20s %-20s %-20s\n", "x", "z", "ppJ_x", "ppJ_z", "pnJ_x", "pnJ_z", "npJ_x", "npJ_z", "nnJ_x", "nnJ_z", "isJ_x", "isJ_z", "ivJ_x", "ivJ_z")
        @inbounds for N in 1:N_map
            r, theta = x_map[1,N], x_map[2,N]
            x = r * sin(theta)
            z = r * cos(theta)
            ppJ_x, pnJ_x, npJ_x, nnJ_x, isJ_x, ivJ_x = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
            ppJ_z, pnJ_z, npJ_z, nnJ_z, isJ_z, ivJ_z = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
            @inbounds for nu_ind in 1:N_phonon
                ppJ_x += pM[nu_ind] * pJ_map[1,nu_ind,N]
                ppJ_z += pM[nu_ind] * pJ_map[3,nu_ind,N]
                pnJ_x += pM[nu_ind] * nJ_map[1,nu_ind,N]
                pnJ_z += pM[nu_ind] * nJ_map[3,nu_ind,N]
                npJ_x += nM[nu_ind] * pJ_map[1,nu_ind,N]
                nnJ_x += nM[nu_ind] * nJ_map[1,nu_ind,N]
                npJ_z += nM[nu_ind] * pJ_map[3,nu_ind,N]
                nnJ_z += nM[nu_ind] * nJ_map[3,nu_ind,N]
                isJ_x += 0.5 * isM[nu_ind] * (pJ_map[1,nu_ind,N] + nJ_map[1,nu_ind,N])
                isJ_z += 0.5 * isM[nu_ind] * (pJ_map[3,nu_ind,N] + nJ_map[3,nu_ind,N])
                ivJ_x += 0.5 * ivM[nu_ind] * (pJ_map[1,nu_ind,N] - nJ_map[1,nu_ind,N])
                ivJ_z += 0.5 * ivM[nu_ind] * (pJ_map[3,nu_ind,N] - nJ_map[3,nu_ind,N])
            end
            # Perform the averaging ...
            @printf(Write_File, "%-20.8f %-20.8f %-20.8f %-20.8f %-20.8f %-20.8f %-20.8f %-20.8f %-20.8f %-20.8f %-20.8f %-20.8f %-20.8f %-20.8f\n", x, z, ppJ_x, ppJ_z, pnJ_x, pnJ_z, npJ_x, npJ_z, nnJ_x, nnJ_z, isJ_x, isJ_z, ivJ_x, ivJ_z)
        end
    end

    display("Export of averaged phonon level transition current maps successfuly completed ...")

    return
end
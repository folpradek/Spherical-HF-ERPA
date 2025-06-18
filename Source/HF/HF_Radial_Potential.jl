function HF_Radial_Potential(Params::Vector{Any},Orb::Vector{NOrb},pU::Matrix{Float64},nU::Matrix{Float64},Vp::Matrix{Float64},Vn::Matrix{Float64})
    # Basic constants ...
    HbarC = 197.326980
    m_n = 939.565346
    m_p = 938.272013

    # Read calculation parameters ...
    HbarOmega = Params[1]
    A = Params[5]
    Z = Params[6]
    N_max = Params[7]
    Output_File = Params[13]

    a_max = div((N_max+1)*(N_max+2),2)

    # Transform the potential from the LHO basis to the HF basis ...
    Vp, Vn = pU' * Vp * pU, nU' * Vn * nU

    println("\nPreparing radial mean-field potentials...")

    # Initialize grid ...
    r1 = 0.0 + 1e-8
    r2 = 2.5 * 1.2 * A^(1/3)
    N_Sampling = 10000
    r_grid = range(r1, stop=r2, length=N_Sampling)
    r_grid = collect(r_grid)

    Vp_rad, Vn_rad = Vector{Float64}(undef,N_Sampling), Vector{Float64}(undef,N_Sampling)

    nu_proton = 0.5 * m_p * HbarOmega / HbarC^2
    nu_neutron = 0.5 * m_n * HbarOmega / HbarC^2

    # Start radial HF potential calculation ...
    @inbounds Threads.@threads  for i in 1:N_Sampling
        r = r_grid[i]
        pSum = 0.0
        nSum = 0.0
        @inbounds for a in 1:a_max
            l_a = Orb[a].l
            j_a = Orb[a].j
            pRad = 0.0
            nRad = 0.0
            @inbounds for k in 1:a_max
                l_k = Orb[k].l
                j_k = Orb[k].j
                n_k = Orb[k].n
                if (j_a == j_k) && (l_a == l_k)
                    pPsi = Psi_rad_LHO(r,n_k,l_k,nu_proton) * pU[k,a]
                    nPsi = Psi_rad_LHO(r,n_k,l_k,nu_neutron) * nU[k,a]
                    pRad += pPsi
                    nRad += nPsi
                end
            end
            pSum += Vp[a,a] * Float64(Orb[a].pO) * (pRad)^2 * Float64(Orb[a].j + 1)
            nSum += Vn[a,a] * Float64(Orb[a].nO) * (nRad)^2 * Float64(Orb[a].j + 1)
        end
        Vp_rad[i] = pSum
        Vn_rad[i] = nSum
    end

    # Radial HF potential export ...
    Density_File_name = "IO/" * Output_File * "/HF/Densities/HF_Radial_Potential.dat"

    open(Density_File_name, "w") do io
        writedlm(io, hcat(r_grid, Vp_rad, Vn_rad), "\t")
    end

    println("\nMean-field potential exported...")
    
    return
end
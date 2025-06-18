function HF_ERPA_Radial_Density(Params::Vector{Any},Orb::Vector{NOrb},pRho::Matrix{Float64},nRho::Matrix{Float64})
    # Read calculation parameters ...
    N_max = Params[4]
    a_max = div((N_max + 1)*(N_max + 2),2)
    Output_File = Params[8]

    CMS = "NoCMS"
    if occursin("CMS1+2B", Output_File)
        CMS = "CMS1+2B"
    elseif occursin("CMS2B", Output_File)
        CMS = "CMS2B"
    end
    
    # Import HF basis transformation matrices ...
    pU_Import = Output_File * "/Bin/pU.bin"
    pU = Matrix{Float64}(undef,a_max,a_max)
    open(pU_Import, "r") do Read_File
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                ME = read(Read_File, Float64)
                pU[a,b] = ME
            end
        end
    end

    nU_Import = Output_File * "/Bin/nU.bin"
    nU = Matrix{Float64}(undef,a_max,a_max)
    open(nU_Import, "r") do Read_File
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                ME = read(Read_File, Float64)
                nU[a,b] = ME
            end
        end
    end

    # Calculation of 1-body & 2-body CMS...
    if CMS == "CMS1+2B" || CMS == "CMS2B"
        R_CMS = HF_ERPA_Radial_CMS(Params,Orb,pU,nU,pRho,nRho)
    else
        R_CMS = [0.0, 0.0, 0.0, 0.0]
    end

    # Evaluate HF densities on a grid ...
    Rho_Grid = HF_ERPA_Radial_Density_Grid(Params,Orb,pU,nU,pRho,nRho)

    # Evaluate HF charged density on a grid & calculte anomalous magnetic moment correction ...
    Rho_Grid, kR2 = HF_ERPA_Radial_ChDensity_Grid(Params,Orb,pU,nU,pRho,nRho,R_CMS,Rho_Grid[1],Rho_Grid[2],Rho_Grid[3])

    # Calculate HF mean-field nuclear radii ...
    pR2, nR2, chR2 = HF_ERPA_Radial_Radii(Rho_Grid[1],Rho_Grid[2],Rho_Grid[3],Rho_Grid[4])

    # Export HF mean-field radial densities ...
    HF_ERPA_Radial_Export(Params,Rho_Grid)

    # Write HF mean-field summary for radii ...
    HF_ERPA_Radial_Summary(Params,pR2,nR2,chR2,R_CMS,kR2)

    return
end

function HF_ERPA_Radial_CMS(Params::Vector{Any},Orb::Vector{NOrb},pU::Matrix{Float64},nU::Matrix{Float64},pRho::Matrix{Float64},nRho::Matrix{Float64})
    # Load parameters ...
    A = Params[1]
    Z = Params[2]
    HbarOmega = Params[3]
    N_max = Params[4]
    a_max = div((N_max + 1)*(N_max + 2),2)

    Radial_r1_LHO = Radial_Moment_Matrix_LHO(1,HbarOmega,Orb)
    Radial_r2_LHO = Radial_Moment_Matrix_LHO(2,HbarOmega,Orb)

    pRadial_r1 = zeros(Float64,a_max,a_max)
    pRadial_r2 = zeros(Float64,a_max,a_max)
    nRadial_r1 = zeros(Float64,a_max,a_max)
    nRadial_r2 = zeros(Float64,a_max,a_max)

    # Transform Radial Moment matrices to the HF-basis ...
    @inbounds for a in 1:a_max
        l_a = Orb[a].l
        j_a = Orb[a].j
        @inbounds for b in 1:a_max
            l_b = Orb[b].l
            j_b = Orb[b].j
            pR1Sum = 0.0
            pR2Sum = 0.0
            nR1Sum = 0.0
            nR2Sum = 0.0
            @inbounds for k in 1:a_max
                l_k = Orb[k].l
                j_k = Orb[k].j
                if (l_k == l_a) && (j_k == j_a)
                    @inbounds for l in 1:a_max
                        l_l = Orb[l].l
                        j_l = Orb[l].j
                        if (l_l == l_b) && (j_l == j_b)
                            pME1 = pU[k,a] * pU[l,b] * Radial_r1_LHO[k,l]
                            pME2 = pU[k,a] * pU[l,b] * Radial_r2_LHO[k,l]
                            pR1Sum += pME1
                            pR2Sum += pME2

                            nME1 = nU[k,a] * nU[l,b] * Radial_r1_LHO[k,l]
                            nME2 = nU[k,a] * nU[l,b] * Radial_r2_LHO[k,l]
                            nR1Sum += nME1
                            nR2Sum += nME2
                        end
                    end
                end
            end
            pRadial_r1[a,b] = pR1Sum
            pRadial_r2[a,b] = pR2Sum
            nRadial_r1[a,b] = nR1Sum
            nRadial_r2[a,b] = nR2Sum
        end
    end

    # Calculate CMS radius corrections ...

    pR2_CMS1 = 0.0
    pR2_CMS2 = 0.0
    nR2_CMS1 = 0.0
    nR2_CMS2 = 0.0

    # Calculate 2-body centre of mass to proton & neutron radii ...
    @inbounds for a in 1:a_max
        l_a = Orb[a].l
        j_a = Orb[a].j
        @inbounds for b in 1:a_max
            l_b = Orb[b].l
            j_b = Orb[b].j
            Amp = fCG(j_a,j_b,2,1,-1,0) * fCG(j_b,j_a,2,1,-1,0)
            pME = Amp * ((pRho[a,a] * pRho[b,b]  * abs(pRadial_r1[a,b])^2) * (2.0 / Float64(Z*A) - 1.0 / Float64(A^2)) - 1.0 / Float64(A^2) * nRho[a,a] * nRho[b,b] * abs(nRadial_r1[a,b])^2)
            nME = Amp * ((nRho[a,a] * nRho[b,b]  * abs(nRadial_r1[a,b])^2) * (2.0 / Float64(Z*A) - 1.0 / Float64(A^2)) - 1.0 / Float64(A^2) * pRho[a,a] * pRho[b,b] * abs(pRadial_r1[a,b])^2)
            pR2_CMS2 += pME
            nR2_CMS2 += nME
        end
        pME = ((1.0 / Float64(A^2) - 2.0 / Float64(Z*A)) * pRho[a,a] * pRadial_r2[a,a] + nRho[a,a] * nRadial_r2[a,a] / Float64(A^2)) * (Float64(j_a) + 1.0)
        nME = ((1.0 / Float64(A^2) - 2.0 / Float64((A-Z)*A)) * nRho[a,a] * nRadial_r2[a,a] + pRho[a,a] * pRadial_r2[a,a] / Float64(A^2)) * (Float64(j_a) + 1.0)
        pR2_CMS1 += pME
        nR2_CMS1 += nME
    end

    return [pR2_CMS1, pR2_CMS2, nR2_CMS1 ,nR2_CMS2]
end

function HF_ERPA_Radial_Density_Grid(Params::Vector{Any},Orb::Vector{NOrb},pU::Matrix{Float64},nU::Matrix{Float64},pRho::Matrix{Float64},nRho::Matrix{Float64})
    # Basic constants ...
    HbarC = 197.326980
    m_n = 939.565346
    m_p = 938.272013

    # Load parameters ...
    A = Params[1]
    HbarOmega = Params[3]
    N_max = Params[4]
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Preallocate grid ...
    r1 = 0.0 + 1e-8
    r2 = 2.5 * 1.2 * A^(1/3)
    N_Sampling = 10000
    r_grid = range(r1, stop = r2, length = N_Sampling)
    r_grid = collect(r_grid)

    pRho_rad = zeros(Float64,N_Sampling)
    nRho_rad = zeros(Float64,N_Sampling)

    nu_proton = 0.5 * m_p * HbarOmega / HbarC^2
    nu_neutron = 0.5 * m_n * HbarOmega / HbarC^2

    # Start grid calculation of density ...
    @inbounds Threads.@threads for i in 1:N_Sampling
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
            pSum += pRho[a,a] * pRad^2 * Float64(j_a + 1)
            nSum += nRho[a,a] * nRad^2 * Float64(j_a + 1)
        end
        pRho_rad[i] = pSum / (4.0 * π)
        nRho_rad[i] = nSum / (4.0 * π)
    end
    Rho_Grid = [r_grid, pRho_rad, nRho_rad]
    return Rho_Grid
end

function HF_ERPA_Radial_ChDensity_Grid(Params::Vector{Any},Orb::Vector{NOrb},pU::Matrix{Float64},nU::Matrix{Float64},pRho::Matrix{Float64},nRho::Matrix{Float64},R_CMS::Vector{Float64},r_grid::Vector{Float64},pRho_rad::Vector{Float64},nRho_rad::Vector{Float64})
    # Basic constants ...
    HbarC = 197.326980
    m_n = 939.565346
    m_p = 938.272013
    k_p = 1.793
    k_n = -1.913

    # Load parameters ...
    A = Params[1]
    Z = Params[2]
    HbarOmega = Params[3]
    N_max = Params[4]
    a_max = div((N_max + 1)*(N_max + 2),2)

    # CMS radius corrections ...
    pR_CMS = R_CMS[1] + R_CMS[2]
    nR_CMS = R_CMS[3] + R_CMS[4]

    # Preallocate charged densitiy grid ...
    N_Sampling = 10000
    chRho_rad = zeros(Float64,N_Sampling)
    kR2 = 0.0

    nu_proton = 0.5 * m_p * HbarOmega / HbarC^2
    nu_neutron = 0.5 * m_n * HbarOmega / HbarC^2
    
    # Evaluate the charged density ...
    @inbounds Threads.@threads for i in 1:N_Sampling
        r = r_grid[i]
        chSum = 0.0
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
            # Spin-Orbit interaction ...
            SO = Float64(l_a * KroneckerDelta(j_a,2*l_a+1) - (l_a+1) * KroneckerDelta(j_a,2*l_a-1))
            chME = (k_n /(m_n)^2 * nRho[a,a] * nRad^2 + k_p /(m_p)^2 * pRho[a,a] * pRad^2) *
                    (HbarC)^2 * SO * (Float64(j_a) + 1.0) / Float64(Z)
            chSum += chME
        end
        chME = Fold_pRho(r,r_grid,pRho_rad,pR_CMS)
        chSum += chME
        chME = Fold_nRho(r,r_grid,nRho_rad,nR_CMS)
        chSum += chME
        chRho_rad[i] += chSum
    end

    # Separately we also evaluate the anomalous magnetic moment contribution to the charged radius ...
    @inbounds for a in 1:a_max
        l_a = Orb[a].l
        j_a = Orb[a].j
        SO = Float64(l_a * KroneckerDelta(j_a,2*l_a+1) - (l_a+1) * KroneckerDelta(j_a,2*l_a-1))
        kME = (HbarC)^2 * (k_p / m_p^2 * pRho[a,a] + k_n / m_n^2 * nRho[a,a]) * (Float64(j_a) + 1.0)^2 * SO / (4.0 * π * Float64(Z))
        kR2 += kME
    end

    Rho_Grid = [r_grid, pRho_rad, nRho_rad, chRho_rad]

    return Rho_Grid, kR2
end

function HF_ERPA_Radial_Radii(r_grid::Vector{Float64},pRho_rad::Vector{Float64},nRho_rad::Vector{Float64},chRho_rad::Vector{Float64})
    pR2 = Integrate_Trap(r_grid, r_grid.^4 .* pRho_rad) / Integrate_Trap(r_grid, r_grid.^2 .* pRho_rad)
    nR2 = Integrate_Trap(r_grid, r_grid.^4 .* nRho_rad) / Integrate_Trap(r_grid, r_grid.^2 .* nRho_rad)
    chR2 = Integrate_Trap(r_grid, r_grid.^4 .* chRho_rad) / Integrate_Trap(r_grid, r_grid.^2 .* chRho_rad)
    return pR2, nR2, chR2
end

function HF_ERPA_Radial_Export(Params::Vector{Any},Rho_Grid::Vector{Vector{Float64}})
    # Read parameters ...
    Output_File = Params[8]

    # Export r_grid radius & pRho, nRho, chRho HF mean-field densities ...
    Density_File_name = Output_File * "/ERPA/Densities/HF_ERPA_Radial_Densities.dat"
    open(Density_File_name, "w") do Export_File
        writedlm(Export_File, hcat(Rho_Grid[1], Rho_Grid[2], Rho_Grid[3], Rho_Grid[4]), "\t")
    end
    return
end

function HF_ERPA_Radial_Summary(Params::Vector{Any},pR2::Float64,nR2::Float64,chR2::Float64,R_CMS::Vector{Float64},kR2::Float64)
    # Read parameters ...
    Orthogon = Params[5]
    Output_Path = Params[8]

    if Orthogon == true
        Output_File = Output_Path * "/ERPA/HF_ERPA_Summary_Ortho.dat"
    else
        Output_File = Output_Path * "/ERPA/HF_ERPA_Summary_Spur.dat"
    end

    println("\nResulting point-proton, charged & point-neutron radii ...")
    println("\nr_p  = " * string(round(sqrt(pR2 + R_CMS[1] + R_CMS[2] + 0.5 * (197.326980 / 938.272013)^2), sigdigits=6)) * " fm \t\t ... \t ERPA ground state point proton radius")
    println("r_ch = " * string(round(sqrt(chR2 + R_CMS[1] + R_CMS[2] + 0.5 * (197.326980 / 938.272013)^2), sigdigits=6)) * " fm \t\t ... \t ERPA ground state charged radius")
    println("r_n  = " * string(round(sqrt(nR2 + R_CMS[3] + R_CMS[4] + 0.5 * (197.326980 / 939.565346)^2), sigdigits=6)) * " fm \t\t ... \t ERPA ground state point neutron radius")
    println("\nFor more details see HF_EPRA_Summary.dat file ...")

    # Export resulsts & contributions to radii ...
    println("\nExporting information on HF_ERPA radii ...")

    Summary_File =  open(Output_File, "a")
        println(Summary_File, "\nInformation on HF-ERPA radii:")
        println(Summary_File, "\nr_p  = " * string(round(sqrt(pR2 + R_CMS[1] + R_CMS[2] + 0.5 * (197.326980 / 938.272013)^2), sigdigits=6)) * " fm \t\t ... \t ERPA ground state point proton radius")
        println(Summary_File, "r_ch = " * string(round(sqrt(chR2 + R_CMS[1] + R_CMS[2] + 0.5 * (197.326980 / 938.272013)^2), sigdigits=6)) * " fm \t\t ... \t ERPA ground state charged radius")
        println(Summary_File, "r_n  = " * string(round(sqrt(nR2 + R_CMS[3] + R_CMS[4] + 0.5 * (197.326980 / 939.565346)^2), sigdigits=6)) * " fm \t\t ... \t ERPA ground state point neutron radius")
        println(Summary_File, "\nr_p  = " * string(round(sqrt(pR2), sigdigits=6)) * "\t fm \t\t ... \t Uncorrected point-proton radius")
        println(Summary_File, "r_ch = " * string(round(sqrt(chR2), sigdigits=6)) * "\t fm \t\t ... \t Uncorrected charged radius")
        println(Summary_File, "r_n  = " * string(round(sqrt(nR2), sigdigits=6)) * "\t fm \t\t ... \t Uncorrected point-neutron radius")
        println(Summary_File, "\nDetails on radial effects:\n")
        println(Summary_File, "r_p2_cms1^2 = " * string(round(R_CMS[1], digits=6)) * "\t fm^2")
        println(Summary_File, "r_p2_cms2^2 = " * string(round(R_CMS[2], digits=6)) * "\t fm^2")
        println(Summary_File, "r_n2_cms1^2 = " * string(round(R_CMS[3], digits=6)) * "\t fm^2")
        println(Summary_File, "r_n2_cms2^2 = " * string(round(R_CMS[4], digits=6)) * "\t fm^2")
        println(Summary_File, "r_ch2_k^2   = " * string(round(kR2, sigdigits=7)) * "\t fm^2")
        println(Summary_File, "r_pDF2^2   = " * string(round(0.5 * (197.326980 / 938.272013)^2, sigdigits=7)) * "\t fm^2")
        println(Summary_File, "r_nDF2^2   = " * string(round(0.5 * (197.326980 / 939.565346)^2, sigdigits=7)) * "\t fm^2")
    close(Summary_File)

    println("\nHF-ERPA radial densities and radii exported ...")
end
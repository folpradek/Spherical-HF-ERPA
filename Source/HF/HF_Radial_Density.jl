function HF_Radial_Density(Params::Vector{Any},Orb::Vector{NOrb},pU::Matrix{Float64},nU::Matrix{Float64})
    CMS = Params[10]

    # Calculation of 1-body & 2-body CMS...
    if CMS == "CMS1+2B" || CMS == "CMS2B"
        R_CMS = HF_Radial_CMS(Params,Orb,pU,nU)
    else
        R_CMS = [0.0, 0.0, 0.0, 0.0]
    end

    # Evaluate HF densities on a grid ...
    Rho_Grid = HF_Radial_Density_Grid(Params,Orb,pU,nU)

    # Evaluate HF charged density on a grid & calculte anomalous magnetic moment correction ...
    Rho_Grid, kR2 = HF_Radial_ChDensity_Grid(Params,Orb,pU,nU,R_CMS,Rho_Grid[1],Rho_Grid[2],Rho_Grid[3])

    # Calculate HF mean-field nuclear radii ...
    pR2, nR2, chR2 = HF_Radial_Radii(Rho_Grid[1],Rho_Grid[2],Rho_Grid[3],Rho_Grid[4])

    # Export HF mean-field radial densities ...
    HF_Radial_Export(Params,Rho_Grid)

    # Write HF mean-field summary for radii ...
    HF_Radial_Summary(Params,pR2,nR2,chR2,R_CMS,kR2)

    return
end

function HF_Radial_CMS(Params::Vector{Any},Orb::Vector{NOrb},pU::Matrix{Float64},nU::Matrix{Float64})
    # Load parameters ...
    HbarOmega = Params[1]
    A = Params[5]
    Z = Params[6]
    N_max = Params[7]
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Initilize arrays ...
    pRho, nRho = HF_Density_Operator(Orb,a_max,pU,nU)
    pRho = pU' * pRho * pU
    nRho = nU' * nRho * nU

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

function HF_Radial_Density_Grid(Params::Vector{Any},Orb::Vector{NOrb},pU::Matrix{Float64},nU::Matrix{Float64})
    # Basic constants ...
    HbarC = 197.326980
    m_n = 939.565346
    m_p = 938.272013

    # Load parameters ...
    HbarOmega = Params[1]
    A = Params[5]
    N_max = Params[7]
    a_max = div((N_max + 1)*(N_max + 2),2)

    pRho, nRho = HF_Density_Operator(Orb,a_max,pU,nU)

    # Transform densities to the HF basis ...
    pRho = pU' * pRho * pU
    nRho = nU' * nRho * nU

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
            pSum += pRho[a,a] * pRad^2 * (Float64(j_a) + 1.0)
            nSum += nRho[a,a] * nRad^2 * (Float64(j_a) + 1.0)
        end
        pRho_rad[i] = pSum / (4.0 * π)
        nRho_rad[i] = nSum / (4.0 * π)
    end
    Rho_Grid = [r_grid, pRho_rad, nRho_rad]
    return Rho_Grid
end

@inline function Fold_pRho(r::Float64,r_grid::Vector{Float64},pRho_rad::Vector{Float64},pR_CMS::Float64)
    N_Sampling = 10000
    Rho = Vector{Float64}(undef,N_Sampling)
    @inbounds for i in 1:N_Sampling
        Rho[i] = pConvolution(r_grid[i],r,pR_CMS) * pRho_rad[i]
    end
    chRho = Integrate_Trap(r_grid,Rho)
    return chRho
end

@inline function Fold_nRho(r::Float64,r_grid::Vector{Float64},nRho_rad::Vector{Float64},nR_CMS::Float64)
    N_Sampling = 10000
    Rho = Vector{Float64}(undef,N_Sampling)
    @inbounds for i in 1:N_Sampling
        Rho[i] = nConvolution(r_grid[i],r,nR_CMS) * nRho_rad[i]
    end
    chRho = Integrate_Trap(r_grid,Rho)
    return chRho
end

@inline function pConvolution(x::Float64,r::Float64,pR_CMS::Float64)
    a_1 = 0.506373
    a_2 = 0.327922
    a_3 = 0.165705
    r_1 = sqrt(abs(0.431566 + pR_CMS + 0.5 * (197.326980 / 938.272013)^2))
    r_2 = sqrt(abs(0.139140 + pR_CMS + 0.5 * (197.326980 / 938.272013)^2))
    r_3 = sqrt(abs(1.525540 + pR_CMS + 0.5 * (197.326980 / 938.272013)^2))
    rho = x / (r * sqrt(pi)) * (a_1 / r_1 * (exp(-((r-x)/r_1)^2) - exp(-((r+x)/r_1)^2)) +
                                a_2 / r_2 * (exp(-((r-x)/r_2)^2) - exp(-((r+x)/r_2)^2)) +
                                a_3 / r_3 * (exp(-((r-x)/r_3)^2) - exp(-((r+x)/r_3)^2)))
    return rho
end

@inline function nConvolution(x::Float64,r::Float64,nR_CMS::Float64)
    r_p = sqrt(abs(0.4828 - 0.038664 + nR_CMS + 0.5 * (197.326980 / 939.565346)^2))
    r_m = sqrt(abs(0.4828 + 0.038664 + nR_CMS + 0.5 * (197.326980 / 939.565346)^2))
    rho = x / (r * sqrt(pi)) * ((exp(-((r-x)/r_p)^2) - exp(-((r+x)/r_p)^2)) / r_p -
                                (exp(-((r-x)/r_m)^2) - exp(-((r+x)/r_m)^2)) / r_m)
    return rho
end

function HF_Radial_ChDensity_Grid(Params::Vector{Any},Orb::Vector{NOrb},pU::Matrix{Float64},nU::Matrix{Float64},R_CMS::Vector{Float64},r_grid::Vector{Float64},pRho_rad::Vector{Float64},nRho_rad::Vector{Float64})
    # Basic constants ...
    HbarC = 197.326980
    m_n = 939.565346
    m_p = 938.272013
    k_p = 1.793
    k_n = -1.913

    # Load parameters ...
    HbarOmega = Params[1]
    A = Params[5]
    Z = Params[6]
    N_max = Params[7]
    a_max = div((N_max + 1)*(N_max + 2),2)

    pRho, nRho = HF_Density_Operator(Orb,a_max,pU,nU)

    # Transform densities to the HF basis ...
    pRho = pU' * pRho * pU
    nRho = nU' * nRho * nU

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
            # Spin-Orbit interaction matrix element ...
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

function HF_Radial_Radii(r_grid::Vector{Float64},pRho_rad::Vector{Float64},nRho_rad::Vector{Float64},chRho_rad::Vector{Float64})
    pR2 = Integrate_Trap(r_grid, r_grid.^4 .* pRho_rad) / Integrate_Trap(r_grid, r_grid.^2 .* pRho_rad)
    nR2 = Integrate_Trap(r_grid, r_grid.^4 .* nRho_rad) / Integrate_Trap(r_grid, r_grid.^2 .* nRho_rad)
    chR2 = Integrate_Trap(r_grid, r_grid.^4 .* chRho_rad) / Integrate_Trap(r_grid, r_grid.^2 .* chRho_rad)
    return pR2, nR2, chR2
end

function HF_Radial_Export(Params::Vector{Any},Rho_Grid::Vector{Vector{Float64}})
    # Read parameters ...
    Output_File = Params[13]

    # Export r_grid radius & pRho, nRho, chRho HF mean-field densities ...
    Density_File_name = "IO/" * Output_File * "/HF/Densities/HF_Radial_Densities.dat"
    open(Density_File_name, "w") do Export_File
        writedlm(Export_File, hcat(Rho_Grid[1], Rho_Grid[2], Rho_Grid[3], Rho_Grid[4]), "\t")
    end
    return
end

function HF_Radial_Summary(Params::Vector{Any},pR2::Float64,nR2::Float64,chR2::Float64,R_CMS::Vector{Float64},kR2::Float64)
    # Read parameters ...
    Output_File = Params[13]

    println("\nResulting point-proton, charged & point-neutron radii ...")
    println("\nr_p  = " * string(round(sqrt(pR2 + R_CMS[1] + R_CMS[2] + 0.5 * (197.326980 / 938.272013)^2), sigdigits=6)) * " fm \t\t ... \t Mean-field ground state point proton radius")
    println("r_ch = " * string(round(sqrt(chR2 + R_CMS[1] + R_CMS[2] + 0.5 * (197.326980 / 938.272013)^2), sigdigits=6)) * " fm \t\t ... \t Mean-field ground state charged radius")
    println("r_n  = " * string(round(sqrt(nR2 + R_CMS[3] + R_CMS[4] + 0.5 * (197.326980 / 939.565346)^2), sigdigits=6)) * " fm \t\t ... \t Mean-field ground state point neutron radius")
    println("\nFor more details see HF_Summary.dat file ...")

    # Export resulsts & contributions to radii ...
    println("\nExporting information on HF mean-field radii ...")
    Summary_File =  open(string("IO/", Output_File, "/HF/HF_Summary.dat"), "a")
        println(Summary_File, "\nr_p  = " * string(round(sqrt(pR2 + R_CMS[1] + R_CMS[2] + 0.5 * (197.326980 / 938.272013)^2), sigdigits=6)) * " fm \t\t ... \t Mean-field ground state point proton radius")
        println(Summary_File, "r_ch = " * string(round(sqrt(chR2 + R_CMS[1] + R_CMS[2] + 0.5 * (197.326980 / 938.272013)^2), sigdigits=6)) * " fm \t\t ... \t Mean-field ground state charged radius")
        println(Summary_File, "r_n  = " * string(round(sqrt(nR2 + R_CMS[3] + R_CMS[4] + 0.5 * (197.326980 / 939.565346)^2), sigdigits=6)) * " fm \t\t ... \t Mean-field ground state point neutron radius")
        println(Summary_File, "\nr_p  = " * string(round(sqrt(pR2), sigdigits=6)) * "\t fm \t\t ... \t Uncorrected point-proton radius")
        println(Summary_File, "r_ch = " * string(round(sqrt(chR2), sigdigits=6)) * "\t fm \t\t ... \t Uncorrected charged radius")
        println(Summary_File, "r_n  = " * string(round(sqrt(nR2), sigdigits=6)) * "\t fm \t\t ... \t Uncorrected point-neutron radius")
        println(Summary_File, "\nDetails on radii corrections:\n")
        println(Summary_File, "r_p2_cms1^2 = " * string(round(R_CMS[1], digits=6)) * "\t fm^2")
        println(Summary_File, "r_p2_cms2^2 = " * string(round(R_CMS[2], digits=6)) * "\t fm^2")
        println(Summary_File, "r_n2_cms1^2 = " * string(round(R_CMS[3], digits=6)) * "\t fm^2")
        println(Summary_File, "r_n2_cms2^2 = " * string(round(R_CMS[4], digits=6)) * "\t fm^2")
        println(Summary_File, "r_ch2_k^2   = " * string(round(kR2, sigdigits=7)) * "\t fm^2")
        println(Summary_File, "r_pDF2^2   = " * string(round(0.5 * (197.326980 / 938.272013)^2, sigdigits=7)) * "\t fm^2")
        println(Summary_File, "r_nDF2^2   = " * string(round(0.5 * (197.326980 / 939.565346)^2, sigdigits=7)) * "\t fm^2")
    close(Summary_File)

    println("\nMean-field radial densities and radii exported ...")
end
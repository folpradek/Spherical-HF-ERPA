function HFMBPT_Radial_Density(Params::Vector{Any},Orb::Vector{NOrb},dpRho::Matrix{Float64},dnRho::Matrix{Float64})
    # Read parameters ...
    N_max = Params[7]
    a_max = div((N_max + 1)*(N_max + 2),2)
    CMS = Params[10]
    Output_File = Params[13]

    # Import HF basis transformation matrices ...
    pU_Import = "IO/" * Output_File * "/Bin/pU.bin"
    pU = Matrix{Float64}(undef,a_max,a_max)
    open(pU_Import, "r") do Read_File
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                ME = read(Read_File, Float64)
                pU[a,b] = ME
            end
        end
    end

    nU_Import = "IO/" * Output_File * "/Bin/nU.bin"
    nU = Matrix{Float64}(undef,a_max,a_max)
    open(nU_Import, "r") do Read_File
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                ME = read(Read_File, Float64)
                nU[a,b] = ME
            end
        end
    end

    println("\nPreparing MBPT(3) corrections to densities & radii ...")

    # Make HF density operator ...
    pRho, nRho = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)

    @inbounds for a in 1:a_max
        pRho[a,a] = Orb[a].pO
        nRho[a,a] = Orb[a].nO
    end

    pRho_HF, nRho_HF = deepcopy(pRho), deepcopy(nRho)

    # Evaluate point-proton & point-neutron radii - non-perturbative HF ...
    Rho_Grid = HF_Radial_Density_Grid(Params,Orb,pU,nU)
    pR2_HF, nR2_HF = HFMBPT_NonPT_Radial_Radii(Rho_Grid[1],Rho_Grid[2],Rho_Grid[3])

    # Include MBPT(3) density perturbation ...
    pRho .= pRho .+ dpRho
    nRho .= nRho .+ dnRho

    # Calculation of 1-body & 2-body CMS & Darwin-Foldy corrections ...
    if CMS == "CMS1+2B" || CMS == "CMS2B"
        R_CMS = HFMBPT_Radial_CMS(Params,Orb,pU,nU,pRho,nRho)
    else
        R_CMS = [0.0, 0.0, 0.0, 0.0]
    end

    # HF-MBPT(2) determinant overlap ...
    HFMBPT_Overlap(Params,Orb,pRho,nRho,pRho_HF,nRho_HF)

    # Evaluate LO MBPT(3) corrections to radii ...
    pR2_MBPT, nR2_MBPT = HFMBPT_Radial_Correction(Params,Orb,pU,nU,dpRho,dnRho)

    # Evaluate radial densities on a grid ...
    Rho_Grid = HFMBPT_Radial_Density_Grid(Params,Orb,pU,nU,pRho,nRho)

    # Evaluate charged radial densitiy on a grid & anomalous magnetic momenent correction ...
    Rho_Grid, kR2 = HFMBPT_Radial_ChDensity_Grid(Params,Orb,pU,nU,pRho,nRho,R_CMS,Rho_Grid[1],Rho_Grid[2],Rho_Grid[3])

    # Evaluate MBPT(3) radial radii ...
    pR2, nR2, chR2 = HFMBPT_Radial_Radii(Rho_Grid[1],Rho_Grid[2],Rho_Grid[3],Rho_Grid[4])

    # Export MBPT(3) radial densities ...
    HFMBPT_Radial_Export(Params,Rho_Grid)

    # Write summary for MBPT(3) results on radii ...
    HFMBPT_Radial_Summary(Params,pR2,nR2,chR2,R_CMS,kR2,pR2_MBPT,nR2_MBPT,pR2_HF,nR2_HF)

    return
end

function HFMBPT_Radial_CMS(Params::Vector{Any},Orb::Vector{NOrb},pU::Matrix{Float64},nU::Matrix{Float64},pRho::Matrix{Float64},nRho::Matrix{Float64})
    # Load parameters ...
    HbarOmega = Params[1]
    A = Params[5]
    Z = Params[6]
    N_max = Params[7]
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Initilize arrays ...
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

    # Calculate 2-body centre of mass & MBPT(3) corrections to proton & neutron radii ...

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

function HFMBPT_Radial_Correction(Params::Vector{Any},Orb::Vector{NOrb},pU::Matrix{Float64},nU::Matrix{Float64},dpRho::Matrix{Float64},dnRho::Matrix{Float64})
    # Read parameters ...
    HbarOmega = Params[1]
    A = Params[5]
    Z = Params[6]
    N_max = Params[7]
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Initilize arrays ...
    Radial_r2_LHO = Radial_Moment_Matrix_LHO(2,HbarOmega,Orb)

    pRadial_r2 = zeros(Float64,a_max,a_max)
    nRadial_r2 = zeros(Float64,a_max,a_max)

    # Transform Radial Moment matrices to the HF-basis ...
    @inbounds for a in 1:a_max
        l_a = Orb[a].l
        j_a = Orb[a].j
        @inbounds for b in 1:a_max
            l_b = Orb[b].l
            j_b = Orb[b].j
            pR2Sum = 0.0
            nR2Sum = 0.0
            @inbounds for k in 1:a_max
                l_k = Orb[k].l
                j_k = Orb[k].j
                if (l_k == l_a) && (j_k == j_a)
                    @inbounds for l in 1:a_max
                        l_l = Orb[l].l
                        j_l = Orb[l].j
                        if (l_l == l_b) && (j_l == j_b)
                            pME2 = pU[k,a] * pU[l,b] * Radial_r2_LHO[k,l]
                            pR2Sum += pME2

                            nME2 = nU[k,a] * nU[l,b] * Radial_r2_LHO[k,l]
                            nR2Sum += nME2
                        end
                    end
                end
            end
            pRadial_r2[a,b] = pR2Sum
            nRadial_r2[a,b] = nR2Sum
        end
    end

    pR2_MBPT = 0.0
    nR2_MBPT = 0.0

    # Evaluate MBPT(3) LO corrections to radii ...
    @inbounds for a in 1:a_max
        l_a = Orb[a].l
        j_a = Orb[a].j
        @inbounds for b in 1:a_max
            l_b = Orb[b].l
            j_b = Orb[b].j
            if (j_a == j_b) && (l_a == l_b)
                pME = pRadial_r2[a,b] * dpRho[a,b] / Float64(Z) * (Float64(j_a) + 1.0)
                nME = nRadial_r2[a,b] * dnRho[a,b] / Float64(A-Z) * (Float64(j_a) + 1.0)
                pR2_MBPT += pME
                nR2_MBPT += nME
            end
        end
    end

    return pR2_MBPT, nR2_MBPT
end

function HFMBPT_Radial_Density_Grid(Params::Vector{Any},Orb::Vector{NOrb},pU::Matrix{Float64},nU::Matrix{Float64},pRho::Matrix{Float64},nRho::Matrix{Float64})
    # Basic constants ...
    HbarC = 197.326980
    m_n = 939.565346
    m_p = 938.272013

    # Load parameters ...
    HbarOmega = Params[1]
    A = Params[5]
    N_max = Params[7]
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
            @inbounds for b in 1:a_max
                l_b = Orb[b].l
                j_b = Orb[b].j
                if (j_a == j_b) && (l_a == l_b)
                    pRad = 0.0
                    nRad = 0.0
                    @inbounds for k in 1:a_max
                        l_k = Orb[k].l
                        j_k = Orb[k].j
                        n_k = Orb[k].n
                        if (j_a == j_k) && (l_a == l_k)
                            @inbounds for l in 1:a_max
                                l_l = Orb[l].l
                                j_l = Orb[l].j
                                n_l = Orb[l].n
                                if (j_b == j_l) && (l_b == l_l)
                                    pPsi = Psi_rad_LHO(r,n_k,l_k,nu_proton) * Psi_rad_LHO(r,n_l,l_l,nu_proton) * pU[k,a] * pU[l,b]
                                    nPsi = Psi_rad_LHO(r,n_k,l_k,nu_neutron) * Psi_rad_LHO(r,n_l,l_l,nu_neutron) * nU[k,a] * nU[l,b]
                                    pRad += pPsi
                                    nRad += nPsi
                                end
                            end
                        end
                    end
                    pSum += pRho[a,b] * pRad * (Float64(j_a) + 1.0)
                    nSum += nRho[a,b] * nRad * (Float64(j_a) + 1.0)
                end
            end
        end
        pRho_rad[i] = pSum / (4.0 * π)
        nRho_rad[i] = nSum / (4.0 * π)
    end
    Rho_Grid = [r_grid, pRho_rad, nRho_rad]
    return Rho_Grid
end

function HFMBPT_Radial_ChDensity_Grid(Params::Vector{Any},Orb::Vector{NOrb},pU::Matrix{Float64},nU::Matrix{Float64},pRho::Matrix{Float64},nRho::Matrix{Float64},R_CMS::Vector{Float64},r_grid::Vector{Float64},pRho_rad::Vector{Float64},nRho_rad::Vector{Float64})
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
            @inbounds for b in 1:a_max
                l_b = Orb[b].l
                j_b = Orb[b].j
                if j_a == j_b && l_a == l_b
                    pRad = 0.0
                    nRad = 0.0
                    @inbounds for k in 1:a_max
                        l_k = Orb[k].l
                        j_k = Orb[k].j
                        n_k = Orb[k].n
                        if (j_a == j_k) && (l_a == l_k)
                            @inbounds for l in 1:a_max
                                l_l = Orb[l].l
                                j_l = Orb[l].j
                                n_l = Orb[l].n
                                if (j_b == j_l) && (l_b == l_l)
                                    pPsi = Psi_rad_LHO(r,n_k,l_k,nu_proton) * Psi_rad_LHO(r,n_l,l_l,nu_proton) * pU[k,a] * pU[l,b]
                                    nPsi = Psi_rad_LHO(r,n_k,l_k,nu_neutron) * Psi_rad_LHO(r,n_l,l_l,nu_neutron) * nU[k,a] * nU[l,b]
                                    pRad += pPsi
                                    nRad += nPsi
                                end
                            end
                        end
                    end
                    # Spin-Orbit interaction matrix element ...
                    SO = Float64(l_a * KroneckerDelta(j_a,2*l_a+1) - (l_a+1) * KroneckerDelta(j_a,2*l_a-1))
                    chME = (k_n /(m_n)^2 * nRho[a,b] * nRad^2 + k_p /(m_p)^2 * pRho[a,b] * pRad^2) *
                            (HbarC)^2 * SO * (Float64(j_a) + 1.0) / Float64(Z)
                    chSum += chME
                end
            end 
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
        @inbounds for b in 1:a_max
            l_b = Orb[b].l
            j_b = Orb[b].j
            if j_a == j_b && l_a == l_b
                SO = Float64(l_a * KroneckerDelta(j_a,2*l_a+1) - (l_a+1) * KroneckerDelta(j_a,2*l_a-1))
                kME = (HbarC)^2 * (k_p / m_p^2 * pRho[a,b] + k_n / m_n^2 * nRho[a,b]) * (Float64(j_a) + 1.0)^2 * SO / (4.0 * π * Float64(Z))
                kR2 += kME
            end
        end
    end

    Rho_Grid = [r_grid, pRho_rad, nRho_rad, chRho_rad]

    return Rho_Grid, kR2
end

function HFMBPT_NonPT_Radial_Radii(r_grid::Vector{Float64},pRho_rad::Vector{Float64},nRho_rad::Vector{Float64})
    pR2 = Integrate_Trap(r_grid, r_grid.^4 .* pRho_rad) / Integrate_Trap(r_grid, r_grid.^2 .* pRho_rad)
    nR2 = Integrate_Trap(r_grid, r_grid.^4 .* nRho_rad) / Integrate_Trap(r_grid, r_grid.^2 .* nRho_rad)
    return pR2, nR2
end

function HFMBPT_Radial_Radii(r_grid::Vector{Float64},pRho_rad::Vector{Float64},nRho_rad::Vector{Float64},chRho_rad::Vector{Float64})
    pR2 = Integrate_Trap(r_grid, r_grid.^4 .* pRho_rad) / Integrate_Trap(r_grid, r_grid.^2 .* pRho_rad)
    nR2 = Integrate_Trap(r_grid, r_grid.^4 .* nRho_rad) / Integrate_Trap(r_grid, r_grid.^2 .* nRho_rad)
    chR2 = Integrate_Trap(r_grid, r_grid.^4 .* chRho_rad) / Integrate_Trap(r_grid, r_grid.^2 .* chRho_rad)
    return pR2, nR2, chR2
end

function HFMBPT_Radial_Export(Params::Vector{Any},Rho_Grid::Vector{Vector{Float64}})
    # Read parameters ...
    Output_File = Params[13]

    # Export r_grid radius & pRho, nRho, chRho HF-MBPT(3) densities ...
    Density_File_name = "IO/" * Output_File * "/HF/Densities/HFMBPT_Radial_Densities.dat"
    open(Density_File_name, "w") do Export_File
        writedlm(Export_File, hcat(Rho_Grid[1], Rho_Grid[2], Rho_Grid[3], Rho_Grid[4]), "\t")
    end
    return
end

function HFMBPT_Radial_Summary(Params::Vector{Any},pR2::Float64,nR2::Float64,chR2::Float64,R_CMS::Vector{Float64},kR2::Float64,pR2_MBPT::Float64,nR2_MBPT::Float64,pR2_HF::Float64,nR2_HF::Float64)
    # Read parameters ...
    Output_File = Params[13]

    println("\nResulting point-proton, charged & point-neutron radii from HF-MBPT(3) ...")
    println("\nr_p  = " * string(round(sqrt(pR2 + R_CMS[1] + R_CMS[2] + 0.5 * (197.326980 / 938.272013)^2), sigdigits=6)) * " fm \t\t ... \t Density calculated MBPT(2) proton radius")
    println("r_ch = " * string(round(sqrt(chR2 + R_CMS[1] + R_CMS[2] + 0.5 * (197.326980 / 938.272013)^2), sigdigits=6)) * " fm \t\t ... \t Density calculated MBPT(2) charged radius")
    println("r_n  = " * string(round(sqrt(nR2 + R_CMS[3] + R_CMS[4] + 0.5 * (197.326980 / 939.565346)^2), sigdigits=6)) * " fm \t\t ... \t Density calculated MBPT(2) neutron radius")
    println("\nr_p  = " * string(round(sqrt(pR2_HF + R_CMS[1] + R_CMS[2] + 0.5 * (197.326980 / 938.272013)^2 + pR2_MBPT), sigdigits=6)) * " fm \t\t ... \t Calculated MBPT(2) proton radius")
    println("r_n  = " * string(round(sqrt(nR2_HF + R_CMS[3] + R_CMS[4] + 0.5 * (197.326980 / 939.565346)^2 + nR2_MBPT), sigdigits=6)) * " fm \t\t ... \t Calculated MBPT(2) neutron radius")
    println("\nFor more details see HF_Summary.dat file ...")

    # Export resulsts & contributions to radii ...
    println("\nExporting information on HF-MBPT(2) radii ...")
    Summary_File =  open(string("IO/", Output_File, "/HF/HF_Summary.dat"), "a")
        println(Summary_File, "\nr_p  = " * string(round(sqrt(pR2 + R_CMS[1] + R_CMS[2] + 0.5 * (197.326980 / 938.272013)^2), sigdigits=6)) * " fm \t\t ... \t Density calculated MBPT(2) proton radius")
        println(Summary_File, "r_ch = " * string(round(sqrt(chR2 + R_CMS[1] + R_CMS[2] + 0.5 * (197.326980 / 938.272013)^2), sigdigits=6)) * " fm \t\t ... \t Density calculated MBPT(2) charged radius")
        println(Summary_File, "r_n  = " * string(round(sqrt(nR2 + R_CMS[3] + R_CMS[4] + 0.5 * (197.326980 / 939.565346)^2), sigdigits=6)) * " fm \t\t ... \t Density calculated MBPT(2) neutron radius")
        println(Summary_File, "\nr_p  = " * string(round(sqrt(pR2_HF + R_CMS[1] + R_CMS[2] + 0.5 * (197.326980 / 938.272013)^2 + pR2_MBPT), sigdigits=6)) * " fm \t\t ... \t Calculated MBPT(2) proton radius")
        println(Summary_File, "r_n  = " * string(round(sqrt(nR2_HF + R_CMS[3] + R_CMS[4] + 0.5 * (197.326980 / 939.565346)^2 + nR2_MBPT), sigdigits=6)) * " fm \t\t ... \t Calculated MBPT(2) neutron radius")
        println(Summary_File, "\nr_p  = " * string(round(sqrt(pR2), sigdigits=6)) * " fm \t\t ... \t Density calculated MBPT(2) uncorrected point-proton radius")
        println(Summary_File, "r_ch = " * string(round(sqrt(chR2), sigdigits=6)) * " fm \t\t ... \t Density calculated MBPT(2) uncorrected charged radius")
        println(Summary_File, "r_n  = " * string(round(sqrt(nR2), sigdigits=6)) * " fm \t\t ... \t Density calculated MBPT(2) uncorrected point-neutron radius")
        println(Summary_File, "\nr_p  = " * string(round(sqrt(pR2_HF + pR2_MBPT), sigdigits=6)) * " fm \t\t ... \t Calculated MBPT(2) uncorrected point-proton radius")
        println(Summary_File, "r_n  = " * string(round(sqrt(nR2_HF + nR2_MBPT), sigdigits=6)) * " fm \t\t ... \t Calculated MBPT(2) uncorrected point-neutron radius")
        println(Summary_File, "\nDetails on radii corrections:\n")
        println(Summary_File, "r_p2_cms1^2 = " * string(round(R_CMS[1], digits=6)) * "\t fm^2")
        println(Summary_File, "r_p2_cms2^2 = " * string(round(R_CMS[2], digits=6)) * "\t fm^2")
        println(Summary_File, "r_n2_cms1^2 = " * string(round(R_CMS[3], digits=6)) * "\t fm^2")
        println(Summary_File, "r_n2_cms2^2 = " * string(round(R_CMS[4], digits=6)) * "\t fm^2")
        println(Summary_File, "r_p2_mbpt^2   = " * string(round(pR2_MBPT, digits=6)) * "\t fm^2")
        println(Summary_File, "r_n2_mbpt^2   = " * string(round(nR2_MBPT, digits=6)) * "\t fm^2")
        println(Summary_File, "r_ch2_k^2   = " * string(round(kR2, sigdigits=7)) * "\t fm^2")
        println(Summary_File, "r_pDF2^2   = " * string(round(0.5 * (197.326980 / 938.272013)^2, sigdigits=7)) * "\t fm^2")
        println(Summary_File, "r_nDF2^2   = " * string(round(0.5 * (197.326980 / 939.565346)^2, sigdigits=7)) * "\t fm^2")
    close(Summary_File)

    println("\nHF-MBPT(2) radial densities and radii exported ...")
end

function HFMBPT_Overlap(Params::Vector{Any},Orb::Vector{NOrb},pRho::Matrix{Float64},nRho::Matrix{Float64},pRho_HF::Matrix{Float64},nRho_HF::Matrix{Float64})
    # Read parameters ...
    A = Params[5]
    Z = Params[6]
    N_max = Params[7]
    a_max = div((N_max + 1)*(N_max + 2),2)
    Output_File = Params[13]
    
    # Initialize traces ...
    Overlap = 0.0

    # Calculate overlap ...

    pProjection, nProjection = pRho * pRho_HF, nRho * nRho_HF

    @inbounds for a in 1:a_max
        Overlap += pProjection[a,a] * Float64(Orb[a].j + 1) + nProjection[a,a] * Float64(Orb[a].j + 1)
    end

    Overlap = Overlap / A

    println("\nHF-MBPT(2) slater determinant overlap with HF determinant = " * string(round(Overlap * 100.0, digits = 5)) * " %")

    println("\nExporting information on HF-MBPT(2) Slater determinant over with HF determinant ...")
    Summary_File =  open(string("IO/", Output_File, "/HF/HF_Summary.dat"), "a")
        println(Summary_File, "\nOverlap = " * string(round(Overlap, digits = 6)) * "\t ... Overlap of HF and HF-MBPT(2) determinants...")
    close(Summary_File)

    return
end
function HF_ERPA_rM(Params::Vector{Any},N_nu::Matrix{Int64},Orb_Phonon::Matrix{Vector{Int64}},Phonon::Vector{PhState},Particle::pnSVector,Hole::pnSVector,X_RPA::Matrix{Matrix{ComplexF64}},Y_RPA::Matrix{Matrix{ComplexF64}},TrOp::TranOper,pRho::Matrix{Float64},nRho::Matrix{Float64})
    # Read calculation parameters ...
    Orthogon = Params[5]

    println("\nCalculating reduced transition matrix elements rM ...")

    # E0
    J = 0
    P = 1
    N_ph = N_nu[J+1,P]
    prME0 = Vector{ComplexF64}(undef,N_ph)
    nrME0 = Vector{ComplexF64}(undef,N_ph)
    @inbounds for nu in 1:N_ph
        pME0Sum = ComplexF64(0.0)
        nME0Sum = ComplexF64(0.0)
        @inbounds for Ind_ph in 1:N_ph
            ph = Orb_Phonon[J+1,P][Ind_ph]
            p,h,t_ph,J_ph,P_ph,a_p,l_p,j_p,E_p,a_h,l_h,j_h,E_h = Phonon_Ind(ph,Phonon,Particle,Hole)

            if t_ph == -1
                ME0 = TrOp.E0.p[a_p,a_h] * (X_RPA[J+1,P][Ind_ph,nu] + Y_RPA[J+1,P][Ind_ph,nu]) * sqrt(pRho[a_h,a_h] - pRho[a_p,a_p])
                pME0Sum += ME0
            end

            if t_ph == 1
                ME0 = TrOp.E0.n[a_p,a_h] * (X_RPA[J+1,P][Ind_ph,nu] + Y_RPA[J+1,P][Ind_ph,nu]) * sqrt(nRho[a_h,a_h] - nRho[a_p,a_p])
                nME0Sum += ME0
            end

        end
        prME0[nu] = pME0Sum
        nrME0[nu] = nME0Sum
    end

    # E1
    J = 1
    P = 2
    N_ph = N_nu[J+1,P]
    prME1 = Vector{ComplexF64}(undef,N_ph)
    nrME1 = Vector{ComplexF64}(undef,N_ph)
    @inbounds for nu in 1:N_ph
        if (Orthogon == true && nu != 1) || (Orthogon == false)
            pME1Sum = ComplexF64(0.0)
            nME1Sum = ComplexF64(0.0)
            @inbounds for Ind_ph in 1:N_ph
                ph = Orb_Phonon[J+1,P][Ind_ph]
                p,h,t_ph,J_ph,P_ph,a_p,l_p,j_p,E_p,a_h,l_h,j_h,E_h = Phonon_Ind(ph,Phonon,Particle,Hole)

                if t_ph == -1
                    ME1 = TrOp.E1.p[a_p,a_h] * (-1.0 * X_RPA[J+1,P][Ind_ph,nu] + Y_RPA[J+1,P][Ind_ph,nu]) * sqrt(pRho[a_h,a_h] - pRho[a_p,a_p])
                    pME1Sum += ME1
                end

                if t_ph == 1
                    ME1 = TrOp.E1.n[a_p,a_h] * (-1.0 * X_RPA[J+1,P][Ind_ph,nu] + Y_RPA[J+1,P][Ind_ph,nu]) * sqrt(nRho[a_h,a_h] - nRho[a_p,a_p])
                    nME1Sum += ME1
                end

            end
            prME1[nu] = pME1Sum
            nrME1[nu] = nME1Sum
        end
    end

    # E2
    J = 2
    P = 1
    N_ph = N_nu[J+1,P]
    prME2 = Vector{ComplexF64}(undef,N_ph)
    nrME2 = Vector{ComplexF64}(undef,N_ph)
    @inbounds for nu in 1:N_ph
        pME2Sum = ComplexF64(0.0)
        nME2Sum = ComplexF64(0.0)
        @inbounds for Ind_ph in 1:N_ph
            ph = Orb_Phonon[J+1,P][Ind_ph]
            p,h,t_ph,J_ph,P_ph,a_p,l_p,j_p,E_p,a_h,l_h,j_h,E_h = Phonon_Ind(ph,Phonon,Particle,Hole)

            if t_ph == -1
                ME2 = TrOp.E2.p[a_p,a_h] * (X_RPA[J+1,P][Ind_ph,nu] + Y_RPA[J+1,P][Ind_ph,nu]) * sqrt(pRho[a_h,a_h] - pRho[a_p,a_p])
                pME2Sum += ME2
            end

            if t_ph == 1
                ME2 = TrOp.E2.n[a_p,a_h] * (X_RPA[J+1,P][Ind_ph,nu] + Y_RPA[J+1,P][Ind_ph,nu]) * sqrt(nRho[a_h,a_h] - nRho[a_p,a_p])
                nME2Sum += ME2
            end

        end
        prME2[nu] = pME2Sum
        nrME2[nu] = nME2Sum
    end

    # E3
    J = 3
    P = 2
    N_ph = N_nu[J+1,P]
    prME3 = Vector{ComplexF64}(undef,N_ph)
    nrME3 = Vector{ComplexF64}(undef,N_ph)
    @inbounds for nu in 1:N_ph
        pME3Sum = ComplexF64(0.0)
        nME3Sum = ComplexF64(0.0)
        @inbounds for Ind_ph in 1:N_ph
            ph = Orb_Phonon[J+1,P][Ind_ph]
            p,h,t_ph,J_ph,P_ph,a_p,l_p,j_p,E_p,a_h,l_h,j_h,E_h = Phonon_Ind(ph,Phonon,Particle,Hole)

            if t_ph == -1
                ME3 = TrOp.E3.p[a_p,a_h] * (-1.0 * X_RPA[J+1,P][Ind_ph,nu] + Y_RPA[J+1,P][Ind_ph,nu]) * sqrt(pRho[a_h,a_h] - pRho[a_p,a_p])
                pME3Sum += ME3
            end

            if t_ph == 1
                ME3 = TrOp.E3.n[a_p,a_h] * (-1.0 * X_RPA[J+1,P][Ind_ph,nu] + Y_RPA[J+1,P][Ind_ph,nu]) * sqrt(nRho[a_h,a_h] - nRho[a_p,a_p])
                nME3Sum += ME3
            end

        end
        prME3[nu] = pME3Sum
        nrME3[nu] = nME3Sum
    end

    rM0 = pnCVector(prME0,nrME0)
    rM1 = pnCVector(prME1,nrME1)
    rM2 = pnCVector(prME2,nrME2)
    rM3 = pnCVector(prME3,nrME3)

    rM = ReducedMultipole(rM0,rM1,rM2,rM3)

    println("\nReduced transition matrix elements rM calculated ...")

    return rM
end

function HF_ERPA_rB(Params::Vector{Any},N_nu::Matrix{Int64},rM::ReducedMultipole)
    # Read calculation parameters ...
    Orthogon = Params[5]

    println("\nCalculating reduced transition intensities rB ...")

    # E0
    J = 0
    P = 1
    N_ph = N_nu[J+1,P]

    rB_phE0 = zeros(Float64,N_ph)
    rB_isE0 = zeros(Float64,N_ph)
    rB_ivE0 = zeros(Float64,N_ph)

    @inbounds for nu in 1:N_ph
        rB_phE0[nu] = abs(rM.E0.p[nu])^2
        rB_isE0[nu] = 0.25 * abs(rM.E0.p[nu] + rM.E0.n[nu])^2
        rB_ivE0[nu] = 0.25 * abs(rM.E0.p[nu] - rM.E0.n[nu])^2
    end

    # E1
    J = 1
    P = 2
    N_ph = N_nu[J+1,P]

    rB_phE1 = zeros(Float64,N_ph)
    rB_isE1 = zeros(Float64,N_ph)
    rB_ivE1 = zeros(Float64,N_ph)

    if Orthogon == true
        @inbounds for nu in 1:N_ph
            rB_phE1[nu] = abs(rM.E1.p[nu])^2
            rB_isE1[nu] = 0.25 * abs(rM.E1.p[nu] + rM.E1.n[nu])^2
            rB_ivE1[nu] = 0.25 * abs(rM.E1.p[nu] - rM.E1.n[nu])^2
        end
    else
        A = Params[1]
        Z = Params[2]
        e_p = Float64(A - Z) / Float64(A)
        e_n = Float64(Z) / Float64(A)
        @inbounds for nu in 1:N_ph
            rB_phE1[nu] = abs(rM.E1.p[nu])^2
            rB_isE1[nu] = 0.25 * abs(rM.E1.p[nu] + rM.E1.n[nu])^2
            rB_ivE1[nu] = abs(e_p * rM.E1.p[nu] - e_n * rM.E1.n[nu])^2
        end
    end

    # E2
    J = 2
    P = 1
    N_ph = N_nu[J+1,P]

    rB_phE2 = zeros(Float64,N_ph)
    rB_isE2 = zeros(Float64,N_ph)
    rB_ivE2 = zeros(Float64,N_ph)

    @inbounds for nu in 1:N_ph
        rB_phE2[nu] = abs(rM.E2.p[nu])^2
        rB_isE2[nu] = 0.25 * abs(rM.E2.p[nu] + rM.E2.n[nu])^2
        rB_ivE2[nu] = 0.25 * abs(rM.E2.p[nu] - rM.E2.n[nu])^2
    end

    # E3
    J = 3
    P = 2
    N_ph = N_nu[J+1,P]

    rB_phE3 = zeros(Float64,N_ph)
    rB_isE3 = zeros(Float64,N_ph)
    rB_ivE3 = zeros(Float64,N_ph)

    @inbounds for nu in 1:N_ph
        rB_phE3[nu] = abs(rM.E3.p[nu])^2
        rB_isE3[nu] = 0.25 * abs(rM.E3.p[nu] + rM.E3.n[nu])^2
        rB_ivE3[nu] = 0.25 * abs(rM.E3.p[nu] - rM.E3.n[nu])^2
    end

    rB_E0 = Transition(rB_phE0,rB_isE0,rB_ivE0)
    rB_E1 = Transition(rB_phE1,rB_isE1,rB_ivE1)
    rB_E2 = Transition(rB_phE2,rB_isE2,rB_ivE2)
    rB_E3 = Transition(rB_phE3,rB_isE3,rB_ivE3)

    rB_ph = ReducedTransition(rB_E0,rB_E1,rB_E2,rB_E3)

    println("\nReduced transition intensities rB evaluated ...")

    return rB_ph
end

function HF_ERPA_Transition_Densities_Export(Params::Vector{Any},Orb::Vector{NOrb},N_nu::Matrix{Int64},Orb_Phonon::Matrix{Vector{Int64}},Phonon::Vector{PhState},Particle::pnSVector,Hole::pnSVector,X_RPA::Matrix{Matrix{ComplexF64}},Y_RPA::Matrix{Matrix{ComplexF64}},TrOp::TranOper)
    # Read parameters ...
    A = Params[1]
    Z = Params[2]
    HbarOmega = Params[3]
    N_max = Params[4]
    N_2max = 2*N_max
    a_max = div((N_max+1)*(N_max+2),2)
    J_max = N_2max + 1
    Orthogon = Params[5]
    Output_File = Params[8]

    # Calc params ...
    J = 1
    P = 2
    N_nu_max = 15
    N_ph = N_nu[J+1,P]

    # Basic constants ...
    HbarC = 197.326980
    m_n = 939.565346
    m_p = 938.272013
    nu_proton = 0.5 * m_p * HbarOmega / HbarC^2
    nu_neutron = 0.5 * m_n * HbarOmega / HbarC^2

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

    # Preallocate grid ...
    r1 = 0.0 + 1e-8
    r2 = 2.5 * 1.2 * A^(1/3)
    N_Sampling = 10000
    r_grid = range(r1, stop = r2, length = N_Sampling)
    r_grid = collect(r_grid)
    pRho_ERPA_rad = Matrix{Float64}(undef,N_Sampling,N_nu_max)
    nRho_ERPA_rad = Matrix{Float64}(undef,N_Sampling,N_nu_max)

    # Evaluate radial densities on grid ...
    @inbounds Threads.@threads for i in 1:N_Sampling
        r = r_grid[i]
        @inbounds for nu in 1:N_nu_max
            pRad, nRad = 0.0, 0.0
            @inbounds for ph in 1:N_ph
                ph_ind = Orb_Phonon[J+1,P][ph]
                p,h,t_ph,J_ph,P_ph,a_p,l_p,j_p,E_p,a_h,l_h,j_h,E_h = Phonon_Ind(ph_ind,Phonon,Particle,Hole)
                pPsi, nPsi = 0.0, 0.0
                @inbounds for a_k in 1:a_max
                    l_k = Orb[a_k].l
                    j_k = Orb[a_k].j
                    n_k = Orb[a_k].n
                    if (j_p == j_k) && (l_p == l_k)
                        @inbounds for a_l in 1:a_max
                            l_l = Orb[a_l].l
                            j_l = Orb[a_l].j
                            n_l = Orb[a_l].n
                            if (j_h == j_l) && (l_h == l_l)
                                if t_ph == -1
                                    Psi = Psi_rad_LHO(r,n_k,l_k,nu_proton) * Psi_rad_LHO(r,n_l,l_l,nu_proton) * pU[a_k,a_p] * pU[a_l,a_h]
                                    pPsi += Psi
                                elseif t_ph == 1
                                    Psi = Psi_rad_LHO(r,n_k,l_k,nu_neutron) * Psi_rad_LHO(r,n_l,l_l,nu_neutron) * nU[a_k,a_p] * nU[a_l,a_h]
                                    nPsi += Psi
                                end
                            end
                        end
                    end
                end
                if t_ph == -1
                    Amp = TrOp.E1.p[a_p,a_h] * (-1.0 * real(X_RPA[J+1,P][ph,nu]) + real(Y_RPA[J+1,P][ph,nu])) / Float64(2*J + 1)
                    pRad += pPsi * Amp
                elseif t_ph == 1
                    Amp = TrOp.E1.n[a_p,a_h] * (-1.0 * real(X_RPA[J+1,P][ph,nu]) + real(Y_RPA[J+1,P][ph,nu])) / Float64(2*J + 1)
                    nRad += nPsi * Amp
                end
            end
            pRho_ERPA_rad[i,nu] = pRad * r^2
            nRho_ERPA_rad[i,nu] = nRad * r^2
        end
    end

    @views r_grid = round.(r_grid,digits = 7)
    @views pRho_ERPA_rad = round.(pRho_ERPA_rad,digits = 7)
    @views nRho_ERPA_rad = round.(nRho_ERPA_rad,digits = 7)

    # Export files ...
    if Orthogon == true
        Output_Path = Output_File * "/ERPA/Densities/HF_ERPA_Radial_Transition_Densities_Ortho.dat"
    else
        Output_Path = Output_File * "/ERPA/Densities/HF_ERPA_Radial_Transition_Densities_Spur.dat"
    end

    open(Output_Path, "w") do Export_File
        @inbounds for i in 1:N_Sampling
            print(Export_File, r_grid[i])
            print(Export_File, "\t")
            @inbounds for nu in 1:N_nu_max
                print(Export_File, pRho_ERPA_rad[i,nu])
                print(Export_File, "\t")
                print(Export_File, nRho_ERPA_rad[i,nu])
                print(Export_File, "\t")
            end
            print(Export_File, "\n")
        end
    end

    return
end
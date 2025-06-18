function HF_RPA_rM(Params::Vector{Any},N_nu::Matrix{Int64},Orb_Phonon::Matrix{Vector{Int64}},Phonon::Vector{PhState},Particle::pnSVector,Hole::pnSVector,X_TDA::Matrix{Matrix{Float64}},X_RPA::Matrix{Matrix{ComplexF64}},Y_RPA::Matrix{Matrix{ComplexF64}},TrOp::TranOper)
    # Read calculation parameters ...
    Orthogon = Params[5]

    println("\nCalculating reduced transition matrix elements rM ...")

    # E0
    J = 0
    P = 1
    N_ph = N_nu[J+1,P]
    prME0_TDA = Vector{ComplexF64}(undef,N_ph)
    prME0_RPA = Vector{ComplexF64}(undef,N_ph)
    nrME0_TDA = Vector{ComplexF64}(undef,N_ph)
    nrME0_RPA = Vector{ComplexF64}(undef,N_ph)
    @inbounds for nu in 1:N_ph
        pME0Sum_TDA = ComplexF64(0.0)
        pME0Sum_RPA = ComplexF64(0.0)
        nME0Sum_TDA = ComplexF64(0.0)
        nME0Sum_RPA = ComplexF64(0.0)
        @inbounds  for Ind_ph in 1:N_ph
            ph = Orb_Phonon[J+1,P][Ind_ph]
            p,h,t_ph,J_ph,P_ph,a_p,l_p,j_p,E_p,a_h,l_h,j_h,E_h = Phonon_Ind(ph,Phonon,Particle,Hole)

            if t_ph == -1
                ME0_TDA = TrOp.E0.p[a_p,a_h] * X_TDA[J+1,P][Ind_ph,nu]
                pME0Sum_TDA += ME0_TDA

                ME0_RPA = TrOp.E0.p[a_p,a_h] * (X_RPA[J+1,P][Ind_ph,nu] + Y_RPA[J+1,P][Ind_ph,nu])
                pME0Sum_RPA += ME0_RPA
            end

            if t_ph == 1
                ME0_TDA = TrOp.E0.n[a_p,a_h] * X_TDA[J+1,P][Ind_ph,nu]
                nME0Sum_TDA += ME0_TDA

                ME0_RPA = TrOp.E0.n[a_p,a_h] * (X_RPA[J+1,P][Ind_ph,nu] + Y_RPA[J+1,P][Ind_ph,nu])
                nME0Sum_RPA += ME0_RPA
            end

        end
        prME0_TDA[nu] = pME0Sum_TDA
        nrME0_TDA[nu] = nME0Sum_TDA
        prME0_RPA[nu] = pME0Sum_RPA
        nrME0_RPA[nu] = nME0Sum_RPA
    end

    # E1
    J = 1
    P = 2
    N_ph = N_nu[J+1,P]
    prME1_TDA = Vector{ComplexF64}(undef,N_ph)
    prME1_RPA = Vector{ComplexF64}(undef,N_ph)
    nrME1_TDA = Vector{ComplexF64}(undef,N_ph)
    nrME1_RPA = Vector{ComplexF64}(undef,N_ph)
    @inbounds  for nu in 1:N_ph
        if (Orthogon == true && nu != 1) || (Orthogon == false)
            pME1Sum_TDA = ComplexF64(0.0)
            pME1Sum_RPA = ComplexF64(0.0)
            nME1Sum_TDA = ComplexF64(0.0)
            nME1Sum_RPA = ComplexF64(0.0)
            @inbounds  for Ind_ph in 1:N_ph
                ph = Orb_Phonon[J+1,P][Ind_ph]
                p,h,t_ph,J_ph,P_ph,a_p,l_p,j_p,E_p,a_h,l_h,j_h,E_h = Phonon_Ind(ph,Phonon,Particle,Hole)

                if t_ph == -1
                    ME1_TDA = TrOp.E1.p[a_p,a_h] * X_TDA[J+1,P][Ind_ph,nu] * -1.0
                    pME1Sum_TDA += ME1_TDA

                    ME1_RPA = TrOp.E1.p[a_p,a_h] * (-1.0 * X_RPA[J+1,P][Ind_ph,nu] + Y_RPA[J+1,P][Ind_ph,nu])
                    pME1Sum_RPA += ME1_RPA
                end

                if t_ph == 1
                    ME1_TDA = TrOp.E1.n[a_p,a_h] * X_TDA[J+1,P][Ind_ph,nu] * -1.0
                    nME1Sum_TDA += ME1_TDA

                    ME1_RPA = TrOp.E1.n[a_p,a_h] * (-1.0 * X_RPA[J+1,P][Ind_ph,nu] + Y_RPA[J+1,P][Ind_ph,nu])
                    nME1Sum_RPA += ME1_RPA
                end

            end
            prME1_TDA[nu] = pME1Sum_TDA
            nrME1_TDA[nu] = nME1Sum_TDA
            prME1_RPA[nu] = pME1Sum_RPA
            nrME1_RPA[nu] = nME1Sum_RPA
        end
    end

    # E2
    J = 2
    P = 1
    N_ph = N_nu[J+1,P]
    prME2_TDA = Vector{ComplexF64}(undef,N_ph)
    prME2_RPA = Vector{ComplexF64}(undef,N_ph)
    nrME2_TDA = Vector{ComplexF64}(undef,N_ph)
    nrME2_RPA = Vector{ComplexF64}(undef,N_ph)
    @inbounds  for nu in 1:N_ph
        pME2Sum_TDA = ComplexF64(0.0)
        pME2Sum_RPA = ComplexF64(0.0)
        nME2Sum_TDA = ComplexF64(0.0)
        nME2Sum_RPA = ComplexF64(0.0)
        @inbounds  for Ind_ph in 1:N_ph
            ph = Orb_Phonon[J+1,P][Ind_ph]
            p,h,t_ph,J_ph,P_ph,a_p,l_p,j_p,E_p,a_h,l_h,j_h,E_h = Phonon_Ind(ph,Phonon,Particle,Hole)

            if t_ph == -1
                ME2_TDA = TrOp.E2.p[a_p,a_h] * X_TDA[J+1,P][Ind_ph,nu]
                pME2Sum_TDA += ME2_TDA

                ME2_RPA = TrOp.E2.p[a_p,a_h] * (X_RPA[J+1,P][Ind_ph,nu] + Y_RPA[J+1,P][Ind_ph,nu])
                pME2Sum_RPA += ME2_RPA
            end

            if t_ph == 1
                ME2_TDA = TrOp.E2.n[a_p,a_h] * X_TDA[J+1,P][Ind_ph,nu]
                nME2Sum_TDA += ME2_TDA

                ME2_RPA = TrOp.E2.n[a_p,a_h] * (X_RPA[J+1,P][Ind_ph,nu] + Y_RPA[J+1,P][Ind_ph,nu])
                nME2Sum_RPA += ME2_RPA
            end

        end
        prME2_TDA[nu] = pME2Sum_TDA
        nrME2_TDA[nu] = nME2Sum_TDA
        prME2_RPA[nu] = pME2Sum_RPA
        nrME2_RPA[nu] = nME2Sum_RPA
    end

    # E3
    J = 3
    P = 2
    N_ph = N_nu[J+1,P]
    prME3_TDA = Vector{ComplexF64}(undef,N_ph)
    prME3_RPA = Vector{ComplexF64}(undef,N_ph)
    nrME3_TDA = Vector{ComplexF64}(undef,N_ph)
    nrME3_RPA = Vector{ComplexF64}(undef,N_ph)
    @inbounds  for nu in 1:N_ph
        pME3Sum_TDA = ComplexF64(0.0)
        pME3Sum_RPA = ComplexF64(0.0)
        nME3Sum_TDA = ComplexF64(0.0)
        nME3Sum_RPA = ComplexF64(0.0)
        @inbounds  for Ind_ph in 1:N_ph
            ph = Orb_Phonon[J+1,P][Ind_ph]
            p,h,t_ph,J_ph,P_ph,a_p,l_p,j_p,E_p,a_h,l_h,j_h,E_h = Phonon_Ind(ph,Phonon,Particle,Hole)

            if t_ph == -1
                ME3_TDA = TrOp.E3.p[a_p,a_h] * X_TDA[J+1,P][Ind_ph,nu] * -1.0
                pME3Sum_TDA += ME3_TDA

                ME3_RPA = TrOp.E3.p[a_p,a_h] * (-1.0 * X_RPA[J+1,P][Ind_ph,nu] + Y_RPA[J+1,P][Ind_ph,nu])
                pME3Sum_RPA += ME3_RPA
            end

            if t_ph == 1
                ME3_TDA = TrOp.E3.n[a_p,a_h] * X_TDA[J+1,P][Ind_ph,nu] * -1.0
                nME3Sum_TDA += ME3_TDA

                ME3_RPA = TrOp.E3.n[a_p,a_h] * (-1.0 * X_RPA[J+1,P][Ind_ph,nu] + Y_RPA[J+1,P][Ind_ph,nu])
                nME3Sum_RPA += ME3_RPA
            end

        end
        prME3_TDA[nu] = pME3Sum_TDA
        nrME3_TDA[nu] = nME3Sum_TDA
        prME3_RPA[nu] = pME3Sum_RPA
        nrME3_RPA[nu] = nME3Sum_RPA
    end

    rM0_TDA = pnCVector(prME0_TDA,nrME0_TDA)
    rM1_TDA = pnCVector(prME1_TDA,nrME1_TDA)
    rM2_TDA = pnCVector(prME2_TDA,nrME2_TDA)
    rM3_TDA = pnCVector(prME3_TDA,nrME3_TDA)

    rM0_RPA = pnCVector(prME0_RPA,nrME0_RPA)
    rM1_RPA = pnCVector(prME1_RPA,nrME1_RPA)
    rM2_RPA = pnCVector(prME2_RPA,nrME2_RPA)
    rM3_RPA = pnCVector(prME3_RPA,nrME3_RPA)

    rM_TDA = ReducedMultipole(rM0_TDA,rM1_TDA,rM2_TDA,rM3_TDA)
    rM_RPA = ReducedMultipole(rM0_RPA,rM1_RPA,rM2_RPA,rM3_RPA)

    println("\nReduced transition matrix elements rM calculated ...")

    return rM_TDA, rM_RPA
end

function HF_RPA_rB(Params::Vector{Any},N_nu::Matrix{Int64},rM::ReducedMultipole)
    # Read parameters ...
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

function HF_RPA_Transition_Densities_Export(Params::Vector{Any},Orb::Vector{NOrb},N_nu::Matrix{Int64},Orb_Phonon::Matrix{Vector{Int64}},Phonon::Vector{PhState},Particle::pnSVector,Hole::pnSVector,X_RPA::Matrix{Matrix{ComplexF64}},Y_RPA::Matrix{Matrix{ComplexF64}},TrOp::TranOper)
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
    pRho_RPA_rad = Matrix{Float64}(undef,N_Sampling,N_nu_max)
    nRho_RPA_rad = Matrix{Float64}(undef,N_Sampling,N_nu_max)

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
            pRho_RPA_rad[i,nu] = pRad * r^2
            nRho_RPA_rad[i,nu] = nRad * r^2
        end
    end

    @views r_grid = round.(r_grid,digits = 7)
    @views pRho_RPA_rad = round.(pRho_RPA_rad,digits = 7)
    @views nRho_RPA_rad = round.(nRho_RPA_rad,digits = 7)

    # Export files ...
    if Orthogon == true
        Output_Path = Output_File * "/RPA/Densities/HF_RPA_Radial_Transition_Densities_Ortho.dat"
    else
        Output_Path = Output_File * "/RPA/Densities/HF_RPA_Radial_Transition_Densities_Spur.dat"
    end

    open(Output_Path, "w") do Export_File
        @inbounds for i in 1:N_Sampling
            print(Export_File, r_grid[i])
            print(Export_File, "\t")
            @inbounds for nu in 1:N_nu_max
                print(Export_File, pRho_RPA_rad[i,nu])
                print(Export_File, "\t")
                print(Export_File, nRho_RPA_rad[i,nu])
                print(Export_File, "\t")
            end
            print(Export_File, "\n")
        end
    end

    return
end
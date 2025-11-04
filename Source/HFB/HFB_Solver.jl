function HFB_Solver(Params::Parameters)
    # Iteration precision parameter ...
    epsilon = 1e-7

    # Make single-particle orbitals - NuHamil ordering ...
    Orb = Make_Orbitals(Params.Calc.A,Params.Calc.Z,Params.Int.Nmax)

    # Load 1-body kinetic operator ...
    T = T1B(Params.Int.Nmax,Orb,Params.Int.hw)

    # 2-body NN interaction & Orbitals ...
    @time VNN, Orb_NN = V2B_Read(Params,Orb)

    # 3-body NNN interaction & Orbitals ...
    @time VNNN, Orb_NNN = V3B_NO2B_Read(Params,Orb)

    # Solve HFB equations ...
    @time Lambda, SQE, U, V, SQE_C, u_C, v_C, C, Rho, Kappa, H, Delta, Iteration = HFB_Solve(Params,Orb,Orb_NN,Orb_NNN,T,VNN,VNNN,epsilon)

    # Calculation of the total HFB mean-field ground-state energy ...
    @time E_HFB = HFB_Energy(Params,Rho,Kappa,Orb,Orb_NN,Orb_NNN,T,VNN,VNNN)

    # Calculate the total HFB kinetic energy ...
    T_HFB = Kinetic_Energy(Params,Rho,Orb,T)

    # Particle number fluctuation calculation ...
    dA = HFB_Particle_Number_Dispersion(Params,Rho,Orb)

    # Calculation summary ...
    HFB_Summary(Params,E_HFB,T_HFB,Lambda,dA,epsilon,Iteration)

    # Evaluate HFB charge radii & radial densities ...
    Summary_File = "IO/" * Params.Calc.Path * "/HFB/HFB_Summary.dat"
    Densities_File = "IO/" * Params.Calc.Path * "/HFB/Densities/HFB_Radial_Densities.dat"
    OBDM_Export(Params,Summary_File,Densities_File,Rho,C,Orb)

    # Export of single-quasiparticle energies, amplitudes U & V & also possibly radial densities ...
    HFB_SQS_Summary(Params,SQE,SQE_C,pnMatrix(C.p' * Rho.p * C.p,C.n' * Rho.n * C.n),Orb)

    return
end

function HFB_Solve(Params::Parameters,Orb::Vector{NOrb},Orb_NN::NNOrb,Orb_NNN::NNNOrb,T::Matrix{Float64},VNN::NNInt,VNNN::Array{Vector{Vector{Float32}},4},epsilon::Float64)
    # Read calculation parameters ...
    A_target, Z_target, N_target = Float64(Params.Calc.A), Float64(Params.Calc.Z), Float64(Params.Calc.A - Params.Calc.Z)
    hw = Params.Int.hw
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1) * (N_max + 2), 2)

    # Setup local iteration variables ...
    Iteration, Iteration_max = 0, 250
    dE, dZ, dN = 0.1, 0.1, 0.1

    # Preallocate arrays ...
        # Matrices for U & V HFB ...
    V = pnMatrix(zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max))
    U = pnMatrix(zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max))
         # Matrices for densities Rho & Kappa ...
    Rho = pnMatrix(zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max))
    Kappa = pnMatrix(zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max))
        # Matrices for the single-particle field H & pairing field Delta ...
    H = pnMatrix(zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max))
    Delta = pnMatrix(zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max))
        # Vectors for single-(quasi)particle energies ...
    SQE = pnVector(zeros(Float64,a_max),zeros(Float64,a_max))
        # Identity matrix ...
    Eye = diagm(ones(Float64,a_max))

    # Setup particle numbers ...
    Z, N = 0.1, 0.1

    # Initial guess on chemical potentials lambda ...
    Lambda = pnFloat(0.5,0.5)

    # Initial guess on densities Rho & Kappa
    Rho, Kappa = HFB_Density_Operator_Initialize(Params,Orb)

    # Preallocate Rho_old, Kappa_old & SQE_old ...
    Rho_old = pnMatrix(Rho.p, Rho.n)
    Kappa_old = pnMatrix(Kappa.p, Kappa.n)
    SQE_old = pnVector(SQE.p, SQE.n)

    # Initialize particle numbers & chemical potentials for secant method ...
    Z_1, Z_2 = 0.0, 0.0
    N_1, N_2 = 0.0, 0.0
    Lambda_1, Lambda_2 = pnFloat(0.75 * Lambda.p, 0.75 * Lambda.n), pnFloat(Lambda.p, Lambda.n)

    # Solve the spherical HFB equations ... by the means of self-consistent iteration ...
    println("\nStarting iteration of HFB equations ...\n")

    while ((dE > epsilon) || (dZ > epsilon) || (dN > epsilon)) && (Iteration < Iteration_max)

        # Perform several secant iterations for the chemical potentil Lambda ...
        @inbounds for L in 1:6
            # Initialization of the secant method ...
            if (L == 1) && (Iteration == 0)
                # Allocate the single-particle fields H_1, H_2 and pairing fields Delta_1, Delta_2 ...
                H_1, Delta_1 = HFB_Allocate(Params,Lambda_1,Rho_old,Kappa_old,Orb,Orb_NN,Orb_NNN,T,VNN,VNNN)
                H_2, Delta_2 = HFB_Allocate(Params,Lambda_2,Rho_old,Kappa_old,Orb,Orb_NN,Orb_NNN,T,VNN,VNNN)

                # Diagonalize the HFB equations ... basis is reordered as needed ...
                SQE_1, U_1, V_1 = HFB_Diagonalize(Params,H_1,Delta_1,Orb)
                SQE_2, U_2, V_2 = HFB_Diagonalize(Params,H_2,Delta_2,Orb)
                
                # Generate temporary densities Rho & Kappa ...
                Rho_1, Kappa_1 = HFB_Density_Operator(Params,U_1,V_1,Orb)
                Rho_2, Kappa_2 = HFB_Density_Operator(Params,U_2,V_2,Orb)

                # Determine new average particle numbers ...
                Z_1, N_1 = HFB_Particle_Number(Params,Rho_1,Orb)
                Z_2, N_2 = HFB_Particle_Number(Params,Rho_2,Orb)

                # Update average particle numbers ...
                Z, N = Z_2, N_2

                # Perform secant iteration to determine new optimal value of Lambda ...
                Lambda = HFB_Lambda_Secant(Lambda_1,pnFloat(Z_target - Z_1, N_target - N_1),Lambda_2,pnFloat(Z_target - Z_2, N_target - N_2))

                # Update chemical potentials ...
                Lambda_1 = pnFloat(Lambda_2.p, Lambda_2.n)
                Lambda_2 = pnFloat(Lambda.p, Lambda.n)

                # Update single-quasiparticle energies ...
                SQE = pnVector(SQE_2.p, SQE_2.n)

                # Update amplitudes U & V ...
                U, V = pnMatrix(U_2.p,U_2.n), pnMatrix(V_2.p,V_2.n)

                # Update fields H & Delta ...
                H, Delta = pnMatrix(H_2.p,H_2.n), pnMatrix(Delta_2.p,Delta_2.n)

            end

            # Allocate the single-particle field H and the pairing field Delta ...
                # Smart re-allocation of H ... only Lambda is tweaked, Delta remains the same ...
            H = pnMatrix(H.p .+ (Lambda_1.p - Lambda.p ) .* Eye, H.n .+ (Lambda_1.n - Lambda.n ) .* Eye)
                # Not so smart allocation ... full iteration ... redundant & slow ...
            #H, Delta = HFB_Allocate(Params,Lambda,Rho_old,Kappa_old,Orb,Orb_NN,Orb_NNN,T,VNN,VNNN)

            # Diagonalize the HFB equations ... basis is reordered as needed ...
            SQE, U, V = HFB_Diagonalize(Params,H,Delta,Orb)
            
            # Generate temporary densities Rho & Kappa ...
            Rho_temp, Kappa_temp = HFB_Density_Operator(Params,U,V,Orb)

            # Determine new average particle numbers ...
            Z, N = HFB_Particle_Number(Params,Rho_temp,Orb)

            # Update average particle numbers ...
            Z_1, N_1 = Z_2, N_2
            Z_2, N_2 = Z, N

            # Perform secant iteration to determine new optimal value of Lambda ...
            Lambda = HFB_Lambda_Secant(Lambda_1,pnFloat(Z_target - Z_1, N_target - N_1),Lambda_2,pnFloat(Z_target - Z_2, N_target - N_2))

            # Update chemical potentials ...
            Lambda_1 = pnFloat(Lambda_2.p, Lambda_2.n)
            Lambda_2 = pnFloat(Lambda.p, Lambda.n)

            # Evaluate iteration of the chemical potential Lambda ...
            dZ, dN = abs(Z_target - Z), abs(N_target - N)

            # Check finite differences for particle numbers ....
            if dZ < epsilon && dN < epsilon
                break
            end
        end

        # Generate new densities Rho & Kappa ...
        Rho, Kappa = HFB_Density_Operator(Params,U,V,Orb)

        # Perform update of densities ...
        Rho, Kappa = HFB_Density_Operator_Update(a_max,Rho,Rho_old,Kappa,Kappa_old)
            #Rho = pnMatrix(0.5 .* Rho_old.p .+ 0.5 .* Rho.p, 0.5 .* Rho_old.n .+ 0.5 .* Rho.n)
            #Kappa = pnMatrix(0.5 .* Kappa_old.p .+ 0.5 .* Kappa.p, 0.5 .* Kappa_old.n .+ 0.5 .* Kappa.n)

        # Allocate the single-particle field H and the pairing field Delta ...
        H, Delta = HFB_Allocate(Params,Lambda,Rho,Kappa,Orb,Orb_NN,Orb_NNN,T,VNN,VNNN)

        # Diagonalize the HFB equations ... basis is reordered as needed ...
        SQE, U, V = HFB_Diagonalize(Params,H,Delta,Orb)

        # Evluate the particle numbers ...
        Z, N = HFB_Particle_Number(Params,Rho,Orb)

        # Evaluate iteration of the single-quasiparticle energies ...
        dE = (sum(abs.(SQE.p .- SQE_old.p )) + sum(abs.(SQE.n .- SQE_old.n))) / Float64(2 * a_max)
        dZ, dN = abs(Z_target - Z), abs(N_target - N)

        Iteration += 1

        println("\n\nHFB iteration number:   " * string(Iteration) * "   Single-quasiparticle energy difference:   " * string(round(dE, sigdigits=8))
                * "   Proton number difference:   " * string(round(dZ, sigdigits=8)) * "   Neutron number difference:   " * string(round(dN, sigdigits=8)))

        # Store old values of SQE, Rho & Kappa ...
        SQE_old = pnVector(SQE.p, SQE.n)
        Rho_old = pnMatrix(Rho.p, Rho.n)
        Kappa_old = pnMatrix(Kappa.p, Kappa.n)

    end

    println("\nHFB iteration with residual NN interaction has converged ...")

    # Construct the canonical basis & evaluate approximate (BCS-like) amplitudes
    # & single-quasiparticle energies u_C, v_C & SQE_C ...
    #   U, V, Rho, Kappa, H, Delta ... remain expressed in the reference LHO basis ...
    SQE_C, C, u_C, v_C = HFB_Canonical_Basis(Params,Rho,H,Delta,Orb)

    return Lambda, SQE, U, V, SQE_C, u_C, v_C, C, Rho, Kappa, H, Delta, Iteration
end

function HFB_Density_Operator_Update(a_max::Int64,Rho::pnMatrix,Rho_old::pnMatrix,Kappa::pnMatrix,Kappa_old::pnMatrix)
    # Basic parameters ...
    q, q_min = 0.95, 1e-8

    # Allocate finite differences for densities ...
    dRho = pnMatrix(Rho.p .- Rho_old.p, Rho.n .- Rho_old.n)
    dKappa = pnMatrix(Kappa.p .- Kappa_old.p, Kappa.n .- Kappa_old.n)

    # Evaluate current finite difference d_current ...
    D = (sum(abs.(dRho.p)) + sum(abs.(dKappa.p)) + sum(abs.(dRho.n)) + sum(abs.(dKappa.n))) / Float64(4*a_max^2)

    # Iterate the quenching of proton densities ...
    while q > q_min
        # Perform trial step ...
        pRho_trial, pKappa_trial = (1.0 - q) * Rho_old.p .+ q * Rho.p, (1.0 - q) * Kappa_old.p .+ q * Kappa.p
        nRho_trial, nKappa_trial = (1.0 - q) * Rho_old.n .+ q * Rho.n, (1.0 - q) * Kappa_old.n .+ q * Kappa.n

        # Symmetrize trial densities ...
        pRho_trial .= 0.5 * (pRho_trial .+ pRho_trial')
        pKappa_trial .= 0.5 * (pKappa_trial .+ pKappa_trial')
        nRho_trial .= 0.5 * (nRho_trial .+ nRho_trial')
        nKappa_trial .= 0.5 * (nKappa_trial .+ nKappa_trial')

        D_trial = (sum(abs.(pRho_trial .- Rho_old.p)) + sum(abs.(pKappa_trial .- Kappa_old.p)) + sum(abs.(nRho_trial .- Rho_old.n)) + sum(abs.(nKappa_trial .- Kappa_old.n))) / Float64(4*a_max^2)

        # Condition on accepting the current trial step ...
        if D_trial / D > 1.1
            q = 0.5 * q
        #elseif D_trial / D < 0.1
        #    q = 1.1 * q
        else
            return pnMatrix(pRho_trial,nRho_trial), pnMatrix(pKappa_trial,nKappa_trial)
        end
    end


    # No improvement due to quenching ... return old densities ...
    return pnMatrix(0.975 * Rho.p, 0.975 * Rho.n), pnMatrix(0.975 * Kappa.p, 0.975 * Kappa.n)
end

function HFB_Density_Operator_Update2(a_max::Int64,Rho::pnMatrix,Rho_old::pnMatrix,Kappa::pnMatrix,Kappa_old::pnMatrix)
    # Basic parameters ...
    pQ, nQ, Q_min = 0.95, 0.95, 1e-9

    pRho, pKappa = 0.975 * Rho_old.p, 0.975 * Kappa_old.p
    nRho, nKappa = 0.975 * Rho_old.n, 0.975 * Kappa_old.n

    # Allocate finite differences for densities ...
    dRho = pnMatrix(Rho.p .- Rho_old.p, Rho.n .- Rho_old.n)
    dKappa = pnMatrix(Kappa.p .- Kappa_old.p, Kappa.n .- Kappa_old.n)

    # Evaluate current finite difference d_current ...
    pD = (sum(abs.(dRho.p)) + sum(abs.(dKappa.p))) / Float64(2*a_max^2)
    nD = (sum(abs.(dRho.n)) + sum(abs.(dKappa.n))) / Float64(2*a_max^2)

    # Iterate the quenching of proton densities ...
    while pQ > Q_min
        # Perform trial step ...
        pRho_trial, pKappa_trial = (1.0 - pQ) * Rho_old.p .+ pQ * Rho.p, (1.0 - pQ) * Kappa_old.p .+ pQ * Kappa.p

        # Symmetrize trial densities ...
        pRho_trial .= 0.5 * (pRho_trial .+ pRho_trial')
        pKappa_trial .= 0.5 * (pKappa_trial .+ pKappa_trial')

        pD_trial = (sum(abs.(pRho_trial .- Rho_old.p)) + sum(abs.(pKappa_trial .- Kappa_old.p))) / Float64(2*a_max^2)

        # Condition on accepting the current trial step ...
        if pD_trial / pD > 1.1
            pQ = 0.5 * pQ
        elseif pD_trial / pD < 0.5
            pQ = 1.5 * pQ
        else
            pRho .= pRho_trial
            pKappa .= pKappa_trial
            break
        end
    end

    # Iterate the quenching of proton densities ...
    while nQ > Q_min
        # Perform trial step ...
        nRho_trial, nKappa_trial = (1.0 - nQ) * Rho_old.n .+ nQ * Rho.n, (1.0 - nQ) * Kappa_old.n .+ nQ * Kappa.n

        # Symmetrize trial densities ...
        nRho_trial .= 0.5 * (nRho_trial .+ nRho_trial')
        nKappa_trial .= 0.5 * (nKappa_trial .+ nKappa_trial')

        nD_trial = (sum(abs.(nRho_trial .- Rho_old.n)) + sum(abs.(nKappa_trial .- Kappa_old.n))) / Float64(2*a_max^2)

        # Condition on accepting the current trial step ...
        if nD_trial / nD > 1.1
            nQ = 0.5 * nQ
        elseif nD_trial / nD < 0.5
            nQ = 1.5 * nQ
        else
            nRho .= nRho_trial
            nKappa .= nKappa_trial
            break
        end
    end


    # No improvement due to quenching ... return old densities ...
    return pnMatrix(pRho,nRho), pnMatrix(pKappa,nKappa)
end

function HFB_Density_Operator_Update3(a_max::Int64,Rho::pnMatrix,Rho_old::pnMatrix,Kappa::pnMatrix,Kappa_old::pnMatrix)
    # Adaptive linear mixing based on residual magnitudes (simple, robust)
    # Compute proton/neutron residual magnitudes (same measure as before)
    pD = (sum(abs.(Rho.p .- Rho_old.p)) + sum(abs.(Kappa.p .- Kappa_old.p))) / Float64(2*a_max^2)
    nD = (sum(abs.(Rho.n .- Rho_old.n)) + sum(abs.(Kappa.n .- Kappa_old.n))) / Float64(2*a_max^2)

    # Adaptive mixing parameter: smaller mixing (more damping) for large residuals,
    # larger mixing (faster update) when residuals small. Clamp to [min_mix, max_mix].
    min_mix, max_mix = 0.08, 0.95
    # simple mapping: mix -> decreases with residual, roughly in (min_mix..max_mix)
    p_mix = clamp( max_mix * exp(-5.0 * pD) , min_mix, max_mix )
    n_mix = clamp( max_mix * exp(-5.0 * nD) , min_mix, max_mix )

    # Form mixed densities and enforce Hermiticity / symmetry
    pRho = (1.0 - p_mix) .* Rho_old.p .+ p_mix .* Rho.p
    pKappa = (1.0 - p_mix) .* Kappa_old.p .+ p_mix .* Kappa.p
    nRho = (1.0 - n_mix) .* Rho_old.n .+ n_mix .* Rho.n
    nKappa = (1.0 - n_mix) .* Kappa_old.n .+ n_mix .* Kappa.n

    # Symmetrize
    pRho .= 0.5 .* (pRho .+ pRho')
    nRho .= 0.5 .* (nRho .+ nRho')
    pKappa .= 0.5 .* (pKappa .+ pKappa')
    nKappa .= 0.5 .* (nKappa .+ nKappa')

    return pnMatrix(pRho,nRho), pnMatrix(pKappa,nKappa)
end

function HFB_Canonical_Basis(Params::Parameters,Rho::pnMatrix,H::pnMatrix,Delta::pnMatrix,Orb::Vector{NOrb})
    # Read calculation parameters ...
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1) * (N_max + 2), 2)

    # Allocate density matrices ...
    pRho, nRho = Rho.p, Rho.n

    # Symmetrize density matrices ...
    pRho .= 0.5 .* (pRho .+ pRho')
    nRho .= 0.5 .* (nRho .+ nRho')

    # Eliminate any possible numerical noise spoiling block-diagonal structure ...
    @inbounds for a in 1:a_max
        @inbounds for b in 1:a_max
            if (Orb[a].j != Orb[b].j) || (Orb[a].l != Orb[b].l)
                pRho[a,b] = 0.0
                nRho[a,b] = 0.0
            end
        end
    end

    # Add a tiny deterministic diagonal splitting to lift accidental degeneracies ...
    @inbounds for i in 1:a_max
        pRho[i,i] += 1e-10 * i
        nRho[i,i] += 1e-10 * i
    end

    # Diagonalize 1-body HFB density matrix Rho ...
    pn_C, pC = eigen(Symmetric(pRho))
    nn_C, nC = eigen(Symmetric(nRho))

    # Clean numerical noise in C ...
    @inbounds for a in 1:a_max
        @inbounds for b in 1:a_max
            pCME, nCME = abs(pC[a,b]), abs(nC[a,b])
            if pCME < 1e-6
                pC[a,b] = 0.0
            end
            if nCME < 1e-6
                nC[a,b] = 0.0
            end
        end
    end

    # Reorder the transformation matrix C & occupation probabilities n_C ...
    pn_C, pC, nn_C, nC = HFB_Canonical_Basis_Particle_Reordering(Params,pnMatrix(pC,nC),pnVector(pn_C,nn_C),Orb)

    # Calculate canonical amplitudes v & u ...
        # Allocate v_C ... from the occupation probabilities
    pv_C, nv_C = abs.(pn_C) .+ 1e-14, abs.(nn_C) .+ 1e-14
        # Calculate u_C ... from the normalization condition ... |u_k|^2 + |v_k|^2 = 1
    pu_C, nu_C = abs.(ones(Float64,a_max) .- pv_C) .+ 1e-14, abs.(ones(Float64,a_max) .- nv_C) .+ 1e-14

    # Proper normalization of u_C & v_C ... square-root ...
    pu_C .= sqrt.(pu_C)
    pv_C .= sqrt.(pv_C)
    nu_C .= sqrt.(nu_C)
    nv_C .= sqrt.(nv_C)

    pu_C, pv_C = diagm(pu_C), diagm(pv_C)
    nu_C, nv_C = diagm(nu_C), diagm(nv_C)

    # Calculate the canonical single-quasiparticle energies (approximate to exact HFB SQE!)...
    SQE_C = HFB_Canonical_Basis_SQE(Params,pnMatrix(pC,nC),H,Delta)

    # Reorder the single-quasiparticle orbitals ... u_C, v_C & SQE_C ...
    SQE_C, u_C, v_C = HFB_Canonical_Basis_Quasiparticle_Reordering(Params,SQE_C,pnMatrix(pu_C,nu_C),pnMatrix(pv_C,nv_C),Orb)

    return SQE_C, pnMatrix(pC,nC), u_C, v_C
end

function HFB_Canonical_Basis_SQE(Params::Parameters,C::pnMatrix,H::pnMatrix,Delta::pnMatrix)
    # Read calculation parameters ...
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1) * (N_max + 2), 2)

    # Preallocate vectors for single-quasiparticle energies in the canonical basis ...
    pSQE_C, nSQE_C = zeros(Float64,a_max), zeros(Float64,a_max)

    # Transform H & Delta into the canonical basis ...
    pH_C, pDelta_C = C.p' * H.p * C.p, C.p' * Delta.p * C.p
    nH_C, nDelta_C = C.n' * H.n * C.n, C.n' * Delta.n * C.n

    # Calculate single-quasiparticle energies in the canonical basis ...
    @inbounds for a in 1:a_max
        pE_C = sqrt(pH_C[a,a]^2 + pDelta_C[a,a]^2)
        nE_C = sqrt(nH_C[a,a]^2 + nDelta_C[a,a]^2)

        pSQE_C[a] = pE_C
        nSQE_C[a] = nE_C
    end

    return pnVector(pSQE_C,nSQE_C)
end

function HFB_Canonical_Basis_Particle_Reordering(Params::Parameters,C::pnMatrix,n_C::pnVector,Orb::Vector{NOrb})
    # Read calculation parameters ...
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1) * (N_max + 2), 2)

    # Read needed arrays ...
    pC, pn_C = C.p, n_C.p
    nC, nn_C = C.n, n_C.n

    # Hungarian algorithm reordering ... first pre-sort
        # Evaluate the basis overlaps ...
    pOverlap = diagm(ones(Float64,a_max)) - abs.(pC)
    nOverlap = diagm(ones(Float64,a_max)) - abs.(nC)
        # Apply the Hungarian algorithm ...
    pOrb_order, Temp = hungarian(pOverlap)
    nOrb_order, Temp = hungarian(nOverlap)
        # Apply the Hungarian reordering ...
    @views pC .= pC[:,pOrb_order]
    @views pn_C .= pn_C[pOrb_order]

    @views nC .= nC[:,nOrb_order]
    @views nn_C .= nn_C[nOrb_order]

    # Continue with reordering in l & j numbers ...

    # Preallocate temporary arrays ...
    pOrb_order, nOrb_order = Vector{Int64}(undef,a_max), Vector{Int64}(undef,a_max)
    pOrb_mask, nOrb_mask = falses(a_max), falses(a_max)

    pl_values, pj_values = zeros(Float64,a_max), zeros(Float64,a_max)
    nl_values, nj_values = zeros(Float64,a_max), zeros(Float64,a_max)

    # Evaluate values of j & l for single-particle orbitals ...
    @inbounds for a in 1:a_max
        pjSum, plSum = 0.0, 0.0
        njSum, nlSum = 0.0, 0.0
        @inbounds for b in 1:a_max
            l_b, j_b = Float64(Orb[b].l), Float64(Orb[b].j)
            pN, nN = pC[b,a]^2, nC[b,a]^2
            pjSum, plSum = pjSum + j_b * pN, plSum + l_b * pN
            njSum, nlSum = njSum + j_b * nN, nlSum + l_b * nN
        end
        pl_values[a], pj_values[a] = plSum, pjSum
        nl_values[a], nj_values[a] = nlSum, njSum
    end

    # Find the ordering for single-particle orbitals ...
    @inbounds for a in 1:a_max
        pl, pj = pl_values[a], pj_values[a]
        nl, nj = nl_values[a], nj_values[a]
        @inbounds for b in 1:a_max
            if (pOrb_mask[b] == false) && (abs(Float64(Orb[b].l) - pl) < 1e-3) && (abs(Float64(Orb[b].j) - pj) < 1e-3)
                pOrb_order[b] = a
                pOrb_mask[b] = true
                break
            end
        end

        @inbounds for b in 1:a_max
            if (nOrb_mask[b] == false) && (abs(Float64(Orb[b].l) - nl) < 1e-3) && (abs(Float64(Orb[b].j) - nj) < 1e-3)
                nOrb_order[b] = a
                nOrb_mask[b] = true
                break
            end
        end
    end

    # Perform reordering of single-particle orbitals ...
    @views pC .= pC[:,pOrb_order]
    @views pn_C .= pn_C[pOrb_order]

    @views nC .= nC[:,nOrb_order]
    @views nn_C .= nn_C[nOrb_order]

    return pn_C, pC, nn_C, nC
end

function HFB_Canonical_Basis_Quasiparticle_Reordering(Params::Parameters,SQE_C::pnVector,u_C::pnMatrix,v_C::pnMatrix,Orb::Vector{NOrb})
    # Read calculation parameters ...
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1) * (N_max + 2), 2)

    # Read needed arrays ...
    pSQE_C, pu_C, pv_C = SQE_C.p, u_C.p, v_C.p
    nSQE_C, nu_C, nv_C = SQE_C.n, u_C.n, v_C.n

    # Reorder basis according to energy ...
    pOrbs_Sort = sortperm(pSQE_C)
    nOrbs_Sort = sortperm(nSQE_C)

    # Proton single-quasiparticle orbitals ...
    @views pSQE_C .= pSQE_C[pOrbs_Sort]
    @views pu_C .= pu_C[:,pOrbs_Sort]
    @views pv_C .= pv_C[:,pOrbs_Sort]

    # Neutron single-quasiparticle orbitals ...
    @views nSQE_C .= nSQE_C[nOrbs_Sort]
    @views nu_C .= nu_C[:,nOrbs_Sort]
    @views nv_C .= nv_C[:,nOrbs_Sort]

    # Next perform reordering according to numbers j & l ...

    # Preallocate temporary arrays ...
    pOrb_order, nOrb_order = Vector{Int64}(undef,a_max), Vector{Int64}(undef,a_max)
    pOrb_mask, nOrb_mask = falses(a_max), falses(a_max)

    pl_values, pj_values = zeros(Float64,a_max), zeros(Float64,a_max)
    nl_values, nj_values = zeros(Float64,a_max), zeros(Float64,a_max)

    # Evaluate values of j & l for single-quasiparticle orbitals ...
    @inbounds for a in 1:a_max
        pjSum, plSum = 0.0, 0.0
        njSum, nlSum = 0.0, 0.0
        @inbounds for b in 1:a_max
            l_b, j_b = Float64(Orb[b].l), Float64(Orb[b].j)
            pN, nN = pu_C[b,a]^2 + pv_C[b,a]^2, nu_C[b,a]^2 + nv_C[b,a]^2
            pjSum, plSum = pjSum + j_b * pN, plSum + l_b * pN
            njSum, nlSum = njSum + j_b * nN, nlSum + l_b * nN
        end
        pl_values[a], pj_values[a] = plSum, pjSum
        nl_values[a], nj_values[a] = nlSum, njSum
    end

    # Find the ordering for single-particle orbitals ...
    @inbounds for a in 1:a_max
        pl, pj = pl_values[a], pj_values[a]
        nl, nj = nl_values[a], nj_values[a]
        @inbounds for b in 1:a_max
            if pOrb_mask[b] == false && abs(Float64(Orb[b].l) - pl) < 1e-7 && abs(Float64(Orb[b].j) - pj) < 1e-7
                pOrb_order[b] = a
                pOrb_mask[b] = true
                break
            end
        end

        @inbounds for b in 1:a_max
            if nOrb_mask[b] == false && abs(Float64(Orb[b].l) - nl) < 1e-7 && abs(Float64(Orb[b].j) - nj) < 1e-7
                nOrb_order[b] = a
                nOrb_mask[b] = true
                break
            end
        end
    end

    # Perform reordering of single-quasiparticle orbitals ...
    @views pSQE_C .= pSQE_C[pOrb_order]
    @views pu_C .= pu_C[:,pOrb_order]
    @views pv_C .= pv_C[:,pOrb_order]

    @views nSQE_C .= nSQE_C[nOrb_order]
    @views nu_C .= nu_C[:,nOrb_order]
    @views nv_C .= nv_C[:,nOrb_order]

    return pnVector(pSQE_C,nSQE_C), pnMatrix(pu_C,nu_C), pnMatrix(pv_C,nv_C)
end

function HFB_Lambda_Secant(Lambda_1::pnFloat,dA_1::pnFloat,Lambda_2::pnFloat,dA_2::pnFloat)
    pLambda, nLambda = 0.0, 0.0

    pSecant = (Lambda_1.p * dA_2.p - Lambda_2.p * dA_1.p) / (dA_2.p - dA_1.p)
    nSecant = (Lambda_1.n * dA_2.n - Lambda_2.n * dA_1.n) / (dA_2.n - dA_1.n)

    if abs(pSecant) < 3.5
        pLambda = pSecant
    else
        pLambda = Lambda_2.p + 0.25 * dA_2.p
    end

    if abs(nSecant) < 3.5
        nLambda = (Lambda_1.n * dA_2.n - Lambda_2.n * dA_1.n) / (dA_2.n - dA_1.n)
    else
        nLambda = Lambda_2.n + 0.25 * dA_2.n
    end

    return pnFloat(pLambda,nLambda)
end

function HFB_Density_Operator_Initialize(Params::Parameters,Orb::Vector{NOrb})
    # Read parameters ...
    Z_target, N_target = Float64(Params.Calc.Z), Float64(Params.Calc.A - Params.Calc.Z)
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1) * (N_max + 2), 2)

    # Initialize density matrices ...
    pRho, nRho = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)
    pKappa, nKappa = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)

    # Initialize particle numbers ...
    Z, N = 0.0, 0.0

    # Fill the proton density matrix ...
    @inbounds for a in 1:a_max
        if (Z - Z_target) < 1e-3
            if ((Z_target - Z) - Float64(Orb[a].j + 1)) > 1e-7
                pRho[a,a] = 1.0
                Z += Float64(Orb[a].j + 1)
            elseif ((Z_target - Z) - Float64(Orb[a].j + 1)) < 1e-7
                pRho[a,a] = abs(Z_target - Z) / Float64(Orb[a].j + 1)
                Z += Float64(Orb[a].j + 1)
            end
        else
            break
        end
    end

    # Fill the neutron density matrix ...
    @inbounds for a in 1:a_max
        if (N - N_target) < 1e-3
            if ((N_target - N) - Float64(Orb[a].j + 1)) > 1e-7
                nRho[a,a] = 1.0
                N += Float64(Orb[a].j + 1)
            elseif ((N_target - N) - Float64(Orb[a].j + 1)) < 1e-7
                nRho[a,a] = abs(N_target - N) / Float64(Orb[a].j + 1)
                N += Float64(Orb[a].j + 1)
            end
        else
            break
        end
    end

    # Initialize pairing tensors Kappa ...
    pKappa, nKappa = 0.1 .* diagm(ones(Float64,a_max)), 0.1 .* diagm(ones(Float64,a_max))

    return pnMatrix(pRho,nRho), pnMatrix(pKappa,nKappa)
end

function HFB_Allocate_Indices(Params::Parameters,Orb::Vector{NOrb})
    # Read # initialize parameters ...
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1)*(N_max + 2),2)
    ab_count = 0

    # Count how many pairs of (ab) are there ...
    @inbounds for a in 1:a_max
        l_a = Orb[a].l
        j_a = Orb[a].j
        @inbounds for b in 1:a
            l_b = Orb[b].l
            j_b = Orb[b].j
            if j_a == j_b && l_a == l_b
                ab_count += 1
            end
        end
    end

    # Initialite the array ab ...
    ab, ab_count = zeros(Int64,3,2*ab_count), 0

    # Allocate the array ab for H field ...
    @inbounds for a in 1:a_max
        l_a = Orb[a].l
        j_a = Orb[a].j
        @inbounds for b in 1:a
            l_b = Orb[b].l
            j_b = Orb[b].j
            if j_a == j_b && l_a == l_b
                ab_count += 1
                ab[1,ab_count] = a
                ab[2,ab_count] = b
                ab[3,ab_count] = 1
            end
        end
    end

    # Allocate the array ab for Delta field ...
    @inbounds for a in 1:a_max
        l_a = Orb[a].l
        j_a = Orb[a].j
        @inbounds for b in 1:a
            l_b = Orb[b].l
            j_b = Orb[b].j
            if j_a == j_b && l_a == l_b
                ab_count += 1
                ab[1,ab_count] = a
                ab[2,ab_count] = b
                ab[3,ab_count] = 2
            end
        end
    end

    return ab, ab_count
end

function HFB_Allocate(Params::Parameters,Lambda::pnFloat,Rho::pnMatrix,Kappa::pnMatrix,Orb::Vector{NOrb},Orb_NN::NNOrb,Orb_NNN::NNNOrb,T::Matrix{Float64},VNN::NNInt,VNNN::Array{Vector{Vector{Float32}},4})
    # Read parameters ...
    A = Params.Calc.A
    N_max = Params.Calc.Nmax
    N_2max = Params.Calc.N2max
    N_3max = Params.Calc.N3max
    CMS = Params.Calc.CMS
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Allocate new HF Hamiltonian matrices ...
    pH, nH = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)
    pDelta, nDelta = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)

    # Allocate  indices for HFB iteration ...
    ad, ad_count = HFB_Allocate_Indices(Params,Orb)

    # Allocate the fields H & Delta ...
    @inbounds Threads.@threads for ad_i in 1:ad_count
        a, d, Type = ad[1,ad_i], ad[2,ad_i], ad[3,ad_i]
        n_a, l_a, j_a = Orb[a].n, Orb[a].l, Orb[a].j
        n_d, l_d, j_d = Orb[d].n, Orb[d].l, Orb[d].j

        j_a_hat = sqrt(Float64(j_a + 1))
        s_j_a_hat = Float64(j_a + 1)
        is_j_a_hat = 1.0 / (Float64(j_a) + 1.0)

        # Single-particle field H ...
        if Type == 1
            pHSum, nHSum = 0.0, 0.0

            @inbounds for b in 1:a_max
                n_b = Orb[b].n
                l_b = Orb[b].l
                if (2*(n_a + n_b) + l_a + l_b) <= N_2max
                    j_b = Orb[b].j
                    @inbounds for e in 1:a_max
                        n_e = Orb[e].n
                        l_e = Orb[e].l
                        j_e = Orb[e].j

                        # Normal Density-dependent part ...
                        if l_b == l_e && j_b == j_e && (2*(n_d + n_e) + l_d + l_e) <= N_2max
                            pRho_be, nRho_be = Rho.p[b,e], Rho.n[b,e]

                            @inbounds for J = div(abs(j_d - j_e),2):div(j_d + j_e,2)
                                
                                # 2-body NN interaction part ...
                                if rem(l_a + l_b, 2) == rem(l_d + l_e, 2)
                                    J_j_a_hat = Float64(2*J + 1) * is_j_a_hat
                
                                    pHSum += J_j_a_hat * V2B(a,b,d,e,J,1,VNN.pp,Orb,Orb_NN) * pRho_be
                                    pHSum += J_j_a_hat * V2B(a,b,d,e,J,0,VNN.pn,Orb,Orb_NN) * nRho_be
                                    nHSum += J_j_a_hat * V2B(a,b,d,e,J,1,VNN.nn,Orb,Orb_NN) * nRho_be
                                    nHSum += J_j_a_hat * V2B(b,a,e,d,J,0,VNN.pn,Orb,Orb_NN) * pRho_be

                                end

                                # 3-body NNN interaction part ...
                                @inbounds for c in 1:a_max
                                    n_c = Orb[c].n
                                    l_c = Orb[c].l
                                    if (2*(n_a + n_b + n_c) + l_a + l_b + l_c) <= N_3max
                                        P = rem(l_a + l_b + l_c, 2) + 1
                                        j_c = Orb[c].j
                                        @inbounds for f in 1:a_max
                                            n_f = Orb[f].n
                                            l_f = Orb[f].l
                                            if (2*(n_d + n_e + n_f) + l_d + l_e + l_f) <= N_3max && l_c == l_f && P == (rem(l_d + l_e + l_f,2) + 1)
                                                j_f = Orb[f].j
                                                if j_c == j_f
                                                    pRho_cf, nRho_cf = Rho.p[c,f], Rho.n[c,f]

                                                    ME001 = V3B_NO2B(a,b,c,0,d,e,f,0,J,1,P,VNNN,Orb,Orb_NNN)
                                                    ME101 = V3B_NO2B(a,b,c,1,d,e,f,0,J,1,P,VNNN,Orb,Orb_NNN)
                                                    ME011 = V3B_NO2B(a,b,c,0,d,e,f,1,J,1,P,VNNN,Orb,Orb_NNN)
                                                    ME111 = V3B_NO2B(a,b,c,1,d,e,f,1,J,1,P,VNNN,Orb,Orb_NNN)
                                                    ME113 = V3B_NO2B(a,b,c,1,d,e,f,1,J,3,P,VNNN,Orb,Orb_NNN)

                                                    pHSum += is_j_a_hat * (0.5*ME113*pRho_be*pRho_cf + 0.25 * (ME001 +
                                                            sqrt(1.0/3.0)*ME101 + sqrt(1.0/3.0)*ME011 + ME111/3.0 + 2.0/3.0*ME113)*nRho_be*nRho_cf +
                                                            1.0/3.0 * (2.0*ME111 + ME113)*pRho_be*nRho_cf)

                                                    nHSum += is_j_a_hat * (0.5*ME113*nRho_be*nRho_cf + 0.25 * (ME001 +
                                                            sqrt(1.0/3.0)*ME101 + sqrt(1.0/3.0)*ME011 + ME111/3.0 + 2.0/3.0*ME113)*pRho_be*pRho_cf +
                                                            1.0/3.0 * (2.0*ME111 + ME113)*nRho_be*pRho_cf)

                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end

                        # Anomal Density-dependent part ... Only 3-body NNN interaction part ...
                        if (2*(n_d + n_e) + l_d + l_e) <= N_2max
                            A_NNN_Amp = 0.25 * sqrt(Float64((j_b + 1) * (j_e + 1))) * is_j_a_hat^2
                            @inbounds for c in 1:a_max
                                n_c = Orb[c].n
                                l_c = Orb[c].l
                                if l_c == l_b && (2*(n_a + n_b + n_c) + l_a + l_b + l_c) <= N_3max
                                    P = rem(l_a + l_b + l_c, 2) + 1
                                    j_c = Orb[c].j
                                    if j_c == j_b
                                        pKappa_cb, nKappa_cb = Kappa.p[c,b], Kappa.n[c,b]
                                        @inbounds for f in 1:a_max
                                            n_f = Orb[f].n
                                            l_f = Orb[f].l
                                            if (2*(n_d + n_e + n_f) + l_d + l_e + l_f) <= N_3max && P == (rem(l_d + l_e + l_f,2) + 1) && l_e == l_f
                                                j_f = Orb[f].j
                                                if j_e == j_f
                                                    pKappa_ef, nKappa_ef = Kappa.p[e,f], Kappa.n[e,f]

                                                    ME111 = V3B_NO2B(b,c,a,1,e,f,d,1,0,1,P,VNNN,Orb,Orb_NNN)
                                                    ME113 = V3B_NO2B(b,c,a,1,e,f,d,1,0,3,P,VNNN,Orb,Orb_NNN)

                                                    pHSum += A_NNN_Amp * (ME113 * pKappa_cb * pKappa_ef +
                                                            1.0 / 3.0 * (2.0 * ME111 + ME113) * nKappa_cb * nKappa_ef)
                                                    nHSum += A_NNN_Amp * (ME113 * nKappa_cb * nKappa_ef +
                                                            1.0 / 3.0 * (2.0 * ME111 + ME113) * pKappa_cb * pKappa_ef)

                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end

                    end
                end
            end

            # Include the 1-body kinetic energy & inclusion of Center-of-Mass motion (CM) correction ...
                # Combined 1- + 2-body kinetic operator with CM correction ...
            if CMS == "CMS1+2B"
                pHSum += T[a,d] * (1.0 - 1.0 / A)
                nHSum += T[a,d] * (1.0 - 1.0 / A)
                # Pure 1-body kinetic operator with no CM correction ...
            elseif CMS != "CMS2B"
                pHSum += T[a,d]
                nHSum += T[a,d]
            end
                # No contribution for pure 2-body kinetic operator with CM correction ...


            # Add 0-body chemical potential Lambda ...
            if a == d
                pHSum -= Lambda.p
                nHSum -= Lambda.n
            end

            # Allocate pH & nH ...
            if a != d
                pH[a,d], pH[d,a] = pHSum, pHSum
                nH[a,d], nH[d,a] = nHSum, nHSum
            elseif a == d
                pH[a,a], nH[a,a] = pHSum, nHSum
            end
        end

        # Single-quasiparticle pairing field Delta ...
        if Type == 2
            if ((2*(n_a + n_d) + l_a + l_d) <= N_2max)
                pDeltaSum, nDeltaSum = 0.0, 0.0

                @inbounds for b in 1:a_max
                    n_b = Orb[b].n
                    l_b = Orb[b].l
                    j_b = Orb[b].j
                    j_b_hat = sqrt(Float64(j_b + 1))
                    NN_Amp = 0.5 * j_b_hat / j_a_hat
                    @inbounds for e in 1:a_max
                        n_e = Orb[e].n
                        l_e = Orb[e].l
                        j_e = Orb[e].j
                        if (l_b == l_e) && (j_b == j_e) && ((2*(n_b + n_e) + l_b + l_e) <= N_2max)
                            pKappa_be, nKappa_be = Kappa.p[b,e], Kappa.n[b,e]

                            # 2-body NN interaction part ...
                            pDeltaSum += NN_Amp * V2B(a,d,b,e,0,1,VNN.pp,Orb,Orb_NN) * pKappa_be
                            nDeltaSum += NN_Amp * V2B(a,d,b,e,0,1,VNN.nn,Orb,Orb_NN) * nKappa_be

                            # 3-body NNN interaction part ...
                            @inbounds for c in 1:a_max
                                n_c = Orb[c].n
                                l_c = Orb[c].l
                                if (2*(n_a + n_d + n_c) + l_a + l_d + l_c) <= N_3max
                                    P = rem(l_a + l_d + l_c, 2) + 1
                                    j_c = Orb[c].j

                                    NNN_Amp = 0.5 * j_b_hat * j_a_hat / Float64(j_c + 1)^2

                                    @inbounds for f in 1:a_max
                                        n_f = Orb[f].n
                                        l_f = Orb[f].l
                                        if (2*(n_b + n_e + n_f) + l_b + l_e + l_f) <= N_3max && l_c == l_f && P == (rem(l_b + l_e + l_f,2) + 1)
                                            j_f = Orb[f].j
                                            if j_e == j_f
                                                pRho_cf, nRho_cf = Rho.p[c,f], Rho.n[c,f]

                                                ME111 = V3B_NO2B(a,d,c,1,b,e,f,1,0,1,P,VNNN,Orb,Orb_NNN)
                                                ME113 = V3B_NO2B(a,d,c,1,b,e,f,1,0,3,P,VNNN,Orb,Orb_NNN)

                                                pDeltaSum += NNN_Amp * (ME113 * pKappa_be * pRho_cf +
                                                            (2.0 * ME111 + ME113) / 3.0 * pKappa_be * nRho_cf)
                                                nDeltaSum += NNN_Amp * (ME113 * nKappa_be * nRho_cf +
                                                            (2.0 * ME111 + ME113) / 3.0 * nKappa_be * pRho_cf)
            
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end

                if a != d
                    pDelta[a,d], pDelta[d,a] = pDeltaSum, pDeltaSum
                    nDelta[a,d], nDelta[d,a] = nDeltaSum, nDeltaSum
                elseif a == d
                    pDelta[a,d] = pDeltaSum
                    nDelta[a,d] = nDeltaSum
                end

            end
        end

    end

    return pnMatrix(pH,nH), pnMatrix(pDelta,nDelta)
end

function HFB_Diagonalize(Params::Parameters,H::pnMatrix,Delta::pnMatrix,Orb::Vector{NOrb})
    # Read parameters ...
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1) * (N_max + 2), 2)

    # Preallocate arrays for solutions ...
    pU, pV, pSQE = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max), zeros(Float64,a_max)
    nU, nV, nSQE = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max), zeros(Float64,a_max)

    # Allocate the HFB eigenvalue system ...
    pHFB = [H.p Delta.p; Delta.p -1.0 * H.p]
    nHFB = [H.n Delta.n; Delta.n -1.0 * H.n]

    # Solve the HFB equations ...
    pE, pC = eigen(Symmetric(pHFB), sortby = +)
    nE, nC = eigen(Symmetric(nHFB), sortby = +)

    # Only positive energy solutions are extracted ...
    @inbounds for a in 1:a_max
        pSQE[a], nSQE[a] = pE[a+a_max], nE[a+a_max]
        @views pU[:,a] .= pC[1:a_max,a+a_max]
        @views pV[:,a] .= pC[a_max+1:2*a_max,a+a_max]
        @views nU[:,a] .= nC[1:a_max,a+a_max]
        @views nV[:,a] .= nC[a_max+1:2*a_max,a+a_max]
    end

    # Perform reordering - to match quantum numbers j & l & ascending in SQE ...
    SQE, U, V = HFB_Orbital_Ordering(Params,pnVector(pSQE,nSQE),pnMatrix(pU,nU),pnMatrix(pV,nV),Orb)

    return SQE, U, V
end

function HFB_Orbital_Ordering(Params::Parameters,SQE::pnVector,U::pnMatrix,V::pnMatrix,Orb::Vector{NOrb})
    # Read parameters ...
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1) * (N_max + 2), 2)

    # Read input arrays ...
    pU, pV, pSQE =  U.p, V.p, SQE.p
    nU, nV, nSQE =  U.n, V.n, SQE.n

    # Preallocate temporary arrays ...
    pOrb_order, nOrb_order = Vector{Int64}(undef,a_max), Vector{Int64}(undef,a_max)
    pOrb_mask, nOrb_mask = falses(a_max), falses(a_max)

    pl_values, pj_values = zeros(Float64,a_max), zeros(Float64,a_max)
    nl_values, nj_values = zeros(Float64,a_max), zeros(Float64,a_max)

    # Evaluate values of j & l for single-quasiparticle orbitals ...
    @inbounds for a in 1:a_max
        pjSum, plSum = 0.0, 0.0
        njSum, nlSum = 0.0, 0.0
        @inbounds for b in 1:a_max
            l_b, j_b = Float64(Orb[b].l), Float64(Orb[b].j)
            pN, nN = pU[b,a]^2 + pV[b,a]^2, nU[b,a]^2 + nV[b,a]^2
            pjSum, plSum = pjSum + j_b * pN, plSum + l_b * pN
            njSum, nlSum = njSum + j_b * nN, nlSum + l_b * nN
        end
        pl_values[a], pj_values[a] = plSum, pjSum
        nl_values[a], nj_values[a] = nlSum, njSum
    end

    # Find the ordering for single-quasiparticle orbitals ...
    @inbounds for a in 1:a_max
        pl, pj = pl_values[a], pj_values[a]
        nl, nj = nl_values[a], nj_values[a]
        @inbounds for b in 1:a_max
            if pOrb_mask[b] == false && abs(Float64(Orb[b].l) - pl) < 1e-7 && abs(Float64(Orb[b].j) - pj) < 1e-7
                pOrb_order[b] = a
                pOrb_mask[b] = true
                break
            end
        end

        @inbounds for b in 1:a_max
            if nOrb_mask[b] == false && abs(Float64(Orb[b].l) - nl) < 1e-7 && abs(Float64(Orb[b].j) - nj) < 1e-7
                nOrb_order[b] = a
                nOrb_mask[b] = true
                break
            end
        end
    end

    # Perform reordering of single-quasiparticle orbitals ...
    @views pSQE .= pSQE[pOrb_order]
    @views pU .= pU[:,pOrb_order]
    @views pV .= pV[:,pOrb_order]

    @views nSQE .= nSQE[nOrb_order]
    @views nU .= nU[:,nOrb_order]
    @views nV .= nV[:,nOrb_order]

    return pnVector(pSQE,nSQE), pnMatrix(pU,nU), pnMatrix(pV,nV)
end

function HFB_Density_Operator(Params::Parameters,U::pnMatrix,V::pnMatrix,Orb::Vector{NOrb})
    # Read parameters ...
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1) * (N_max + 2), 2)

    # Preallocate matrices for Rho & Kappa density operators ...
    pRho, nRho = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)
    pKappa, nKappa = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)

    # Allocate density operators Rho & Kappa ...
    @inbounds for a in 1:a_max
        @inbounds for b in 1:a_max
            if (Orb[a].j == Orb[b].j) && (Orb[a].l == Orb[b].l)
                pRhoSum, pKappaSum = 0.0, 0.0
                nRhoSum, nKappaSum = 0.0, 0.0
                @inbounds for c in 1:a_max
                    if (Orb[a].j == Orb[c].j) && (Orb[a].l == Orb[c].l)
                        pMERho, pMEKappa = V.p[a,c] * V.p[b,c], V.p[a,c] * U.p[b,c]
                        nMERho, nMEKappa = V.n[a,c] * V.n[b,c], V.n[a,c] * U.n[b,c]
                        pRhoSum, pKappaSum = pRhoSum + pMERho, pKappaSum + pMEKappa
                        nRhoSum, nKappaSum = nRhoSum + nMERho, nKappaSum + nMEKappa
                    end
                end
                pRho[a,b], pKappa[a,b] = pRhoSum, pKappaSum
                nRho[a,b], nKappa[a,b] = nRhoSum, nKappaSum
            end
        end
    end

    return pnMatrix(pRho,nRho), pnMatrix(pKappa,nKappa)
end

function HFB_Particle_Number(Params::Parameters,Rho::pnMatrix,Orb::Vector{NOrb})
    # Read parameters ...
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1) * (N_max + 2), 2)

    # Initialize particle numbers ...
    Z, N = 0,0, 0,0

    # Evaluate the angular-momentum weighted trace of density matrices ...
    @inbounds for a in 1:a_max
        Z += Rho.p[a,a] * Float64(Orb[a].j + 1)
        N += Rho.n[a,a] * Float64(Orb[a].j + 1)
    end

    return Z, N
end

function HFB_Particle_Number_Dispersion(Params::Parameters,Rho::pnMatrix,Orb::Vector{NOrb})
    # Read parameters ...
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1)*(N_max + 2),2)
    dZ, dN = 0.0, 0.0

    # Calculate the particle number dispersion ...
    println("\nCalculating the HFB dispersion of proton & neutron particle numbers ...")

    @inbounds for a in 1:a_max
        j_a, l_a = Orb[a].j, Orb[a].l
        j_a_hat = Float64(j_a + 1)
        dZ += 2.0 * j_a_hat * Rho.p[a,a] * (1.0 - Rho.p[a,a])
        dN += 2.0 * j_a_hat * Rho.n[a,a] * (1.0 - Rho.n[a,a])
    end

    # Calculate square roots of dispersion numbers ...
    dZ, dN = sqrt(dZ), sqrt(dN)

    println("\nHFB dispersion of nucleons numbers are ...")
    println("dZ = " * string(round(dZ,digits=5)))
    println("dN = " * string(round(dN,digits=5)))

    return pnFloat(dZ,dN)
end

function HFB_Energy(Params::Parameters,Rho::pnMatrix,Kappa::pnMatrix,Orb::Vector{NOrb},Orb_NN::NNOrb,Orb_NNN::NNNOrb,T::Matrix{Float64},VNN::NNInt,VNNN::Array{Vector{Vector{Float32}},4})
    # Read parameters ...
    A = Params.Calc.A
    N_max = Params.Calc.Nmax
    N_2max = Params.Calc.N2max
    N_3max = Params.Calc.N3max
    CMS = Params.Calc.CMS
    a_max = div((N_max + 1)*(N_max + 2),2)

    #Calculate the HF energy ...

    E_HFB = 0.0

    println("\nCalculating total HFB ground-state energy ...")

    E_h_partial = Threads.Atomic{Float64}[Threads.Atomic{Float64}(0.0) for _ in 1:Threads.nthreads()]
    E_Delta_partial = Threads.Atomic{Float64}[Threads.Atomic{Float64}(0.0) for _ in 1:Threads.nthreads()]

    # Single-particle mean-field energy ...
    @inbounds Threads.@threads for a = 1:a_max
        thread_id = Threads.threadid()
        n_a = Orb[a].n
        l_a = Orb[a].l
        j_a = Orb[a].j
        @inbounds for b = 1:a_max
            n_b = Orb[b].n
            l_b = Orb[b].l
            if (2*(n_a + n_b) + l_a + l_b) <= N_2max
                j_b = Orb[b].j
                @inbounds for d = 1:a_max
                    l_d = Orb[d].l
                    j_d = Orb[d].j
                    if (l_a == l_d) && (j_a == j_d)
                        n_d = Orb[d].n
                        pRho_ad = Rho.p[a,d]
                        nRho_ad = Rho.n[a,d]
                        @inbounds for e = 1:a_max
                            n_e = Orb[e].n
                            l_e = Orb[e].l
                            j_e = Orb[e].j
                            if (l_b == l_e && j_b == j_e) && ((2*(n_d + n_e) + l_d + l_e) <= N_2max)
                                pRho_be = Rho.p[b,e]
                                nRho_be = Rho.n[b,e]
                                @inbounds for J = div(abs(j_a - j_b),2):div((j_a + j_b),2)

                                    # 2-body NN interaction ...
                                    if (rem(l_a + l_b, 2) == rem(l_d + l_e, 2))

                                        Hat =  Float64(2*J + 1)
                                        @views E_h_partial[thread_id][] += 0.5 * Hat * pRho_ad * pRho_be * V2B(a,b,d,e,J,1,VNN.pp,Orb,Orb_NN)
                                        @views E_h_partial[thread_id][] += 0.5 * Hat * nRho_ad * nRho_be * V2B(a,b,d,e,J,1,VNN.nn,Orb,Orb_NN)
                                        @views E_h_partial[thread_id][] += Hat * pRho_ad * nRho_be * V2B(a,b,d,e,J,0,VNN.pn,Orb,Orb_NN)

                                    end
                                    
                                    # 3-body NNN interaction ...
                                    @inbounds for c = 1:a_max
                                        n_c = Orb[c].n
                                        l_c = Orb[c].l
                                        if ((2*(n_a + n_b + n_c) + l_a + l_b + l_c) <= N_3max)
                                            j_c = Orb[c].j
                                            P = rem(l_a + l_b + l_c, 2) + 1
                                            @inbounds for f = 1:a_max
                                                n_f = Orb[f].n
                                                l_f = Orb[f].l
                                                j_f = Orb[f].j
                                                if (2*(n_d + n_e + n_f) + l_d + l_e + l_f) <= N_3max && l_c == l_f && j_c == j_f && P == (rem(l_d + l_e + l_f, 2)+1)
                                                    pRho_cf, nRho_cf = Rho.p[c,f], Rho.n[c,f]

                                                    me1 = V3B_NO2B(a,b,c,0,d,e,f,0,J,1,P,VNNN,Orb,Orb_NNN)
                                                    me2 = V3B_NO2B(a,b,c,1,d,e,f,0,J,1,P,VNNN,Orb,Orb_NNN)
                                                    me3 = V3B_NO2B(a,b,c,0,d,e,f,1,J,1,P,VNNN,Orb,Orb_NNN)
                                                    me4 = V3B_NO2B(a,b,c,1,d,e,f,1,J,1,P,VNNN,Orb,Orb_NNN)
                                                    me5 = V3B_NO2B(a,b,c,1,d,e,f,1,J,3,P,VNNN,Orb,Orb_NNN)

                                                    @views E_h_partial[thread_id][] += 1.0/6.0 * (me5 * pRho_ad * pRho_be * pRho_cf + (2.0 * me4 + me5) *
                                                    pRho_ad * pRho_be * nRho_cf + (1.5 * me1 + sqrt(3.0/4.0) * me2 +
                                                    sqrt(3.0/4.0) * me3 + 0.5 * me4 + me5)  * pRho_ad * nRho_be * nRho_cf +
                                                    me5 * nRho_ad * nRho_be * nRho_cf)
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            if CMS == "CMS1+2B"
                @views E_h_partial[thread_id][] += (Rho.p[a,b] + Rho.n[a,b]) * T[a,b] * Float64(j_a + 1) * (1.0 - 1.0 / Float64(A))
            elseif CMS == "CMS2B"
                @views E_h_partial[thread_id][] += 0
            else
                @views E_h_partial[thread_id][] += (Rho.p[a,b] + Rho.n[a,b]) * T[a,b] * Float64(j_a + 1)
            end
        end
    end

    # Pairing field energy ...
    @inbounds Threads.@threads for a = 1:a_max
        thread_id = Threads.threadid()
        n_a = Orb[a].n
        l_a = Orb[a].l
        j_a = Orb[a].j
        j_a_hat = sqrt(Float64(Orb[a].j + 1)) 
        @inbounds for b = 1:a_max
            n_b = Orb[b].n
            l_b = Orb[b].l
            if (2*(n_a + n_b) + l_a + l_b) <= N_2max && l_a == l_b
                j_b = Orb[b].j
                if j_a == j_b
                    pKappa_ba, nKappa_ba = Kappa.p[b,a], Kappa.n[b,a]
                    @inbounds for d = 1:a_max
                        n_d = Orb[d].n
                        l_d = Orb[d].l
                        j_d = Orb[d].j
                        j_d_hat = sqrt(Float64(Orb[d].j + 1)) 
                        @inbounds for e = 1:a_max
                            n_e = Orb[e].n
                            l_e = Orb[e].l
                            if ((2*(n_d + n_e) + l_d + l_e) <= N_2max) && l_d == l_e
                                j_e = Orb[e].j
                                if j_d == j_e
                                    pKappa_de, nKappa_de = Kappa.p[d,e], Kappa.n[d,e]

                                    # 2-body NN interaction ...
                                    @views E_Delta_partial[thread_id][] += 0.25 * j_a_hat * j_d_hat * pKappa_ba * pKappa_de * V2B(a,b,d,e,0,1,VNN.pp,Orb,Orb_NN)
                                    @views E_Delta_partial[thread_id][] += 0.25 * j_a_hat * j_d_hat * nKappa_ba * nKappa_de * V2B(a,b,d,e,0,1,VNN.nn,Orb,Orb_NN)
                                
                                    # 3-body NNN interaction ...
                                    @inbounds for c = 1:a_max
                                        n_c = Orb[c].n
                                        l_c = Orb[c].l
                                        if ((2*(n_a + n_b + n_c) + l_a + l_b + l_c) <= N_3max)
                                            j_c = Orb[c].j
                                            j_c_hat = Float64(Orb[c].j + 1)
                                            P = rem(l_a + l_b + l_c, 2) + 1
                                            @inbounds for f = 1:a_max
                                                n_f = Orb[f].n
                                                l_f = Orb[f].l
                                                j_f = Orb[f].j
                                                if (2*(n_d + n_e + n_f) + l_d + l_e + l_f) <= N_3max && l_c == l_f && j_c == j_f && P == (rem(l_d + l_e + l_f, 2)+1)
                                                    pRho_cf, nRho_cf = Rho.p[c,f], Rho.n[c,f]

                                                    ME111 = V3B_NO2B(a,b,c,1,d,e,f,1,0,1,P,VNNN,Orb,Orb_NNN)
                                                    ME113 = V3B_NO2B(a,b,c,1,d,e,f,1,0,3,P,VNNN,Orb,Orb_NNN)

                                                    @views E_Delta_partial[thread_id][] += 0.25 * j_a_hat * j_d_hat / j_c_hat * (pKappa_ba * pKappa_de *
                                                                                            pRho_cf * ME113 + pKappa_ba * pKappa_de * nRho_cf * 1.0 / 3.0 *
                                                                                            (2.0 * ME111 + ME113) + nKappa_ba * nKappa_de * nRho_cf * ME113 +
                                                                                            nKappa_ba * nKappa_de * pRho_cf * 1.0 / 3.0 *(2.0 * ME111 + ME113))
                                                end
                                            end
                                        end
                                    end

                                end
                            end
                        end
                    end

                end
            end
        end
    end

    E_HFB += sum(x[] for x in E_h_partial)
    E_HFB += sum(x[] for x in E_Delta_partial)

    println("\nHFB energy    ...   E_HFB = " * string(E_HFB) * " MeV")

    return E_HFB
end
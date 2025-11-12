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
    if Params.Calc.HFB.Broyden == true
        @time Lambda, SQE, U, V, SQE_C, u_C, v_C, C, Rho, Kappa, H, Delta, Iteration = HFB_Solve_Broyden(Params,Orb,Orb_NN,Orb_NNN,T,VNN,VNNN,epsilon)
    else
        @time Lambda, SQE, U, V, SQE_C, u_C, v_C, C, Rho, Kappa, H, Delta, Iteration = HFB_Solve(Params,Orb,Orb_NN,Orb_NNN,T,VNN,VNNN,epsilon)
    end

    # Calculation of the total HFB mean-field ground-state energy ...
    @time E_HFB = HFB_Energy(Params,Rho,Kappa,Orb,Orb_NN,Orb_NNN,T,VNN,VNNN)

    # Calculate the total HFB ground-state kinetic energy ...
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
    Iteration, Iteration_max = 0, 150
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
        @inbounds for L in 1:Iteration_max
            # Initialization of the secant method ...
            if (L == 1) && (Iteration == 0)
                # Allocate the single-particle fields H_1, H_2 and pairing fields Delta_1, Delta_2 ...
                H_1, Delta_1 = HFB_Allocate(Params,Lambda_1,Rho_old,Kappa_old,Orb,Orb_NN,Orb_NNN,T,VNN,VNNN)
                H_2, Delta_2 = HFB_Allocate(Params,Lambda_2,Rho_old,Kappa_old,Orb,Orb_NN,Orb_NNN,T,VNN,VNNN)

                # Diagonalize the HFB equations ... basis is reordered as needed ...
                SQE_1, U_1, V_1 = HFB_Diagonalize(Params,H_1,Delta_1,U,V,Orb)
                SQE_2, U_2, V_2 = HFB_Diagonalize(Params,H_2,Delta_2,U_1,V_1,Orb)
                
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

            # Diagonalize the HFB equations ... basis is reordered as needed ...
            SQE, U, V = HFB_Diagonalize(Params,H,Delta,U,V,Orb)
            
            # Generate temporary densities Rho & Kappa ...
            Rho_temp, Kappa_temp = HFB_Density_Operator(Params,U,V,Orb)

            # Determine new average particle numbers ...
            Z, N = HFB_Particle_Number(Params,Rho_temp,Orb)

            # Update average particle numbers ...
            Z_1, N_1 = Z_2, N_2
            Z_2, N_2 = Z, N

            # Perform secant iterations to determine new optimal value of Lambda ...
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
        Rho, Kappa = HFB_Mixing_Update(a_max,Rho,Rho_old,Kappa,Kappa_old)

        # Allocate the single-particle field H and the pairing field Delta ...
        H, Delta = HFB_Allocate(Params,Lambda,Rho,Kappa,Orb,Orb_NN,Orb_NNN,T,VNN,VNNN)

        # Diagonalize the HFB equations ... basis is reordered as needed ...
        SQE, U, V = HFB_Diagonalize(Params,H,Delta,U,V,Orb)

        # Evaluate particle numbers ...
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

    # Perform final evaluation of resulting densities ...
    Rho, Kappa = HFB_Density_Operator(Params,U,V,Orb)

    println("\nHFB iteration with residual NN interaction has converged ...")

    # Construct the canonical basis & evaluate approximate (BCS-like) amplitudes
    # & single-quasiparticle energies u_C, v_C & SQE_C ...
    #   U, V, Rho, Kappa, H, Delta ... remain expressed in the reference LHO basis ...
    SQE_C, C, u_C, v_C = HFB_Canonical_Basis(Params,Rho,H,Delta,Orb)

    return Lambda, SQE, U, V, SQE_C, u_C, v_C, C, Rho, Kappa, H, Delta, Iteration
end

function HFB_Solve_Broyden(Params::Parameters,Orb::Vector{NOrb},Orb_NN::NNOrb,Orb_NNN::NNNOrb,T::Matrix{Float64},VNN::NNInt,VNNN::Array{Vector{Vector{Float32}},4},epsilon::Float64)
    # Read calculation parameters ...
    A_target, Z_target, N_target = Float64(Params.Calc.A), Float64(Params.Calc.Z), Float64(Params.Calc.A - Params.Calc.Z)
    hw = Params.Int.hw
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1) * (N_max + 2), 2)

    # Setup local iteration variables ...
    Iteration, Iteration_max = 0, 150
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

    # Initialite the Broyden's method ...
    Broyden, BroyVec = HFB_Broyden_Initialize(Params,Rho,Kappa,Orb)

    # Preallocate SQE_old ...
    SQE_old = pnVector(SQE.p, SQE.n)

    # Initialize particle numbers & chemical potentials for secant method ...
    Z_1, Z_2 = 0.0, 0.0
    N_1, N_2 = 0.0, 0.0
    Lambda_1, Lambda_2 = pnFloat(0.75 * Lambda.p, 0.75 * Lambda.n), pnFloat(Lambda.p, Lambda.n)

    # Solve the spherical HFB equations ... by the means of self-consistent iteration ...
    println("\nStarting iteration of HFB equations ...\n")

    while ((dE > epsilon) || (dZ > epsilon) || (dN > epsilon)) && (Iteration < Iteration_max)
        # Perform several secant iterations for the chemical potentil Lambda ...
        @inbounds for L in 1:Iteration_max
            # Initialization of the secant method ...
            if (L == 1) && (Iteration == 0)
                # Allocate the single-particle fields H_1, H_2 and pairing fields Delta_1, Delta_2 ...
                H_1, Delta_1 = HFB_Allocate(Params,Lambda_1,Rho,Kappa,Orb,Orb_NN,Orb_NNN,T,VNN,VNNN)
                H_2, Delta_2 = HFB_Allocate(Params,Lambda_2,Rho,Kappa,Orb,Orb_NN,Orb_NNN,T,VNN,VNNN)

                # Diagonalize the HFB equations ... basis is reordered as needed ...
                SQE_1, U_1, V_1 = HFB_Diagonalize(Params,H_1,Delta_1,U,V,Orb)
                SQE_2, U_2, V_2 = HFB_Diagonalize(Params,H_2,Delta_2,U_1,V_1,Orb)
                
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

            # Diagonalize the HFB equations ... basis is reordered as needed ...
            SQE, U, V = HFB_Diagonalize(Params,H,Delta,U,V,Orb)
            
            # Generate temporary densities Rho & Kappa ...
            Rho_temp, Kappa_temp = HFB_Density_Operator(Params,U,V,Orb)

            # Determine new average particle numbers ...
            Z, N = HFB_Particle_Number(Params,Rho_temp,Orb)

            # Update average particle numbers ...
            Z_1, N_1 = Z_2, N_2
            Z_2, N_2 = Z, N

            # Perform secant iterations to determine new optimal value of Lambda ...
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
        Rho, Kappa, BroyVec = HFB_Broyden_Update(Params,Iteration,Rho,Kappa,Broyden,BroyVec,Orb)

        # Allocate the single-particle field H and the pairing field Delta ...
        H, Delta = HFB_Allocate(Params,Lambda,Rho,Kappa,Orb,Orb_NN,Orb_NNN,T,VNN,VNNN)

        # Diagonalize the HFB equations ... basis is reordered as needed ...
        SQE, U, V = HFB_Diagonalize(Params,H,Delta,U,V,Orb)

        # Evaluate particle numbers ...
        Z, N = HFB_Particle_Number(Params,Rho,Orb)

        # Evaluate iteration of the single-quasiparticle energies ...
        dE = (sum(abs.(SQE.p .- SQE_old.p )) + sum(abs.(SQE.n .- SQE_old.n))) / Float64(2 * a_max)
        dZ, dN = abs(Z_target - Z), abs(N_target - N)

        Iteration += 1

        println("\n\nHFB iteration number:   " * string(Iteration) * "   Single-quasiparticle energy difference:   " * string(round(dE, sigdigits=8))
                * "   Proton number difference:   " * string(round(dZ, sigdigits=8)) * "   Neutron number difference:   " * string(round(dN, sigdigits=8)))

        # Store old values of SQE ...
        SQE_old = pnVector(SQE.p, SQE.n)

    end

    # Perform final evaluation of resulting densities ...
    Rho, Kappa = HFB_Density_Operator(Params,U,V,Orb)

    println("\nHFB iteration with residual NN interaction has converged ...")

    # Construct the canonical basis & evaluate approximate (BCS-like) amplitudes
    # & single-quasiparticle energies u_C, v_C & SQE_C ...
    #   U, V, Rho, Kappa, H, Delta ... remain expressed in the reference LHO basis ...
    SQE_C, C, u_C, v_C = HFB_Canonical_Basis(Params,Rho,H,Delta,Orb)

    return Lambda, SQE, U, V, SQE_C, u_C, v_C, C, Rho, Kappa, H, Delta, Iteration
end

function HFB_Diagonalize(Params::Parameters,H::pnMatrix,Delta::pnMatrix,U::pnMatrix,V::pnMatrix,Orb::Vector{NOrb})
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
    SQE, U_new, V_new = HFB_Orbital_Ordering(Params,pnVector(pSQE,nSQE),pnMatrix(pU,nU),pnMatrix(pV,nV),Orb)

    # Set phases of U & V to match with previous iteration ...
    @inbounds for a in 1:a_max
        pU_Overlap = dot(U_new.p[:,a],U.p[:,a])
        pV_Overlap = dot(V_new.p[:,a],V.p[:,a])
        nU_Overlap = dot(U_new.n[:,a],U.n[:,a])
        nV_Overlap = dot(V_new.n[:,a],V.n[:,a])
        if pU_Overlap < -1e-10; @views U_new.p[:,a] .= -1.0 .* U_new.p[:,a]; end
        if pV_Overlap < -1e-10; @views V_new.p[:,a] .= -1.0 .* V_new.p[:,a]; end
        if nU_Overlap < -1e-10; @views U_new.n[:,a] .= -1.0 .* U_new.n[:,a]; end
        if nV_Overlap < -1e-10; @views V_new.n[:,a] .= -1.0 .* V_new.n[:,a]; end
        #if pU_Overlap + pV_Overlap < -1e-8; @views U_new.p[:,a] .= -1.0 .* U_new.p[:,a]; @views V_new.p[:,a] .= -1.0 .* V_new.p[:,a]; end
        #if nU_Overlap + nV_Overlap < -1e-8; @views U_new.n[:,a] .= -1.0 .* U_new.n[:,a]; @views V_new.n[:,a] .= -1.0 .* V_new.n[:,a]; end
    end

    return SQE, U_new, V_new
end
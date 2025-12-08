function HFB_Solver(Params::Parameters)
    # Iteration precision parameter ...
    epsilon = 1e-7

    # Make single-particle orbitals - NuHamil ordering ...
    Orb = Make_Orbitals(Params.Calc.A,Params.Calc.Z,Params.Int.Nmax)

    # Load 1-body kinetic operator ...
    T = T1B(Params.Int.Nmax,Orb,Params.Int.hw)

    # 2-body NN interaction & Orbitals ...
    @time V_NN, Orb_NN = V2B_Read(Params,Orb)

    # 3-body NNN interaction & Orbitals ...
    @time V_NNN, Orb_NNN = V3B_NO2B_Read(Params,Orb)

    # Solve HFB equations ...
    @time Lambda, SQE, U, V, SQE_C, u_C, v_C, C, Rho, Kappa, H, Delta, Iteration = HFB_Solve(Params,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN,epsilon)

    # Calculation of the total HFB mean-field ground-state energy ...
    @time E_HFB = HFB_Energy(Params,Rho,Kappa,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN)

    # Calculate the total HFB ground-state kinetic energy ...
    @time T_HFB = Kinetic_Energy(Params,Rho,Orb,T)

    # Particle number fluctuation calculation ...
    @time dA = HFB_Particle_Number_Dispersion(Params,Rho,Orb)

    # Calculation summary ...
    @time HFB_Summary(Params,E_HFB,T_HFB,Lambda,dA,epsilon,Iteration)

    # Evaluate HFB charge radii & radial densities ...
    Summary_File = "IO/" * Params.Calc.Path * "/HFB/HFB_Summary.dat"
    Densities_File = "IO/" * Params.Calc.Path * "/HFB/Densities/HFB_Radial_Densities.dat"
    @time OBDM_Export(Params,Summary_File,Densities_File,Rho,C,Orb)

    # Export of single-quasiparticle energies, amplitudes U & V & also possibly radial densities ...
    @time HFB_SQS_Summary(Params,SQE,SQE_C,pnMatrix(C.p' * Rho.p * C.p,C.n' * Rho.n * C.n),Orb)

    return

    # Allocate H1B ... 1-body HFB Hamiltonian in the quasiparticle basis ...
    @time H_N = HFB_allocate_H1b(Params,SQE)

    # Make density dependent residual 2-body interaction in the LHO basis ...
    @time V_NN = H2b_res_no2b(Params,Orb,Orb_NN,Orb_NNN,V_NN,V_NNN,Rho)

    # Deallocate V_NNN 3-body interaction ...
    V_NNN = nothing

    # Perform the Garbage Collection ...
    GC.gc()

    # Perform transformation of U and V to the canonical basis ...
    U, V = pnMatrix(C.p' * U.p, C.n' * U.n), pnMatrix(C.p' * V.p, C.n' * V.n)

    # Allocate H2B ... residual interaction Hamiltonian in the quasiparticle basis ...
    @time H_NN = qpH2b(Params,Orb,Orb_NN,V_NN,C,U,V)

    # Perform final export of HFB solution into binary files ...
    @time HFB_Export(Params,Orb_NN,C,U,V,H_N,H_NN)

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
    dE, d2E, dZ, dN = 0.1, 0.1, 0.1, 0.1
    Degeneracy = false

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
        # Arrays for degenerate solutions ...
    Rho1 = pnMatrix(zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max))
    Kappa1 = pnMatrix(zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max))
    Lambda1 = pnFloat(0.0,0.0)
    Rho2 = pnMatrix(zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max))
    Kappa2 = pnMatrix(zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max))
    Lambda2 = pnFloat(0.0,0.0)

    # Setup particle numbers ...
    Z, N = 0.1, 0.1

    # Initial guess on chemical potentials lambda ...
        # By default set to -5 MeV ... can be adjusted in Pairing() parameters ...
    Lambda = pnFloat(Params.Calc.Pairing.pL0, Params.Calc.Pairing.nL0)

    # Initial guess on densities Rho & Kappa
    Rho, Kappa = HFB_Density_Operator_Initialize(Params,Orb)

    # Initialite the Broyden's method ...
    Broyden, BroyVec = HFB_Broyden_Initialize(Params,Rho,Kappa,Orb)

    # Preallocate SQE_old ...
    SQE_old = pnVector(SQE.p,SQE.n)

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

        #display(diag(Kappa.n)[1:10])

        # Allocate the single-particle field H and the pairing field Delta ...
        H, Delta = HFB_Allocate(Params,Lambda,Rho,Kappa,Orb,Orb_NN,Orb_NNN,T,VNN,VNNN)

        # Diagonalize the HFB equations ... basis is reordered as needed ...
        SQE, U, V = HFB_Diagonalize(Params,H,Delta,U,V,Orb)

        # Evaluate particle numbers ...
        Z, N = HFB_Particle_Number(Params,Rho,Orb)

        # Evaluate iteration of the single-quasiparticle energies ...
        d2E = abs(dE - (sum(abs.(SQE.p .- SQE_old.p )) + sum(abs.(SQE.n .- SQE_old.n))) / Float64(2 * a_max))
        dE = (sum(abs.(SQE.p .- SQE_old.p )) + sum(abs.(SQE.n .- SQE_old.n))) / Float64(2 * a_max)
        dZ, dN = abs(Z_target - Z), abs(N_target - N)
        Iteration += 1

        println("\n\nHFB iteration number:   " * string(Iteration) * "   dE = " * string(round(dE, sigdigits=8))
                * " MeV,   dZ = " * string(round(dZ, sigdigits=8)) * ",   dN = " * string(round(dN, sigdigits=8)))

        # Store old values of SQE ...
        SQE_old = pnVector(SQE.p,SQE.n)

        # Resolve the problem of degenerate HFB solutions if they appear ...
        if Degeneracy == true
            Rho2, Kappa2 = pnMatrix(Rho.p,Rho.n), pnMatrix(Kappa.p,Kappa.n)
            Lambda2 = pnFloat(Lambda.p,Lambda.n)
            # Calculate the energies of degenerate HFB solutions ...
            println("\nCalculating the energies of degenerate HFB solutions ...")
            E1_HFB = HFB_Energy(Params,Rho1,Kappa1,Orb,Orb_NN,Orb_NNN,T,VNN,VNNN)
            E2_HFB = HFB_Energy(Params,Rho2,Kappa2,Orb,Orb_NN,Orb_NNN,T,VNN,VNNN)

            # Select the better degenerate solution ... by energy ...
            if E1_HFB < E2_HFB
                Rho, Kappa, Lambda = pnMatrix(Rho1.p,Rho1.n), pnMatrix(Kappa1.p,Kappa1.n), pnFloat(Lambda1.p,Lambda1.n)
            else
                Rho, Kappa, Lambda = pnMatrix(Rho2.p,Rho2.n), pnMatrix(Kappa2.p,Kappa2.n), pnFloat(Lambda2.p,Lambda2.n)
            end

            # Allocate the resulting fields H & Delta ...
            H, Delta = HFB_Allocate(Params,Lambda,Rho,Kappa,Orb,Orb_NN,Orb_NNN,T,VNN,VNNN)

            # Determine the resulting SQE & amplitudes U & V ...
            SQE, U, V = HFB_Diagonalize(Params,H,Delta,U,V,Orb)

            println("\nHFB iteration terminated due to degeneracy ...")
            println("\tThe solution with lower energy was selected ...\n")
            break
        end

        # Chek for degenerate solutions ...
        if d2E < epsilon * 1e-1 && Degeneracy == false
            println("\nHFB iteration stucked at a degenerate solution ... Degenerate solutions will be analyzed ...\n")
            Degeneracy = true

            Rho1, Kappa1 = pnMatrix(Rho.p,Rho.n), pnMatrix(Kappa.p,Kappa.n)
            Lambda1 = pnFloat(Lambda.p,Lambda.n)
        end

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
        pInd, pMax = 0, 0.0
        nInd, nMax = 0, 0.0

        # Find the most dominant pair of amplitudes U & V ...
        @inbounds for b in 1:a_max
            if (abs(U_new.p[b,a])^2 + abs(V_new.p[b,a])^2) > pMax
                pInd = b
                pMax = (abs(U_new.p[b,a])^2 + abs(V_new.p[b,a])^2)
            end

            if (abs(U_new.n[b,a])^2 + abs(V_new.n[b,a])^2) > nMax
                nInd = b
                nMax = (abs(U_new.n[b,a])^2 + abs(V_new.n[b,a])^2)
            end
        end


        # Check phase change of U & V ...
        if U_new.p[pInd,a] < 0.0
            U_new.p[:,a] .*= -1.0
        end
        if V_new.p[pInd,a] < 0.0
            V_new.p[:,a] .*= -1.0
        end

        if U_new.n[nInd,a] < 0.0
            U_new.n[:,a] .*= -1.0
        end
        if V_new.n[nInd,a] < 0.0
            V_new.n[:,a] .*= -1.0
        end

    end

    return SQE, U_new, V_new
end
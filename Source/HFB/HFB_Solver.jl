function HFB_solver(Params::Parameters)
    # Make single-particle orbitals - NuHamil ordering ...
    Orb = orbitals_make(Params)

    # 1-body kinetic operator ...
    T = T1b(Params,Orb)

    # 2-body NN interaction & Orbitals ...
    @time V_NN, Orb_NN = V2b_read(Params,Orb)

    # 3-body NNN interaction & Orbitals ...
    @time V_NNN, Orb_NNN = V3b_no2b_read(Params,Orb)

    # Solve HFB equations ...
    @time Lambda, SQE, U, V, SQE_C, u_C, v_C, C, Rho, Kappa, H, Delta, Convergence, Iteration = HFB_solve(Params,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN)

    # Calculation of the total HFB mean-field ground-state energy ...
    @time E_HFB = HFB_energy(Params,Rho,Kappa,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN)

    # Calculate the total HFB ground-state kinetic energy ...
    @time T_HFB = T1b_energy(Params,Rho,Orb,T)

    # Particle number fluctuation calculation ...
    @time dA = HFB_particle_number_dispersion(Params,Rho,Orb)

    # Calculation summary ...
    @time HFB_summary(Params,E_HFB,T_HFB,Lambda,dA,Convergence,Iteration)

    # Evaluate HFB charge radii & radial densities ...
    Summary_File = "IO/" * Params.Calc.Path * "/HFB/HFB_Summary.dat"
    Densities_File = "IO/" * Params.Calc.Path * "/HFB/Densities/HFB_Radial_Densities.dat"
    @time OBDM_export(Params,Orb,Summary_File,Densities_File,Rho,C)

    # Export of single-quasiparticle energies, amplitudes U & V & also possibly radial densities ...
    @time HFB_summary_SQS(Params,SQE,SQE_C,O1B(C.p' * Rho.p * C.p,C.n' * Rho.n * C.n),Orb)

    if Params.Calc.HFB.BMF == true
        # Allocate H1B ... 1-body HFB Hamiltonian in the quasiparticle basis ...
        @time H_N = HFB_allocate_H1b(Params,SQE)

        # Make density dependent residual 2-body interaction in the LHO basis ...
        @time V_NN = V2b_residual_no2b(Params,Orb,Orb_NN,Orb_NNN,Rho,V_NN,V_NNN)

        # Perform transformation of U and V to the canonical basis ...
            # I think there should be no transpose when transforming U matrix ,,,
        U, V = O1B(C.p * U.p, C.n * U.n), O1B(C.p' * V.p, C.n' * V.n)

        # Allocate H_NN ... residual interaction 2-body Hamiltonian in the quasiparticle basis ...
        @time H_NN = qpO2b(Params,Orb,Orb_NN,V_NN,C,U,V)

        # Perform final export of HFB solution into binary files ...
        @time HFB_export(Params,Orb,Orb_NN,C,U,V,H_N,H_NN)

        # Perform the HFB-BMBPT(2) calculation of correlation energy ...
        @time HFB_BMBPT_energy(Params,Orb,Orb_NN,H_N,H_NN)
    end

    # Deallocate V_NN 2-body & V_NNN 3-body interaction ...
    V_NN = nothing
    V_NNN = nothing

    # Perform the Garbage Collection ...
    GC.gc()

    return
end

function HFB_solve(Params::Parameters,Orb::Vector{Orb1B},Orb_NN::Orb2B,Orb_NNN::Orb3B,T::O1B,V_NN::O2B,V_NNN::Array{Vector{Vector{Float32}},4})
    # Read calculation parameters ...
    Z_target, N_target = Float64(Params.Calc.Z), Float64(Params.Calc.A - Params.Calc.Z)
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1) * (N_max + 2), 2)
    epsilon = Params.Calc.HFB.Tol

    # Setup local iteration variables ...
    Iteration, Iteration_max, Convergence, Degeneracy = 0, Params.Calc.HFB.IMax, false, false
    dE, d2E, dZ, dN = 0.1, 0.1, 0.1, 0.1

    # Preallocate arrays ...
        # Matrices for U & V HFB ...
    V = O1B(zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max))
    U = O1B(zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max))
         # Matrices for densities Rho & Kappa ...
    Rho = O1B(zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max))
    Kappa = O1B(zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max))
        # Matrices for the single-particle field H & pairing field Delta ...
    H = O1B(zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max))
    Delta = O1B(zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max))
        # Vectors for single-(quasi)particle energies ...
    SQE = pnVector(zeros(Float64,a_max),zeros(Float64,a_max))
        # Identity matrix ...
    Eye = diagm(ones(Float64,a_max))
        # Arrays for degenerate solutions ...
        E_HFB_1, E_HFB_2 = 0.0, 0.0
        Lambda_1, Lambda_2 = pnFloat(0.0,0.0), pnFloat(0.0,0.0)
        U_1, U_2 = O1B(zeros(Float64,a_max,a_max),zeros(Float64,a_max,a_max)), O1B(zeros(Float64,a_max,a_max),zeros(Float64,a_max,a_max))
        V_1, V_2 = O1B(zeros(Float64,a_max,a_max),zeros(Float64,a_max,a_max)), O1B(zeros(Float64,a_max,a_max),zeros(Float64,a_max,a_max))
        Rho_1, Rho_2 = O1B(zeros(Float64,a_max,a_max),zeros(Float64,a_max,a_max)), O1B(zeros(Float64,a_max,a_max),zeros(Float64,a_max,a_max))
        Kappa_1, Kappa_2 = O1B(zeros(Float64,a_max,a_max),zeros(Float64,a_max,a_max)), O1B(zeros(Float64,a_max,a_max),zeros(Float64,a_max,a_max))
        H_1, H_2 = O1B(zeros(Float64,a_max,a_max),zeros(Float64,a_max,a_max)), O1B(zeros(Float64,a_max,a_max),zeros(Float64,a_max,a_max))
        Delta_1, Delta_2 = O1B(zeros(Float64,a_max,a_max),zeros(Float64,a_max,a_max)), O1B(zeros(Float64,a_max,a_max),zeros(Float64,a_max,a_max))
        SQE_1, SQE_2 = O1B(zeros(Float64,a_max,a_max),zeros(Float64,a_max,a_max)), O1B(zeros(Float64,a_max,a_max),zeros(Float64,a_max,a_max))


    # Setup particle numbers ...
    Z, N = 0.1, 0.1

    # Initial guess on chemical potentials lambda ...
    Lambda = pnFloat(Params.Calc.HFB.pL0, Params.Calc.HFB.nL0)

    # Initial guess on densities Rho & Kappa ...
    Rho, Kappa = HFB_density_operator_initialize(Params,Orb)

    # Initialite the Broyden's method ...
    Broyden, BroyVec = HFB_Broyden_initialize(Params,Rho,Kappa,Orb)

    # Preallocate SQE_old ...
    SQE_old = pnVector(SQE.p,SQE.n)

    # Initialize particle numbers & chemical potentials for secant method ...
    Z_1, Z_2 = 0.0, 0.0
    N_1, N_2 = 0.0, 0.0
    Lambda_1, Lambda_2 = pnFloat(0.75 * Lambda.p, 0.75 * Lambda.n), pnFloat(Lambda.p, Lambda.n)

    # Solve the spherical HFB equations ... by the means of self-consistent iteration ...
    println("\nStarting iteration of HFB equations ...")

    while ((dE > epsilon) || (dZ > epsilon) || (dN > epsilon)) && (Iteration < Iteration_max)
        # Perform several secant iterations for the chemical potentil Lambda ...
        @inbounds for L in 1:Iteration_max
            # Initialization of the secant method ...
            if (L == 1) && (Iteration == 0)
                # Allocate the single-particle fields H_1, H_2 and pairing fields Delta_1, Delta_2 ...
                H_1, Delta_1 = HFB_allocate(Params,Lambda_1,Rho,Kappa,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN)
                H_2, Delta_2 = HFB_allocate(Params,Lambda_2,Rho,Kappa,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN)

                # Diagonalize the HFB equations ... basis is reordered as needed ...
                SQE_1, U_1, V_1 = HFB_diagonalize(Params,H_1,Delta_1,Orb)
                SQE_2, U_2, V_2 = HFB_diagonalize(Params,H_2,Delta_2,Orb)
                
                # Generate temporary densities Rho & Kappa ...
                Rho_1, Kappa_1 = HFB_density_operator(Params,U_1,V_1,Orb)
                Rho_2, Kappa_2 = HFB_density_operator(Params,U_2,V_2,Orb)

                # Determine new average particle numbers ...
                Z_1, N_1 = HFB_particle_number(Params,Rho_1,Orb)
                Z_2, N_2 = HFB_particle_number(Params,Rho_2,Orb)

                # Update average particle numbers ...
                Z, N = Z_2, N_2

                # Perform secant iteration to determine new optimal value of Lambda ...
                Lambda = HFB_Lambda_secant(Params,Lambda_1,pnFloat(Z_target - Z_1, N_target - N_1),Lambda_2,pnFloat(Z_target - Z_2, N_target - N_2))

                # Update chemical potentials ...
                Lambda_1 = pnFloat(Lambda_2.p,Lambda_2.n)
                Lambda_2 = pnFloat(Lambda.p,Lambda.n)

                # Update single-quasiparticle energies ...
                SQE = pnVector(SQE_2.p, SQE_2.n)

                # Update amplitudes U & V ...
                U, V = O1B(U_2.p,U_2.n), O1B(V_2.p,V_2.n)

                # Update fields H & Delta ...
                H, Delta = O1B(H_2.p,H_2.n), O1B(Delta_2.p,Delta_2.n)

            end

            # Allocate the single-particle field H and the pairing field Delta ...
                # Smart re-allocation of H ... only Lambda is tweaked, Delta remains the same ...
            H = O1B(H.p .+ (Lambda_1.p - Lambda.p ) .* Eye, H.n .+ (Lambda_1.n - Lambda.n ) .* Eye)

            # Diagonalize the HFB equations ... basis is reordered as needed ...
            SQE, U, V = HFB_diagonalize(Params,H,Delta,Orb)
            
            # Generate temporary densities Rho & Kappa ...
            Rho_Temp, Kappa_Temp = HFB_density_operator(Params,U,V,Orb)

            # Determine new average particle numbers ...
            Z, N = HFB_particle_number(Params,Rho_Temp,Orb)

            # Update average particle numbers ...
            Z_1, N_1 = Z_2, N_2
            Z_2, N_2 = Z, N

            # Determine new chemical potential Lambda ...
            Lambda = HFB_Lambda_secant(Params,Lambda_1,pnFloat(Z_target - Z_1, N_target - N_1),Lambda_2,pnFloat(Z_target - Z_2, N_target - N_2))

            # Update chemical potentials ...
            Lambda_1 = pnFloat(Lambda_2.p,Lambda_2.n)
            Lambda_2 = pnFloat(Lambda.p,Lambda.n)

            # Evaluate iteration of the chemical potential Lambda ...
            dZ, dN = abs(Z_target - Z), abs(N_target - N)

            # Check finite differences for particle numbers ....
            if dZ < epsilon && dN < epsilon
                break
            end
        end

        # Generate new densities Rho & Kappa ...
        Rho, Kappa = HFB_density_operator(Params,U,V,Orb)

        # Perform update of densities ...
        Rho, Kappa, BroyVec = HFB_Broyden_update(Params,Iteration,Rho,Kappa,Broyden,BroyVec,Orb)

        # Allocate the single-particle field H and the pairing field Delta ...
        H, Delta = HFB_allocate(Params,Lambda,Rho,Kappa,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN)

        # Diagonalize the HFB equations ... basis is reordered as needed ...
        SQE, U, V = HFB_diagonalize(Params,H,Delta,Orb)

        # Evaluate particle numbers ...
        Z, N = HFB_particle_number(Params,Rho,Orb)

        # Evaluate iteration of the single-quasiparticle energies ...
        d2E = abs(dE - (sum(abs.(SQE.p .- SQE_old.p )) + sum(abs.(SQE.n .- SQE_old.n))) / Float64(2 * a_max))
        dE = (sum(abs.(SQE.p .- SQE_old.p )) + sum(abs.(SQE.n .- SQE_old.n))) / Float64(2 * a_max)
        dZ, dN = abs(Z_target - Z), abs(N_target - N)

        if Degeneracy == false
            Iteration += 1
            @printf("\tHFB iteration number: %4d   dE = %12.9f MeV,   dZ = %12.9f,   dN = %12.9f\n", Iteration, dE, dZ, dN)
        end

        # Store old values of SQE ...
        SQE_old = pnVector(SQE.p,SQE.n)

        # Resolve degenerate solutions ...
        if Degeneracy == true
            H_2, Delta_2, Lambda_2 = O1B(H.p,H.n), O1B(Delta.p,Delta.n), pnFloat(Lambda.p,Lambda.n)
            SQE_2, U_2, V_2 = HFB_diagonalize(Params,H_2,Delta_2,Orb)
            Rho_2, Kappa_2 = HFB_density_operator(Params,U_2,V_2,Orb)
            Out = "/dev/null"
            if Sys.iswindows()
                Out = "NUL"
            end
            open(Out, "w") do devnull_io
                redirect_stdout(devnull_io) do
                    redirect_stderr(devnull_io) do
                        E_HFB_2 = HFB_energy(Params,Rho_2,Kappa_2,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN)
                        return E_HFB_2
                    end
                end
            end
            if sum(E_HFB_1) < sum(E_HFB_2)
                U, V, Lambda = O1B(U_1.p,U_1.n), O1B(V_1.p,V_1.n), pnFloat(Lambda_1.p,Lambda_1.n)
                H, Delta, SQE = O1B(H_1.p,H_1.n), O1B(Delta_1.p,Delta_1.n), pnVector(SQE_1.p,SQE_1.n)
                Rho, Kappa = O1B(Rho_1.p,Rho_1.n), O1B(Kappa_1.p,Kappa_1.n)
            else
                U, V, Lambda = O1B(U_2.p,U_2.n), O1B(V_2.p,V_2.n), pnFloat(Lambda_2.p,Lambda_2.n)
                H, Delta, SQE = O1B(H_2.p,H_2.n), O1B(Delta_2.p,Delta_2.n), pnVector(SQE_2.p,SQE_2.n)
                Rho, Kappa = O1B(Rho_2.p,Rho_2.n), O1B(Kappa_2.p,Kappa_2.n)
            end
            println("\nThe degenerate solution has been chosen ... HFB iteration terminates\n")
            Convergence = true
            break
        end

        # If convergence is reached & terminate the loop ...
        if dE < epsilon && dZ < epsilon && dN < epsilon
            println("\nHFB iteration with residual NN interaction has converged ...")
            Convergence = true
            break
        end

        # Check for degenerate solutions ...
        if d2E < dE * 1e-4 && dE < 1e-1 && dZ < epsilon && dN < epsilon
            println("\nHFB iteration stucked at a degenerate solution ... Degenerate solutions will be analyzed ...")
            Degeneracy = true
            H_1, Delta_1, Lambda_1 = O1B(H.p,H.n), O1B(Delta.p,Delta.n), pnFloat(Lambda.p,Lambda.n)
            SQE_1, U_1, V_1 = HFB_diagonalize(Params,H_1,Delta_1,Orb)
            Rho_1, Kappa_1 = HFB_density_operator(Params,U_1,V_1,Orb)
            Out = "/dev/null"
            if Sys.iswindows()
                Out = "NUL"
            end
            open(Out, "w") do devnull_io
                redirect_stdout(devnull_io) do
                    redirect_stderr(devnull_io) do
                        E_HFB_1 = HFB_energy(Params,Rho_1,Kappa_1,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN)
                        return E_HFB_1
                    end
                end
            end
        end

        # Terminate the loop if maximum number of iterations is reached ...
        if Iteration == Iteration_max
            println("\nMaximum number of HFB iterations reached ... Iteration terminates & convergence was not reached ...\n")
            Convergence = false
            break
        end

    end

    # Perform final evaluation of resulting densities ...
    #Rho, Kappa = HFB_density_operator(Params,U,V,Orb)

    # self-consistency check
    #=
    # Generate new densities Rho & Kappa ...
    Rho, Kappa = HFB_density_operator(Params,U,V,Orb)

    # Allocate the single-particle field H and the pairing field Delta ...
    H, Delta = HFB_allocate(Params,Lambda,Rho,Kappa,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN)

    # Diagonalize the HFB equations ... basis is reordered as needed ...
    SQE1, U, V = HFB_diagonalize(Params,H,Delta,Orb)

    dE = (sum(abs.(SQE.p .- SQE1.p )) + sum(abs.(SQE.n .- SQE1.n))) / Float64(2 * a_max)

    println("dE = " * string(dE))
    =#



    # Construct the canonical basis & evaluate approximate (BCS-like) amplitudes
    # & single-quasiparticle energies u_C, v_C & SQE_C ...
    #   U, V, Rho, Kappa, H, Delta ... remain expressed in the reference LHO basis ...
    SQE_C, C, u_C, v_C = HFB_canonical_basis(Params,Rho,H,Delta,Orb)

    return Lambda, SQE, U, V, SQE_C, u_C, v_C, C, Rho, Kappa, H, Delta, Convergence, Iteration
end

function HFB_diagonalize(Params::Parameters,H::O1B,Delta::O1B,Orb::Vector{Orb1B})
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

    # Perform reordering - to match quantum numbers j & l ascending in E ...
    SQE, U_New, V_New = HFB_orbital_ordering(Params,pnVector(pSQE,nSQE),O1B(pU,nU),O1B(pV,nV),Orb)

    # Set phases of U & V to match with previous iteration ...

    # Old algorithm ...
    #=
    @inbounds for a in 1:a_max
        # Setup local variables ...
        pInd, pMax = 0, 0.0
        nInd, nMax = 0, 0.0

        # Find the most dominant pair of amplitudes U & V ...
        @inbounds for b in 1:a_max
            if (abs(U_New.p[b,a])^2 + abs(V_New.p[b,a])^2) > pMax
                pInd = b
                pMax = (abs(U_New.p[b,a])^2 + abs(V_New.p[b,a])^2)
            end

            if (abs(U_New.n[b,a])^2 + abs(V_New.n[b,a])^2) > nMax
                nInd = b
                nMax = (abs(U_New.n[b,a])^2 + abs(V_New.n[b,a])^2)
            end
        end

        # Check phase change of U & V ...
        if U_New.p[pInd,a] < 0.0
            @views U_New.p[:,a] .*= -1.0
        end
        if V_New.p[pInd,a] < 0.0
            @views V_New.p[:,a] .*= -1.0
        end

        if U_New.n[nInd,a] < 0.0
            @views U_New.n[:,a] .*= -1.0
        end
        if V_New.n[nInd,a] < 0.0
            @views V_New.n[:,a] .*= -1.0
        end

    end
    =#

    #=
    # Modified phase algorithm ...
    @inbounds for a in 1:a_max
        # Setup local variables ...
        pInd, pMax, pType = 1, -Inf, :U
        nInd, nMax, nType = 1, -Inf, :U

        # Find the most dominant pair of amplitudes U & V ...
        @inbounds for b in 1:a_max
            if (abs(U_New.p[b,a])^2 + abs(V_New.p[b,a])^2) > pMax
                pInd = b
                if abs(U_New.p[b,a]) < 1e-8
                    pType = :V
                end
                pMax = abs(U_New.p[b,a])^2 + abs(V_New.p[b,a])^2
            end

            if (abs(U_New.n[b,a])^2 + abs(V_New.n[b,a])^2) > nMax
                nInd = b
                if abs(U_New.n[b,a]) < 1e-8
                    nType = :V
                end
                nMax = abs(U_New.n[b,a])^2 + abs(V_New.n[b,a])^2
            end
        end

        # Check phase change of U & V ...
        if (pType == :U && U_New.p[pInd,a] < 0.0) || (pType == :V && V_New.p[pInd,a] < 0.0)
            @views U_New.p[:,a] .*= -1.0
            @views V_New.p[:,a] .*= -1.0
        end

        if (nType == :U && U_New.n[nInd,a] < 0.0) || (nType == :V && V_New.n[nInd,a] < 0.0)
            @views U_New.n[:,a] .*= -1.0
            @views V_New.n[:,a] .*= -1.0
        end

    end
    =#
 
    return SQE, U_New, V_New
end
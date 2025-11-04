struct HFB_Broyden
    Map::Matrix{Int64}
    Key::Vector{Matrix{Int64}}
    M::Int64
end




# Cheap Broyden HFB ... HFB iteration using simplified & cheap block Broyden method ... doesnt work yet ...
function HFB_Solve_CheapBroyden(Params::Parameters,Orb::Vector{NOrb},Orb_NN::NNOrb,Orb_NNN::NNNOrb,T::Matrix{Float64},VNN::NNInt,VNNN::Array{Vector{Vector{Float32}},4},epsilon::Float64)
    # Read calculation parameters ...
    A_target, Z_target, N_target = Float64(Params.Calc.A), Float64(Params.Calc.Z), Float64(Params.Calc.A - Params.Calc.Z)
    hw = Params.Int.hw
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1) * (N_max + 2), 2)

    # Setup local iteration variables ...
    Iteration, Iteration_max = 0, 250
    dE, dZ, dN = 1.0, 1.0, 1.0

    # Preallocate arrays ...
        # Matrices for U & V HFB ...
    V = pnMatrix(zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max))
    U = pnMatrix(zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max))
         # Matrices for densities Rho & Kappa ...
    Rho, Kappa = HFB_Density_Operator_Initialize(Params,Orb)
    Rho_old, Kappa_old = pnMatrix(Rho.p,Rho.n), pnMatrix(Kappa.p,Kappa.n)
        # Matrices for the single-particle field H & pairing field Delta ...
    H = pnMatrix(zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max))
    H_old = pnMatrix(zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max))
    Delta = pnMatrix(zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max))
    Delta_old = pnMatrix(zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max))
        # Vectors for single-(quasi)particle energies ...
    SQE = pnVector(zeros(Float64,a_max),zeros(Float64,a_max))
    SQE_old = pnVector(zeros(Float64,a_max),zeros(Float64,a_max))

    # Setup particle numbers ...
    Z, Z_1, Z_2 = Z_target, 0.0, 0.0
    N, N_1, N_2 = N_target, 0.0, 0.0

    # Initial guess on chemical potentials lambda ...
    Lambda, Lambda_2, Lambda_1 = pnFloat(0.5,0.5), pnFloat(0.5,0.5), pnFloat(1.0,1.0)

    # Initialize the Broyden Jacobian B ...
    B = diagm(ones(Float64,4))

    # Solve the spherical HFB equations ... by the means of self-consistent cheap Broyden iteration ...
    println("\nStarting iteration of HFB equations ...\n")

    while ((dE > epsilon) || (dZ > epsilon) || (dN > epsilon)) && (Iteration < Iteration_max)

        # Perform a few secant iterations for the chemical potentials Lambda ...
        #println("\nStarting to iterate the chemical potential Lambda ...")
        @inbounds for L in 1:6
            # Initialization of the secant method ...
            if (Iteration == 0) && (L == 1)
                # Allocate the single-particle fields H_1, H_2 and pairing fields Delta_1, Delta_2 ...
                H_1, Delta_1 = HFB_Allocate(Params,Lambda_1,Rho,Kappa,Orb,Orb_NN,Orb_NNN,T,VNN,VNNN)
                H_2, Delta_2 = HFB_Allocate(Params,Lambda_2,Rho,Kappa,Orb,Orb_NN,Orb_NNN,T,VNN,VNNN)

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
                Lambda_1, Lambda_2 = pnFloat(Lambda_2.p, Lambda_2.n), pnFloat(Lambda.p, Lambda.n)

                # Update single-quasiparticle energies ...
                SQE = pnVector(SQE_2.p, SQE_2.n)

                # Update amplitudes U & V ...
                U, V = pnMatrix(U_2.p,U_2.n), pnMatrix(V_2.p,V_2.n)

                # Update the matrices H, Delta, Rho, Kappa
                H, H_old = pnMatrix(H_2.p,H_2.n), pnMatrix(H_1.p,H_1.n)
                Delta, Delta_old = pnMatrix(Delta_2.p,Delta_2.n), pnMatrix(Delta_1.p,Delta_1.n)
                Rho, Rho_old = pnMatrix(Rho_2.p,Rho_2.n), pnMatrix(Rho_1.p,Rho_1.n)
                Kappa, Kappa_old = pnMatrix(Kappa_2.p,Kappa_2.n), pnMatrix(Kappa_1.p,Kappa_1.n)

            end

            # Allocate the single-particle field H and the pairing field Delta ...
            H, Delta = HFB_Allocate(Params,Lambda,Rho,Kappa,Orb,Orb_NN,Orb_NNN,T,VNN,VNNN)

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

            if dZ < epsilon && dN < epsilon
                break
            end

        end

        # Allocate the single-particle field H and the pairing field Delta ...
        H, Delta = HFB_Allocate(Params,Lambda,Rho,Kappa,Orb,Orb_NN,Orb_NNN,T,VNN,VNNN)

        # Diagonalize the HFB equations ... basis is reordered as needed ...
        SQE, U, V = HFB_Diagonalize(Params,H,Delta,Orb)

        # Generate new densities Rho & Kappa ...
        Rho, Kappa = HFB_Density_Operator(Params,U,V,Orb)

        # Perform a cheap Broyden iteration - Update of Rho, Kappa & Lambda ...
        Rho, Kappa, B = HFB_CheapBroyden(Params,Rho,Kappa,H,Delta,Rho_old,Kappa_old,H_old,Delta_old,B,Orb)

        # Evaluate the particle numbers ...
        Z_1, N_1 = Z, N
        Z, N = HFB_Particle_Number(Params,Rho,Orb)
        Z_2, N_2 = Z, N

        # Evaluate iteration of the single-quasiparticle energies ...
        Iteration += 1

        dE = (sum(abs.(SQE.p .- SQE_old.p )) + sum(abs.(SQE.n .- SQE_old.n))) / Float64(2 * a_max)
        dZ, dN = abs(Z_target - Z), abs(N_target - N)

        println("\n\nHFB iteration number:   " * string(Iteration) * "   Single-quasiparticle energy difference:   " * string(round(dE, sigdigits=8))
                * "   Proton number difference:   " * string(round(dZ, sigdigits=8)) * "   Neutron number difference:   " * string(round(dN, sigdigits=8)))

        # Store old values of SQE, Rho, Kappa, Lambda, H, Delta, N & B ...
        SQE_old, Lambda_old = pnVector(SQE.p,SQE.n), pnFloat(Lambda.p, Lambda.n)
        Rho_old, Kappa_old = pnMatrix(Rho.p,Rho.n), pnMatrix(Kappa.p,Kappa.n)
        H_old, Delta_old = pnMatrix(H.p,H.n), pnMatrix(Delta.p,Delta.n)

        println("pLambda = " * string(Lambda.p))
        println("nLambda = " * string(Lambda.n))
        println("Z = " * string(Z))
        println("N = " * string(N))

    end

    println("\nHFB iteration with residual NN interaction has converged ...")

    # Perform transformation of U & V into the canonical basis ...
    #   Rho, Kappa, H, Delta ... remain expressed in the reference LHO basis ...
    SQE_C, C, U_C, V_C = HFB_Canonical_Basis(Params,U,V,Rho,H,Delta,Orb)

    return Lambda, SQE, U, V, SQE_C, U_C, V_C, C, Rho, Kappa, H, Delta, Iteration
end

function HFB_CheapBroyden(Params::Parameters,Rho::pnMatrix,Kappa::pnMatrix,H::pnMatrix,Delta::pnMatrix,Rho_old::pnMatrix,Kappa_old::pnMatrix,H_old::pnMatrix,Delta_old::pnMatrix,B_old::Matrix{Float64},Orb::Vector{NOrb})
    # Read calculation parameters ...
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1) * (N_max + 2), 2)

    # Initialize needed arrays ...
    dx, df, B = zeros(Float64,4), zeros(Float64,4), zeros(Float64,4,4)

    # Evaluate finite differences ...
    d_pRho, d_nRho = Rho.p .- Rho_old.p, Rho.n .- Rho_old.n
    d_pKappa, d_nKappa = Kappa.p .- Kappa_old.p, Kappa.n .- Kappa_old.n
    d_pH, d_nH = H.p .- H_old.p, H.n .- H_old.n
    d_pDelta, d_nDelta = Delta.p .- Delta_old.p, Delta.n .- Delta_old.n

    # Allocate finite difference vectors dx & df ...
    dx[1], dx[2] = sum(d_pRho), sum(d_pKappa)
    dx[3], dx[4] = sum(d_nRho), sum(d_nKappa)

    df[1], df[2] = sum(d_pH), sum(d_pDelta)
    df[3], df[4] = sum(d_nH), sum(d_nDelta)

    # Renormalize dx & df ...
    dx .= dx ./ Float64(a_max^2)
    df .= df ./ Float64(a_max^2)

    # Update the Broyden Jacobian B ...
    D = dot(dx,dx)
    if D < 1e-14
        D += 1e-14
    end
    B .= B_old + ((df - B_old * dx) * dx') / D

    # Invert the Broyden Jacobian B ...
    eta = 1e-12 * maximum(abs.(diag(B))) + 1e-14
    iB = inv(B + eta * diagm(ones(Float64,4)))

    #display(iB)
    #display(H.n)
    #display(H.n .- H_old.n)

    # Evaluate the Broyden finite difference steps ...
    dRho = pnMatrix(iB[1,1] * (H.p .- H_old.p) .+ iB[1,2] * (Delta.p .- Delta_old.p) .+ iB[1,3] * (H.n .- H_old.n) .+ iB[1,4] * (Delta.n .- Delta_old.n),
                    iB[3,1] * (H.p .- H_old.p) .+ iB[3,2] * (Delta.p .- Delta_old.p) .+ iB[3,3] * (H.n .- H_old.n) .+ iB[3,4] * (Delta.n .- Delta_old.n))

    dKappa = pnMatrix(iB[2,1] * (H.p .- H_old.p) .+ iB[2,2] * (Delta.p .- Delta_old.p) .+ iB[2,3] * (H.n .- H_old.n) .+ iB[2,4] * (Delta.n .- Delta_old.n),
                      iB[4,1] * (H.p .- H_old.p) .+ iB[4,2] * (Delta.p .- Delta_old.p) .+ iB[4,3] * (H.n .- H_old.n) .+ iB[4,4] * (Delta.n .- Delta_old.n))

    Z, N = HFB_Particle_Number(Params,Rho,Orb)
    dZ, dN = HFB_Particle_Number(Params,dRho,Orb)

    if dZ / Z > 0.05
        dRho = pnMatrix(0.05 * Z / dZ .* dRho.p, dRho.n)
        dKappa = pnMatrix(0.05 * Z / dZ .* dKappa.p, dKappa.n)
    end
    if dN / N > 0.05
        dRho = pnMatrix(dRho.p, 0.05 * N / dN .* dRho.n)
        dKappa = pnMatrix(dKappa.p, 0.05 * N / dN .* dKappa.n)
    end

    #display(Rho.n)
    #display(Kappa.n)

    # Perform the Broyden update of densities Rho & Kappa ...
    Rho, Kappa = HFB_Broyden_Update(a_max,Rho,dRho,Kappa,dKappa)

    #display(Rho.n)
    #display(Kappa.n)
    #throw("Stop here")

    return Rho, Kappa, B
end

function HFB_Broyden_Update(a_max::Int64,Rho::pnMatrix,dRho::pnMatrix,Kappa::pnMatrix,dKappa::pnMatrix)
    # Basic parameters ...
    q, q_min,g, c = 1.0, 1e-4, 0.5, 1e-12

    # Evaluate current finite difference d_current ...
    d_current = (sum(abs.(dRho.p)) + sum(abs.(dKappa.p)) + sum(abs.(dRho.n)) + sum(abs.(dKappa.n))) / Float64(4*a_max^2)

    # Iterate on the quenching factor q ...
    while q >= q_min
        # Perform trial step ...
        pRho_trial, pKappa_trial = (1.0 - q) * Rho.p .- q * dRho.p, (1.0 - q) * Kappa.p .- q * dKappa.p
        nRho_trial, nKappa_trial = (1.0 - q) * Rho.p .- q * dRho.n, (1.0 - q) * Kappa.n .- q * dKappa.n

        # Symmetrize trial densities ...
        pRho_trial .= 0.5 * (pRho_trial .+ pRho_trial')
        nRho_trial .= 0.5 * (nRho_trial .+ nRho_trial')
        pKappa_trial .= 0.5 * (pKappa_trial .+ pKappa_trial')
        nKappa_trial .= 0.5 * (nKappa_trial .+ nKappa_trial')

        # Evaluate trial finite difference d_trial ...
        d_trial = (sum(abs.(pRho_trial .- Rho.p)) + sum(abs.(nRho_trial .- Rho.n)) + sum(abs.(pKappa_trial .- Kappa.p)) + sum(abs.(nKappa_trial .- Kappa.n))) / Float64(4*a_max^2)

        # Armijo-like condition on accepting the current trial step ...
        if (d_trial + c * q) < d_current
            return pnMatrix(pRho_trial, nRho_trial), pnMatrix(pKappa_trial, nKappa_trial)
        end

        q = g * q
    end

    # No improvement due to quenching ... return old densities ...
    return Rho, Kappa
end

function HFB_Broyden_Update_Difference(pRho::Matrix{Float64},nRho::Matrix{Float64},pKappa::Matrix{Float64},nKappa::Matrix{Float64},Rho::pnMatrix,Kappa::pnMatrix)
    d = sum(abs.(pRho .- Rho.p)) + sum(abs.(nRho .- Rho.n)) + sum(abs.(pKappa .- Kappa.p)) + sum(abs.(nKappa .- Kappa.n))
    return d
end
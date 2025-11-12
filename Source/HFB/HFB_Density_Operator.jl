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
    pKappa, nKappa = 0.5 .* diagm(ones(Float64,a_max)), 0.5 .* diagm(ones(Float64,a_max))

    return pnMatrix(pRho,nRho), pnMatrix(pKappa,nKappa)
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

function HFB_Mixing_Update(a_max::Int64,Rho::pnMatrix,Rho_old::pnMatrix,Kappa::pnMatrix,Kappa_old::pnMatrix)
    # Basic parameters ...
    q, q_min = 0.9, 1e-8

    # Allocate finite differences for densities ...
    dRho = pnMatrix(Rho.p .- Rho_old.p, Rho.n .- Rho_old.n)
    dKappa = pnMatrix(Kappa.p .- Kappa_old.p, Kappa.n .- Kappa_old.n)

    # Evaluate current finite difference D ...
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
        else
            return pnMatrix(pRho_trial,nRho_trial), pnMatrix(pKappa_trial,nKappa_trial)
        end
    end

    # No improvement due to quenching ... return old densities ...
    return pnMatrix(0.975 * Rho.p, 0.975 * Rho.n), pnMatrix(0.975 * Kappa.p, 0.975 * Kappa.n)
end
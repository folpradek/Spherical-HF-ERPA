function HFB_Density_Operator_Initialize(Params::Parameters,Orb::Vector{NOrb})
    # Read parameters ...
    Z_target, N_target = Float64(Params.Calc.Z), Float64(Params.Calc.A - Params.Calc.Z)
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1) * (N_max + 2), 2)

    # Initialize the valence space shell number N ...
    pN_valence, nN_valence = 0, 0

    # Determine the valence space shell number N ...
    @inbounds for N in 0:N_max
        if div((N  + 1) * (N + 2),2) >= Z_target
            pN_valence = N - 1
        end
        if div((N  + 1) * (N + 2),2) >= N_target
            nN_valence = N - 1
        end
    end

    # Initialize particle numbers ...
    Z, N = 0.0, 0.0

    # Initialize density matrices ...
    pRho, nRho = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)
    pKappa, nKappa = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)

    # Allocate the density matrix Rho ...
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

    # Allocate the pairing tensor Kappa ...
    @inbounds for a in 1:a_max
        n_a, l_a = Orb[a].n, Orb[a].l
        N_a = 2 * n_a + l_a
        pK0, pdN0 = Params.Calc.Pairing.pK0, Params.Calc.Pairing.pdN0
        nK0, ndN0 = Params.Calc.Pairing.nK0, Params.Calc.Pairing.ndN0
        if abs(pK0) > (0.5 - 1e-8)
            pK0 = 0.5 * pK0 / abs(pK0)
        end
        if abs(nK0) > (0.5 - 1e-8)  
            nK0 = 0.5 * nK0 / abs(nK0)
        end
        pKappa[a,a] = pK0 * exp(-Float64((N_a - pN_valence)^2) / abs(pdN0)^2)
        nKappa[a,a] = nK0 * exp(-Float64((N_a - nN_valence)^2) / abs(ndN0)^2)
    end

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
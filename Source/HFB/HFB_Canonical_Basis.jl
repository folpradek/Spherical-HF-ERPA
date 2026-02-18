function HFB_canonical_basis(Params::Parameters,Rho::O1B,Orb::Vector{Orb1B})
    # Read calculation parameters ...
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1) * (N_max + 2), 2)

    # Allocate density matrices ...
    pRho, nRho = deepcopy(Rho.p), deepcopy(Rho.n)

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
        pRho[i,i] += 1e-11 * Float64(a_max-i)
        nRho[i,i] += 1e-11 * Float64(a_max-i)
    end

    # Diagonalize 1-body HFB density matrix Rho ...
    pn_C, pC = eigen(Symmetric(pRho),sortby=-)
    nn_C, nC = eigen(Symmetric(nRho),sortby=-)

    # Clean numerical noise in C ...
    @inbounds for a in 1:a_max
        @inbounds for b in 1:a_max
            pCME, nCME = abs(pC[a,b]), abs(nC[a,b])
            if pCME < 1e-13
                pC[a,b] = 0.0
            end
            if nCME < 1e-13
                nC[a,b] = 0.0
            end
        end
    end

    # Reorder the transformation matrix C ...
    C = HFB_canonical_basis_particle_reordering(Params,O1B(pC,nC),Orb)

    return C
end

function HFB_canonical_basis_BCS(Params::Parameters,H::O1B,Orb::Vector{Orb1B})
    # Read calculation parameters ...
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1) * (N_max + 2), 2)

    # Allocate density matrices ...
    pH, nH = deepcopy(H.p), deepcopy(H.n)

    # Symmetrize density matrices ...
    pH .= 0.5 .* (pH .+ pH')
    nH .= 0.5 .* (nH .+ nH')

    # Eliminate any possible numerical noise spoiling block-diagonal structure ...
    @inbounds for a in 1:a_max
        @inbounds for b in 1:a_max
            if (Orb[a].j != Orb[b].j) || (Orb[a].l != Orb[b].l)
                pH[a,b] = 0.0
                nH[a,b] = 0.0
            end
        end
    end

    # Diagonalize 1-body HFB density matrix Rho ...
    pe, pC = eigen(Symmetric(pH),sortby=+)
    ne, nC = eigen(Symmetric(nH),sortby=+)

    # Clean numerical noise in C ...
    @inbounds for a in 1:a_max
        @inbounds for b in 1:a_max
            pCME, nCME = abs(pC[a,b]), abs(nC[a,b])
            if pCME < 1e-13
                pC[a,b] = 0.0
            end
            if nCME < 1e-13
                nC[a,b] = 0.0
            end
        end
    end

    # Reorder the transformation matrix C ...
    C = HFB_canonical_basis_particle_reordering(Params,O1B(pC,nC),Orb)

    return C
end
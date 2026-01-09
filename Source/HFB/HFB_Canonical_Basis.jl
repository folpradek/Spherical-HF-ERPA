function HFB_canonical_basis(Params::Parameters,Rho::O1B,H::O1B,Delta::O1B,Orb::Vector{Orb1B})
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
            if pCME < 1e-11
                pC[a,b] = 0.0
            end
            if nCME < 1e-11
                nC[a,b] = 0.0
            end
        end
    end

    # Reorder the transformation matrix C & occupation probabilities n_C ...
    pn_C, pC, nn_C, nC = HFB_canonical_basis_particle_reordering(Params,pnVector(pn_C,nn_C),O1B(pC,nC),Orb)

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
    SQE_C = HFB_canonical_basis_SQE(Params,O1B(pC,nC),H,Delta)

    # Reorder the single-quasiparticle orbitals ... u_C, v_C & SQE_C ...
    SQE_C, u_C, v_C = HFB_canonical_basis_quasiparticle_reordering(Params,SQE_C,O1B(pu_C,nu_C),O1B(pv_C,nv_C),Orb)

    return SQE_C, O1B(pC,nC), u_C, v_C
end

function HFB_canonical_basis_SQE(Params::Parameters,C::O1B,H::O1B,Delta::O1B)
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
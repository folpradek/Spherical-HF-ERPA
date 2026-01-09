function HFB_orbital_ordering(Params::Parameters,SQE::pnVector,U::O1B,V::O1B,Orb::Vector{Orb1B})
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

    return pnVector(pSQE,nSQE), O1B(pU,nU), O1B(pV,nV)
end

function HFB_canonical_basis_particle_reordering(Params::Parameters,n_C::pnVector,C::O1B,Orb::Vector{Orb1B})
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

function HFB_canonical_basis_quasiparticle_reordering(Params::Parameters,SQE_C::pnVector,u_C::O1B,v_C::O1B,Orb::Vector{Orb1B})
    # Read calculation parameters ...
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1) * (N_max + 2), 2)

    # Read needed arrays ...
    pSQE_C, pu_C, pv_C = SQE_C.p, u_C.p, v_C.p
    nSQE_C, nu_C, nv_C = SQE_C.n, u_C.n, v_C.n

    # Reorder basis according to energy ...
    pOrbs_sort = sortperm(pSQE_C)
    nOrbs_sort = sortperm(nSQE_C)

    # Proton single-quasiparticle orbitals ...
    @views pSQE_C .= pSQE_C[pOrbs_sort]
    @views pu_C .= pu_C[:,pOrbs_sort]
    @views pv_C .= pv_C[:,pOrbs_sort]

    # Neutron single-quasiparticle orbitals ...
    @views nSQE_C .= nSQE_C[nOrbs_sort]
    @views nu_C .= nu_C[:,nOrbs_sort]
    @views nv_C .= nv_C[:,nOrbs_sort]

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

    return pnVector(pSQE_C,nSQE_C), O1B(pu_C,nu_C), O1B(pv_C,nv_C)
end
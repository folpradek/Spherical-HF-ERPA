function HFB_orbital_ordering(Params::Parameters,SQE::pnVector,U::O1B,V::O1B,Orb::Vector{Orb1B};Final_Ordering::Bool=false,C::O1B=O1B(zeros(Float64,1,1),zeros(Float64,1,1)))
    # Read parameters ...
    Tol = 1e-3
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1) * (N_max + 2), 2)

    # Read input arrays ...
    pU, pV, pSQE =  deepcopy(U.p), deepcopy(V.p), deepcopy(SQE.p)
    nU, nV, nSQE =  deepcopy(U.n), deepcopy(V.n), deepcopy(SQE.n)

    # Final reordering by occupation numbers ...
    if Final_Ordering == true && Params.Calc.HFB.Pairing != "BCS"
        # Reorder ... descending in U ...
        pU_norm, nU_norm = zeros(Float64,a_max), zeros(Float64,a_max)

        @inbounds for a in 1:a_max
            pU_norm[a] = @views norm(pU[:,a])^2
            nU_norm[a] = @views norm(nU[:,a])^2
        end

        pOrbs_sort = sortperm(pU_norm)
        nOrbs_sort = sortperm(nU_norm)

        # Proton single-quasiparticle orbitals ...
        @views pSQE .= pSQE[pOrbs_sort]
        @views pU .= pU[:,pOrbs_sort]
        @views pV .= pV[:,pOrbs_sort]

        # Neutron single-quasiparticle orbitals ...
        @views nSQE .= nSQE[nOrbs_sort]
        @views nU .= nU[:,nOrbs_sort]
        @views nV .= nV[:,nOrbs_sort]

            # To ensure the Residual Hamiltonian is correctly built ...
        #=
            pOrbs_sort = sortperm(pSQE)
            nOrbs_sort = sortperm(nSQE)

            # Proton single-quasiparticle orbitals ...
            @views pSQE .= pSQE[pOrbs_sort]
            @views pU .= pU[:,pOrbs_sort]
            @views pV .= pV[:,pOrbs_sort]

            # Neutron single-quasiparticle orbitals ...
            @views nSQE .= nSQE[nOrbs_sort]
            @views nU .= nU[:,nOrbs_sort]
            @views nV .= nV[:,nOrbs_sort]
        =#

    end

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
            if pOrb_mask[b] == false && abs(Float64(Orb[b].l) - pl) < Tol && abs(Float64(Orb[b].j) - pj) < Tol
                pOrb_order[b] = a
                pOrb_mask[b] = true
                break
            end
        end

        @inbounds for b in 1:a_max
            if nOrb_mask[b] == false && abs(Float64(Orb[b].l) - nl) < Tol && abs(Float64(Orb[b].j) - nj) < Tol
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

function HFB_BMZ_orbital_ordering(Params::Parameters,C::O1B,D::O1B,Orb::Vector{Orb1B})
    # Read parameters ...
    Tol = 1e-3
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1) * (N_max + 2), 2)

    # Read input arrays ...
    pC, pD =  deepcopy(C.p), deepcopy(D.p)
    nC, nD =  deepcopy(C.n), deepcopy(D.n)

    # Preallocate temporary arrays ...
    pC_order, nC_order = Vector{Int64}(undef,a_max), Vector{Int64}(undef,a_max)
    pC_mask, nC_mask = falses(a_max), falses(a_max)

    pD_order, nD_order = Vector{Int64}(undef,a_max), Vector{Int64}(undef,a_max)
    pD_mask, nD_mask = falses(a_max), falses(a_max)

    pC_l_values, pC_j_values = zeros(Float64,a_max), zeros(Float64,a_max)
    nC_l_values, nC_j_values = zeros(Float64,a_max), zeros(Float64,a_max)

    pD_l_values, pD_j_values = zeros(Float64,a_max), zeros(Float64,a_max)
    nD_l_values, nD_j_values = zeros(Float64,a_max), zeros(Float64,a_max)

    # Evaluate values of j & l for single-quasiparticle orbitals ...
    @inbounds for a in 1:a_max
        pCjSum, pClSum = 0.0, 0.0
        nCjSum, nClSum = 0.0, 0.0
        pDjSum, pDlSum = 0.0, 0.0
        nDjSum, nDlSum = 0.0, 0.0
        @inbounds for b in 1:a_max
            l_b, j_b = Float64(Orb[b].l), Float64(Orb[b].j)
            pCN, nCN = pC[b,a]^2, nC[b,a]^2
            pCjSum, pClSum = pCjSum + j_b * pCN, pClSum + l_b * pCN
            nCjSum, nClSum = nCjSum + j_b * nCN, nClSum + l_b * nCN
            pDN, nDN = pD[b,a]^2, nD[b,a]^2
            pDjSum, pDlSum = pDjSum + j_b * pDN, pDlSum + l_b * pDN
            nDjSum, nDlSum = nDjSum + j_b * nDN, nDlSum + l_b * nDN
        end
        pC_l_values[a], pC_j_values[a] = pClSum, pCjSum
        nC_l_values[a], nC_j_values[a] = nClSum, nCjSum
        pD_l_values[a], pD_j_values[a] = pDlSum, pDjSum
        nD_l_values[a], nD_j_values[a] = nDlSum, nDjSum
    end

    # Find the ordering for single-quasiparticle orbitals ...
    @inbounds for a in 1:a_max
        pCl, pCj = pC_l_values[a], pC_j_values[a]
        nCl, nCj = nC_l_values[a], nC_j_values[a]
        @inbounds for b in 1:a_max
            if pC_mask[b] == false && abs(Float64(Orb[b].l) - pCl) < Tol && abs(Float64(Orb[b].j) - pCj) < Tol
                pC_order[b] = a
                pC_mask[b] = true
                break
            end
        end

        @inbounds for b in 1:a_max
            if nC_mask[b] == false && abs(Float64(Orb[b].l) - nCl) < Tol && abs(Float64(Orb[b].j) - nCj) < Tol
                nC_order[b] = a
                nC_mask[b] = true
                break
            end
        end

        pDl, pDj = pD_l_values[a], pD_j_values[a]
        nDl, nDj = nD_l_values[a], nD_j_values[a]
        @inbounds for b in 1:a_max
            if pD_mask[b] == false && abs(Float64(Orb[b].l) - pDl) < Tol && abs(Float64(Orb[b].j) - pDj) < Tol
                pD_order[b] = a
                pD_mask[b] = true
                break
            end
        end

        @inbounds for b in 1:a_max
            if nD_mask[b] == false && abs(Float64(Orb[b].l) - nDl) < Tol && abs(Float64(Orb[b].j) - nDj) < Tol
                nD_order[b] = a
                nD_mask[b] = true
                break
            end
        end
    end

    # Perform reordering of single-quasiparticle orbitals ...
    @views pC .= pC[:,pC_order]
    @views pD .= pD[:,pD_order]

    @views nC .= nC[:,nC_order]
    @views nD .= nD[:,nD_order]

    return pC, nC, pD, nD
end

function HFB_canonical_basis_particle_reordering(Params::Parameters,C::O1B,Orb::Vector{Orb1B})
    # Read calculation parameters ...
    Tol = 1e-3
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1) * (N_max + 2), 2)

    # Read needed arrays ...
    pC, nC = deepcopy(C.p), deepcopy(C.n)

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
            if (pOrb_mask[b] == false) && (abs(Float64(Orb[b].l) - pl) < Tol) && (abs(Float64(Orb[b].j) - pj) < Tol)
                pOrb_order[b] = a
                pOrb_mask[b] = true
                break
            end
        end

        @inbounds for b in 1:a_max
            if (nOrb_mask[b] == false) && (abs(Float64(Orb[b].l) - nl) < Tol) && (abs(Float64(Orb[b].j) - nj) < Tol)
                nOrb_order[b] = a
                nOrb_mask[b] = true
                break
            end
        end
    end

    # Perform reordering of single-particle orbitals ...
    @views pC .= pC[:,pOrb_order]
    @views nC .= nC[:,nOrb_order]

    return O1B(pC,nC)
end

function HF_Orbital_Ordering(Orb::Vector{NOrb},a_max::Int64,ppU::Matrix{Float64},nnU::Matrix{Float64},pSPEnergies::Vector{Float64},nSPEnergies::Vector{Float64})
    pOrb_order = Vector{Int}(undef,a_max)
    nOrb_order = Vector{Int}(undef,a_max)
    pOrb_mask = falses(a_max)
    nOrb_mask = falses(a_max)

    pRhos = @views [ppU[:,a] .* ppU[:,a]' for a in 1:a_max]
    nRhos = @views [nnU[:,a] .* nnU[:,a]' for a in 1:a_max]

    lp_values = [sum(pRhos[a][i,i] * Orb[i].l for i in 1:a_max) for a in 1:a_max]
    jp_values = [sum(pRhos[a][i,i] * Orb[i].j for i in 1:a_max) for a in 1:a_max]
    ln_values = [sum(nRhos[a][i,i] * Orb[i].l for i in 1:a_max) for a in 1:a_max]
    jn_values = [sum(nRhos[a][i,i] * Orb[i].j for i in 1:a_max) for a in 1:a_max]

    @inbounds for a in 1:a_max
        lp, jp = lp_values[a], jp_values[a]
        ln, jn = ln_values[a], jn_values[a]

        @inbounds for c in 1:a_max
            if !pOrb_mask[c] && abs(Orb[c].l - lp) < 1e-8 && abs(Orb[c].j - jp) < 1e-8
                pOrb_order[c] = a
                pOrb_mask[c] = true
                break
            end
        end

        @inbounds for c in 1:a_max
            if !nOrb_mask[c] && abs(Orb[c].l - ln) < 1e-8 && abs(Orb[c].j - jn) < 1e-8
                nOrb_order[c] = a
                nOrb_mask[c] = true
                break
            end
        end
    end

    pSPEnergies .= pSPEnergies[pOrb_order]
    nSPEnergies .= nSPEnergies[nOrb_order]
    ppU .= @views ppU[:, pOrb_order]
    nnU .= @views nnU[:, nOrb_order]

    return ppU, pSPEnergies, nnU, nSPEnergies
end
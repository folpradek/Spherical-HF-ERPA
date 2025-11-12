function HF_Orbital_Ordering(Orb::Vector{NOrb},a_max::Int64,U::pnMatrix,SPEnergies::pnVector)
    pU, nU =  U.p, U.n
    pSPEnergies, nSPEnergies = SPEnergies.p, SPEnergies.n

    pOrb_order, nOrb_order = Vector{Int}(undef,a_max), Vector{Int}(undef,a_max)
    pOrb_mask, nOrb_mask = falses(a_max), falses(a_max)

    pRhos, nRhos = @views [pU[:,a] .* pU[:,a]' for a in 1:a_max], @views [nU[:,a] .* nU[:,a]' for a in 1:a_max]

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
    pU .= @views pU[:, pOrb_order]
    nU .= @views nU[:, nOrb_order]

    return pnMatrix(pU,nU), pnVector(pSPEnergies,nSPEnergies)
end
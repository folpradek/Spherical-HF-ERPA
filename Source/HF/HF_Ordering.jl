function HF_orbital_ordering(Orb::Vector{Orb1B},a_max::Int64,C::O1B,SPE::pnVector)
    # Read HF orbitals & energies ...
    pC, nC =  C.p, C.n
    pSPE, nSPE = SPE.p, SPE.n

    # Initialize arrays for ordering & masks ...
    pOrb_order, nOrb_order = Vector{Int}(undef,a_max), Vector{Int}(undef,a_max)
    pOrb_mask, nOrb_mask = falses(a_max), falses(a_max)

    # Make temporary HF densities ...
    pRho, nRho = @views [pC[:,a] .* pC[:,a]' for a in 1:a_max], @views [nC[:,a] .* nC[:,a]' for a in 1:a_max]

    # Allocate arrays with corresponding values of l & j numbers ...
    lp_values = [sum(pRho[a][i,i] * Orb[i].l for i in 1:a_max) for a in 1:a_max]
    jp_values = [sum(pRho[a][i,i] * Orb[i].j for i in 1:a_max) for a in 1:a_max]
    ln_values = [sum(nRho[a][i,i] * Orb[i].l for i in 1:a_max) for a in 1:a_max]
    jn_values = [sum(nRho[a][i,i] * Orb[i].j for i in 1:a_max) for a in 1:a_max]

    # Find ordering & performing masking ...
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

    # Reshuffle SPE & C ...
    pSPE .= pSPE[pOrb_order]
    nSPE .= nSPE[nOrb_order]
    pC .= @views pC[:,pOrb_order]
    nC .= @views nC[:,nOrb_order]

    return O1B(pC,nC), pnVector(pSPE,nSPE)
end
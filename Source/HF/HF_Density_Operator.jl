function HF_Density_Operator(Orb::Vector{NOrb},a_max::Int64,U::pnMatrix)
    pRho, nRho = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)
    pU, nU = U.p, U.n
    @inbounds for k in 1:a_max
        @inbounds for l in 1:a_max
            if Orb[k].l == Orb[l].l && Orb[k].j == Orb[l].j
                pSum = 0.0
                nSum = 0.0
                @inbounds for m in 1:a_max
                    if (Orb[m].pO == 1) && Orb[m].l == Orb[k].l && Orb[m].j == Orb[k].j
                        pSum += pU[k,m] * pU[l,m]
                    end
                    if (Orb[m].nO == 1) && Orb[m].l == Orb[k].l && Orb[m].j == Orb[k].j
                        nSum += nU[k,m] * nU[l,m]
                    end
                end
                pRho[k,l] = pSum
                nRho[k,l] = nSum
            end
        end
    end
    return pnMatrix(pRho,nRho)
end
function HF_density_operator(a_max::Int64,C::O1B,Orb::Vector{Orb1B})
    # Initialize density Rho ...
    pRho, nRho = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)

    # Allocate density Rho ...
    @inbounds for k in 1:a_max
        @inbounds for l in 1:a_max
            if Orb[k].l == Orb[l].l && Orb[k].j == Orb[l].j
                pSum, nSum = 0.0, 0.0
                @inbounds for m in 1:a_max
                    if (Orb[m].pO == 1) && Orb[m].l == Orb[k].l && Orb[m].j == Orb[k].j
                        pSum += C.p[k,m] * C.p[l,m]
                    end
                    if (Orb[m].nO == 1) && Orb[m].l == Orb[k].l && Orb[m].j == Orb[k].j
                        nSum += C.n[k,m] * C.n[l,m]
                    end
                end
                pRho[k,l] = pSum
                nRho[k,l] = nSum
            end
        end
    end
    
    return O1B(pRho,nRho)
end
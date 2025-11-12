function BCS_Density_Operator(a_max::Int64,V::pnVector)
    # Initialize density matrices ...
    pRho, nRho = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)

    # Compute the 1-body density matrix ... diagonal in canonical basis ...
    @inbounds for a in 1:a_max
        pRho[a,a] = V.p[a]^2
        nRho[a,a] = V.n[a]^2
    end

    return pnMatrix(pRho,nRho)
end

function BCS_Pairing_Operator(a_max::Int64,U::pnVector,V::pnVector)
    # Initialize pairing tensors ...
    pKappa, nKappa = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)

    # Compute the pairing tensor ... diagonal in canonical basis ...
    @inbounds for a in 1:a_max
        pKappa[a,a] = V.p[a] * U.p[a]
        nKappa[a,a] = V.n[a] * U.n[a]
    end

    return pnMatrix(pKappa,nKappa)
end
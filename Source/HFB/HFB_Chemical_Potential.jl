function HFB_Lambda_Secant(Lambda_1::pnFloat,dA_1::pnFloat,Lambda_2::pnFloat,dA_2::pnFloat)
    ###=
    # Basic constants ...
        # Denominator epsilon
    Denom_eps = 1e-7
        # Maximum allowed Lambda step ...
    d_Lambda_max = 3.0
    
    # Initialize Lambda ...
    pLambda, d_pLambda = Lambda_2.p, 0.0
    nLambda, d_nLambda = Lambda_2.n, 0.0

    # Evaluate secant iteration ...
        # Proton secant ...
    pDenom = (dA_2.p - dA_1.p)
    if abs(pDenom) > Denom_eps
        d_pLambda = (Lambda_1.p * dA_2.p - Lambda_2.p * dA_1.p) / pDenom - pLambda
    elseif isfinite(dA_2.p)
        d_pLambda = 0.1 * dA_2.p
    else
        d_pLambda = 0.05 * pLambda
    end
        # Neutron secant ...
    nDenom = (dA_2.n - dA_1.n)
    if abs(nDenom) > Denom_eps
        d_nLambda = (Lambda_1.n * dA_2.n - Lambda_2.n * dA_1.n) / nDenom - nLambda
    elseif isfinite(dA_2.n)
        d_nLambda = 0.1 * dA_2.n
    else
        d_nLambda = 0.05 * nLambda
    end

    # Check if steps are not too large ...
        # Proton step ...
    if abs(d_pLambda) > d_Lambda_max
        d_pLambda = sign(d_pLambda) * d_Lambda_max
    end
        # Neutron step ...
    if abs(d_nLambda) > d_Lambda_max
        d_nLambda = sign(d_nLambda) * d_Lambda_max
    end

    # Update Lambda ...
    pLambda += d_pLambda
    nLambda += d_nLambda

    #=
    pSecant = (Lambda_1.p * dA_2.p - Lambda_2.p * dA_1.p) / (dA_2.p - dA_1.p)
    nSecant = (Lambda_1.n * dA_2.n - Lambda_2.n * dA_1.n) / (dA_2.n - dA_1.n)

    if abs(pSecant) < 3.5
        pLambda = pSecant
    else
        pLambda = Lambda_2.p + 0.25 * dA_2.p
    end

    if abs(nSecant) < 3.5
        nLambda = (Lambda_1.n * dA_2.n - Lambda_2.n * dA_1.n) / (dA_2.n - dA_1.n)
    else
        nLambda = Lambda_2.n + 0.25 * dA_2.n
    end
    =#

    return pnFloat(pLambda,nLambda)
end
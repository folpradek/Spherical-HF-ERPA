function QRPA_OBDM(Params::Parameters,Orb::Vector{Orb1B},Orb_2qp::qpOrb2B,U::O1B,V::O1B,Y_QRPA::Matrix{Matrix{ComplexF64}})
    # Read calculation parameters ..
    N_max = Params.Calc.Nmax
    N_2max = 2*N_max
    J_max = N_2max + 1
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Make the list of values of J & P for iteration ...
    JP_list = JP_initialize(J_max)

    # Evaluate the reference mean-field OBDM ...
    pRho, nRho = V.p * V.p', V.n * V.n'

    pN, nN = Vector{Float64}(undef,a_max), Vector{Float64}(undef,a_max)

    # Calculate the QRPA corrections to the reference mean-field OBDM ...
    @inbounds for JP in JP_list
        J, P = JP[1], JP[2]
        N_qp = Orb_2qp.N[P,J+1]

        J_hat = Float64(2*J+ 1)

        @inbounds for qp in 1:N_qp
            a, b, T_ab = Orb_2qp.i[P,J+1][qp].a, Orb_2qp.i[P,J+1][qp].b, Orb_2qp.i[P,J+1][qp].T
            j_a, j_b = Orb[a].j, Orb[b].j
            j_a_hat, j_b_hat = Float64(j_a + 1), Float64(j_b + 1)

            aAmp = J_hat / j_a_hat
            bAmp = J_hat / j_b_hat

            aSum, bSum = 0.0, 0.0

            @inbounds for nu in 1:N_qp
                aME = aAmp * abs2(Y_QRPA[J+1,P][qp,nu])
                bME = bAmp * abs2(Y_QRPA[J+1,P][qp,nu])

                aSum += aME
                bSum += bME
            end

            if T_ab == -1
                pN[a] += aSum
                pN[b] += bSum
            elseif T_ab == 1
                nN[a] += aSum
                nN[b] += bSum
            end

        end
    end

    # Update the OBDM elements using the calculated QRPA 1-quasiparticle occupations ...

    return O1B(pRho,nRho)
end
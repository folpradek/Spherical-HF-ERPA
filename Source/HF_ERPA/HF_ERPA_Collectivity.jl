function HF_ERPA_Collectivity(Params::Vector{Any},N_nu::Matrix{Int64},X_RPA::Matrix{Matrix{ComplexF64}},Y_RPA::Matrix{Matrix{ComplexF64}})
    # Read calculation parameters ...
    N_max = Params[4]
    N_2max = 2*N_max
    J_max = N_2max + 1

    println("\nEvaluating collectivity of RPA & TDA phonon states ...")

    s = Matrix{Vector{Float64}}(undef,J_max+1,2)

    # "Phonon-states entropy" - meassure of collectivity?
    @inbounds for J in 0:J_max
        @inbounds for P in 1:2
            N_ph = N_nu[J+1,P]
            s_JP = Vector{Float64}(undef,N_ph)
            X_RPA_JP = X_RPA[J+1,P]
            Y_RPA_JP = Y_RPA[J+1,P]

            @inbounds for ph in 1:N_ph
                sSum = 0.0
                @inbounds for qg in 1:N_ph
                    x_RPA = abs(X_RPA_JP[qg, ph])^2
                    y_RPA = abs(Y_RPA_JP[qg, ph])^2
                    if (x_RPA - y_RPA) > 1e-5
                        ME = - (x_RPA - y_RPA) * log(x_RPA - y_RPA)
                        sSum += ME
                    end
                end
                s_JP[ph] = sSum
            end
            s[J+1,P] = s_JP / N_ph
        end
    end

    return s
end
function HF_RPA_Collectivity(Params::Vector{Any},N_nu::Matrix{Int64},X_TDA::Matrix{Matrix{Float64}},X_RPA::Matrix{Matrix{ComplexF64}},Y_RPA::Matrix{Matrix{ComplexF64}})
    # Read calculation parameters ...
    N_max = Params[4]
    N_2max = 2*N_max
    J_max = N_2max + 1

    println("\nEvaluating collectivity of RPA & TDA phonon states ...")

    s_TDA = Matrix{Vector{Float64}}(undef,J_max+1,2)
    s_RPA = Matrix{Vector{Float64}}(undef,J_max+1,2)

    # "Phonon-states entropy" - meassure of collectivity?
    @inbounds for J in 0:J_max
        @inbounds for P in 1:2
            N_ph = N_nu[J+1,P]
            s_TDA_JP = Vector{Float64}(undef,N_ph)
            s_RPA_JP = Vector{Float64}(undef,N_ph)
            X_TDA_JP = X_TDA[J+1,P]
            X_RPA_JP = X_RPA[J+1,P]
            Y_RPA_JP = Y_RPA[J+1,P]


            @inbounds for ph in 1:N_ph
                sSum_TDA = 0.0
                sSum_RPA = 0.0
                @inbounds for qg in 1:N_ph
                    x_TDA = abs(X_TDA_JP[qg, ph])^2
                    if x_TDA > 1e-5
                        ME = - x_TDA * log(x_TDA)
                        sSum_TDA += ME
                    end
                    x_RPA = abs(X_RPA_JP[qg, ph])^2
                    y_RPA = abs(Y_RPA_JP[qg, ph])^2
                    if (x_RPA - y_RPA) > 1e-5
                        ME = - (x_RPA - y_RPA) * log(x_RPA - y_RPA)
                        sSum_RPA += ME
                    end
                end
                s_TDA_JP[ph] = sSum_TDA
                s_RPA_JP[ph] = sSum_RPA
            end
            s_TDA[J+1,P] = s_TDA_JP / N_ph
            s_RPA[J+1,P] = s_RPA_JP / N_ph
        end
    end

    return s_TDA, s_RPA
end
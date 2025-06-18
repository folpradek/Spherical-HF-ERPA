function HF_Export(Params::Vector{Any},pU::Matrix{Float64},nU::Matrix{Float64},pSPEnergies::Vector{Float64},nSPEnergies::Vector{Float64})
    N_max = Params[7]
    Output_File = Params[13]

    a_max = div((N_max + 1)*(N_max + 2),2)

    pU_Export = "IO/" * Output_File * "/Bin/pU.bin"
    open(pU_Export, "w") do Export_File
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                ME = @views pU[a,b]
                write(Export_File, Float64(ME))
            end
        end
    end

    nU_Export = "IO/" * Output_File * "/Bin/nU.bin"
    open(nU_Export, "w") do Export_File
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                ME = @views nU[a,b]
                write(Export_File, Float64(ME))
            end
        end
    end

    pE_Export = "IO/" * Output_File * "/Bin/pE.bin"
    open(pE_Export, "w") do Export_File
        @inbounds for a in 1:a_max
            ME = @views pSPEnergies[a]
            write(Export_File, Float64(ME))
        end
    end

    nE_Export = "IO/" * Output_File * "/Bin/nE.bin"
    open(nE_Export, "w") do Export_File
        @inbounds for a in 1:a_max
            ME = @views nSPEnergies[a]
            write(Export_File, Float64(ME))
        end
    end

    return
end
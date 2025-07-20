function HF_Export(Params::Parameters,U::pnMatrix,SPEnergies::pnVector)
    # Read parameters
    N_max = Params.Calc.Nmax
    Output_File = Params.Calc.Path

    a_max = div((N_max + 1)*(N_max + 2),2)

    # Export densities ...
    pU_Export = "IO/" * Output_File * "/Bin/pU_HF.bin"
    open(pU_Export, "w") do Export_File
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                ME = @views U.p[a,b]
                write(Export_File, Float64(ME))
            end
        end
    end

    nU_Export = "IO/" * Output_File * "/Bin/nU_HF.bin"
    open(nU_Export, "w") do Export_File
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                ME = @views U.n[a,b]
                write(Export_File, Float64(ME))
            end
        end
    end

    pE_Export = "IO/" * Output_File * "/Bin/pE_HF.bin"
    open(pE_Export, "w") do Export_File
        @inbounds for a in 1:a_max
            ME = @views SPEnergies.p[a]
            write(Export_File, Float64(ME))
        end
    end

    nE_Export = "IO/" * Output_File * "/Bin/nE_HF.bin"
    open(nE_Export, "w") do Export_File
        @inbounds for a in 1:a_max
            ME = @views SPEnergies.n[a]
            write(Export_File, Float64(ME))
        end
    end

    return
end
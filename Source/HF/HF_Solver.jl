function HF_Solver(Int_Params::Vector{Any}, NN_File::String, NNN_File::String, Calc_Params::Vector{Any})
    # Calculation parameters
    Params = Any[]
    for k in 1:length(Int_Params)
        Params = push!(Params, Int_Params[k])
    end
    for k in 1:length(Calc_Params)
        Params = push!(Params,Calc_Params[k])
    end
    Params = push!(Params, "A" * string(Params[5]) * "_Z" * string(Params[6]) *
                    "_hw" * string(Params[1]) * "_Nmax" * string(Params[7]) *
                    "_N2max" * string(Params[8]) * "_N3max" * string(Params[9]) *
                    "_" * string(Params[10]))

    println("Starting calculations with NO2B NN+NNN interaction")
    println("\nInteraction data:")
    println("HbarOmega = " * string(Params[1]) * " , N_max = " * string(Params[2]) * " , N_2max = " * string(Params[3]) *
            " , N_3max = " * string(Params[4]) * ", J-scheme basis size = " * string(div((Params[2]+1)*(Params[2]+2),2)) * ", M-scheme basis size = " *
            string(div((Params[2]+1)*(Params[2]+2)*(Params[2]+3),6)))
    println("\nCalculation data:")
    println("A = " * string(Params[5]) * " , Z = " * string(Params[6]) * " , N_max = " * string(Params[7]) *
            " , N_2max = " * string(Params[8]) * " , N_3max = " * string(Params[9]) * ", J-scheme basis size = " *
            string(div((Params[7]+1)*(Params[7]+2),2)) * ", M-scheme basis size = " * string(div((Params[7]+1)*(Params[7]+2)*(Params[7]+3),6)))
    println("Center of mass correction option is set to:    " * string(Params[10]))
    if Params[10] == "CMS1+2B"
        println("\nCombined 1-body + 2-body center of mass motion correction is included ...")
    elseif Params[10] == "CMS2B"
        println("\nOnly pure 2-body center of mass motion correction is included ...")
    else
        println("\nNo center of mass motion correction is included ...")
    end

    if Threads.nthreads() > 1
        println("\n_________________________________________________________")
        println("Multiple active threads detected ...")
        println("Parallelization report:    Number of active threads = " * string(Threads.nthreads()))
        println("_________________________________________________________")
    end

    # Make new directory for results
    if isdir("IO/" * Params[13])
        rm("IO/" * Params[13], recursive = true)
    end
    mkdir("IO/" * Params[13])
    mkdir("IO/" * Params[13] * "/Bin")
    mkdir("IO/" * Params[13] * "/HF")
    mkdir("IO/" * Params[13] * "/HF/Densities")

    # Make orbitals - NuHamil ordering
    Orb = Make_Orbitals(Params[5],Params[6],Params[2])

    # Start calculations
    println("\nStarting spherical Hartree-Fock calculation ...")

    @time VNN, Orb_NN = HF_NO2B(NN_File,NNN_File,Params,Orb)

    if Params[11] == true

        @time VNN, Orb_NN = V2B_Res(Params,Orb,Orb_NN,VNN)

        @time HF_MBPT(Params,Orb,Orb_NN,VNN)

        @time V2B_Res_Export(Params,Orb_NN,VNN)

        if Params[12] == "HRBin" || Params[12] == "HR"
            @time Orbitals_Export(Params,Orb)
        end

    end
    
    println("\nAll calculations have finished ...\n")

    return
    
end
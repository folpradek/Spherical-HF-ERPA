function HF_RPA_Solver(Input_File::String,Calc_Params::Vector{Any})
    # Calculation parameters
    Params = Any[]
    for k in 1:length(Calc_Params)
        Params = push!(Params,Calc_Params[k])
    end
    Params = push!(Params, Input_File)

    println("Starting RPA (TDA) calculations with residual NO2B NN+NNN interaction")
    println("\nCalculation data:")
    println("A = " * string(Params[1]) * " , Z = " * string(Params[2]) * " , HbarOmega = " * string(Params[3]) *
            " MeV , N_max = " * string(Params[4]) * " , J-scheme LHO basis size = " * string(div((Params[4]+1)*(Params[4]+2),2)) *
            " , M-scheme LHO basis size = " * string(div((Params[4]+1)*(Params[4]+2)*(Params[4]+3),6)))

    if Threads.nthreads() > 1
        println("\n_________________________________________________________")
        println("Multiple active threads detected ...")
        println("Parallelization report:    Number of active threads = " * string(Threads.nthreads()))
        println("_________________________________________________________")
    end

    # Make new directories for results ...
    if !(isdir(Params[8]))
        println("\nError! ... No precomputed HF solution is available in given Input_File path ... run HF solver first ...")
        return
    end
    if (isdir(Params[8] * "/RPA"))
        rm(Params[8] * "/RPA", recursive = true)
    end
    mkdir(Params[8] * "/RPA")
    if !(isdir(Params[8] * "/RPA/Densities"))
        mkdir(Params[8] * "/RPA/Densities")
    end
    if !(isdir(Params[8] * "/RPA/Transitions"))
        mkdir(Params[8] * "/RPA/Transitions")
        mkdir(Params[8] * "/RPA/Transitions/E0")
        mkdir(Params[8] * "/RPA/Transitions/E1")
        mkdir(Params[8] * "/RPA/Transitions/E2")
        mkdir(Params[8] * "/RPA/Transitions/E3")
    end
    if !(isdir(Params[8] * "/RPA/Spectra"))
        mkdir(Params[8] * "/RPA/Spectra")
    end
    if !(isdir(Params[8] * "/RPA/Amplitudes"))
        mkdir(Params[8] * "/RPA/Amplitudes")
    end
    # Make s.p. orbitals - NuHamil ordering convention ...
    Orb = Make_Orbitals(Params[1],Params[2],Params[4])

    # Import Residual 2-body Interaction ...
    println("\nImporting Residual 2-body interaction ...")
    @time VNN, Orb_NN = V2B_Res_Import(Params[4],Params[8],Orb)
    
    println("\nStarting RPA & TDA calculations ...")
    # Start RPA & TDA calculations ...
    @time HF_RPA(Params,Orb,Orb_NN,VNN)

    
    println("\nAll calculations have finished ...\n")

    return
end
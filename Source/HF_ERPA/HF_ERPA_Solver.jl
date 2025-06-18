function HF_ERPA_Solver(Int_Params::Vector{Any},NN_File::String,NNN_File::String,Input_File::String,Calc_Params::Vector{Any})
    # Calculation parameters
    Params = Any[]
    for k in 1:length(Calc_Params)
        Params = push!(Params,Calc_Params[k])
    end
    Params = push!(Params, Input_File)

    println("Starting Extended-RPA calculations with residual NO2B NN+NNN interaction")
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
    if (isdir(Params[8] * "/ERPA"))
        rm(Params[8] * "/ERPA", recursive = true)
    end
    mkdir(Params[8] * "/ERPA")
    if !(isdir(Params[8] * "/ERPA/Densities"))
        mkdir(Params[8] * "/ERPA/Densities")
    end
    if !(isdir(Params[8] * "/ERPA/Transitions"))
        mkdir(Params[8] * "/ERPA/Transitions")
        mkdir(Params[8] * "/ERPA/Transitions/E0")
        mkdir(Params[8] * "/ERPA/Transitions/E1")
        mkdir(Params[8] * "/ERPA/Transitions/E2")
        mkdir(Params[8] * "/ERPA/Transitions/E3")
    end
    if !(isdir(Params[8] * "/ERPA/Spectra"))
        mkdir(Params[8] * "/ERPA/Spectra")
    end
    if !(isdir(Params[8] * "/ERPA/Amplitudes"))
        mkdir(Params[8] * "/ERPA/Amplitudes")
    end

    # Make s.p. orbitals - NuHamil ordering convention ...
    Orb = Make_Orbitals(Params[1],Params[2],Int_Params[2])

    # Import Residual 2-body Interaction ...
    println("\nImporting Residual 2-body interaction ...")
    @time VNN, Orb_NN = V2B_Res_Import(Params[4],Params[8],Orb)
    
    println("\nStarting ERPA calculations ...")
    # Start RPA & TDA calculations ...
    @time HF_ERPA(Int_Params,NN_File,NNN_File,Params,Orb,Orb_NN,VNN)

    println("\nAll calculations have finished ...\n")

    return
end
function HF_RPA_Solver(Params::Parameters)
    # Initialize Wigner symbols ...
    wigner_init_float(Params.Int.N2max + 2, "Jmax", 6)

    # Calculation parameters
    println("Starting RPA (TDA) calculations with residual NO2B NN+NNN interaction")
    println("\nCalculation data:")
    println("A = " * string(Params.Calc.A) * " , Z = " * string(Params.Calc.Z) * " , HbarOmega = " * string(Params.Calc.Z) *
            " MeV , N_max = " * string(Params.Calc.Nmax) * " , J-scheme LHO basis size = " * string(div((Params.Calc.Nmax+1)*(Params.Calc.Nmax+2),2)) *
            " , M-scheme LHO basis size = " * string(div((Params.Calc.Nmax+1)*(Params.Calc.Nmax+2)*(Params.Calc.Nmax+3),6)))

    if Threads.nthreads() > 1
        println("\n_________________________________________________________")
        println("Multiple active threads detected ...")
        println("Parallelization report:    Number of active threads = " * string(Threads.nthreads()))
        println("_________________________________________________________")
    end

    # Make new directories for results ...
    println(Params.Calc.Path)
    if !(isdir("IO/" * Params.Calc.Path))
        println("\nError! ... No precomputed HF solution is available in given Input_File path ... run HF solver first ...")
        return
    end
    if (isdir("IO/" * Params.Calc.Path * "/RPA"))
        rm("IO/" * Params.Calc.Path * "/RPA", recursive = true)
    end
    mkdir("IO/" * Params.Calc.Path * "/RPA")
    if !(isdir("IO/" * Params.Calc.Path * "/RPA/Densities"))
        mkdir("IO/" * Params.Calc.Path * "/RPA/Densities")
    end
    if !(isdir("IO/" * Params.Calc.Path * "/RPA/Transitions"))
        mkdir("IO/" * Params.Calc.Path * "/RPA/Transitions")
        mkdir("IO/" * Params.Calc.Path * "/RPA/Transitions/E0")
        mkdir("IO/" * Params.Calc.Path * "/RPA/Transitions/E1")
        mkdir("IO/" * Params.Calc.Path * "/RPA/Transitions/E2")
        mkdir("IO/" * Params.Calc.Path * "/RPA/Transitions/E3")
    end
    if !(isdir("IO/" * Params.Calc.Path * "/RPA/Spectra"))
        mkdir("IO/" * Params.Calc.Path * "/RPA/Spectra")
    end
    if !(isdir("IO/" * Params.Calc.Path * "/RPA/Amplitudes"))
        mkdir("IO/" * Params.Calc.Path * "/RPA/Amplitudes")
    end
    
    # Make s.p. orbitals - NuHamil ordering convention ...
    Orb = Make_Orbitals(Params.Calc.A,Params.Calc.Z,Params.Int.Nmax)

    # Import Residual 2-body Interaction ...
    println("\nImporting Residual 2-body interaction ...")
    @time VNN, Orb_NN = V2B_Res_Import(Params,Orb)
    
    println("\nStarting RPA & TDA calculation ...")
    # Start RPA & TDA calculation ...
    @time HF_RPA(Params,Orb,Orb_NN,VNN)

    
    println("\nAll calculations have finished ...\n")

    return
end
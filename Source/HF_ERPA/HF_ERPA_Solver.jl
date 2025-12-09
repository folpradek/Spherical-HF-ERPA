function HF_ERPA_Solver(Params::Parameters)
    # Initialize Wigner symbols ...
    wigner_init_float(Params.Int.N2max + 2, "Jmax", 6)

    # Calculation parameters
    println("Starting Extended-RPA calculations with residual NO2B NN+NNN interaction")
    println("\nCalculation data:")
    println("A = " * string(Params.Calc.A) * " , Z = " * string(Params.Calc.Z) * " , HbarOmega = " * string(Params.Calc.hw) *
            " MeV , N_max = " * string(Params.Calc.Nmax) * " , J-scheme LHO basis size = " * string(div((Params.Calc.Nmax+1)*(Params.Calc.Nmax+2),2)) *
            " , M-scheme LHO basis size = " * string(div((Params.Calc.Nmax+1)*(Params.Calc.Nmax+2)*(Params.Calc.Nmax+3),6)))

    if Threads.nthreads() > 1
        println("\n_________________________________________________________")
        println("Multiple active threads detected ...")
        println("Parallelization report:    Number of active threads = " * string(Threads.nthreads()))
        println("_________________________________________________________")
    end

   # Make new directories for results ...
    println("\nInput/Output path to the calculation reads ... ''IO/" * Params.Calc.Path * "''")
    if !(isdir("IO/" * Params.Calc.Path))
        println("\nError! ... No precomputed HF solution is available in given Input_File path ... run HF solver first ...")
        return
    end
    if (isdir("IO/" * Params.Calc.Path * "/ERPA"))
        rm("IO/" * Params.Calc.Path * "/ERPA", recursive = true)
    end
    mkdir("IO/" * Params.Calc.Path * "/ERPA")
    if !(isdir("IO/" * Params.Calc.Path * "/ERPA/Densities"))
        mkdir("IO/" * Params.Calc.Path * "/ERPA/Densities")
    end
    if !(isdir("IO/" * Params.Calc.Path * "/ERPA/Transitions"))
        mkdir("IO/" * Params.Calc.Path * "/ERPA/Transitions")
        mkdir("IO/" * Params.Calc.Path * "/ERPA/Transitions/E0")
        mkdir("IO/" * Params.Calc.Path * "/ERPA/Transitions/E1")
        mkdir("IO/" * Params.Calc.Path * "/ERPA/Transitions/E2")
        mkdir("IO/" * Params.Calc.Path * "/ERPA/Transitions/E3")
    end
    if !(isdir("IO/" * Params.Calc.Path * "/ERPA/Spectra"))
        mkdir("IO/" * Params.Calc.Path * "/ERPA/Spectra")
    end
    if !(isdir("IO/" * Params.Calc.Path * "/ERPA/Amplitudes"))
        mkdir("IO/" * Params.Calc.Path * "/ERPA/Amplitudes")
    end

    # Make s.p. orbitals - NuHamil ordering convention ... Interaction basis size ...
    Orb = Make_Orbitals(Params.Calc.A,Params.Calc.Z,Params.Int.Nmax)

    # Make 1-body kinetic operator ...
    @time T = T1B(Params.Int.Nmax,Orb,Params.Int.hw)

    # 2-body NN interaction & Orbitals ...
    @time VNN_bare, Orb_NN_bare = V2B_Read(Params,Orb)

    # 3-body NNN interaction & Orbitals ...
    @time VNNN_bare, Orb_NNN_bare = V3B_NO2B_Read(Params,Orb)

    # Import Residual 2-body Interaction ...
    #println("\nImporting Residual 2-body interaction ...")
    #@time VNN_res, Orb_NN_res = V2B_Res_Import(Params,Orb)
    
    # Start RPA & TDA calculation ...
    println("\nStarting ERPA calculation ...")
    @time HF_ERPA(Params,Orb,Orb_NN_bare,Orb_NNN_bare,T,VNN_bare,VNNN_bare)

    println("\nAll calculations have finished ...\n")

    return
end
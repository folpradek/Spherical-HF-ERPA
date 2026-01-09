function QTDA(Params::Parameters)
    # Initialize Wigner symbols ...
    wigner_init_float(Params.Int.N2max + 2, "Jmax", 6)

    # Calculation parameters
    println("Starting QTDA calculation with residual NO2B NN+NNN interaction")
    println("\nCalculation data:")
    println("A = " * string(Params.Calc.A) * " , Z = " * string(Params.Calc.Z) * " , hw = " * string(Params.Calc.Z) *
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
        println("\nError! ... No precomputed reference mean-field (HFB, BCS) solution is available in given Input_File path ... run HFB/BCS solver first ...")
        return
    end
    if (isdir("IO/" * Params.Calc.Path * "/QTDA"))
        rm("IO/" * Params.Calc.Path * "/QTDA", recursive = true)
    end
    mkdir("IO/" * Params.Calc.Path * "/QTDA")
    if !(isdir("IO/" * Params.Calc.Path * "/QTDA/Densities"))
        mkdir("IO/" * Params.Calc.Path * "/QTDA/Densities")
    end
    if !(isdir("IO/" * Params.Calc.Path * "/QTDA/Transitions"))
        mkdir("IO/" * Params.Calc.Path * "/QTDA/Transitions")
        mkdir("IO/" * Params.Calc.Path * "/QTDA/Transitions/E0")
        mkdir("IO/" * Params.Calc.Path * "/QTDA/Transitions/E1")
        mkdir("IO/" * Params.Calc.Path * "/QTDA/Transitions/E2")
        mkdir("IO/" * Params.Calc.Path * "/QTDA/Transitions/E3")
    end
    if !(isdir("IO/" * Params.Calc.Path * "/QTDA/Spectra"))
        mkdir("IO/" * Params.Calc.Path * "/QTDA/Spectra")
    end
    if !(isdir("IO/" * Params.Calc.Path * "/QTDA/Amplitudes"))
        mkdir("IO/" * Params.Calc.Path * "/QTDA/Amplitudes")
    end
    
    println("\nStarting QTDA calculation ...")
    # Start QTDA calculation ...
    @time QTDA_solver(Params)
    
    println("\nAll calculations have finished ...\n")

    return 
end



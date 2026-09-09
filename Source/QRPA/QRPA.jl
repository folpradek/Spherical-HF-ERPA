function QRPA(Params::Parameters)
    # Initialize Wigner symbols ...
    wigner_init_float(75,"Jmax",9)

    # Calculation parameters ...
    println("Starting spherical QRPA calculation with residual NO2B NN+NNN interaction based on available reference quasiparticle mean-field calculation")
    println("\nCalculation data:")
    println("A = " * string(Params.Calc.A) * " , Z = " * string(Params.Calc.Z) * " , hw = " * string(Params.Calc.hw) *
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
        println("\nError! ... No precomputed reference quasiparticle mean-field (HFB, BCS) solution is available in given Input_File path ... run HFB/BCS solver first ...")
        return
    end
    if (isdir("IO/" * Params.Calc.Path * "/QRPA"))
        rm("IO/" * Params.Calc.Path * "/QRPA", recursive = true)
    end
    mkdir("IO/" * Params.Calc.Path * "/QRPA")
    if !(isdir("IO/" * Params.Calc.Path * "/QRPA/Densities"))
        mkdir("IO/" * Params.Calc.Path * "/QRPA/Densities")
    end
    if !(isdir("IO/" * Params.Calc.Path * "/QRPA/Transitions"))
        mkdir("IO/" * Params.Calc.Path * "/QRPA/Transitions")
        mkdir("IO/" * Params.Calc.Path * "/QRPA/Transitions/E0")
        mkdir("IO/" * Params.Calc.Path * "/QRPA/Transitions/E1")
        mkdir("IO/" * Params.Calc.Path * "/QRPA/Transitions/E2")
        mkdir("IO/" * Params.Calc.Path * "/QRPA/Transitions/E3")
    end
    if !(isdir("IO/" * Params.Calc.Path * "/QRPA/Spectra"))
        mkdir("IO/" * Params.Calc.Path * "/QRPA/Spectra")
    end
    if !(isdir("IO/" * Params.Calc.Path * "/QRPA/Amplitudes"))
        mkdir("IO/" * Params.Calc.Path * "/QRPA/Amplitudes")
    end
    
    println("\nStarting QRPA calculation ...")

    # Start QRPA calculation ...
    @time QRPA_solver(Params)
    
    println("\nSpherical QRPA solver run completed ...\n")

    return 
end


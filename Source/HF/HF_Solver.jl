function HF_Solver(Params::Parameters)
    # Initialize Wigner symbols ...
    wigner_init_float(Params.Int.N2max + 2, "Jmax", 6)

    # Calculation parameters
    println("Starting calculations with NO2B NN+NNN interaction")

    println("\nInteraction data:")
    println("HbarOmega = " * string(Params.Int.hw) * " , N_max = " * string(Params.Int.Nmax) * " , N_2max = " * string(Params.Int.N2max) *
            " , N_3max = " * string(Params.Int.N3max) * ", J-scheme basis size = " * string(div((Params.Int.Nmax+1)*(Params.Int.Nmax+2),2)) *
            ", M-scheme basis size = " * string(div((Params.Int.Nmax+1)*(Params.Int.Nmax+2)*(Params.Int.Nmax+3),6)))

    println("\nCalculation data:")
    println("A = " * string(Params.Calc.A) * " , Z = " * string(Params.Calc.Z) * " , N_max = " * string(Params.Calc.Nmax) *
            " , N_2max = " * string(Params.Calc.N2max) * " , N_3max = " * string(Params.Calc.N3max) * ", J-scheme basis size = " *
            string(div((Params.Calc.Nmax+1)*(Params.Calc.Nmax+2),2)) * ", M-scheme basis size = " *
            string(div((Params.Calc.Nmax+1)*(Params.Calc.Nmax+2)*(Params.Calc.Nmax+3),6)))

    println("Center of mass correction option is set to:    " * Params.Calc.CMS)
    if Params.Calc.CMS == "CMS1+2B"
        println("\nCombined 1-body + 2-body center of mass motion correction is included ...")
    elseif Params.Calc.CMS == "CMS2B"
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
    if isdir("IO/" * Params.Calc.Path)
        rm("IO/" * Params.Calc.Path, recursive = true)
    end
    mkdir("IO/" * Params.Calc.Path)
    mkdir("IO/" * Params.Calc.Path * "/Bin")
    mkdir("IO/" * Params.Calc.Path * "/HF")
    mkdir("IO/" * Params.Calc.Path * "/HF/Densities")

    # Make orbitals - NuHamil ordering
    Orb = Make_Orbitals(Params.Calc.A,Params.Calc.Z,Params.Int.Nmax)

    # Start calculations
    println("\nStarting spherical Hartree-Fock calculation ...")

    @time VNN, Orb_NN = HF_NO2B(Params,Orb)

    if Params.Calc.BMF == true

        U = pnMatrix(Read_Transformation_Matrix(Params,"IO/" * Params.Calc.Path * "/Bin/pU_HF.bin"),
                     Read_Transformation_Matrix(Params,"IO/" * Params.Calc.Path * "/Bin/nU_HF.bin"))

        @time VNN, Orb_NN = V2B_Res(Params,Orb,Orb_NN,VNN,U)

        @time HF_MBPT(Params,Orb,Orb_NN,VNN)

        @time V2B_Res_Export(Params,Orb_NN,VNN)

        if Params.Calc.Format == "HRBin" || Params.Calc.Format == "HR"
            @time Orbitals_Export(Params,Orb)
        end

    end
    
    println("\nAll calculations have finished ...\n")

    return
end
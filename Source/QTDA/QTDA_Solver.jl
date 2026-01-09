function QTDA_solver(Params::Parameters)
    # Make s.p. orbitals - NuHamil ordering convention ...
    @time Orb = orbitals_make(Params)

    # Make 2qp orbitals  ...
    Orb_2qp = orbitals_2qp_make(Params,Orb)

    # Import LHO -> (HFB) quasiparticle mean-field transformation matrices ...
    @time C = O1b_import(Params,Orb,"IO/" * Params.Calc.Path * "/Bin/C.bin")

    # Import U & V quasiparticle amplitudes ...
    @time U = O1b_import(Params,Orb,"IO/" * Params.Calc.Path * "/Bin/U.bin")
    @time V = O1b_import(Params,Orb,"IO/" * Params.Calc.Path * "/Bin/V.bin")

    # Import H_N ... 1-body quasiparticle Hamiltonian ...
    @time H_N = qpO1B_import(Params,Orb,"IO/" * Params.Calc.Path * "/Bin/qpH1B.bin")

    # Import residual 2-body NN interaction quasiparticle Hamiltonian ...
    @time H_NN, Orb_NN = qpO2b_import(Params,Orb,"IO/" * Params.Calc.Path * "/Bin/qpH2B.bin")

    # Initialize the quasiparticle transition operators ...
        # THIS ONE NEEDS TO BE ADJUSTED TO QP FORMALISM ... QP 1-BODY TRANSITION OPERATORS ... TO BE DONE ...
    #@time qpTrOp = Tr1b_initialize(Params,Orb)

    # Allocate the QTDA matrix A ...
    @time A = QTDA_allocate(Params,Orb,Orb_NN,Orb_2qp,H_N,H_NN)

    # Solve QTDA eigenvalue problem ...
    @time E_QTDA, X_QTDA= QTDA_diagonalize(Params,Orb,Orb_2qp,A)

    # Export of QTDA solutions ...
    @time QTDA_export(Params,Orb,Orb_2qp,E_QTDA,X_QTDA)

    println("\nQTDA solver executed properly ...\n")

    return
end
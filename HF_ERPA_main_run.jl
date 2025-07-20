include("Source/Import_Module.jl")
using .HF_RPA_Module

function run_main()

    # Interaction parameters structure ...
    IntParams = Interaction_Parameters(
                NN_File = "IO/NN.bin",          #   Path to the NN interaction binary file
                NNN_File = "IO/NNN.bin",        #   Path to the NNN interaction binary file
                hw = 16.0,
                Nmax = 3,
                N2max = 6,
                N3max = 9
                )

    # Calculation parameters structure ...
    CalcParams = Calculation_Parameters(
                A = 16,
                Z = 8,
                hw = 16.0,
                Nmax = 3,
                N2max = 6,
                N3max = 9,
                Ortho = true,                                                               # Optional parameter ... Orthogonalization of spurious 1- state - true/false, by default set to true
                CMS = "CMS1+2B",                                                            # Options of CMS corrections ...  "CMS1+2B",  "CMS2B" or "NoCMS"
                BMF = true,                                                                 # HF-MBPT calculation & export of residual interaction ... true/false
                #Path = "DN2LO_GO_394_A4_Z2_hw16.0_Nmax3_N2max6_N3max9_CMS1+2B",            # Optional parameter ... Calculation data path ... under "IO/"
                Format = "Bin",                                                             # Optional parameter ... Export of residual interaction & Orbitals - "Bin", "HRBin" & "HR"
                ERPA = Parameters_ERPA(OBDM = "Full", ScOBH = true, Sc3N = true),           # Optional Parameters of ERPA ... OBDM iteration mode ... "Diagonal" / "Full", Selfconsistent 1-body Hamiltonian & NNN interaction true/false
                #cV_res = 0.01                                                              # Optional parameter ... Rescalling parameter for the residual interaction ...
                )

    # HF calculation call ...
    HF_Solver(Parameters(IntParams,CalcParams))

    # HF-RPA calculation call ...
    HF_RPA_Solver(Parameters(IntParams,CalcParams))

    # HF-ERPA calculation call ...
    HF_ERPA_Solver(Parameters(IntParams,CalcParams))

end

@time run_main()
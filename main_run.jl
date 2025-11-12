include("Source/Import_NuJuliet_Module.jl")
using .NuJuliet

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
                A = 18,
                Z = 8,
                hw = 16.0,
                Nmax = 3,
                N2max = 6,
                N3max = 9,
                Ortho = true,
                CMS = "CMS1+2B",
                BMF = true,
                #Path = "TestRun_DN2LOGO394_A16_Z8_hw16.0_Nmax3_N2max6_N3max9_CMS1+2B",
                #Format = "Bin",
                #ERPA = Parameters_ERPA(OBDM = "Full", ScOBH = true, Sc3N = true),
                #cV_res = 0.01,
                Pairing = Parameters_Pairing(ScBCS = false),
                HFB = Parameters_HFB(Broyden = true)
                )

    # HF calculation call ...
    #HF(Parameters(IntParams,CalcParams))

    # HF-RPA calculation call ...
    #pHF_RPA_Solver(Parameters(IntParams,CalcParams))

    # HF-ERPA calculation call ...
    #HF_ERPA_Solver(Parameters(IntParams,CalcParams))

    # HFB calculation call ...
    #BCS(Parameters(IntParams,CalcParams))

    # HFB calculation call ...
    HFB(Parameters(IntParams,CalcParams))

end

@time run_main()
include("../Source/Import_NuJuliet_Module.jl")
using .NuJuliet

function test_run()

    # Interaction parameters structure ...
    IntParams = Interaction_Parameters(
                NN_File = "../../../../storage/folprecht/ChiPots/NN/NN-DN2LOGO450_vbare2.0_hw16_emax14_e2max28.me2j.bin",
                NNN_File = "../../../../storage/folprecht/ChiPots/NNN/NNN_DN2LOGO450_hw16_Nmax14_N2max28_N3max26_no2b.stream.bin",
                hw = 16.0,
                Nmax = 14,
                N2max = 28,
                N3max = 26
                )

    # Calculation parameters structure ...
    CalcParams = Calculation_Parameters(
                A = 16,
                Z = 8,
                hw = 16.0,
                Nmax = 10,
                N2max = 20,
                N3max = 18,
                CMS = "CMS1+2B",
                Path = "DN2LOGO450_A16_Z8_hw16.0_Nmax10_N2max20_N3max18",
                HF = HF_Parameters(Tol = 1e-7, IMax = 100, BMF = true, HRF = false),
                RPA = RPA_Parameters(Ortho = true),
                RRPA = RRPA_Parameters(Ortho = true, Tol = 1e-7, IMax = 100, ScV3N = true, ScOBH = true, dOBDM = false)
                )

    # HF calculation call ...
    HF(Parameters(IntParams,CalcParams))

    # HF-RRPA calculation call ...
    #HF_RRPA(Parameters(IntParams,CalcParams))

end

test_run()

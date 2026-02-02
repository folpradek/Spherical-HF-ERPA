include("../Source/Import_NuJuliet_Module.jl")
using .NuJuliet

# Note all the commented parameter options are optional and
# can be safely omitted for solver calls - default values
# will be used instead ...
function sample_script_run()

    # Interaction parameters structure ...
    IntParams = Interaction_Parameters(
                NN_File = "IO/NN.bin",
                NNN_File = "IO/NNN.bin",
                hw = 16.0,
                Nmax = 3,
                N2max = 6,
                N3max = 9,
                #cRes = 1.0,
                #cP2N = 1.0,
                #cP3N = 1.0
                )

    # Calculation parameters structure ...
    CalcParams = Calculation_Parameters(
                A = 16,
                Z = 8,
                hw = 16.0,
                Nmax = 3,
                N2max = 6,
                N3max = 9,
                #CMS = "CMS1+2B",
                #Path = "A16_Z8_hw16.0_Nmax3_N2max6_N3max9",
                #HF = HF_Parameters(Tol = 1e-7, Imax = 100, BMF = true, HRF = false, EFA = HF_EFA_Parameters(A = 0, Z = 0, BMF = false)),
                #RPA = RPA_Parameters(Ortho = true),
                #RRPA = RRPA_Parameters(Ortho = true, Tol = 1e-7, Imax = 100, ScV3N = true, ScOBH = true, dOBDM = false)
                )

    # HF calculation call ...
    HF(Parameters(IntParams,CalcParams))

    # HF-RPA calculation call ...
    HF_RPA(Parameters(IntParams,CalcParams))

    # HF-RRPA calculation call ...
    HF_RRPA(Parameters(IntParams,CalcParams))

    CalcParams = Calculation_Parameters(
                A = 18,
                Z = 8,
                hw = 16.0,
                Nmax = 3,
                N2max = 6,
                N3max = 9,
                #CMS = "CMS1+2B",
                #Path = "ScBCS_A18_Z8_hw16.0_Nmax3_N2max6_N3max9",
                #HF = HF_Parameters(Tol = 1e-7, IMax = 100),
                #BCS = BCS_Parameters(ScBCS = true, Tol = 1e-7, Imax = 500, q = 0.05, pD0 = 0.5, nD0 = 0.5, BMF = false, HRF = false)
                )

    # HF-BCS calculation call ...
    BCS(Parameters(IntParams,CalcParams))

    CalcParams = Calculation_Parameters(
                A = 18,
                Z = 8,
                hw = 16.0,
                Nmax = 3,
                N2max = 6,
                N3max = 9,
                #CMS = "CMS1+2B",
                #Path = "HFB_A18_Z8_hw16.0_Nmax3_N2max6_N3max9",
                #HFB = HFB_Parameters(Tol = 1e-8, Imax = 150, Pairing = "Full", dLmax = 0.5, BMF = true, HRF = false),
                #QTDA = QTDA_Parameters(Ortho = true)
                )

    # HFB calculation call ...
    HFB(Parameters(IntParams,CalcParams))

    # QTDA calculation call ...
    QTDA(Parameters(IntParams,CalcParams))

end

sample_script_run()

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
                #cP3N = 1.0,
                #pG2N = 0.0,
                #nG2N = 0.0
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
    
    # Function for export of transition currents ... works only for TDA phonons so far ... no spin current yet ...
        # JP specifies spin & partiy, M specifies projection, I suggest using M = J, nu_list the list of phonon to plot,
        # note phonons are ordered by energy, so nu = 1 is the lowest phonon of given JP, nu = 2 the second lowest ...,
        # File_Name is the name of the output files in .dat format
    HF_RPA_transition_currents(Parameters(IntParams,CalcParams); JP = "1-", M = 1, nu_list=[6,7,8,9,10], File_Name ="Phonon_Currents")


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
                #HFB = HFB_Parameters(Tol = 1e-8, Imax = 150, Pairing = "Full", dLmax = 0.5, BMF = true, LNT = false, LNRes = false, HRF = false),
                #QTDA = QTDA_Parameters(Ortho = true)
                #QRPA = QRPA_Parameters(Ortho = true)
                )

    # HFB calculation call ...
    HFB(Parameters(IntParams,CalcParams))

    # QTDA calculation call ...
    QTDA(Parameters(IntParams,CalcParams))
    
    # QRPA calculation call ...
    QRPA(Parameters(IntParams,CalcParams))

end

sample_script_run()

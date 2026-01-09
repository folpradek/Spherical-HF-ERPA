include("../Source/Import_NuJuliet_Module.jl")
using .NuJuliet

function test_run()

    p_list = [1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0,3.25,3.5,0.75]

    for p in p_list

        # Interaction parameters structure ...
        IntParams = Interaction_Parameters(
                    NN_File = "../../../../storage/folprecht/ChiPots/NN/NN-DN2LOGO450_vbare2.0_hw16_emax14_e2max28.me2j.bin",
                    NNN_File = "../../../../storage/folprecht/ChiPots/NNN/NNN_DN2LOGO450_hw16_Nmax14_N2max28_N3max26_no2b.stream.bin",
                    hw = 16.0,
                    Nmax = 14,
                    N2max = 28,
                    N3max = 26,
                    cP2N = p,
                    cP3N = p
                    )

        # Calculation parameters structure ...
        CalcParams = Calculation_Parameters(
                    A = 18,
                    Z = 8,
                    hw = 16.0,
                    Nmax = 10,
                    N2max = 20,
                    N3max = 18,
                    CMS = "CMS1+2B",
                    Path = "HFB_A18_Z8_DN2LOGO450_hw16.0_Nmax10_N2max20_N3max18_c" * string(p),
                    HF = HF_Parameters(Tol = 1e-7, IMax = 100),
                    BCS = BCS_Parameters(ScBCS = true, Tol = 1e-7, IMax = 500),
                    HFB = HFB_Parameters(Tol = 1e-7, IMax = 150)
                    )

        # HFB calculation call ...
        HFB(Parameters(IntParams,CalcParams))

        CalcParams = Calculation_Parameters(
                    A = 38,
                    Z = 20,
                    hw = 16.0,
                    Nmax = 10,
                    N2max = 20,
                    N3max = 18,
                    CMS = "CMS1+2B",
                    Path = "HFB_A38_Z20_DN2LOGO450_hw16.0_Nmax10_N2max20_N3max18_c" * string(p),
                    HF = HF_Parameters(Tol = 1e-7, IMax = 100),
                    BCS = BCS_Parameters(ScBCS = true, Tol = 1e-7, IMax = 500),
                    HFB = HFB_Parameters(Tol = 1e-7, IMax = 150)
                    )

        # HFB calculation call ...
        HFB(Parameters(IntParams,CalcParams))

    end

end

test_run()

#______________________________________________________________________________________________________________________________________________________
# Parameter menu - HF_Solver()
#
# Int_Params      - Interaction parameters given in array [HbarOmega(MeV), N_max, N_2max N_3max],                                                     
#                   these have to be same for both NN & NNN interaction.
# NN_File         - path to the NN interaction file - uses NuHamil binary format.
# NNN_File        - path to the NNN interaction file - uses NuHamil binary format. For no NNN interaction modify V3B_NO2B_Read.
# Calc_Params     - Array of individual calculation parameters, given as follows [A, Z, N_max, N_2max, N_3max, CMS, MBPT+V_res, Format],
#                   where CMS takes values:   "CMS1+2B", "CMS2B" - for inclusion of 1-body + 2-body Center of Mass Correction,
#                   or to use pure 2-body Center of Mass Motion Correction, MBPT+Beyond takes values: true/false for calculation of
#                   Many-Body Perturbation Theory + export of Residual 2-body Interaction, Format defines output format for
#                   residual 2-body interaction & export of s.p. orbitals, has 3 options:
#                       (1) Format = "Bin"   - uses internal binary format for IO, no orbitals are exported.
#                       (2) Format = "HR"    - uses human readable format for IO files - can be read with notepad.
#                       (3) Format = "HRBin" - uses human readable format for IO file, binary interaction file with no header.
#
# Parameter menu  - HF_RPA_Solver()
#
# RPA_Input_Path  - Path to HF calculation output files - HF transformation matrices, single-particle energies and residual 2-body interaction.
#                   Example: "Output/A16_Z8_hw16.0_Nmax3_N2max6_N3max9_CMS1+2B"
# RPA_Calc_Params - Array of individual calculation parameters, given as follows [A, Z, HbarOmega(MeV), N_max, Orthogonalization, [E_min, E_max, Delta], [E_min, E_max, Delta]],
#                   where N_max is the number of oscillatory quanta used in HF calculation, and Orthogonalization enables/disables treatment of spurious
#                   1- states from RPA & TDA spectra by the means of Orthogonalization Procedure. The last two vectors define energy range and smearing width Delta
#                   for S0 strenght functions and total photoabsorbtion cross section in units of MeV.
#
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

include("Source/Import_Module.jl")
using .HF_RPA_Module

function HF_ERPA_dev_main()

    Int_Params = Any[16.0, 3, 6, 9]
    NN_File = "IO/NN.bin"
    NNN_File = "IO/NNN.bin"
    Calc_Params = Any[16, 8, 3, 6, 9, "CMS1+2B", true, "Bin"]
    HF_Solver(Int_Params, NN_File, NNN_File, Calc_Params)

    RPA_Input_Path = "IO/A16_Z8_hw16.0_Nmax3_N2max6_N3max9_CMS1+2B"
    RPA_Calc_Params = Any[16, 8, 16.0, 3, true, [0.0, 50.0, 0.5], [0.0, 50.0, 3.0]]
    HF_RPA_Solver(RPA_Input_Path, RPA_Calc_Params)

    RPA_Input_Path = "IO/A16_Z8_hw16.0_Nmax3_N2max6_N3max9_CMS1+2B"
    RPA_Calc_Params = Any[16, 8, 16.0, 3, true, [0.0, 50.0, 0.5], [0.0, 50.0, 3.0]]
    HF_ERPA_Solver(Int_Params, NN_File, NNN_File, RPA_Input_Path, RPA_Calc_Params)

end

HF_ERPA_dev_main()
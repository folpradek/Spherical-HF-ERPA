function HF_NO2B(NN_File::String,NNN_File::String,Params::Vector{Any},Orb::Vector{NOrb})
    # Read calculation parameters ...
    HbarOmega = Params[1]
    N_max = Params[7]
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Load 1-body kinetic operator ...
    T = T1B(N_max,Orb,HbarOmega)

    # 2-body NN interaction & Orbitals - Note that 2-body kinetic operator is computed following read ...
    @time VNN, Orb_NN = V2B_Read(NN_File,Params,Orb)

    # 3-body NNN interaction & Orbitals ...
    @time VNNN, Orb_NNN = V3B_NO2B_Read(NNN_File,Params,Orb)

    # Preallocate some needed arrays ...
    pSPEnergies = zeros(Float64,a_max,1)
    nSPEnergies = zeros(Float64,a_max,1)
    pSPEnergies_Old = zeros(Float64,a_max,1)
    nSPEnergies_Old = zeros(Float64,a_max,1)
    pH = zeros(Float64,a_max,a_max)
    nH = zeros(Float64,a_max,a_max)

    # Initial guess on density operator ...
    pU = diagm(ones(Float64,a_max))
    nU = diagm(ones(Float64,a_max))
    pRho, nRho = HF_Density_Operator(Orb,a_max,pU,nU)

    # Iteration technical parameters - Epsilon precision & Interation_max number may be changed ...
    Iteraction_max = 100
    Iteration = 0
    Epsilon = 1e-6
    Delta = 1.0

    # Hartree-Fock (HF) eqs. iteration ...
    println("\nStarting iteration of HF eqs. with NO2B NN+NNN interaction ...\n")
    @time while (Iteration < Iteraction_max) && (Delta > Epsilon)
        pH .= 0.0
        nH .= 0.0
        pH, nH = HF_Iteration(Params,pH,nH,pRho,nRho,Orb,Orb_NN,Orb_NNN,T,VNN,VNNN)

        # Hamiltonian diagonalization
        pSPEnergies, pU = eigen(Symmetric(pH), sortby=+)
        nSPEnergies, nU = eigen(Symmetric(nH), sortby=+)

        # Single-particle orbitals ordering according to NuHamil orbitals
        pU, pSPEnergies, nU, nSPEnergies = HF_Orbital_Ordering(Orb, a_max, pU, nU, pSPEnergies, nSPEnergies)

        # Generation of density operator
        pRho, nRho = HF_Density_Operator(Orb, a_max, pU, nU)
    
        # Convergence critterion conditon
        Delta = (sum(abs.(pSPEnergies .- pSPEnergies_Old )) + sum(abs.(nSPEnergies .- nSPEnergies_Old ))) / (2*a_max)

        pSPEnergies_Old = deepcopy(pSPEnergies)
        nSPEnergies_Old = deepcopy(nSPEnergies)
        Iteration += 1

        println("Iteration number:   " * string(Iteration) * "   Energy difference:   " * string(round(Delta, sigdigits=8)) * " MeV")
    end
    println("\nHF eqs. with NO2B NN+NNN interaction iteration converged ...")

    # HF energy calculation ...
    E_HF = HF_Energy(Params,pRho,nRho,Orb,Orb_NN,Orb_NNN,T,VNN,VNNN)

    # Kinetic energy effects extraction ...
    T_HF, pT, nT = HF_Kinetic_Energy(Params,pRho,nRho,Orb,T)

    # Calculation summary ...
    HF_Summary(Params,E_HF,T_HF,Epsilon,Iteration)

    # Radial Density & Mean-field potential evaluation & proton radius calculation ...
    HF_Radial_Density(Params,Orb,pU,nU)
    HF_Radial_Potential(Params,Orb,pU,nU,(pH .- pT),(nH .- nT))

    # Single-Particle States export ...
    HF_SPS_Summary(Params,pSPEnergies,nSPEnergies,Orb)

    # Include NO2B NNN interaction to VNN ...
    V2B_Res_Density(Params,Orb,Orb_NN,Orb_NNN,VNN,VNNN,pRho,nRho)

    # Export of single-particle energies & densities ...
    HF_Export(Params,pU,nU,pSPEnergies,nSPEnergies)

    return VNN, Orb_NN

end

function HF_Iteration(Params::Vector{Any},pH::Matrix{Float64},nH::Matrix{Float64},pRho::Matrix{Float64},nRho::Matrix{Float64},Orb::Vector{NOrb},Orb_NN::NNOrb,Orb_NNN::NNNOrb,T::Matrix{Float64},VNN::NNInt,VNNN::Array{Vector{Vector{Float32}},4})
    A = Params[5]
    N_max = Params[7]
    N_2max = Params[8]
    N_3max = Params[9]
    CMS = Params[10]

    a_max = div((N_max + 1)*(N_max + 2),2)

    @inbounds Threads.@threads for a = 1:a_max
        n_a = Orb[a].n
        l_a = Orb[a].l
        j_a = Orb[a].j
        Hat = 1.0/(Float64(j_a) + 1.0)
        @inbounds for d = 1:a
            l_d = Orb[d].l
            j_d = Orb[d].j
            if l_a == l_d && j_a == j_d
                n_d = Orb[d].n
                pSum = 0.0
                nSum = 0.0
                @inbounds for b = 1:a_max
                    n_b = Orb[b].n
                    l_b = Orb[b].l
                    if (2*(n_a + n_b) + l_a + l_b) <= N_2max
                        j_b = Orb[b].j
                        @inbounds for e = 1:a_max
                            n_e = Orb[e].n
                            l_e = Orb[e].l
                            j_e = Orb[e].j
                            if l_b == l_e && j_b == j_e && (2*(n_d + n_e) + l_d + l_e) <= N_2max
                                pRho_be = pRho[b,e]
                                nRho_be = nRho[b,e]

                                @inbounds for J = div(abs(j_d - j_e),2):div(j_d + j_e,2)

                                    if rem(l_a + l_b, 2) == rem(l_d + l_e, 2)

                                        Hat_2b = (Float64(2*J) + 1.0) * Hat
                                        @views pSum += Hat_2b * V2B(a,b,d,e,J,1,VNN.pp,Orb,Orb_NN) * pRho_be
                                        @views pSum += Hat_2b * V2B(a,b,d,e,J,0,VNN.pn,Orb,Orb_NN) * nRho_be
                                        @views nSum += Hat_2b * V2B(a,b,d,e,J,1,VNN.nn,Orb,Orb_NN) * nRho_be
                                        @views nSum += Hat_2b * V2B(b,a,e,d,J,0,VNN.pn,Orb,Orb_NN) * pRho_be

                                    end

                                    # 3-body NNN interaction
                                    @inbounds for c = 1:a_max
                                        n_c = Orb[c].n
                                        l_c = Orb[c].l
                                        if (2*(n_a + n_b + n_c) + l_a + l_b + l_c) <= N_3max
                                            P = rem(l_a + l_b + l_c, 2) + 1
                                            j_c = Orb[c].j
                                            @inbounds for f = 1:a_max
                                                n_f = Orb[f].n
                                                l_f = Orb[f].l
                                                if (2*(n_d + n_e + n_f) + l_d + l_e + l_f) <= N_3max && l_c == l_f && P == (rem(l_d + l_e + l_f,2) + 1)
                                                    j_f = Orb[f].j
                                                    if j_c == j_f
                                                        pRho_cf = pRho[c,f]
                                                        nRho_cf = nRho[c,f]

                                                        ME001 = V3B_NO2B(a,b,c,0,d,e,f,0,J,1,P,VNNN,Orb,Orb_NNN)
                                                        ME101 = V3B_NO2B(a,b,c,1,d,e,f,0,J,1,P,VNNN,Orb,Orb_NNN)
                                                        ME011 = V3B_NO2B(a,b,c,0,d,e,f,1,J,1,P,VNNN,Orb,Orb_NNN)
                                                        ME111 = V3B_NO2B(a,b,c,1,d,e,f,1,J,1,P,VNNN,Orb,Orb_NNN)
                                                        ME113 = V3B_NO2B(a,b,c,1,d,e,f,1,J,3,P,VNNN,Orb,Orb_NNN)

                                                        pSum += Hat * (0.5*ME113*pRho_be*pRho_cf + 0.25 * (ME001 +
                                                                sqrt(1.0/3.0)*ME101 + sqrt(1.0/3.0)*ME011 + 1.0/3.0*ME111 + 2.0/3.0*ME113)*nRho_be*nRho_cf +
                                                                1.0/3.0 * (2*ME111 + ME113)*pRho_be*nRho_cf)

                                                        nSum += Hat * (0.5*ME113*nRho_be*nRho_cf + 0.25 * (ME001 +
                                                                sqrt(1.0/3.0)*ME101 + sqrt(1.0/3.0)*ME011 + 1.0/3.0*ME111 + 2.0/3.0*ME113)*pRho_be*pRho_cf +
                                                                1.0/3.0 * (2*ME111 + ME113)*nRho_be*pRho_cf)

                                                    end
                                                end
                                            end
                                        end
                                    end
                                end

                            end
                        end
                    end
                end

                # 1-body kinetic operator & CMS correction
                if CMS == "CMS1+2B"
                    pH[a,d] = pSum + T[a,d] * (1.0 - 1.0/A)
                    nH[a,d] = nSum + T[a,d] * (1.0 - 1.0/A)
                    pH[d,a] = pH[a,d]
                    nH[d,a] = nH[a,d]
                elseif CMS == "CMS2B"
                    pH[a,d] = pSum
                    nH[a,d] = nSum
                    pH[d,a] = pH[a,d]
                    nH[d,a] = nH[a,d]
                else
                    pH[a,d] = pSum + T[a,d]
                    nH[a,d] = nSum + T[a,d]
                    pH[d,a] = pH[a,d]
                    nH[d,a] = nH[a,d]
                end
            end
        end
    end

    return pH, nH
end


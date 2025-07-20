function HF_NO2B(Params::Parameters,Orb::Vector{NOrb})
    # Read calculation parameters ...
    hw = Params.Int.hw
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1) * (N_max + 2), 2)

    # Load 1-body kinetic operator ...
    T = T1B(N_max,Orb,hw)

    # 2-body NN interaction & Orbitals - Note that 2-body kinetic operator is computed following read ...
    @time VNN, Orb_NN = V2B_Read(Params,Orb)

    # 3-body NNN interaction & Orbitals ...
    @time VNNN, Orb_NNN = V3B_NO2B_Read(Params,Orb)

    # Preallocate some needed arrays ...
    SPEnergies = pnVector(zeros(Float64,a_max), zeros(Float64,a_max))
    SPEnergies_Old = pnVector(zeros(Float64,a_max), zeros(Float64,a_max))

    H = pnMatrix(zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max))

    # Initial guess on density operator ...
    U = pnMatrix(diagm(ones(Float64,a_max)), diagm(ones(Float64,a_max)))
    Rho = HF_Density_Operator(Orb,a_max,U)

    # Iteration technical parameters - Epsilon precision & Interation_max number may be changed ...
    Iteraction_max = 100
    Iteration = 0
    Epsilon = 1e-6
    Delta = 1.0

    # Hartree-Fock (HF) eqs. iteration ...
    println("\nStarting iteration of HF eqs. with NO2B NN+NNN interaction ...\n")
    @time while (Iteration < Iteraction_max) && (Delta > Epsilon)
        H = pnMatrix(zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max))

        H = HF_Iteration(Params,H,Rho,Orb,Orb_NN,Orb_NNN,T,VNN,VNNN)

        # Hamiltonian diagonalization
        pSPEnergies, pU = eigen(Symmetric(H.p), sortby=+)
        nSPEnergies, nU = eigen(Symmetric(H.n), sortby=+)

        SPEnergies = pnVector(pSPEnergies,nSPEnergies)
        U = pnMatrix(pU,nU)

        # Single-particle orbitals ordering according to NuHamil orbitals
        U, SPEnergies = HF_Orbital_Ordering(Orb,a_max,U,SPEnergies)

        # Generation of density operator
        Rho = HF_Density_Operator(Orb,a_max,U)
    
        # Convergence critterion conditon
        Delta = (sum(abs.(SPEnergies.p .- SPEnergies_Old.p )) + sum(abs.(SPEnergies.n .- SPEnergies_Old.n))) / Float64(2*a_max)

        SPEnergies_Old  = pnVector(deepcopy(SPEnergies.p), deepcopy(SPEnergies.n))
        Iteration += 1

        println("Iteration number:   " * string(Iteration) * "   Energy difference:   " * string(round(Delta, sigdigits=8)) * " MeV")
    end

    println("\nHF eqs. with NO2B NN+NNN interaction iteration converged ...")

    # HF energy calculation ...
    E_HF = HF_Energy(Params,Rho,Orb,Orb_NN,Orb_NNN,T,VNN,VNNN)

    # Kinetic energy effects extraction ...
    T_HF, T = HF_Kinetic_Energy(Params,Rho,Orb,T)

    # Calculation summary ...
    HF_Summary(Params,E_HF,T_HF,Epsilon,Iteration)

    # Radial Density & Mean-field potential evaluation & proton radius calculation ...
    HF_Radial_Density(Params,Orb,U)
    HF_Radial_Potential(Params,Orb,U,pnMatrix((H.p .- T.p),(H.n .- T.n)))

    # Single-Particle States export ...
    HF_SPS_Summary(Params,SPEnergies,Orb)

    # Include NO2B NNN interaction to VNN ...
    V2B_Res_Density(Params,Orb,Orb_NN,Orb_NNN,VNN,VNNN,Rho)

    # Export of single-particle energies & densities ...
    HF_Export(Params,U,SPEnergies)

    return VNN, Orb_NN
end

function HF_Iteration(Params::Parameters,H::pnMatrix,Rho::pnMatrix,Orb::Vector{NOrb},Orb_NN::NNOrb,Orb_NNN::NNNOrb,T::Matrix{Float64},VNN::NNInt,VNNN::Array{Vector{Vector{Float32}},4})
    # Read parameters ...
    A = Params.Calc.A
    N_max = Params.Calc.Nmax
    N_2max = Params.Calc.N2max
    N_3max = Params.Calc.N3max
    CMS = Params.Calc.CMS

    a_max = div((N_max + 1)*(N_max + 2),2)

    pH, nH = H.p, H.n
    pRho, nRho = Rho.p, Rho.n

    # Fill HF Hamiltonian ...
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

    return pnMatrix(pH,nH)
end


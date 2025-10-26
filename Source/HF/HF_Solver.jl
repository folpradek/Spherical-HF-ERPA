function HF_Solver(Params::Parameters)
    # Make single-particle orbitals - NuHamil ordering
    Orb = Make_Orbitals(Params.Calc.A,Params.Calc.Z,Params.Int.Nmax)

    # Load 1-body kinetic operator ...
    T = T1B(Params.Calc.Nmax,Orb,Params.Int.hw)

    # 2-body bare NN interaction & Orbitals ...
    @time VNN, Orb_NN = V2B_Read(Params,Orb)

    # 3-body bare NNN interaction & Orbitals ...
    @time VNNN, Orb_NNN = V3B_NO2B_Read(Params,Orb)

    # Solve HF equations ...
    @time SPE, C, Rho, h, Iteration = HF_Solve(Params,Orb,Orb_NN,Orb_NNN,T,VNN,VNNN,1e-6)

    # HF energy calculation ...
    E_HF = HF_Energy(Params,Rho,Orb,Orb_NN,Orb_NNN,T,VNN,VNNN)

    # Kinetic energy expectation values T_HF & effective HF kinetic operator t ...
    T_HF, t = HF_Kinetic_Energy(Params,Rho,Orb,T)

    # Calculation summary ...
    HF_Summary(Params,E_HF,T_HF,1e-6,Iteration)

    # Radial Density & Mean-field potential evaluation & proton radius calculation ...
    HF_Radial_Density(Params,Orb,C)
    HF_Radial_Potential(Params,Orb,C,pnMatrix((h.p .- t.p),(h.n .- t.n)))

    # Single-Particle States export ...
    HF_SPS_Summary(Params,SPE,Orb)

    # Export of single-particle energies & densities ...
    HF_Export(Params,C,SPE)

    # Perform Beyond-HF MBPT calculations & export residual interation ...
    if Params.Calc.BMF == true

        # Include NO2B NNN interaction to NN component ...
        @time VNN = V2B_Res_Density(Params,Orb,Orb_NN,Orb_NNN,VNN,VNNN,Rho)

        # Create residual NN (density-dependent) interaction ... in HF basis ...
        @time VNN_res, Orb_NN_res = V2B_Res(Params,Orb,Orb_NN,VNN,C)

        # Perform HF-MBPT calculations ...
        @time HF_MBPT(Params,Orb,Orb_NN_res,VNN_res)

        # Export residual 2-body interaction ...
        @time V2B_Res_Export(Params,Orb_NN_res,VNN_res)

        # Additional export of single-particle orbitals ...
        if Params.Calc.Format == "HRBin" || Params.Calc.Format == "HR"
            @time Orbitals_Export(Params,Orb)
        end
    end

    return
end

function HF_Solve(Params::Parameters,Orb::Vector{NOrb},Orb_NN::NNOrb,Orb_NNN::NNNOrb,T::Matrix{Float64},VNN::NNInt,VNNN::Array{Vector{Vector{Float32}},4},epsilon::Float64)
    # Read calculation parameters ...
    hw = Params.Int.hw
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1) * (N_max + 2), 2)

    # Set some parameters for iteration ...
    Iteration, Iteraction_max = 0, 100
    delta = 1.0

    # Preallocate arrays ...

        # Vectors for single-particle energies ...
    SPE = pnVector(zeros(Float64,a_max),zeros(Float64,a_max))
    SPE_old = pnVector(zeros(Float64,a_max),zeros(Float64,a_max))

        # HF single-particle Hamiltonian matrix ...
    h = pnMatrix(zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max))

        # Initial guess on single-particle states & 1-body density matrix ... LHO orbitals ...
    C = pnMatrix(diagm(ones(Float64,a_max)), diagm(ones(Float64,a_max)))
    Rho = HF_Density_Operator(Orb,a_max,C)

    # Iteration of spherical HF equations ...
    @time while (Iteration < Iteraction_max) && (delta > epsilon)

        # Perform iteration of HF equations ... fills the HF Hamiltonian ...
        h = HF_Allocate(Params,Rho,Orb,Orb_NN,Orb_NNN,T,VNN,VNNN)

        # Diagonalize the HF Hamiltonians ...
        pSPE, pC = eigen(Symmetric(h.p), sortby=+)
        nSPE, nC = eigen(Symmetric(h.n), sortby=+)

        # Allocate single-particle energies & HF orbitals ...
        SPE = pnVector(pSPE,nSPE)
        C = pnMatrix(pC,nC)

        # Reorder HF orbitals & SPEs ...
        C, SPE = HF_Orbital_Ordering(Orb,a_max,C,SPE)

        # Generate new HF density matrix ...
        Rho = HF_Density_Operator(Orb,a_max,C)
    
        # Check on convergence of HF SPEs ...
        delta = (sum(abs.(SPE.p .- SPE_old.p )) + sum(abs.(SPE.n .- SPE_old.n))) / Float64(2 * a_max)

        SPE_old  = pnVector(deepcopy(SPE.p), deepcopy(SPE.n))
        Iteration += 1

        println("Iteration number:   " * string(Iteration) * "   Energy difference:   " * string(round(delta, sigdigits=8)) * " MeV")
    end

    println("\nIteration of HF eqs. with NO2B NN+NNN interaction finished ...")

    return SPE, C, Rho, h, Iteration
end

function HF_Allocate(Params::Parameters,Rho::pnMatrix,Orb::Vector{NOrb},Orb_NN::NNOrb,Orb_NNN::NNNOrb,T::Matrix{Float64},VNN::NNInt,VNNN::Array{Vector{Vector{Float32}},4})
    # Read parameters ...
    A = Params.Calc.A
    N_max = Params.Calc.Nmax
    N_2max = Params.Calc.N2max
    N_3max = Params.Calc.N3max
    CMS = Params.Calc.CMS
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Allocate new HF Hamiltonian matrices ...
    pH, nH = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)

    # Fill HF Hamiltonian ...
    @inbounds Threads.@threads for a = 1:a_max
        n_a = Orb[a].n
        l_a = Orb[a].l
        j_a = Orb[a].j
        ja_hat = 1.0/(Float64(j_a) + 1.0)
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
                                pRho_be = Rho.p[b,e]
                                nRho_be = Rho.n[b,e]
                                @inbounds for J = div(abs(j_d - j_e),2):div(j_d + j_e,2)

                                    # 2-body NN interaction part ...
                                    if rem(l_a + l_b, 2) == rem(l_d + l_e, 2)
                                        J_ja_hat = Float64(2*J + 1) * ja_hat
                                        @views pSum += J_ja_hat * V2B(a,b,d,e,J,1,VNN.pp,Orb,Orb_NN) * pRho_be
                                        @views pSum += J_ja_hat * V2B(a,b,d,e,J,0,VNN.pn,Orb,Orb_NN) * nRho_be
                                        @views nSum += J_ja_hat * V2B(a,b,d,e,J,1,VNN.nn,Orb,Orb_NN) * nRho_be
                                        @views nSum += J_ja_hat * V2B(b,a,e,d,J,0,VNN.pn,Orb,Orb_NN) * pRho_be

                                    end

                                    # 3-body NNN interaction part ...
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
                                                        pRho_cf = Rho.p[c,f]
                                                        nRho_cf = Rho.n[c,f]

                                                        ME001 = V3B_NO2B(a,b,c,0,d,e,f,0,J,1,P,VNNN,Orb,Orb_NNN)
                                                        ME101 = V3B_NO2B(a,b,c,1,d,e,f,0,J,1,P,VNNN,Orb,Orb_NNN)
                                                        ME011 = V3B_NO2B(a,b,c,0,d,e,f,1,J,1,P,VNNN,Orb,Orb_NNN)
                                                        ME111 = V3B_NO2B(a,b,c,1,d,e,f,1,J,1,P,VNNN,Orb,Orb_NNN)
                                                        ME113 = V3B_NO2B(a,b,c,1,d,e,f,1,J,3,P,VNNN,Orb,Orb_NNN)

                                                        pSum += ja_hat * (0.5*ME113*pRho_be*pRho_cf + 0.25 * (ME001 +
                                                                sqrt(1.0/3.0)*ME101 + sqrt(1.0/3.0)*ME011 + 1.0/3.0*ME111 + 2.0/3.0*ME113)*nRho_be*nRho_cf +
                                                                1.0/3.0 * (2*ME111 + ME113)*pRho_be*nRho_cf)

                                                        nSum += ja_hat * (0.5*ME113*nRho_be*nRho_cf + 0.25 * (ME001 +
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

                # 1-body kinetic energy & inclusion of Center-of-Mass motion (CM) correction ...
                if CMS == "CMS1+2B"
                    pH[a,d] = pSum + T[a,d] * (1.0 - 1.0 / A)
                    nH[a,d] = nSum + T[a,d] * (1.0 - 1.0 / A)
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
struct HFB_Broyden_Key
    a::Int64
    b::Int64
end

struct HFB_Broyden
    Map::Matrix{Int64}
    Key::Vector{HFB_Broyden_Key}
    M::Int64
    m::Int64
end

mutable struct HFB_Broyden_Vector
    x_Rho::pnVector
    y_Rho::pnVector
    X_Rho::pnMatrix
    r_Rho::pnVector
    s_Rho::pnVector
    R_Rho::pnMatrix
    x_Kappa::pnVector
    y_Kappa::pnVector
    X_Kappa::pnMatrix
    r_Kappa::pnVector
    s_Kappa::pnVector
    R_Kappa::pnMatrix
end

function HFB_Broyden_Initialize_Mapping(Params::Parameters,m::Int64,Orb::Vector{NOrb})
    # Read calculation parameters ...
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1) * (N_max + 2), 2)

    # Initialize the counter of linearized indices & Broyden mapping ...
    M_count, Broyden_Map = 0, zeros(Int64,a_max,a_max)

    # Enumerate M_count & allocate Broyden mapping ...
    @inbounds for a in 1:a_max
        l_a, j_a = Orb[a].l, Orb[a].j
        @inbounds for b in 1:a_max
            l_b, j_b = Orb[b].l, Orb[b].j
            if (l_a == l_b) && (j_a == j_b)
                M_count += 1
                Broyden_Map[a,b] = M_count 
                #Broyden_Map[b,a] = M_count
            end
        end
    end

    # Initialize the total number of linearized indices & Broyden key ...
    M, Broyden_Key = 0, Vector{HFB_Broyden_Key}(undef,M_count)

    # Enumerate M & allocate Broyden key ...
    @inbounds for a in 1:a_max
        l_a, j_a = Orb[a].l, Orb[a].j
        @inbounds for b in 1:a_max
            l_b, j_b = Orb[b].l, Orb[b].j
            if (l_a == l_b) && (j_a == j_b)
                M += 1
                Broyden_Key[M] = HFB_Broyden_Key(a,b) 
            end
        end
    end

    return HFB_Broyden(Broyden_Map,Broyden_Key,M,m), M
end

function HFB_Broyden_Initialize(Params::Parameters,Rho::pnMatrix,Kappa::pnMatrix,Orb::Vector{NOrb})
    # Read calculation parameters ...
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1) * (N_max + 2), 2)

    # Number of history vectors ...
    m = 8

    # Prepare Broyden ... Map, Key & calculate the number of linearized indices M ...
    Broyden, M = HFB_Broyden_Initialize_Mapping(Params,m,Orb)

    # Initialize arrays x & r ...
    x_Rho, x_Kappa = HFB_Broyden_Linearize(Rho,Kappa,Broyden;Initialize = true)
    y_Rho, y_Kappa = pnVector(zeros(Float64,M),zeros(Float64,M)), pnVector(zeros(Float64,M),zeros(Float64,M))
    r_Rho, r_Kappa = pnVector(zeros(Float64,M),zeros(Float64,M)), pnVector(zeros(Float64,M),zeros(Float64,M))
    s_Rho, s_Kappa = pnVector(zeros(Float64,M),zeros(Float64,M)), pnVector(zeros(Float64,M),zeros(Float64,M))

    # Initialize history vectors X & R ...
    X_Rho, X_Kappa = pnMatrix(zeros(Float64,M,m),zeros(Float64,M,m)), pnMatrix(zeros(Float64,M,m),zeros(Float64,M,m))
    R_Rho, R_Kappa = pnMatrix(zeros(Float64,M,m),zeros(Float64,M,m)), pnMatrix(zeros(Float64,M,m),zeros(Float64,M,m))

    # Initialize Broyden vector ...
    Broyden_Vector = HFB_Broyden_Vector(x_Rho,y_Rho,X_Rho,r_Rho,s_Rho,R_Rho,x_Kappa,y_Kappa,X_Kappa,r_Kappa,s_Kappa,R_Kappa)

    return Broyden, Broyden_Vector
end

function HFB_Broyden_Linearize(Rho::pnMatrix,Kappa::pnMatrix,Broyden::HFB_Broyden;Initialize::Bool=false)
    # Initialize linearized Rho & Kappa ...
    pRho, pKappa = zeros(Float64,Broyden.M), zeros(Float64,Broyden.M)
    nRho, nKappa = zeros(Float64,Broyden.M), zeros(Float64,Broyden.M)

    # Allocate linearized Rho & Kappa ...
        # Paralelize using threads (???) ...
    @inbounds for i in 1:Broyden.M
        a, b = Broyden.Key[i].a, Broyden.Key[i].b
        pRho[i], pKappa[i] = Rho.p[a,b], Kappa.p[a,b]
        nRho[i], nKappa[i] = Rho.n[a,b], Kappa.n[a,b]
    end

    if Initialize == false
        return pRho, pKappa, nRho, nKappa
    elseif Initialize == true
        return pnVector(pRho,nRho), pnVector(pKappa,nKappa)
    end
end

function HFB_Broyden_Reconstruct(Params::Parameters,Rho::pnVector,Kappa::pnVector,Broyden::HFB_Broyden,Orb::Vector{NOrb})
    # Read calculation parameters ...
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1) * (N_max + 2), 2)

    # Initialize density matrices Rho & Kappa ...
    pRho, pKappa = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)
    nRho, nKappa = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)

    # Allocate density matrices Rho & Kappa ...
        # Paralelize using threads (???) ...
    @inbounds for a in 1:a_max
        @inbounds for b in 1:a_max
            i = Broyden.Map[a,b]
            if i != 0
                pRho[a,b], pKappa[a,b] = Rho.p[i], Kappa.p[i]
                nRho[a,b], nKappa[a,b] = Rho.n[i], Kappa.n[i]
            end
        end
    end

    # Clean numerical noise in Rho & Kappa ...
    @inbounds for a in 1:a_max
        @inbounds for b in 1:a_max
            if abs(pRho[a,b]) < 1e-14
                pRho[a,b] = 0.0
            end
            if abs(pKappa[a,b]) < 1e-14
                pKappa[a,b] = 0.0
            end
           if abs(nRho[a,b]) < 1e-14
                nRho[a,b] = 0.0
            end
            if abs(nKappa[a,b]) < 1e-14
                nKappa[a,b] = 0.0
            end
        end
    end

    return pnMatrix(pRho,nRho), pnMatrix(pKappa,nKappa)
end

function HFB_Broyden_History_Update(Broyden::HFB_Broyden,BroyVec::HFB_Broyden_Vector)
    # Define the residual overlap matrix N = R' * R ...
    pRho_N, pKappa_N = BroyVec.R_Rho.p' * BroyVec.R_Rho.p, BroyVec.R_Kappa.p' * BroyVec.R_Kappa.p
    nRho_N, nKappa_N = BroyVec.R_Rho.n' * BroyVec.R_Rho.n, BroyVec.R_Kappa.n' * BroyVec.R_Kappa.n

    # Calculate the regularizing factors epsilon ...
    pRho_Eps = max(1e-4 * opnorm(pRho_N,2), 1e-7) .* diagm(ones(Float64,Broyden.m))
    nRho_Eps = max(1e-4 * opnorm(nRho_N,2), 1e-7) .* diagm(ones(Float64,Broyden.m))
    pKappa_Eps = max(1e-4 * opnorm(pKappa_N,2), 1e-7) .* diagm(ones(Float64,Broyden.m))
    nKappa_Eps = max(1e-4 * opnorm(nKappa_N,2), 1e-7) .* diagm(ones(Float64,Broyden.m))

    # Calculate & allocate beta ...
    pRho_beta = (pRho_N .+ pRho_Eps) \ (BroyVec.R_Rho.p' * BroyVec.r_Rho.p)
    nRho_beta = (nRho_N .+ nRho_Eps) \ (BroyVec.R_Rho.n' * BroyVec.r_Rho.n)
    pKappa_beta = (pKappa_N .+ pKappa_Eps) \ (BroyVec.R_Kappa.p' * BroyVec.r_Kappa.p)
    nKappa_beta = (nKappa_N .+ nKappa_Eps) \ (BroyVec.R_Kappa.n' * BroyVec.r_Kappa.n)

    return pRho_beta, pKappa_beta, nRho_beta, nKappa_beta
end

function HFB_Broyden_Update(Params::Parameters,Iteration::Int64,Rho::pnMatrix,Kappa::pnMatrix,Broyden::HFB_Broyden,BroyVec::HFB_Broyden_Vector,Orb::Vector{NOrb})
    # Initialize quenching flags ...
    Q_pRho, Q_pKappa = false, false
    Q_nRho, Q_nKappa = false, false

    # History vector index m ...
    m = rem(Iteration,Broyden.m) + 1

    # Quenching parameters alpha ...
    alpha_pRho, alpha_pKappa = 0.9, 0.9
    alpha_nRho, alpha_nKappa = 0.9, 0.9

    # Linearize Rho & Kappa ...
    pRho, pKappa, nRho, nKappa = HFB_Broyden_Linearize(Rho,Kappa,Broyden)

        # Remove this later ... allocate as slightly quenched previous iteration ...
    pRho_new, pKappa_new, nRho_new, nKappa_new = zeros(Float64,Broyden.M), zeros(Float64,Broyden.M), zeros(Float64,Broyden.M), zeros(Float64,Broyden.M)

    # Update & allocate vectors r & s ...
        # Case of s = r^(n-1) ...
    BroyVec.s_Rho.p .= BroyVec.r_Rho.p
    BroyVec.s_Rho.n .= BroyVec.r_Rho.n
    BroyVec.s_Kappa.p .= BroyVec.r_Kappa.p
    BroyVec.s_Kappa.n .= BroyVec.r_Kappa.n
        # Case of r = r^(n) ...
    BroyVec.r_Rho.p .= pRho .- BroyVec.x_Rho.p
    BroyVec.r_Rho.n .= nRho .- BroyVec.x_Rho.n
    BroyVec.r_Kappa.p .= pKappa .- BroyVec.x_Kappa.p
    BroyVec.r_Kappa.n .= nKappa .- BroyVec.x_Kappa.n

    # Update history vectors X & R ...
        # Case of X = X^(n) ...
    @views BroyVec.X_Rho.p[:,m] = BroyVec.x_Rho.p[:] .- BroyVec.y_Rho.p[:]
    @views BroyVec.X_Rho.n[:,m] = BroyVec.x_Rho.n[:] .- BroyVec.y_Rho.n[:]
    @views BroyVec.X_Kappa.p[:,m] = BroyVec.x_Kappa.p[:] .- BroyVec.y_Kappa.p[:]
    @views BroyVec.X_Kappa.n[:,m] = BroyVec.x_Kappa.n[:] .- BroyVec.y_Kappa.n[:]
        # Case of R = R^(n) ...
    @views BroyVec.R_Rho.p[:,m] = BroyVec.r_Rho.p[:] .- BroyVec.s_Rho.p[:]
    @views BroyVec.R_Rho.n[:,m] = BroyVec.r_Rho.n[:] .- BroyVec.s_Rho.n[:]
    @views BroyVec.R_Kappa.p[:,m] = BroyVec.r_Kappa.p[:] .- BroyVec.s_Kappa.p[:]
    @views BroyVec.R_Kappa.n[:,m] = BroyVec.r_Kappa.n[:] .- BroyVec.s_Kappa.n[:]

    # Allocate finite difference D for previous residual r ...
    D_pRho, D_pKappa = norm(BroyVec.r_Rho.p), norm(BroyVec.r_Kappa.p)
    D_nRho, D_nKappa = norm(BroyVec.r_Rho.n), norm(BroyVec.r_Kappa.n)

    # Initialize the Broyden's history fit matrix beta ...
    pRho_beta, pKappa_beta = zeros(Float64,Broyden.m), zeros(Float64,Broyden.m)
    nRho_beta, nKappa_beta = zeros(Float64,Broyden.m), zeros(Float64,Broyden.m)
        # If there is enough history, calculate the history update matrix beta ...
    if Iteration > (2*Broyden.m + 1)
        pRho_beta, pKappa_beta, nRho_beta, nKappa_beta = HFB_Broyden_History_Update(Broyden,BroyVec)
    end

    # Perform several quenching iterations for densities ...
    @inbounds for Quench in 1:3

        # Quench pRho ...
        if Q_pRho == false
            alpha_r = alpha_pRho .* BroyVec.r_Rho.p
            X_beta = BroyVec.X_Rho.p * pRho_beta
            dX_beta_max = 0.15 * norm(alpha_r)
            if dX_beta_max > 1e-8 && norm(X_beta) > dX_beta_max
                X_beta .= X_beta .* (dX_beta_max / (norm(X_beta) + eps()))
            end

            #pRho_trial = BroyVec.x_Rho.p .+ alpha_pRho .* BroyVec.r_Rho.p .- BroyVec.X_Rho.p * pRho_beta
            pRho_trial = BroyVec.x_Rho.p .+ alpha_r .- X_beta
            D_pRho_trial = norm(pRho .- pRho_trial)

            if ((D_pRho_trial / D_pRho) < 0.9) || (Quench == 3)
                pRho_new .= pRho_trial
                Q_pRho = true
            else
                alpha_pRho = 0.5 * alpha_pRho
            end
        end

        # Quench pKappa ...
        if Q_pKappa == false
            alpha_r = alpha_pKappa .* BroyVec.r_Kappa.p
            X_beta = BroyVec.X_Kappa.p * pKappa_beta
            dX_beta_max = 0.15 * norm(alpha_r)
            if dX_beta_max > 1e-8 && norm(X_beta) > dX_beta_max
                X_beta .= X_beta .* (dX_beta_max / (norm(X_beta) + eps()))
            end

            #pKappa_trial = BroyVec.x_Kappa.p .+ alpha_pKappa .* BroyVec.r_Kappa.p .- BroyVec.X_Kappa.p * pKappa_beta
            pKappa_trial = BroyVec.x_Kappa.p .+ alpha_r .- X_beta

            D_pKappa_trial = norm(pKappa .- pKappa_trial)

            if ((D_pKappa_trial / D_pKappa) < 0.9) || (Quench == 3)
                pKappa_new .= pKappa_trial
                Q_pKappa = true
            else 
                alpha_pKappa = 0.5 * alpha_pKappa
            end
        end

        # Quench nRho ...
        if Q_nRho == false
            alpha_r = alpha_nRho .* BroyVec.r_Rho.n
            X_beta = BroyVec.X_Rho.n * nRho_beta
            dX_beta_max = 0.15 * norm(alpha_r)
            if dX_beta_max > 1e-8 && norm(X_beta) > dX_beta_max
                X_beta .= X_beta .* (dX_beta_max / (norm(X_beta) + eps()))
            end

            #nRho_trial = BroyVec.x_Rho.n .+ alpha_nRho .* BroyVec.r_Rho.n .- BroyVec.X_Rho.n * nRho_beta
            nRho_trial = BroyVec.x_Rho.n .+ alpha_r .- X_beta

            D_nRho_trial = norm(nRho .- nRho_trial)

            if ((D_nRho_trial / D_nRho) < 0.9) || (Quench == 3)
                nRho_new .= nRho_trial
                Q_nRho = true
            else
                alpha_nRho = 0.5 * alpha_nRho
            end
        end

        # Quench nKappa ...
        if Q_nKappa == false
            alpha_r = alpha_nKappa .* BroyVec.r_Kappa.n
            X_beta = BroyVec.X_Kappa.n * nKappa_beta
            dX_beta_max = 0.15 * norm(alpha_r)
            if dX_beta_max > 1e-8 && norm(X_beta) > dX_beta_max
                X_beta .= X_beta .* (dX_beta_max / (norm(X_beta) + eps()))
            end

            #nKappa_trial = BroyVec.x_Kappa.n .+ alpha_nKappa .* BroyVec.r_Kappa.n .- BroyVec.X_Kappa.n * nKappa_beta
            nKappa_trial = BroyVec.x_Kappa.n .+ alpha_r .- X_beta

            D_nKappa_trial = norm(nKappa .- nKappa_trial)

            if ((D_nKappa_trial / D_nKappa) < 0.9) || (Quench == 3)
                nKappa_new .= nKappa_trial
                Q_nKappa = true
            else
                alpha_nKappa = 0.5 * alpha_nKappa
            end
        end

        # Break loop if update is done ...
        if Q_pRho == true && Q_pKappa == true && Q_nRho == true && Q_nKappa == true 
            break
        end

    end

    # Update previous step vector y ...
    BroyVec.y_Rho.p .= BroyVec.x_Rho.p
    BroyVec.y_Rho.n .= BroyVec.x_Rho.n
    BroyVec.y_Kappa.p .= BroyVec.x_Kappa.p
    BroyVec.y_Kappa.n .= BroyVec.x_Kappa.n

    # Update current step vector x ...
    BroyVec.x_Rho.p .= pRho_new
    BroyVec.x_Rho.n .= nRho_new
    BroyVec.x_Kappa.p .= pKappa_new
    BroyVec.x_Kappa.n .= nKappa_new

    # Update densities Rho & Kappa ...
    Rho, Kappa = HFB_Broyden_Reconstruct(Params,pnVector(pRho_new,nRho_new),pnVector(pKappa_new,nKappa_new),Broyden,Orb)

    return Rho, Kappa, BroyVec
end
function Make_Orbitals(A::Int64,Z::Int64,N_max::Int64)
    println("\nPreparing oscillator basis single-particle orbitals ...")
    a_max = div((N_max + 1)*(N_max + 2),2)
    Orb = Vector{NOrb}(undef, a_max)
    Orb_arr = zeros(Int64, a_max, 6)
    a = 0
    @inbounds for N = 0:N_max
        @inbounds for l in 0:N
            n = Int64(floor((N - l)/2))
            if (2*n + l == N)
                @inbounds for j in Int64(abs(2*l - 1)):2:Int64(2*l + 1)
                    a += 1
                    Orb_arr[a,1] = a
                    Orb_arr[a,2] = n
                    Orb_arr[a,3] = l
                    Orb_arr[a,4] = j
                end
            end
        end
    end
    if A != 0 && Z != 0
        N_proton = 0
        N_neutron = 0
        Neutron_Filling = true
        Proton_Filling = true
        a_max =  Int64((N_max + 1)*(N_max + 2)/2)
        while Neutron_Filling
            @inbounds for N = 0:N_max
                @inbounds for n in 0:div(N, 2)
                    l = N - 2*n
                    @inbounds for j = (2*l + 1):-2:abs(2*l-1)
                        if N_neutron + (j + 1) <= (A - Z)
                            N_neutron += j + 1
                            @inbounds for a = 1:a_max
                                if (Orb_arr[a,2] == n) && (Orb_arr[a,3] == l) && (Orb_arr[a,4] == j)
                                    Orb_arr[a,6] = 1
                                end
                            end
                        else
                            Neutron_Filling = false
                        end
                    end
                end
            end
        end
        while Proton_Filling
            @inbounds for N in 0:N_max
                @inbounds for n = 0:div(N, 2)
                    l = N - 2*n
                    @inbounds for j = (2*l + 1):-2:abs(2*l-1)
                        if N_proton + (j + 1) <= Z
                            N_proton += j + 1
                            @inbounds for a = 1:a_max
                                if (Orb_arr[a,2] == n) && (Orb_arr[a,3] == l) && (Orb_arr[a,4] == j)
                                    Orb_arr[a,5] = 1
                                end
                            end
                        else
                            Proton_Filling = false
                        end
                    end
                end
            end
        end
    end
    @inbounds for b = 1:a_max
        Orb[b] = NOrb(Orb_arr[b,1],Orb_arr[b,2],Orb_arr[b,3],Orb_arr[b,4],Orb_arr[b,5],Orb_arr[b,6])
    end
    println("\nOrbitals initialized ...")
    return Orb
end

function Make_ParticleHole_Orbitals(N_max::Int64,Input_File::String,Orb::Vector{NOrb})
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Import single-particle energies ...
    pE_Import = Input_File * "/Bin/pE.bin"
    pSPE = Vector{Float64}(undef,a_max)
    open(pE_Import, "r") do Read_File
        @inbounds for a in 1:a_max
            ME = read(Read_File, Float64)
            pSPE[a] = ME
        end
    end

    nE_Import = Input_File * "/Bin/nE.bin"
    nSPE = Vector{Float64}(undef,a_max)
    open(nE_Import, "r") do Read_File
        @inbounds for a in 1:a_max
            ME = read(Read_File, Float64)
            nSPE[a] = ME
        end
    end

    N_pParticle = 0
    N_nParticle = 0
    N_pHole = 0
    N_nHole = 0

    @inbounds for a in 1:a_max
        if Orb[a].pO == 1
            N_pHole += 1
        end
        if Orb[a].pO == 0
            N_pParticle += 1
        end
        if Orb[a].nO == 1
            N_nHole += 1
        end
        if Orb[a].nO == 0
            N_nParticle += 1
        end
    end

    pParticle_vec = Vector{SPState}(undef, N_pParticle)
    nParticle_vec = Vector{SPState}(undef, N_nParticle)
    pHole_vec = Vector{SPState}(undef, N_pHole)
    nHole_vec = Vector{SPState}(undef, N_nHole)

    C_pParticle = 1
    C_nParticle = 1
    C_pHole = 1
    C_nHole = 1

    @inbounds for a in 1:a_max
        if Orb[a].pO == 0
            pParticle_vec[C_pParticle] = SPState(Orb[a].a,Orb[a].l,Orb[a].j,pSPE[a])
            C_pParticle += 1
        end

        if Orb[a].pO == 1
            pHole_vec[C_pHole] = SPState(Orb[a].a,Orb[a].l,Orb[a].j,pSPE[a])
            C_pHole += 1
        end

        if Orb[a].nO == 0
            nParticle_vec[C_nParticle] = SPState(Orb[a].a,Orb[a].l,Orb[a].j,nSPE[a])
            C_nParticle += 1
        end

        if Orb[a].nO == 1
            nHole_vec[C_nHole] = SPState(Orb[a].a,Orb[a].l,Orb[a].j,nSPE[a])
            C_nHole += 1
        end

    end

    N_particle = pnInteger(N_pParticle,N_nParticle)
    N_hole = pnInteger(N_pHole,N_nHole)

    Particle = pnSVector(pParticle_vec,nParticle_vec)
    Hole = pnSVector(pHole_vec,nHole_vec)

    return N_particle, Particle, N_hole, Hole

end

function Make_Phonon_Orbitals(N_Particle::pnInteger,Particle::pnSVector,N_Hole::pnInteger,Hole::pnSVector)
    N_Phonon = 0

    for p in 1:N_Particle.p
        for h in 1:N_Hole.p
            for J in div(abs(Particle.p[p].j - Hole.p[h].j),2):div(Particle.p[p].j + Hole.p[h].j,2)
                N_Phonon += 1
            end
        end
    end
    
    for p in 1:N_Particle.n
        for h in 1:N_Hole.n
            for J in div(abs(Particle.n[p].j - Hole.n[h].j),2):div(Particle.n[p].j + Hole.n[h].j,2)
                N_Phonon += 1
            end
        end
    end

    Phonon = Vector{PhState}(undef,N_Phonon)

    N_ph = 0
    for p in 1:N_Particle.p
        for h in 1:N_Hole.p
            for J in div(abs(Particle.p[p].j - Hole.p[h].j),2):div(Particle.p[p].j + Hole.p[h].j,2)
                N_ph += 1
                Phonon[N_ph] = PhState(p,h,-1,J,rem(Particle.p[p].l + Hole.p[h].l,2)+1)
            end
        end
    end

    for p in 1:N_Particle.n
        for h in 1:N_Hole.n
            for J in div(abs(Particle.n[p].j - Hole.n[h].j),2):div(Particle.n[p].j + Hole.n[h].j,2)
                N_ph += 1
                Phonon[N_ph] = PhState(p,h,1,J,rem(Particle.n[p].l + Hole.n[h].l,2)+1)
            end
        end
    end

    return N_Phonon, Phonon
end

function Phonon_Ind(ph::Int64,Phonon::Vector{PhState},Particle::pnSVector,Hole::pnSVector)
    p = Phonon[ph].p
    h = Phonon[ph].h
    t_ph = Phonon[ph].tz
    J = Phonon[ph].J
    P = Phonon[ph].P
    if t_ph == -1
        a_p = Particle.p[p].a
        l_p = Particle.p[p].l
        j_p = Particle.p[p].j
        E_p = Particle.p[p].E
        a_h = Hole.p[h].a
        l_h = Hole.p[h].l
        j_h = Hole.p[h].j
        E_h = Hole.p[h].E
    elseif t_ph == 1
        a_p = Particle.n[p].a
        l_p = Particle.n[p].l
        j_p = Particle.n[p].j
        E_p = Particle.n[p].E
        a_h = Hole.n[h].a
        l_h = Hole.n[h].l
        j_h = Hole.n[h].j
        E_h = Hole.n[h].E
    end
    return p, h, t_ph, J, P, a_p, l_p, j_p, E_p, a_h, l_h ,j_h, E_h
end

function Orbitals_Export(Params::Vector{Any},Orb::Vector{NOrb})
    N_max = Params[7]
    Output_File = Params[13]
    a_max = div((N_max+1)*(N_max+2),2)

    Orbitals_LHO_Export = "IO/" * Output_File * "/Orbitals_LHO.dat"
    Orbitals_HF_Export = "IO/" * Output_File * "/Orbitals_HF.dat"

    # Import HF s.p. energies ...
    pE_Import = "IO/" * Output_File * "/Bin/pE.bin"
    pSPE = Vector{Float64}(undef,a_max)
    open(pE_Import, "r") do Read_File
        @inbounds for a in 1:a_max
            ME = read(Read_File, Float64)
            pSPE[a] = ME
        end
    end

    nE_Import = "IO/" * Output_File * "/Bin/nE.bin"
    nSPE = Vector{Float64}(undef,a_max)
    open(nE_Import, "r") do Read_File
        @inbounds for a in 1:a_max
            ME = read(Read_File, Float64)
            nSPE[a] = ME
        end
    end

    # Export Single-Particle State Orbitals in the LHO basis ...
    println("\nExporting LHO basis s.p. orbitals in human-readable format ...")

    open(Orbitals_LHO_Export, "w") do Export_File
        println(Export_File, "a\tn\tl\t2j\tpOc\tnOc")
        @inbounds for a in 1:a_max
            Row = string(Orb[a].a) * "\t" * string(Orb[a].n) *
                  "\t" * string(Orb[a].l) * "\t" * string(Orb[a].j) *
                  "\t " * string(Orb[a].pO) * "\t " * string(Orb[a].nO)
            println(Export_File, Row)
        end
    
    end
    
    println("\nExported LHO basis s.p. orbitals in human-readable format ...")

    # Export Single-Particle State Orbitals in the HF basis ...
    println("\nExporting HF basis s.p. orbitals in human-readable format ...")

    open(Orbitals_HF_Export, "w") do Export_File
        println(Export_File, "a\tl\t2j\tpOc\tnOc\tE_p\t\t\tE_n")
        @inbounds for a in 1:a_max
            Row = string(Orb[a].a) * "\t" * string(Orb[a].l) * "\t" * string(Orb[a].j) *
                  "\t " * string(Orb[a].pO) * "\t " * string(Orb[a].nO) *
                  "\t" * string(pSPE[a]) * "\t" * string(nSPE[a])
            println(Export_File, Row)
        end
    end
    
    println("\nExported HF basis s.p. orbitals in human-readable format ...")

    return
end

function Make_ParticleHole_List(J_max::Int64,N_Particle::pnInteger,Particle::pnSVector,N_Hole::pnInteger,Hole::pnSVector)
    N_Phonon = zeros(Int64,J_max+1,2)

    pParticleHole = zeros(Int64,J_max+1,2,N_Particle.p,N_Hole.p)
    nParticleHole = zeros(Int64,J_max+1,2,N_Particle.n,N_Hole.n)

    @inbounds for p in 1:N_Particle.p
        @inbounds for h in 1:N_Hole.p
            P = rem(Particle.p[p].l + Hole.p[h].l,2) + 1
            @inbounds for J in div(abs(Particle.p[p].j - Hole.p[h].j),2):div(Particle.p[p].j + Hole.p[h].j,2)
                N_Phonon[J+1,P] += 1
                pParticleHole[J+1,P,p,h] = N_Phonon[J+1,P]
            end
        end
    end
    
    @inbounds for p in 1:N_Particle.n
        @inbounds for h in 1:N_Hole.n
            P = rem(Particle.n[p].l + Hole.n[h].l,2) + 1
            @inbounds for J in div(abs(Particle.n[p].j - Hole.n[h].j),2):div(Particle.n[p].j + Hole.n[h].j,2)
                N_Phonon[J+1,P] += 1
                nParticleHole[J+1,P,p,h] = N_Phonon[J+1,P]
            end
        end
    end

    ParticleHole = pnArray(pParticleHole,nParticleHole)

    return ParticleHole
end
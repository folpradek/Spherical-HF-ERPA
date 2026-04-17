function HF_RRPA_import_binary()

    return
end

function HF_RRPA_transition_density(Params::Parameters;J::Int64;P::Int64;Phonon_list::Vector{Int64},Weightened::Bool=false)
    # Read parameters ...

    # Define basic constants ...
    N_grid = 2^13 + 1 # 8192 + 1 grid points ...

    # Make s.p. orbitals ...
    Orb = orbitals_make(Params)

    # Prepare Particle & Hole orbitals ...
    N_Particle, Particle, N_Hole, Hole = orbitals_ph_make(Params,Orb)

    # Prepare 1p-1h phonon states ...
    N_Phonon, Phonon = orbitals_one_phonon_make(N_Particle,Particle,N_Hole,Hole)

    # Count & pre-index all phonon states in JP subspaces ...
    N_nu, Orb_Phonon = HF_RRPA_phonon_count(Params,N_Phonon,Phonon)

    # Load RRPA solutions ...

    # Calculate radial representation of single-particle states ...

    # If set to true, calculate the radial representation of 1-body transition operator ...
    if Weightened == true
        println("cool")
    end

    # For given phonon nu evaluate the radial transition density ...
    Rho_nu = Vector{Float64}(undef,N_grid)

    @inbounds for i in 1:N_grid
        @inbounds for ph in 1:N_ph

        end
    end

    return
end
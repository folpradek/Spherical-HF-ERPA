function O1b_export(Params::Parameters,Orb::Vector{Orb1B},A::O1B,A_Export_Path::String)
    # Read parameters ...
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Export 1-body operator A under given path ...
    open(A_Export_Path, "w") do Export_File
        # Proton entries ...
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                if Orb[a].l == Orb[b].l &&Orb[a].j == Orb[b].j
                    ME = @views A.p[a,b]
                    write(Export_File, Float64(ME))
                end
            end
        end

        # Neutron entries ...
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                if Orb[a].l == Orb[b].l &&Orb[a].j == Orb[b].j
                    ME = @views A.n[a,b]
                    write(Export_File, Float64(ME))
                end
            end
        end
    end
    
    return
end

function O1b_import(Params::Parameters,Orb::Vector{Orb1B},Import_Path::String)
    # Read parameters ...
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Read the 1-body operator A matrix from given file ...
    pA, nA = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)
    open(Import_Path, "r") do Read_File
        # Proton entries ...
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                if Orb[a].l == Orb[b].l && Orb[a].j == Orb[b].j
                    ME = read(Read_File,Float64)
                    pA[a,b] = ME
                end
            end
        end

        # Neutron entries ...
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                if Orb[a].l == Orb[b].l && Orb[a].j == Orb[b].j
                    ME = read(Read_File,Float64)
                    nA[a,b] = ME
                end
            end
        end
    end
    
    return O1B(pA,nA)
end

function qpO1B_export(Params::Parameters,Orb::Vector{Orb1B},A::qpO1B,A_Export_Path::String)
    # Read parameters ...
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Export 1-body operator A under given path ...
    open(A_Export_Path, "w") do Export_File
            # Export of qp11 part ...
        # Proton entries ...
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                if Orb[a].l == Orb[b].l && Orb[a].j == Orb[b].j
                    ME = A.qp11.p[a,b]
                    write(Export_File, Float64(ME))
                end
            end
        end

        # Neutron entries ...
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                if Orb[a].l == Orb[b].l &&Orb[a].j == Orb[b].j
                    ME = A.qp11.n[a,b]
                    write(Export_File, Float64(ME))
                end
            end
        end

            # Export of qp20 part ...
        # Proton entries ...
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                if Orb[a].l == Orb[b].l &&Orb[a].j == Orb[b].j
                    ME = A.qp20.p[a,b]
                    write(Export_File, Float64(ME))
                end
            end
        end

        # Neutron entries ...
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                if Orb[a].l == Orb[b].l &&Orb[a].j == Orb[b].j
                    ME = A.qp20.n[a,b]
                    write(Export_File, Float64(ME))
                end
            end
        end

    end
    
    return
end

function qpO1B_import(Params::Parameters,Orb::Vector{Orb1B},Import_Path::String)
    # Read parameters ...
    N_max = Params.Calc.Nmax
    a_max = div((N_max + 1)*(N_max + 2),2)

    # Read the 1-body operator A matrix from given file ...
    pA11, pA20 = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)
    nA11, nA20 = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)
    open(Import_Path, "r") do Read_File
            # Read qp11 part ...
        # Proton entries ...
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                if Orb[a].l == Orb[b].l && Orb[a].j == Orb[b].j
                    ME = read(Read_File,Float64)
                    pA11[a,b] = ME
                end
            end
        end

        # Neutron entries ...
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                if Orb[a].l == Orb[b].l && Orb[a].j == Orb[b].j
                    ME = read(Read_File,Float64)
                    nA11[a,b] = ME
                end
            end
        end

            # Read qp20 part ...
        # Proton entries ...
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                if Orb[a].l == Orb[b].l && Orb[a].j == Orb[b].j
                    ME = read(Read_File,Float64)
                    pA20[a,b] = ME
                end
            end
        end

        # Neutron entries ...
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                if Orb[a].l == Orb[b].l && Orb[a].j == Orb[b].j
                    ME = read(Read_File,Float64)
                    nA20[a,b] = ME
                end
            end
        end
    end

    return qpO1B(O1B(pA11,nA11),O1B(pA20,nA20))
end
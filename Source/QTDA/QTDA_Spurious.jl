function QTDA_spurious_PN_initialize(Orb::Vector{Orb1B},Orb_2qp::qpOrb2B,qpN::qpO1B)
    # Allocate the 0+ spurious vector ...
    J, P, N_qp = 0, 1, Orb_2qp.N[1,1]
    Spur_PN = Vector{Float64}(undef,N_qp)
    @inbounds for qp in 1:N_qp
        a, b, T_ab = Orb_2qp.i[P,J+1][qp].a, Orb_2qp.i[P,J+1][qp].b, Orb_2qp.i[P,J+1][qp].T
        j_a = Orb[a].j
        j_a_hat = sqrt(Float64(j_a + 1))
        if T_ab == -1
            ME = j_a_hat * qpN.qp20.p[a,b]
            Spur_PN[qp] = ME
        elseif T_ab == 1
            ME = j_a_hat * qpN.qp20.n[a,b]
            Spur_PN[qp] = ME
        end
    end

    # Renormalize the spurious vector ...
    Spur_PN .= Spur_PN ./ norm(Spur_PN)

    return Spur_PN
end

function QTDA_spurious_CM_initialize(Orb_2qp::qpOrb2B,qpTrOp::qpTr1B)
    # Allocate the 1- spurious vector ...
    J, P, N_qp = 1, 2, Orb_2qp.N[2,2]
    Spur_CM = Vector{Float64}(undef,N_qp)
    @inbounds for qp in 1:N_qp
        a, b, T_ab = Orb_2qp.i[P,J+1][qp].a, Orb_2qp.i[P,J+1][qp].b, Orb_2qp.i[P,J+1][qp].T
        if T_ab == -1
            ME = qpTrOp.E1.qp20.p[a,b]
            Spur_CM[qp] = ME
        elseif T_ab == 1
            ME = qpTrOp.E1.qp20.n[a,b]
            Spur_CM[qp] = ME
        end
    end

    # Renormalize the spurious vector ...
    Spur_CM .= Spur_CM ./ norm(Spur_CM)

    return Spur_CM
end

function QTDA_spurious_orthogonalize(N_qp::Int64,Spur::Vector{Float64})
    # Orthogonalization done by Householder reflection algorithm
    # For explanation see ...
    #
    # https://blogs.mathworks.com/cleve/2016/07/25/compare-gram-schmidt-and-householder-orthogonalization-algorithms/?s_tid=answers_rc2-1_p4_BOTH
    #
    # "More numerically stable than modified Gramm-Schmidt with same numerical cost"

    # Initialize the spurious index ...
    Spur_ind = 0

    # Initialize identity matrix ... 2qp basis ...
    I = diagm(ones(Float64,N_qp))

    # Project-out the spurious states from the given 2qp basis ...
    @inbounds for a in 1:N_qp
        Projection = 0.0
        @inbounds for b in 1:N_qp
            ME = Spur[b] * I[b,a]
            Projection += ME
        end
        @inbounds for b in 1:N_qp
            ME = I[b,a] - Projection * Spur[b]
            I[b,a] = ME
        end
    end

    # Renormalize the basis ...
    @inbounds for a in 1:N_qp
        b = @views I[:,a] ./ norm(I[:,a]) 
        @views I[:,a] = b
    end

    # Allocate temporary matrix U ...
    U = zeros(Float64,N_qp,N_qp)

    # Determine the index of spurious state ...
    @views U[:,1] = I[:,1]

    @inbounds for nu in 2:N_qp
        u = @views I[:,nu]
        @inbounds for mu in 1:(nu - 1)
            u -= @views (U[:, mu]' * u) * U[:, mu]
        end
        if abs(norm(u)) < 1e-4
            Spur_ind = nu
        end
        @views U[:, nu] = u ./ norm(u)
    end
    @views U[:,Spur_ind] = Spur

        # Check orthonormality ... precision in order of 1e-15 or smaller is thanks to the Householder algorithm ... very precise
        #println("Orthogonality check (close to identity?) ...")
        #orthogonality_check = U' * U
        #display(orthogonality_check)

    return U, Spur_ind
end
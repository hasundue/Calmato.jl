struct EquilibResult
    sys::System
    X::Vector{Float64} # molar fractions of components
    Y::Vector{Float64} # molar amount of phases
    x::Array{Float64,2} # molar fractions of components in each phase
    y::Array{Float64,3} # molar fractions of constituents in each sublattice
end

function Base.display(res::EquilibResult)
    I = length(res.sys.elem)
    K = length(res.sys.phas)
    @assert size(res.x) == (K,I)
    for k in 1:K
        phas = res.sys.phas[k]
        res.Y[k] < 1e-5 && continue
        @printf "%s: %.5f mol\n" phas.name res.Y[k]
        S = length(phas.cons)
        for s in 1:S
            phas.cons[s] == [] && continue
            println("\tsublattice $s")
            for j in 1:I
                consname = res.sys.elem[j].name
                @printf "\t\t%s: %.4f\n" consname res.y[k,s,j]
            end
        end
    end
end

function equilib(sys::System, X = equiatom(sys), T = 298.15; eps = 2e-8)
    # JuMP model
    model = sys.model

    I = sys.nelem # number of components
    J = sys.ncons # number of constituents
    K = sys.nphas # number of phases
    S = sys.nlatt # maximum number of sites in a sublattice
    n = sys.nsite # number of sites in each sublattice

    # Variables
    _T = model[:_T] # temperature
    _X = model[:_X] # molar amount of components
    Y = model[:Y] # moalr amount of phases
    x = model[:x] # molar fraction of components in each phase
    y = model[:y] # site fraction in each sublattice

    # Fix temperature, T
    # TODO: Using fix() results in not obtaining a solution for some reason
    set_lower_bound(_T, T)
    set_upper_bound(_T, T)

    # Fix _X[i]
    for i in 1:I
        fix(_X[i], X[i], force = true)
    end

    # Set lower and lower bounds for variables
    Y_max = [ sum( X[i] for i in 1:I ) / sum( n[k,s] for s in 1:S ) for k in 1:K ]
    for k in 1:K
        set_lower_bound(Y[k], eps)
        set_upper_bound(Y[k], Y_max[k])
        for s in 1:S, j in 1:J
            is_fixed(y[k,s,j]) && continue
            set_lower_bound(y[k,s,j], eps)
        end
    end

    optimize!(model)

    return EquilibResult(sys, X, [ value(Y[k]) for k in 1:K ],
                                 [ value(x[k,i]) for k in 1:K, i in 1:I ],
                                 [ value(y[k,s,j]) for k in 1:K, s in 1:S, j in 1:J ])
end

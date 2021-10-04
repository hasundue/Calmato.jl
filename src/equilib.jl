struct EquilibResult
    sys::System
    T::Real # temperature
    X::Vector{Float64} # molar fractions of components
    Y::Vector{Float64} # molar amount of phases
    x::Array{Float64,2} # molar fractions of components in each phase
    y::Array{Float64,3} # molar fractions of constituents in each sublattice
end

function Base.display(res::EquilibResult)
    elems = res.sys.elems
    conss = res.sys.conss
    phass = res.sys.phass
    I = length(elems)
    J = length(conss)
    K = length(phass)
    @assert size(res.x) == (K,I)
    for k in 1:K
        phas = phass[k]
        res.Y[k] < 1e-5 && continue
        @printf "%s; %s: %.4f\n" phas.name constitutionstring(phas) res.Y[k] / sum(res.Y)
        !any(latt -> length(latt) > 1, phas.cons) && continue
        for i in 1:I
            res.x[k,i] < 1e-5 && continue
            @printf "\t%s: %.4f\n" elems[i].name res.x[k,i]
        end
        cons = filter(latt -> !isempty(latt), phas.cons)
        S = length(cons)
        S < 2 && length(cons[1]) == length(elems) && continue
        for s in 1:S
            length(phas.cons[s]) < 2 && continue
            if length(filter(l -> l ≠ [], phas.cons)) > 1
                println("\t" * constitutionstring(phas, s))
            end
            for j in 1:J
                consname = conss[j]
                if consname in cons[s]
                    res.y[k,s,j] < 1e-4 && continue
                    @printf "\t\t%s: %.4f\n" subscript(consname) res.y[k,s,j]
                end
            end
        end
    end
end

Base.print(res::EquilibResult) = Base.display(res)

function equilib(sys::System, X = equiatom(sys), T = 298.15)
    # JuMP model
    model = sys.model

    I = sys.nelem # number of components
    J = sys.ncons # number of constituents
    K = sys.nphas # number of phases
    S = sys.nlatt # maximum number of sites in a sublattice
    n = sys.nsite # number of sites in each sublattice

    if length(X) ≠ I
        error("Invalid length of molar amount vector of components")
    end

    # Variables
    _T = model[:T] # temperature
    _X = model[:X] # molar amount of components
    Y = model[:Y] # moalr amount of phases
    x = model[:x] # molar fraction of components in each phase
    y = model[:y] # site fraction in each sublattice

    # Fix temperature, T
    # TODO: Using fix() results in calling f(T) at T = 0 for some reason
    Tl, Tu = sys.temp
    if T < Tl || T > Tu
        error("Temperature out of range")
    end
    set_lower_bound(_T, T)
    set_upper_bound(_T, T)

    # Fix X[i]
    # TODO: Using fix() results in X[i] = 0 for some reason
    for i in 1:I
        set_lower_bound(_X[i], X[i])
        set_upper_bound(_X[i], X[i])
    end

    # Determine maximum values of Y[k]
    Y_max = [ sum( X[i] for i in 1:I ) / sum( n[k,s] for s in 1:S ) for k in 1:K ]
    for k in 1:K
        set_upper_bound(Y[k], Y_max[k])
        for s in 1:S, j in 1:J
            is_fixed(y[k,s,j]) && continue
        end
    end

    @debug model

    optimize!(model)

    return EquilibResult(sys, T, X, [ value(Y[k]) for k in 1:K ],
                                    [ value(x[k,i]) for k in 1:K, i in 1:I ],
                                    [ value(y[k,s,j]) for k in 1:K, s in 1:S, j in 1:J ])
end

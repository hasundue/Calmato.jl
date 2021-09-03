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

function equilib(sys::System, X = equiatom(sys), T = 298.15; eps = 1e-7)
    elem = sys.elem
    phas = sys.phas

    # 
    # TO DO:
    # Consitituents are identical to components in the current version of Calmato.
    # 
    cons = elem

    I = length(elem) # number of components
    J = length(cons) # number of constituents
    K = length(phas) # number of phases

    # Maximum number of sublattices in a phase
    S = maximum([ length(phas[k].cons) for k in 1:K ])

    # Indexing constituent
    cons_id = Dict{AbstractString,Int}()
    for j in 1:J
        push!(cons_id, cons[j].name => j)
    end

    # Put zero or empty elements in each sublattice for later use
    for k in 1:K
        @assert length(phas[k].equi) == length(phas[k].cons)
        while length(phas[k].equi) < S
            push!(phas[k].equi, zero(Int))
            push!(phas[k].cons, Vector{Int}[])
        end
    end

    # Reconstuction of constitution vector into an array of Int
    constitution = [ map.(name -> get(cons_id, name, 0), phas[k].cons) for k in 1:K ]

    # Number of sites on sublattice s in phase f, n[f,s]
    n = [ phas[k].equi[s] for k in 1:K, s in 1:S ]
    
    #
    # EAGO Optimizer
    #
    model = Model(EAGO.Optimizer)
    set_optimizer_attributes(model, "verbosity" => 1,
                                    "output_iterations" => 1,
                                    "iteration_limit" => 1000,
                                    "absolute_tolerance" => 1e-3,
                                    "relative_tolerance" => 1e-3,
                                    "obbt_tolerance" => eps,
                                    "dbbt_tolerance" => eps,
                                    "absolute_constraint_feas_tolerance" => eps)
    
    # We have to use a special function xlogx() provided by EAGO.
    EAGO.register_eago_operators!(model)

    #
    # Definition of variables.
    # 
    # We define sparse arrays of variables and add constraints for zero elements,
    # instead of defining as many variables as needed, for the sake of programmability
    # and readability of the code, assuming the optimizer is very efficient.
    # 
    # Y[k]: molar amount of each phase
    # x[k,i]: molar fraction of component i in phase f
    # y[k,s,j]: site fraction of a constituent j in sublattice s of phase f
    #
    @variable(model, _T)

    Y_max = [ sum( X[i] for i in 1:I ) / sum( n[k,s] for s in 1:S ) for k in 1:K ]
    @variable(model, 0 <= Y[k=1:K] <= Y_max[k])

    x_max = [ sum( n[k,s] for s in 1:S if i in constitution[k][s] ) / sum( n[k,s] for s in 1:S ) for k in 1:K, i in 1:I ]
    @variable(model, 0 <= x[k=1:K,i=1:I] <= x_max[k,i])

    @variable(model, eps <= y[k=1:K,s=1:S,j=1:J] <= 1)

    # Relationship between molar amount of components and phases
    for i in 1:I
        @constraint(model, sum( sum( n[k,s] for s in 1:S ) * Y[k] * x[k,i]  for k in 1:K ) == X[i])
    end

    # Sum of x[k,i] in each phase equals to unity
    for k in 1:K
        @constraint(model, sum( x[k,i] for i in 1:I ) == 1)
    end

    # Fix site fractions of zero elements to zero
    for k in 1:K, j in 1:J, s in 1:S
        if constitution[k][s] == [] || !in(j, constitution[k][s])
            fix(y[k,s,j], 0, force = true)
        end
    end

    # Sum of y[k,j,s] for all constituents in a sublattice equals to unity.
    for k in 1:K, s in 1:S
        if constitution[k][s] ≠ []
            @constraint(model, sum(y[k,s,j] for j in 1:J) == 1)
        end
    end

    # Relation between mole fractions of components and site fractions
    for k in 1:K, i in 1:I
        @constraint(model, x[k,i] == sum( n[k,s] * y[k,s,i] for s in 1:S if constitution[k][s] ≠ [] )
                                        / sum( n[k,s] for s in 1:S if constitution[k][s] ≠ [] ))
    end

    # Gibbs energy contribution from each parameter
    G_phase = JuMP.NonlinearExpression[]
    for k in 1:K
        G_param = JuMP.NonlinearExpression[]

        m = 0
        for para in phas[k].para
            m += 1
            l = para.order

            # Replace element names with indexes
            comb = map.(name -> get(cons_id, name, 0), para.comb)

            # Value of the parameter
            funcname = getfuncname(para.name)
            funcsym = Symbol(funcname)
            func = getfield(Calmato, funcsym)
            register(model, funcsym, 1, func, autodiff = true)
            @eval push!($G_param, @NLexpression($model, ($funcsym)($_T)))

            # Solo constituents
            s_solo = findall(x -> length(x) == 1, comb)
            solo = [ (s, comb[s][1]) for s in s_solo ]
            for (s,j) in solo
                if length(constitution[k][s]) > 1
                    G_param[m] = @NLexpression(model, G_param[m] * y[k,s,j])
                end
            end

            # Paired constituents
            s_duo = findall(latt -> length(latt) == 2, comb)
            @assert length(s_duo) ≤ 1
            if length(s_duo) == 1
                s = s_duo[1]
                i = comb[s][1]
                j = comb[s][2]
                G_param[m] = @NLexpression(model, G_param[m] * y[k,s,i] * y[k,s,j])
                for ν in 1:l 
                    G_param[m] = @NLexpression(model, G_param[m] * (y[k,s,i] - y[k,s,j]))
                end
            end
        end

        push!(G_phase, @NLexpression(model, sum(G_param[i] for i in 1:m)))

        G_phase[k] = @NLexpression(model, G_phase[k] + R*_T*sum( n[k,s]*xlogx(y[k,s,j])
                                            for s in 1:S, j in 1:J if j in constitution[k][s] && length(constitution[k][s]) > 1 ))
    end

    @NLobjective(model, Min, sum(Y[k]*G_phase[k] for k in 1:K))

    set_lower_bound(_T, T)
    set_upper_bound(_T, T)

    optimize!(model)

    return EquilibResult(sys, X, [ value(Y[k]) for k in 1:K ],
                                 [ value(x[k,i]) for k in 1:K, i in 1:I ],
                                 [ value(y[k,s,j]) for k in 1:K, s in 1:S, j in 1:J ])
end

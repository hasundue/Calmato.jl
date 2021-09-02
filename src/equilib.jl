struct EquilibResult
    sys::System
    x::Vector{Float64} # mole fraction of each phase
    y::Array{Float64,3}
end

function Base.display(res::EquilibResult)
    N = length(res.sys.elem)
    M = length(res.sys.phas)
    @assert length(res.x) == M
    for j in 1:M
        phas = res.sys.phas[j]
        println("$(phas.name): x = $(res.x[j])")
        L = length(phas.cons[1])
        for k in 1:L
            println("\tsublattice $k")
            for i in 1:N
                elname = res.sys.elem[i].name
                println("\t\t$elname: y = $(res.y[i,j,k])")
            end
        end
    end
end

function equilib(sys::System, X = equiatom(sys), T = 298.15, P = 1.0)
    N = length(sys.elem)
    M = length(sys.phas)

    # EAGO
    model = Model(EAGO.Optimizer)
    set_optimizer_attributes(model, "verbosity" => 1,
                                    "output_iterations" => 1,
                                    "iteration_limit" => 1000,
                                    "absolute_tolerance" => 1e-3,
                                    "relative_tolerance" => 1.0)
    
    EAGO.register_eago_operators!(model)
    
    # Define molar fraction of each phase
    @variable(model, 0 <= x[j=1:M] <= 1)
    @constraint(model, sum(x[j] for j in 1:M) == 1)

    # Maximum number of sublattices in a phase
    L = maximum([length(sys.phas[j].cons[1]) for j in 1:M])

    # Molar quantity of each phase
    # mol = [ sum(sys.phas[j].equi) for j in 1:M ]

    # Transform an valence vector to a 3D-array for convenience
    equi = zeros(N,M,L)
    for i in 1:N, j in 1:M
        l = length(sys.phas[j].cons[1])
        for k in 1:l
            equi[i,j,k] = sys.phas[j].equi[k]
        end
    end

    # Define site fraction variable y[i,j,k] where indices i, j, and k
    # correspond to element, phase, and sublattice, respectively.
    @variable(model, 1e-7 <= y[i=1:N,j=1:M,k=1:L] <= 1)
    for j in 1:M
        cons = sys.phas[j].cons[1]
        l = length(cons)
        for k in 1:l
            @constraint(model, sum( y[i,j,k] for i in 1:N ) == 1)
        end
    end

    # Fix blank variables as zeros
    for i in 1:N, j in 1:M, k in 1:L
        cons = sys.phas[j].cons[1]
        if k > length(cons) || !(sys.elem[i].name in cons[k])
            fix(y[i,j,k], 0, force = true)
        end
    end

    # Molar constraint on each element
    for i in 1:N
        nonzero = [ (j,k) for j in 1:M, k in 1:L if !is_fixed(y[i,j,k]) ]
        @constraint(model, sum( x[j] * equi[i,j,k] * y[i,j,k] for (j,k) in nonzero ) == X[i])
    end

    # Gibbs energy contribution from each parameter
    G_expr = JuMP.NonlinearExpression[]
    for j in 1:M
        for para in sys.phas[j].para
            comb = para.comb
            ord = para.order

            funcname = getfuncname(para.name)
            funcsym = Symbol(funcname)
            func = getfield(Calmato, funcsym)
            funcval = func(T,P)

            sololatt = findall(latt -> length(latt) == 1, comb)
            solo = Tuple{Int,Int}[]
            for k in sololatt
                i = findfirst(el -> el.name == comb[k][1], sys.elem)
                if !isnothing(i) && !is_fixed(y[i,j,k])
                    push!(solo, (i,k))
                end
            end

            S = length(solo)
            k = findfirst(latt -> length(latt) == 2, comb)

            # 
            # FIX THIS:
            # The reason of hardcoding with if-elseif-else is that calling prod() for an empty
            # collection causes Inf as a return value of the objective function for some reason.
            # 
            if isnothing(k)
                if S > 2
                    push!(G_expr,
                          @NLexpression(model, funcval * x[j] * prod( y[i,j,k] for (i,k) in solo )))
                elseif S == 1
                    i, k = solo[1]
                    push!(G_expr, @NLexpression(model, funcval * x[j] * y[i,j,k]))
                else
                    @assert S == 0
                    push!(G_expr, @NLexpression(model, funcval * x[j]))
                end
            else
                i₁, i₂ = [ findfirst(el -> el.name == elem, sys.elem) for elem in comb[k] ]
                if S > 2
                    push!(G_expr, @NLexpression(model, funcval * x[j] * prod( y[i,j,k] for (i,k) in solo )
                                                       * y[i₁,j,k] * y[i₂,j,k] * ( y[i₁,j,k] - y[i₂,j,k] )^ord))
                elseif S == 1
                    iₛ, kₛ = solo[1]
                    push!(G_expr, @NLexpression(model, funcval * x[j] * y[iₛ,j,kₛ]
                                                       * y[i₁,j,k] * y[i₂,j,k] * ( y[i₁,j,k] - y[i₂,j,k] )^ord))
                else
                    @assert S == 0
                    push!(G_expr, @NLexpression(model, funcval * x[j] * y[i₁,j,k] * y[i₂,j,k] * ( y[i₁,j,k] - y[i₂,j,k] )^ord))
                end
            end
        end
    end

    # Ideal entropy of mixing
    for i in 1:N, j in 1:M, k in 1:L
        if !is_fixed(y[i,j,k])
            push!(G_expr, @NLexpression(model, x[j] * equi[i,j,k] * R*T * xlogx(y[i,j,k])))
        end
    end

    @NLobjective(model, Min, sum(G_expr[i] for i in 1:length(G_expr)))

    println(model)

    optimize!(model)

    return EquilibResult(sys, [ value(x[j]) for j in 1:M ],
                              [ value(y[i,j,k]) for i in 1:N, j in 1:M, k in 1:L ])
end

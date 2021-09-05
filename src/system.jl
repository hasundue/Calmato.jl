"""
Notations follow from:

Sundman, B., Chen, Q. & Du, Y. A Review of Calphad Modeling of Ordered Phases. J. Phase Equilib. Diffus. 39, 678–693 (2018).  https://doi.org/10.1007/s11669-018-0671-y

"""

struct System
    elem::Vector{Element}
    phas::Vector{Phase}
    nelem::Int # number of elements
    ncons::Int # number of constituents
    nphas::Int # number of phases
    nlatt::Int # maximum number of sublattices in a phase
    nsite::Matrix{Int} # number of sites in each sublattice
    model::Model # JuMP model
end

function init_system(db::Database, elem::Vector{Element}, phas::Vector{<:Phase}; eps = 2e-8)
    # We do this because we modify Phase structs destructively
    phas = deepcopy(phas)

    # Evaluate parsed CALPHAD functions and define them as Julia functions
    for fn in db.func
        eval(Meta.parse(fn.func))
    end
    for ph in phas
        for pr in ph.para
            eval(Meta.parse(pr.func))
        end
    end

    # TODO: Consitituents are identical to components in the current version of Calmato.
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
    # GLPK as the lower solver in EAGO
    #
    glpk = GLPK.Optimizer()
    MOI.set(glpk, MOI.RawParameter("meth"), 2)
    MOI.set(glpk, MOI.RawParameter("tol_bnd"), 1e-3)
    MOI.set(glpk, MOI.RawParameter("tol_dj"), 1e-3)
    MOI.set(glpk, MOI.RawParameter("tol_piv"), 1e-5)

    @debug begin
        MOI.set(glpk, MOI.RawParameter("msg_lev"), 4) # GLP_MSG_ALL
        MOI.set(glpk, MOI.RawParameter("out_frq"), 1) # GLP_MSG_ALL
        "Enabled GLPK debugging"
    end

    #
    # Ipopt as the upper solver in EAGO
    #
    ipopt = Ipopt.Optimizer()

    MOI.set(ipopt, MOI.RawParameter("tol"), 1e-3)
    MOI.set(ipopt, MOI.RawParameter("dual_inf_tol"), 1e+2)
    MOI.set(ipopt, MOI.RawParameter("constr_viol_tol"), eps)
    MOI.set(ipopt, MOI.RawParameter("max_iter"), 100)
    MOI.set(ipopt, MOI.RawParameter("print_level"), 0)

    @debug begin
        MOI.set(ipopt, MOI.RawParameter("print_level"), 5)
        "Enabled Ipopt debugging"
    end
    
    #
    # EAGO Optimizer
    #
    model = Model(optimizer_with_attributes(EAGO.Optimizer,
        "relaxed_optimizer" => glpk,
        "upper_optimizer" => ipopt,
        "verbosity" => 0,
        "output_iterations" => 1,
        "iteration_limit" => 3,
        "obbt_depth" => 0,
        "dbbt_tolerance" => eps,
        "absolute_constraint_feas_tolerance" => eps,
    ))

    @debug begin
        set_optimizer_attributes(model, "verbosity" => 5)
        "Enabled EAGO debugging"
    end

    # We have to use a special function xlogx() provided by EAGO.
    EAGO.register_eago_operators!(model)

    #
    # Definition of variables.
    # 
    # We define sparse arrays of variables and add constraints for zero elements,
    # instead of defining as many variables as needed, for the sake of programmability
    # and readability of the code, assuming the optimizer is very efficient.
    # 
    # _X[i]: molar amount of components
    # Y[k]: molar amount of each phase
    # x[k,i]: molar fraction of component i in phase f
    # y[k,s,j]: site fraction of a constituent j in sublattice s of phase f
    #
    @variable(model, _T)
    @variable(model, _X[i=1:I] >= 0)

    @variable(model, Y[k=1:K] >= eps)

    x_max = [ sum( n[k,s] for s in 1:S if i in constitution[k][s] ) / sum( n[k,s] for s in 1:S ) for k in 1:K, i in 1:I ]
    @variable(model, 0 <= x[k=1:K,i=1:I] <= x_max[k,i])

    @variable(model, eps <= y[k=1:K,s=1:S,j=1:J] <= 1)

    # Relationship between molar amount of components and phases
    for i in 1:I
        @constraint(model, sum( sum( n[k,s] for s in 1:S ) * Y[k] * x[k,i]  for k in 1:K ) == _X[i])
    end

    # Sum of x[k,i] in each phase equals to unity
    for k in 1:K
        @constraint(model, sum( x[k,i] for i in 1:I ) == 1)
    end

    # Fix site fractions of zero elements to zero
    for k in 1:K, s in 1:S, j in 1:J
        if !in(j, constitution[k][s])
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

        G_phase[k] = @NLexpression(model, G_phase[k] + R*_T*sum( n[k,s] * xlogx(y[k,s,j])
                                            for s in 1:S, j in 1:J if j in constitution[k][s] && length(constitution[k][s]) > 1 ))
    end

    @NLobjective(model, Min, sum(Y[k]*G_phase[k] for k in 1:K))

    @debug model

    return System(elem, phas, I, J, K, S, n, model)
end

function init_system(db::Database; eps = 2e-8)
    elem = Vector{Element}()
    phas = Vector{Phase}()

    for ph in db.phas
        if 'O' in ph.model
            @warn """
            Order/disorder transition is not supported currently.
            Phase "$(ph.name)" is not included in the system.
            """
            continue
        end
        push!(phas, ph)
    end

    for el in db.elem
        el.name in ["/-", "VA"] && continue
        push!(elem, el)
    end

    return init_system(db, elem, phas, eps = eps)
end

function equiatom(sys::System)
    N = length(sys.elem)
    return fill(1//N, N)
end

function Base.display(sys::System)
    println("System:")
    nelem = length(sys.elem)
    println("\tElements: $nelem")
    nphas = length(sys.phas)
    println("\tPhases: $nphas")
end
struct System
    elems::Vector{Element}
    phass::Vector{Phase}
    nelem::Int # number of elements
    ncons::Int # number of constituents
    nphas::Int # number of phases
    nlatt::Int # maximum number of sublattices in a phase
    nsite::Matrix{Real} # number of sites in each sublattice
    temp::Tuple{Real,Real} # temprature range where all the functions are valid
    model::Model # JuMP model
end

function init_system(db::Database, elems::Vector{Element}, phass::Vector{<:Phase})
    # We do this because we modify Phase structs destructively
    phass = deepcopy(phass)

    # Evaluate parsed CALPHAD functions and define them as Julia functions and
    # determine minimum and maximum temperature simultaneously.
    for func in db.funcs
        eval(Meta.parse(func.funcstr))
    end
    Tl, Tu = -Inf, Inf
    for phas in phass
        for param in phas.params
            eval(Meta.parse(param.funcstr))
            Tl = param.temp[1] > Tl ? param.temp[1] : Tl
            Tu = param.temp[2] < Tu ? param.temp[2] : Tu
        end
    end
    temp = (Tl, Tu)

    # The elements always exists as consitituents
    I = length(elems) # number of components
    conss = [ elems[i].name for i in 1:I ]

    # Register non-element constituents
    for phas in phass
        for param in phas.params
            for latt in param.comb
                for cons in latt
                    !(cons in conss) && push!(conss, cons)
                end
            end
        end
    end

    J = length(conss) # number of constituents
    K = length(phass) # number of phases

    # Maximum number of sublattices in a phase
    S = maximum([ length(phass[k].cons) for k in 1:K ])

    # Indexing constituent
    cons_ids = Dict{AbstractString,Int}()
    for j in 1:J
        push!(cons_ids, conss[j] => j)
    end

    # Put zero or empty elements in each sublattice for later use
    for k in 1:K
        @assert length(phass[k].sites) == length(phass[k].cons)
        while length(phass[k].sites) < S
            push!(phass[k].sites, zero(Int))
            push!(phass[k].cons, Vector{Int}[])
        end
    end

    # Reconstuction of constitution vector into an array of Int
    constitution = [ map.(name -> get(cons_ids, name, 0), phass[k].cons) for k in 1:K ]

    # Chemical stoichiometry of constituents
    comps = [ stoichiometry(conss[j]) for j in 1:J ]

    # Stoichiometry matrix, a[k,s,j,i]
    a = [ j in constitution[k][s] ? get(comps[j], elems[i].name, 0) : 0 
          for k in 1:K, s in 1:S, j in 1:J, i in 1:I ]

    # Number of sites on sublattice s in phase k, n[k,s]
    n = [ phass[k].sites[s] for k in 1:K, s in 1:S ]

    #
    # GLPK as the lower solver in EAGO
    #
    glpk = GLPK.Optimizer()
    MOI.set(glpk, MOI.RawParameter("meth"), 2)
    MOI.set(glpk, MOI.RawParameter("msg_lev"), 0)

    @debug begin
        MOI.set(glpk, MOI.RawParameter("msg_lev"), 4) # GLP_MSG_ALL
        MOI.set(glpk, MOI.RawParameter("out_frq"), 1)
        "Enabled GLPK debugging"
    end

    #
    # Ipopt as the upper solver in EAGO
    #
    ipopt = Ipopt.Optimizer()

    eps = 2e-8 # TODO: parameterize this

    MOI.set(ipopt, MOI.RawParameter("tol"), 1e-2)
    MOI.set(ipopt, MOI.RawParameter("dual_inf_tol"), 1e-6)
    MOI.set(ipopt, MOI.RawParameter("constr_viol_tol"), eps)
    MOI.set(ipopt, MOI.RawParameter("max_iter"), 1000)
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
    # T: temperature
    # X[i]: molar amount of components
    # Y[k]: molar amount of each phase
    # x[k,i]: molar fraction of component i in phase k
    # y[k,s,j]: site fraction of a constituent j in sublattice s of phase k
    #
    @variable(model, Tl <= T <= Tu)
    @variable(model, X[i=1:I] >= 0)
    @variable(model, Y[k=1:K] >= eps)
    # @variable(model, eps <= x[k=1:K,i=1:I] <= 1)
    @variable(model, eps <= y[k=1:K,s=1:S,j=1:J] <= 1)

    # Sum of y[k,j,s] for all constituents in a sublattice equals to unity.
    for k in 1:K, s in 1:S
        if constitution[k][s] ≠ []
            @constraint(model, sum(y[k,s,j] for j in 1:J) == 1)
        end
    end

    # Fix site fractions of zero elements to zero
    for k in 1:K, s in 1:S, j in 1:J
        if !in(j, constitution[k][s])
            fix(y[k,s,j], 0, force=true)
        end
    end

    # Relationship between x and y
    # for k in 1:K, i in 1:I
        @NLexpression(model, x[k=1:K,i=1:I],
                      sum(n[k,s] * sum(a[k,s,j,i] * y[k,s,j] 
                                       for j in 1:J if a[k,s,j,i] ≠ 0)
                          for s in 1:S if constitution[k][s] ≠ []) /
                      sum(n[k,s] * sum(a[k,s,j,i′] * y[k,s,j]
                                       for j in 1:J, i′ in 1:I if j ≠ J && a[k,s,j,i′] ≠ 0)
                          for s in 1:S if constitution[k][s] ≠ []))
    # end
    
    for k in 1:K, i in 1:I
        @NLconstraint(model, eps <= x[k,i] <= 1)
    end

    # Sum of x[k,i] in each phase equals to unity
    for k in 1:K
        @constraint(model, sum(x[k,i] for i in 1:I) == 1)
    end

    # Determine maximum values of x[k,i]
    xmax = zeros(K,I)
    for k in 1:K, i in 1:I
        @show k, i
        @NLobjective(model, Max, x[k,i])
        optimize!(model)
        xmax[k,i] = @show value(x[k,i])
    end
    return nothing

    # 
    # x[k,i]: molar fraction of component i in phase k
    #
    # TODO: We Assume that amounts of molecular like constituents are small and
    # there's no vacancy. This is because Using NLexpression in NLconstraint
    # results in not obtaining solutions for upper problems.
    # 
    @expression(model, x[k=1:K,i=1:I],
                sum(n[k,s] * sum(a[k,s,j,i] * y[k,s,j] 
                                 for j in 1:J
                                 if a[k,s,j,i] ≠ 0)
                    for s in 1:S if constitution[k][s] ≠ []) /
                sum(n[k,s] for s in 1:S))
                # TODO: Use this as the denominator
                # sum(n[k,s] * sum(a[k,s,j,i′] * y[k,s,j]
                #                  for j in 1:J, i′ in 1:I
                #                  if conss[j] ≠ "Va" && a[k,s,j,i′] ≠ 0)
                #     for s in 1:S if constitution[k][s] ≠ []))

    # Relationship between molar amount of components and phases
    for i in 1:I
        @constraint(model, sum(sum(n[k,s] for s in 1:S) * Y[k] * x[k,i] for k in 1:K) == X[i])
    end

    # Register all the functions for parameters
    for phas in phass
        for param in phas.params
            funcname = getfuncname(param.name)
            funcsym = Symbol(funcname)
            func = getfield(Calmato, funcsym)
            try
                register(model, funcsym, 1, func, autodiff=true)
            catch
                @warn "Duplicated definition of $funcname in the database"
            end
        end
    end

    function G_phas(k::Int; k_param::Int=k, disordered::Bool=false)
        Gs_param = JuMP.NonlinearExpression[]
        m = 0
        for param in phass[k_param].params
            m += 1
            l = param.order

            # Replace element names with indexes
            comb = map.(name -> get(cons_ids, name, 0), param.comb)

            # Value of the parameter
            funcname = getfuncname(param.name)
            funcsym = Symbol(funcname)
            @eval push!($Gs_param, @NLexpression($model, ($funcsym)($T)))

            # Solo constituents
            s_solo = findall(x -> length(x) == 1, comb)
            solo = [ (s, comb[s][1]) for s in s_solo ]
            for (s, j) in solo
                length(constitution[k][s]) < 2 && continue
                Gs_param[m] = @NLexpression(model, Gs_param[m] * ( disordered ? x[k,j] : y[k,s,j] ))
            end

            # Paired constituents
            s_duo = findall(latt -> length(latt) == 2, comb)
            @assert length(s_duo) ≤ 1
            length(s_duo) == 0 && continue
            s = s_duo[1]
            i = comb[s][1]
            j = comb[s][2]
            Gs_param[m] = disordered ? @NLexpression(model, Gs_param[m] * x[k,i] * x[k,j]) :
                                       @NLexpression(model, Gs_param[m] * y[k,s,i] * y[k,s,j])
            for ν in 1:l 
                Gs_param[m] = disordered ? @NLexpression(model, Gs_param[m] * ( x[k,i] - x[k,j] )) :
                                           @NLexpression(model, Gs_param[m] * ( y[k,s,i] - y[k,s,j] ))
            end
        end

        G = @NLexpression(model, sum(Gs_param[i] for i in 1:m))

        # Ideal entropy of configuration
        if !disordered
            G = @NLexpression(model, G + R*T*sum(n[k,s] * xlogx(y[k,s,j])
                                                 for s in 1:S, j in 1:J
                                                 if j in constitution[k][s] &&
                                                 length(constitution[k][s]) > 1 ))
        end

        return G
    end

    # Gibbs energy contribution from each phase
    Gs_phas = JuMP.NonlinearExpression[]
    for k in 1:K
        phas = phass[k]

        k_do = disorder_id(db, phas)
        if k_do ≠ nothing # ordered phase
            push!(Gs_phas, G_phas(k, k_param=k_do, disordered=true)) # disordered part
            # TODO: + 1.0 is an adhock and unphysical parameter of "penalty" on the ordered phase.
            @eval $Gs_phas[$k] = @NLexpression($model, $Gs_phas[$k] + $(G_phas(k)) + 1.0)
            @eval $Gs_phas[$k] = @NLexpression($model, $Gs_phas[$k] - $(G_phas(k, disordered=true)))
        else
            @eval push!($Gs_phas, $(G_phas(k)))
        end
    end

    @NLobjective(model, Min, sum(Y[k] * Gs_phas[k] for k in 1:K))

    @debug model

    return System(elems, phass, I, J, K, S, n, temp, model)
end

function disorder_id(db::Database, phas::Phase)
    code = phas.type
    for char in code
        i = findfirst(type -> type.code == char, db.types)
        @assert i ≠ nothing
        tdef = db.types[i]
        args = split(tdef.args)
        length(args) < 3 && continue
        if args[2] == phas.name && args[3] == "DISORDER_PART"
            @assert args[1] == "AMEND_PHASE_DESCRIPTION"
            phasname = args[4]
            j = findfirst(phas -> phas.name == phasname, db.phass)
            @assert j ≠ nothing
            return j
        end
    end
    return nothing
end

function init_system(db::Database)
    elems = Vector{Element}()
    phass = Vector{Phase}()

    for ph in db.phass
        push!(phass, ph)
    end

    for el in db.elems
        el.name in ["/-", "Va"] && continue
        push!(elems, el)
    end

    return init_system(db, elems, phass)
end

function equiatom(sys::System)
    N = length(sys.elems)
    return fill(1 // N, N)
end

function Base.display(sys::System)
    println("System:")
    nelem = length(sys.elems)
    println("\tElements: $nelem")
    nphas = length(sys.phass)
    println("\tPhases: $nphas")
end

# TODO: Detailed output
Base.print(sys::System) = Base.display(sys)

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

function init_system(db::Database, elems::Vector{Element}, phass::Vector{<:Phase}; eps = 2e-8)
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
    temp = (Tl,Tu)

    # TODO: Consitituents are identical to components in the current version of Calmato.
    conss = elems

    I = length(elems) # number of components
    J = length(conss) # number of constituents
    K = length(phass) # number of phases

    # Maximum number of sublattices in a phase
    S = maximum([ length(phass[k].cons) for k in 1:K ])

    # Indexing constituent
    cons_ids = Dict{AbstractString,Int}()
    for j in 1:J
        push!(cons_ids, conss[j].name => j)
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

    # Number of sites on sublattice s in phase f, n[f,s]
    n = [ phass[k].sites[s] for k in 1:K, s in 1:S ]

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

    tulip = Tulip.Optimizer()
    MOI.set(tulip, MOI.RawParameter("OutputLevel"), 0)

    @debug begin
        MOI.set(tulip, MOI.RawParameter("OutputLevel"), 1)
        "Enabled Tulip debugging"
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
        "relaxed_optimizer" => tulip,
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
    @variable(model, Tl <= _T <= Tu)
    @variable(model, _X[i=1:I] >= 0)

    @variable(model, Y[k=1:K] >= eps)

    x_max = [ sum( n[k,s] for s in 1:S if i in constitution[k][s] ) / sum( n[k,s] for s in 1:S ) for k in 1:K, i in 1:I ]
    @variable(model, eps <= x[k=1:K,i=1:I] <= x_max[k,i])

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

    # Register all the functions for parameters
    for phas in phass
        for param in phas.params
            funcname = getfuncname(param.name)
            funcsym = Symbol(funcname)
            func = getfield(Calmato, funcsym)
            register(model, funcsym, 1, func, autodiff = true)
        end
    end

    function G_phas(k::Int; k_param::Int = k, disordered::Bool = false)
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
            @eval push!($Gs_param, @NLexpression($model, ($funcsym)($_T)))

            # Solo constituents
            s_solo = findall(x -> length(x) == 1, comb)
            solo = [ (s, comb[s][1]) for s in s_solo ]
            for (s,j) in solo
                if length(constitution[k][s]) > 1
                    Gs_param[m] = @NLexpression(model, Gs_param[m] * ( disordered ? x[k,j] : y[k,s,j] ))
                end
            end

            # Paired constituents
            s_duo = findall(latt -> length(latt) == 2, comb)
            @assert length(s_duo) ≤ 1
            if length(s_duo) == 1
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
        end

        if disordered
            return @NLexpression(model, sum(Gs_param[i] for i in 1:m))
        else
            return @NLexpression(model, sum(Gs_param[i] for i in 1:m) + R*_T*sum( n[k,s] * xlogx( y[k,s,j] )
                for s in 1:S, j in 1:J if j in constitution[k][s] && length(constitution[k][s]) > 1 ))
        end
    end

    # Gibbs energy contribution from each phase
    Gs_phas = JuMP.NonlinearExpression[]
    for k in 1:K
        phas = phass[k]

        k_do = disorder_id(db, phas)
        if k_do ≠ nothing # ordered phase
            push!(Gs_phas, G_phas(k, k_param = k_do, disordered = true)) # disordered part
            # TODO: +0.05 is an adhock parameter of "penalty" on the ordered phase.
            @eval $Gs_phas[$k] = @NLexpression($model, $Gs_phas[$k] + $(G_phas(k)) + 0.05)
            @eval $Gs_phas[$k] = @NLexpression($model, $Gs_phas[$k] - $(G_phas(k, disordered = true)))
        else
            @eval push!($Gs_phas, $(G_phas(k)))
        end
    end

    @NLobjective(model, Min, sum(Y[k]*Gs_phas[k] for k in 1:K))

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

function init_system(db::Database; eps = 2e-8)
    elems = Vector{Element}()
    phass = Vector{Phase}()

    for ph in db.phass
        push!(phass, ph)
    end

    for el in db.elems
        el.name in ["/-", "VA"] && continue
        push!(elems, el)
    end

    return init_system(db, elems, phass, eps = eps)
end

function equiatom(sys::System)
    N = length(sys.elems)
    return fill(1//N, N)
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
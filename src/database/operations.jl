function select(db::Database, elnames::Vector{<:AbstractString})
    elems = Element[]
    funcs = GFunction[]
    phass = Phase[]

    for elem in db.elems
        if elem.name in vcat(elnames, ["/-", "Va"])
            push!(elems, elem)
        end
    end
    for phas in db.phass
        phas_selected = select(phas, elnames)
        isnothing(phas_selected) && continue
        push!(phass, phas_selected)
    end

    # Functions used in parameters
    for func in db.funcs
        if any(phas -> isused(func, phas), phass)
            push!(funcs, func)
        end
    end

    # Functions used in functions
    function push_func!(func)
        func in funcs && return false
        if any(_func -> isused(func, _func), funcs)
            push!(funcs, func)
            return true
        else
            return false
        end
    end
    found = true
    while found
        found = any(push_func!, db.funcs)
    end
    return Database(elems, funcs, phass, db.types)
end

function select(db::Database, elnames::AbstractString)
    select(db, split(elnames))
end

function select(phas::Phase, elnames::Vector{<:AbstractString})
    cons, nonempty = select(phas.cons, elnames)
    [] in cons && return nothing

    sites = [ phas.sites[s] for s in nonempty ]

    params = Parameter[]
    for param in phas.params
        if iscomposedof(param.comb, vcat(elnames, "Va"))
            comb = [ param.comb[s] for s in nonempty ]
            push!(params, Parameter(param.name, 
                                    param.symbol,
                                    comb,
                                    param.order,
                                    param.temp,
                                    param.funcstr))
        end
    end

    return Phase(phas.name, phas.state, phas.type, sites, cons, params)
end

function select(cons::Constitution, elnames::Vector{<:AbstractString})
    latts = Sublattice[]
    nonempty = Int[]
    S = length(cons)
    for s in 1:S
        latt = filter(spec -> iscomposedof(spec, vcat(elnames, "Va")), cons[s])
        latt == [] || latt == ["Va"] && continue
        push!(latts, latt)
        push!(nonempty, s)
    end
    return latts, nonempty
end

function Base.merge!(db_master::Database, db_source::Database,
                     args::Union{Pair{Int,Int},Int,Function}...)
    # Elements
    elemnames_master = map(elem -> getfield(elem, :name), db_master.elems)
    for elem in db_source.elems
        elem.name in elemnames_master && continue
        push!(db_master.elems, elem)
    end

    # Phases
    merged = Int[]
    for arg in args
        if typeof(arg) == Pair{Int,Int}
            ks = arg.first
            km = arg.second
            merge!(db_master.phass[km], db_source.phass[ks])
            push!(merged, ks)
        elseif typeof(arg) == Int
            if arg < 0
                push!(merged, -arg)
                continue
            end
            arg in merged && continue
            push!(db_master, db_source.phass[arg])
            push!(merged, arg)
        elseif arg == *
            K = length(db_source.phass)
            for k in 1:K
                k in merged && continue
                push!(db_master, db_source.phass[k])
                push!(merged, k)
            end
        else
            @error "unsupported type of argument: $arg"
        end
    end

    # Functions used in parameters
    for func in db_source.funcs
        if any(phas -> isused(func, phas), db_master.phass)
            push!(db_master.funcs, func)
        end
    end

    # Functions used in functions
    function push_func!(func)
        func in db_master.funcs && return false
        if any(_func -> isused(func, _func), db_master.funcs)
            push!(db_master.funcs, func)
            return true
        else
            return false
        end
    end
    found = true
    while found
        found = any(push_func!, db_source.funcs)
    end

    return nothing
end

function Base.merge!(phas_master::Phase, phas_source::Phase)
    if phas_master.sites ≠ phas_source.sites
        @error "Stoichiometry of both phases must match ($(phas_master.name))"
    end

    if length(phas_master.cons) ≠ length(phas_source.cons)
        @error "Constitution of both phases must match ($(phas_master.name))"
    end

    # Constitution
    S = length(phas_master.cons)
    for s in 1:S
        latt = phas_source.cons[s]
        for spec in latt
            if !(spec in phas_master.cons[s])
                push!(phas_master.cons[s], spec)
            end
        end
    end

    # Parameters
    paramnames_master = map(localname, phas_master.params)
    for param in phas_source.params
        name = localname(param)
        if name in paramnames_master
            @info "$(param.symbol)($name) already exists in $(phas_master.name)"
            continue
        end
        push!(phas_master.params, param)
    end
end

function Base.push!(db::Database, phas::Phase)
    push!(db.phass, phas)
end

function Base.merge(db_master::Database, db_source::Database,
                    args::Union{Pair{Int,Int},Int,Function}...)
    db = deepcopy(db_master)
    merge!(db, db_source, args...)
    return db
end

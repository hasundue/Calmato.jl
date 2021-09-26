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

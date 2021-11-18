function read_db(toml::AbstractString)
    dict_db = TOML.parsefile(toml)

    elems = Element[]
    dicts = get(dict_db, "element", Dict())
    for dict in dicts
        name = get(dict, "name", "")
        refstate = get(dict, "std", "")
        mass = get(dict, "mass", 0.0)
        H298 = get(dict, "H298", 0.0)
        S298 = get(dict, "S298", 0.0)
        elem = Element(name, refstate, mass, H298, S298)
        push!(elems, elem)
    end

    funcs = GFunction[]
    dicts = get(dict_db, "function", Dict())
    for dict in dicts
        name = get(dict, "name", "")
        expr_union = get(dict, "expr", nothing)
        if typeof(expr_union) <: AbstractString
            tmin = get(dict, "tmin", -Inf)
            tmax = get(dict, "tmax", Inf)
            temps = [tmin, tmax]
            func = GFunction(name, temps, [expr_union])
        elseif typeof(expr_union) <: Vector
            temps = Real[]
            exprs = AbstractString[]
            for expr_dict in expr_union
                tmin = get(expr_dict, "tmin", -Inf)
                push!(temps, tmin)

                str = get(expr_dict, "expr", "")
                push!(exprs, str)
            end
            tmax = get(expr_union[end], "tmax", Inf)
            push!(temps, tmax)
            func = GFunction(name, temps, exprs)
        else
            @error "unsupported format for a function $name"
        end
        push!(funcs, func)
    end

    phass = Phase[]
    dicts = get(dict_db, "phase", Dict[])
    for dict in dicts
        phasname = get(dict, "name", "")
        state = get(dict, "state", "")[1]

        str = get(dict, "constitution", "")
        rms = collect(eachmatch(r"\(([^\(\)]*)\)(\d{0,})", str))
        caps = map(rm -> rm.captures, rms)
        sites = [ cap[2] == "" ? 1 : parse(Int, cap[2]) for cap in caps ]
        cons = [ split(cap[1], ",") for cap in caps ]

        params = Parameter[]
        params_dict = get(dict, "param", Dict[])
        for dict in params_dict
            form = get(dict, "form", "")

            funcname = getfuncname(form)
            name = ""
            name *= funcname[1]
            name *= phasname
            name *= funcname[2:end]

            symbol = form[1]

            str = split(form[3:end-1], ';') # Cu,Zn:Zn;0
            @assert length(str) ≤ 2
            comb = parse_combination(str[1])
            order = length(str) == 1 ? 0 : parse(Int, str[2])

            tmin = get(dict, "tmin", -Inf)
            tmax = get(dict, "tmax", Inf)
            temps = [tmin, tmax]

            expr = get(dict, "expr", "")

            param = Parameter(form, symbol, comb, order, temps, [expr], name)
            push!(params, param)
        end

        phas = Phase(phasname, state, sites, cons, params)
        push!(phass, phas)
    end

    return Database(elems, funcs, phass)
end

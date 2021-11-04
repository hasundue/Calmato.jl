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
    println(elems)

    funcs = GFunction[]
    dicts = get(dict_db, "function", Dict())
    for dict in dicts
        name = get(dict, "name", "")
        expr = get(dict, "expr", nothing)
        if typeof(expr) <: AbstractString
            tmin = get(dict, "tmin", -Inf)
            tmax = get(dict, "tmax", Inf)
            temp = (tmin, tmax)
            func = GFunction(name, temp, expr)
        elseif typeof(expr) <: Vector
            temps = Real[]
            exprs = AbstractString[]
            for ex in expr
                tmin = get(ex, "tmin", -Inf)
                push!(temps, tmin)
                str = get(ex, "expr", "")
                push!(exprs, str)
            end
            tmax = get(expr[end], "tmax", Inf)
            push!(temps, tmax)
            func = GFunction(name, temps, exprs)
        else
            @error "unsupported format for a function $name"
        end
        push!(funcs, func)
    end
    println(funcs)

    phass = Phase[]
    dicts = get(dict_db, "phase", Dict())
    for dict in dicts
        name = get(dict, "formula", "")
        state = get(dict, "state", "")[1]
        rms = collect(eachmatch())
    end
end

function parse_function_yaml(str::AbstractString)
end

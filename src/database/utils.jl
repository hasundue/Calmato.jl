function iscomposedof(cons::Constitution, elems::Vector{<:AbstractString})
    for latt in cons
        for spec in latt
            !iscomposedof(spec, elems) && return false
        end
    end
    return true
end

function isused(func::GFunction, phas::Phase)
    for param in phas.params
        occursin(func.name, param.funcstr) && return true
    end
    return false
end

function isused(func1::GFunction, func2::GFunction)
    return occursin(func1.name, func2.funcstr)
end

function constitutionstring(phas::Phase, s::Int)
    cons = phas.cons
    str = ""
    if length(cons[s]) > 1
        str *= '('
    end
    for spec in cons[s]
        str *= subscript(spec)
        if spec ≠ cons[s][end]
            str *= ','
        end
    end
    if length(cons[s]) > 1
        str *= ')'
    end
    n = phas.sites[s]
    if n ≠ one(n)
        str *= subscript(n)
    end
    return str
end

function constitutionstring(phas::Phase)
    S = length(phas.cons)
    str = ""
    for s in 1:S
        if phas.cons[s] ≠ []
            str *= constitutionstring(phas, s)
        end
    end
    return str
end

function print_func(db::Database)
    for func in db.funcs
        println(func.funcstr)
    end
end

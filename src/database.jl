struct Element
    name::AbstractString
    refstate::AbstractString # reference state
    mass::Float64
    H298::Float64
    S298::Float64
end

struct GFunction{T<:Real}
    name::AbstractString
    temp::Tuple{T,T}
    funcstr::AbstractString
end

struct Parameter{S<:AbstractString, T<:Real}
    name::AbstractString
    symbol::Char
    comb::Vector{Vector{S}}
    order::Int
    temp::Tuple{T,T}
    funcstr::AbstractString
end

const Constitution = Vector{Vector{AbstractString}}

mutable struct Phase{T<:Real}
    name::AbstractString
    state::Char
    type::AbstractString
    sites::Vector{T}
    cons::Constitution
    params::Vector{Parameter}

    function Phase(name::AbstractString,
                   state::Char,
                   type::AbstractString,
                   sites::Vector{T}) where T <: Real
        phas = new{T}()
        phas.name = name
        phas.state = state
        phas.type = type
        phas.sites = sites
        phas.params = Parameter[]
        return phas
    end
end

struct Database
    elems::Vector{Element}
    funcs::Vector{GFunction}
    phass::Vector{Phase}
end

function read_tdb(tdb::AbstractString)
    elems = Vector{Element}()
    funcs = Vector{GFunction}()
    phass = Vector{Phase}()

    open(tdb, "r") do io
        while !eof(io)
            block = readuntil(io, '!')

            for keyword in [" ELEMENT ", "Function", "Phase", "Constituent", "Parameter"]
                range = findfirst(keyword, block)
                if range ≠ nothing
                    text = block[first(range):end]
                    if keyword == " ELEMENT "
                        elem = parse_element(text)
                        push!(elems, elem)
                    elseif keyword == "Function"
                        func = parse_gfunction(text, funcs)
                        push!(funcs, func)
                    elseif keyword == "Phase"
                        phas = parse_phase(text)
                        push!(phass, phas)
                    elseif keyword == "Constituent"
                        phasname, cons = parse_constituent(text)
                        phas = filter(phas -> phas.name == phasname, phass)
                        @assert length(phas) == 1
                        phas[1].cons = cons
                    elseif keyword == "Parameter"
                        phasname, param = parse_parameter(text, funcs)
                        phas = filter(phas -> phas.name == phasname, phass)
                        @assert length(phas) == 1
                        push!(phas[1].params, param)
                    end
                end
            end
        end
    end

    return Database(elems, funcs, phass)
end

function parse_element(text::AbstractString)
    strs = split(text)
    @assert length(strs) == 6
    @assert strs[1] == "ELEMENT"
    name = strs[2]
    refstate = strs[3]
    mass, H298, S298 = parse.(Float64, strs[4:6])
    return Element(name, refstate, mass, H298, S298)
end

function parse_gfunction(text::AbstractString, funcs::Vector{GFunction})
    strs = split(text)
    @assert strs[1] == "Function"
    name = strs[2]
    temp, funcstr = parse_function(name, strs[3:end], funcs)
    return GFunction(name, temp, funcstr)
end

function parse_function(name::AbstractString, strs::Vector{SubString{String}}, funcs::Vector{GFunction})
    Ts = Vector{AbstractString}()
    funcstrs = Vector{AbstractString}()

    push!(Ts, strs[1])
    k = 2

    while true
        push!(funcstrs, strs[k])
        k += 1
        while isnothing(findfirst(';', funcstrs[end]))
            funcstrs[end] = funcstrs[end] * strs[k]
            k += 1
        end
        push!(Ts, strs[k])
        k += 1
        if strs[k] == "N"
            k += 1
            break
        elseif strs[k] == "Y"
            k += 1
            continue
        else
            error("unexpected format in the tdb file")
        end
    end

    N = length(funcstrs)
    str = "function $(name)(T)::typeof(T)\n"
    for i in 1:N
        for pair in ["T*LN(T)" => "McCormick.xlogx(T)", "LN(" => "log(", "**" => "^", ".+" => ".0+", ".-" => ".0-", ".*" => ".0*", "./" => ".0/"]
            funcstrs[i] = replace(funcstrs[i], pair)
        end
        funcstrs[i] = strip(funcstrs[i], ';')

        str *= i == 1 ? "\tif " : "\telseif "
        str *= "$(Ts[i]) ≤ T ≤ $(Ts[i+1])\n"
        for func in funcs
            funcstrs[i] = replace(funcstrs[i], func.name => "$(func.name)(T)")
        end
        str *= "\t\t" * funcstrs[i] * "\n"
    end
    str *= "\telseif T ≤ $(Ts[1])\n\t\t$(funcstrs[1])\n"
    str *= "\telseif T ≥ $(Ts[end])\n\t\t$(funcstrs[end])\n"
    str *= "\tend\nend\n"

    Tl, Tu = tryparse.(Int, [Ts[1], Ts[end]])
    if nothing in [Tl, Tu]
        Tl, Tu = parse.(Float64, [Ts[1], Ts[end]])
    end

    return (Tl, Tu), str
end

function parse_phase(text::AbstractString)
    vect = split(text)
    @assert vect[1] == "Phase"
    model = vect[3] # "%", "%&", "%Q+", etc.
    N = parse(Int, vect[4])
    equi = tryparse.(Int, vect[5:end])
    if nothing in equi
        equi = parse.(Float64, vect[5:end])
    end
    @assert length(equi) == N
    vect = split(vect[2], ':') # "Liquid:L", "Bcc", "Fcc", etc.
    name = vect[1]
    state = length(vect) > 1 ? vect[2][1] : 'S'
    return Phase(name, state, model, equi)
end

function parse_constituent(text::AbstractString)
    strs = split(text)
    @assert strs[1] == "Constituent"
    phas = split(strs[2], ':')[1] # "Liquid:L", "Bcc", "Fcc", etc.
    cons_text = strs[3] # ":Cu,Zn", ":Cu,Zn:Cu,Zn:", etc.
    cons_strs = split(cons_text, ':')
    cons_strs = filter(e -> e ≠ "", cons_strs)
    cons = map(text -> split(text, ','), cons_strs)
    return phas, cons
end

function parse_parameter(text::AbstractString, funcs::Vector{GFunction})
    # Parameters
    strs = split(text)
    @assert strs[1] == "Parameter"

    # L(Liquid,Cu,Zn;0), etc.
    text = strs[2]
    name = text
    funcname = getfuncname(text)
    symbol = text[1]
    k = findfirst('(', text)
    l = findfirst(')', text)
    text = text[k+1:l-1]
    k = findfirst(',', text)
    phas = text[1:k-1]
    text = text[k+1:end]
    comb_text, order_text = split(text, ';')
    comb = map(text -> split(text, ','), split(comb_text, ':'))
    order = parse(Int, order_text)
    temp, funcstr = parse_function(funcname, strs[3:end], funcs)
    return phas, Parameter(name, symbol, comb, order, temp, funcstr)
end

function getfuncname(str::AbstractString)
    filter(c -> !in(c, ['(', ')', ',', ':', ';']), str)
end

function print_func(db::Database)
    for func in db.funcs
        println(func.funcstr)
    end
end

function Base.display(db::Database)
    println("Database Summary:")
    nelem = length(db.elems)
    println("\tElements: $nelem")
    nfunc = length(db.funcs)
    println("\tFunctions: $nfunc")
    nphas = length(db.phass)
    println("\tPhases: $nphas")
end

function Base.print(db::Database)
    println("Database:")

    nelem = length(db.elems)
    println("\tElements: $nelem")
    for elem in db.elems
        println("\t\t$(elem.name)")
    end

    nfunc = length(db.funcs)
    println("\tFunctions: $nfunc")
    for func in db.funcs
        println("\t\t$(func.name)")
    end

    nphas = length(db.phass)
    println("\tPhases: $nphas")
    for k in 1:nphas
        phas = db.phass[k]
        print("\t\t$k: $(phas.name); ")
        cons = phas.cons
        nlatt = length(cons)
        for i in 1:nlatt
            print('(')
            for elem in cons[i]
                print(elem)
                if elem ≠ cons[i][end]
                    print(',')
                end
            end
            print(')')
            print(phas.sites[i])
        end
        print('\n')
        for param in phas.params
            print("\t\t\t$(param.symbol)")
            print('(')
            nlatt = length(param.comb)
            for i in 1:nlatt
            latt = param.comb[i]
                for cons in latt
                    print(cons)
                    if cons ≠ latt[end]
                        print(',')
                    end
                end
                if i ≠ nlatt
                    print(':')
                end
            end
            print(";$(param.order)")
            print(")\n")
        end
    end
end

Base.println(db::Database) = Base.print(db)
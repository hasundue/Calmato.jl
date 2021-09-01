struct Element
    name::AbstractString
    phase::AbstractString
    param::Vector{Float64}
end

struct GFunction
    name::AbstractString
    func::AbstractString
end

struct Parameter
    name::AbstractString
    symbol::Char
    comb::Vector{Vector{AbstractString}}
    order::Int
    func::AbstractString
end

const Constituent = Vector{Vector{AbstractString}}

struct Phase{T<:Real}
    name::AbstractString
    state::Char
    model::AbstractString
    equi::Vector{T}
    cons::Vector{Constituent}
    para::Vector{Parameter}
end

struct Database
    elem::Vector{Element}
    func::Vector{GFunction}
    phas::Vector{Phase}
end

function read_tdb(tdb::AbstractString)
    elem = Vector{Element}()
    func = Vector{GFunction}()
    phas = Vector{Phase}()

    open(tdb, "r") do io
        while !eof(io)
            block = readuntil(io, '!')

            for type in [" ELEMENT ", "Function", "Phase", "Constituent", "Parameter"]
                range = findfirst(type, block)
                if range ≠ nothing
                    text = block[first(range):end]
                    if type == " ELEMENT "
                        el = parse_element(text)
                        push!(elem, el)
                    elseif type == "Function"
                        fn = parse_gfunction(text, func)
                        push!(func, fn)
                    elseif type == "Phase"
                        ph = parse_phase(text)
                        push!(phas, ph)
                    elseif type == "Constituent"
                        phasname, cn = parse_constituent(text)
                        ph = filter(ph -> ph.name == phasname, phas)
                        @assert length(ph) == 1
                        push!(ph[1].cons, cn)
                    elseif type == "Parameter"
                        phasname, pr = parse_parameter(text, func)
                        ph = filter(ph -> ph.name == phasname, phas)
                        push!(ph[1].para, pr)
                    end
                end
            end
        end
    end

    return Database(elem, func, phas)
end

function parse_element(text::AbstractString)
    strs = split(text)
    @assert length(strs) == 6
    @assert strs[1] == "ELEMENT"
    name = strs[2]
    phase = strs[3]
    param = parse.(Float64, strs[4:6])
    return Element(name, phase, param)
end

function parse_gfunction(text::AbstractString, funcs::Vector{GFunction})
    strs = split(text)
    @assert strs[1] == "Function"
    name = strs[2]
    func = parse_function(name, strs[3:end], funcs)
    return GFunction(name, func)
end

function parse_function(name::AbstractString, strs::Vector{SubString{String}}, funcs::Vector{GFunction})
    T = Vector{AbstractString}()
    funcstr = Vector{AbstractString}()

    push!(T, strs[1])
    k = 2

    while true
        push!(funcstr, strs[k])
        k += 1
        while isnothing(findfirst(';', funcstr[end]))
            funcstr[end] = funcstr[end] * strs[k]
            k += 1
        end
        push!(T, strs[k])
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

    N = length(funcstr)
    str = "function $(name)(T,P)\n"
    for i in 1:N
        for pair in ["LN(" => "log(", "**" => "^", ".+" => ".0+", ".-" => ".0-", ".*" => ".0*", "./" => ".0/"]
            funcstr[i] = replace(funcstr[i], pair)
        end
        funcstr[i] = strip(funcstr[i], ';')

        str *= i == 1 ? "\tif " : "\telseif "
        str *= "$(T[i]) ≤ T ≤ $(T[i+1])\n"
        for func in funcs
            funcstr[i] = replace(funcstr[i], func.name => "$(func.name)(T,P)")
        end
        str *= "\t\t" * funcstr[i] * "\n"
    end
    str *= "\tend\nend\n"
    return str
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
    return Phase(name, state, model, equi, Vector{Constituent}(), Vector{Parameter}())
end

function parse_constituent(text::AbstractString)
    vect = split(text)
    @assert vect[1] == "Constituent"
    phas = split(vect[2], ':')[1] # "Liquid:L", "Bcc", "Fcc", etc.
    cons_text = vect[3] # ":Cu,Zn", ":Cu,Zn:Cu,Zn:", etc.
    cons_vect = split(cons_text, ':')
    cons_vect = filter(e -> e ≠ "", cons_vect)
    cons = map(text -> split(text, ','), cons_vect)
    return phas, cons
end

function parse_parameter(text::AbstractString, funcs::Vector{GFunction})
    # Parameters
    vect = split(text)
    @assert vect[1] == "Parameter"

    # L(Liquid,Cu,Zn;0), etc.
    text = vect[2]
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
    func = parse_function(funcname, vect[3:end], funcs)
    return phas, Parameter(name, symbol, comb, order, func)
end

function getfuncname(str::AbstractString)
    filter(c -> !in(c, ['(', ')', ',', ':', ';']), str)
end

function print_func(db::Database)
    for fn in db.func
        println(fn.func)
    end
end

function Base.display(db::Database)
    println("Database Summary:")
    nelem = length(db.elem)
    println("\tElements: $nelem")
    nfunc = length(db.func)
    println("\tFunctions: $nfunc")
    nphas = length(db.phas)
    println("\tPhases: $nphas")
end

function Base.print(db::Database)
    println("Database:")

    nelem = length(db.elem)
    println("\tElements: $nelem")
    for elem in db.elem
        println("\t\t$(elem.name)")
    end

    nfunc = length(db.func)
    println("\tFunctions: $nfunc")
    for func in db.func
        println("\t\t$(func.name)")
    end

    nphas = length(db.phas)
    println("\tPhases: $nphas")
    for phas in db.phas
        print("\t\t$(phas.name); ")
        cons = phas.cons[1]
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
            print(phas.equi[i])
        end
        print('\n')
        for para in phas.para
            print("\t\t\t$(para.symbol)")
            print('(')
            nlatt = length(para.comb)
            for i in 1:nlatt
            latt = para.comb[i]
                for el in latt
                    print(el)
                    if el ≠ latt[end]
                        print(',')
                    end
                end
                if i ≠ nlatt
                    print(':')
                end
            end
            print(";$(para.order)")
            print(")\n")
        end
    end
end

Base.println(db::Database) = Base.print(db)
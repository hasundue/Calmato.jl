struct Element
    name::AbstractString
    refstate::AbstractString # reference state
    mass::Float64
    H298::Float64
    S298::Float64
end

abstract type AbstractGFunction end

mutable struct GFunction{T<:Real} <: AbstractGFunction
    name::AbstractString
    temp::Tuple{T,T}
    funcstr::AbstractString
end

mutable struct Parameter{S<:AbstractString, T<:Real} <: AbstractGFunction
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
    type::AbstractString
    sites::Vector{T}
    cons::Constitution
    params::Vector{Parameter}

    function Phase(name::AbstractString,
                   type::AbstractString,
                   sites::Vector{T}) where T <: Real
        phas = new{T}()
        phas.name = name
        phas.type = type
        phas.sites = sites
        phas.params = Parameter[]
        return phas
    end
end

function constitutionstring(phas::Phase, s::Int)
    cons = phas.cons
    str = '('
    for elname in cons[s]
        str *= elname
        if elname ≠ cons[s][end]
            str *= ','
        end
    end
    str *= ')'
    n = phas.sites[s]
    if n ≠ one(n)
        str *= string(n)
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

struct TypeDefinition
    code::Char
    key::AbstractString # SEQ, GES, etc.
    args::AbstractString
end

struct Database
    elems::Vector{Element}
    funcs::Vector{GFunction}
    phass::Vector{Phase}
    types::Vector{TypeDefinition}
end

const keywordlist = [
    "ELEMENT",
    "FUNCTION",
    "FUNCT",
    "TYPE_DEFINITION",
    "PHASE",
    "CONSTITUENT",
    "PARAMETER",
]

function read_tdb(io::IO)
    elems = Vector{Element}()
    funcs = Vector{GFunction}()
    phass = Vector{Phase}()
    types = Vector{TypeDefinition}()

    # Parse each elements
    while !eof(io)
        block = readline(io)

        block = replace(block, '\t' => "") # remove tabs

        isempty(block) && continue # blank line
        block[1] == '$' && continue # comment line

        i = findfirst('!', block)
        if isnothing(i)
            block *= readuntil(io, '!')
        else
            block = block[1:i-1] # remove '!'
        end

        strs = split(block)
        isempty(strs) && continue

        keyword = uppercase(strs[1])
        
        !(keyword in keywordlist) && continue
        try
            if keyword == "ELEMENT"
                elem = parse_element(block)
                push!(elems, elem)
            elseif keyword in ["FUNCTION", "FUNCT"]
                func = parse_gfunction(block)
                push!(funcs, func)
            elseif keyword == "TYPE_DEFINITION"
                type = parse_typedefinition(block)
                push!(types, type)
            elseif keyword == "PHASE"
                phas = parse_phase(block)
                push!(phass, phas)
            elseif keyword == "CONSTITUENT"
                phasname, cons = parse_constituent(block)
                phas = filter(phas -> uppercase(phas.name) == uppercase(phasname), phass)
                @assert length(phas) == 1
                phas[1].cons = cons
            elseif keyword == "PARAMETER"
                phasname, param = parse_parameter(block)
                phas = filter(phas -> uppercase(phas.name) == uppercase(phasname), phass)
                @assert length(phas) == 1
                push!(phas[1].params, param)
            end
        catch
            error("Invalid description in the tdb file: $block")
        end
    end

    # "Julialize" the functions
    for func in funcs
        julialize_funcstr!(func, funcs)
    end
    for phas in phass
        for param in phas.params
            julialize_funcstr!(param, funcs)
        end
    end

    return Database(elems, funcs, phass, types)
end

function julialize_funcstr!(arg::AbstractGFunction, funcs::Vector{GFunction})
    for pair in ["T*LN(T)" => "xlogx(T)", "LN(" => "log(", "**" => "^", ".+" => ".0+", ".-" => ".0-", ".*" => ".0*", "./" => ".0/"]
        arg.funcstr = replace(arg.funcstr, pair)
    end
    for func in funcs
        func.name == arg.name && continue
        arg.funcstr = replace(arg.funcstr, func.name => "$(func.name)(T)")
    end
end

function read_tdb(tdb::AbstractString)
    open(tdb, "r") do io
        read_tdb(io)
    end
end

function parse_element(text::AbstractString)
    strs = split(text)
    @assert length(strs) == 6
    @assert match(r"ELEMENT"i, strs[1]) ≠ ""
    name = strs[2]
    refstate = strs[3]
    mass, H298, S298 = parse.(Float64, strs[4:6])
    return Element(name, refstate, mass, H298, S298)
end

function parse_gfunction(text::AbstractString)
    strs = split(text)
    @assert match(r"FUNCTION"i, strs[1]) ≠ ""
    name = strs[2]
    temp, funcstr = parse_function(name, strs[3:end])
    return GFunction(name, temp, funcstr)
end

function parse_function(name::AbstractString, strs::Vector{SubString{String}})
    Ts = Vector{AbstractString}()
    funcstrs = Vector{AbstractString}()

    push!(Ts, strs[1]) # ex. 298.15
    k = 2

    while true
        push!(funcstrs, replace(strs[k], "#" => ""))
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
        funcstrs[i] = strip(funcstrs[i], ';')

        str *= i == 1 ? "\tif " : "\telseif "
        str *= "$(Ts[i]) ≤ T ≤ $(Ts[i+1])\n"
        str *= "\t\t" * funcstrs[i] * "\n"
    end
    str *= "\tend\nend\n"

    Tl, Tu = tryparse.(Int, [Ts[1], Ts[end]])
    if nothing in [Tl, Tu]
        Tl, Tu = parse.(Float64, [Ts[1], Ts[end]])
    end

    return (Tl, Tu), str
end

function parse_typedefinition(text::AbstractString)
    strs = split(text)
    @assert match(r"TYPE_DEFINITION"i, strs[1]) ≠ ""
    code = strs[2]
    @assert length(code) == 1
    code = code[1]
    key = strs[3]
    if !(key in ["SEQ", "GES"])
        @warn "Unsupported keyword, $key, for TYPE_DEFINITION"
    end
    args = join(strs[4:end], ' ')
    return TypeDefinition(code, key, args)
end

function parse_phase(text::AbstractString)
    strs = split(text)
    @assert match(r"PHASE"i, strs[1]) ≠ ""
    model = strs[3] # "%", "%&", "%Q+", etc.
    N = parse(Int, strs[4])
    sites = tryparse.(Int, strs[5:end])
    if nothing in sites
        sites = parse.(Float64, strs[5:end])
    end
    @assert length(sites) == N
    name = replace(strs[2], " " => "")
    name = split(name, ':')[1]
    return Phase(name, model, sites)
end

function parse_constituent(text::AbstractString)
    strs = split(text)
    @assert match(r"CONSTITUENT"i, strs[1]) ≠ ""
    phas = split(strs[2], ':')[1] # "Liquid", "Bcc", "Fcc", etc.
    cons_text = join(strs[3:end]) # ":Cu,Zn", ":Cu,Zn:Cu,Zn:", etc.
    cons_strs = split(cons_text, ':')
    cons_strs = filter(e -> e ≠ "", cons_strs)
    cons = map(text -> split(text, ','), cons_strs)
    return phas, cons
end

function parse_parameter(text::AbstractString)
    # Parameters
    strs = split(text)
    @assert match(r"PARAMETER"i, strs[1]) ≠ ""

    text = join(strs[2:end], ' ') # ex. G(BCC2,...

    k = findfirst('(', text)
    l = findfirst(')', text)
    name = replace(text[1:l], " " => "") # G(BCC_B2,Cu:Cu,Zn;0)
    funcname = getfuncname(name) # GBCC_B2CuCuZn0
    symbol = name[1] # G

    strs = split(text[l+1:end]) # 298.15...

    text = text[k+1:l-1] # BCC_B2,Cu,Zn:Cu,Zn;0

    k = findfirst(',', text)
    phas = replace(text[1:k-1], " " => "") # BCC_B2

    text = text[k+1:end] # Cu:Cu,Zn;0
    comb_text, order_text = split(text, ';') # Cu:Cu,Zn, 0
    comb_strs = split(comb_text, ':') # ["Cu", "Cu,Zn"]
    comb = map(text -> split(text, ','), comb_strs) # [[Cu], [Cu,Zn]]
    order = parse(Int, order_text)

    temp, funcstr = parse_function(funcname, strs)

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
        print(constitutionstring(phas))
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
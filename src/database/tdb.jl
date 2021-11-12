const keywordlist = [
    "ELEMENT",
    "FUNCTION",
    "FUNCT",
    "TYPE_DEFINITION",
    "PHASE",
    "CONSTITUENT",
    "CONST",
    "PARAMETER",
    "PAR",
]

function read_tdb(io::IO)
    elems = Vector{Element}()
    funcs = Vector{GFunction}()
    phass = Vector{Phase}()
    types = Vector{TypeDefinition}()

    # Parse each elements
    while !eof(io)
        block = readline(io, keep = true)

        block = replace(block, '\t' => "") # remove tabs

        isempty(block) && continue # blank line
        !occursin(r"\S+", block) && continue # blank line
        occursin(r"\$", block) && continue # comment line

        i = findfirst('!', block)
        if isnothing(i)
            block *= readuntil(io, '!')
        else
            block = block[1:i-1] # remove '!'
        end

        block = replace(block, '\r' => "")
        block = replace(block, '\n' => ' ')
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
            elseif keyword in ["CONSTITUENT", "CONST"]
                phasname, cons = parse_constituent(block)
                phas = filter(phas -> uppercase(phas.name) == uppercase(phasname), phass)
                @assert length(phas) == 1
                phas[1].cons = cons
            elseif keyword in ["PARAMETER", "PAR"]
                phasname, param = parse_parameter(block)
                isnothing(param) && continue
                phas = filter(phas -> uppercase(phas.name) == uppercase(phasname), phass)
                @assert length(phas) == 1
                push!(phas[1].params, param)
            end
        catch
            error("Unrecognized description in the tdb file:\n$block")
        end
    end

    return Database(elems, funcs, phass, types)
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
    name = format_element(strs[2])
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
    K = parse(Int, strs[4])
    i = 5 + K - 1
    sites = tryparse.(Int, strs[5:i])
    if nothing in sites
        sites = parse.(Float64, strs[5:i])
    end
    text = replace(strs[2], " " => "")
    strs = split(text, ':')
    @assert length(strs) ≤ 2
    name = strs[1]
    state = length(strs) == 1 ? 'S' : strs[2][1]
    return Phase(name, state[1], model, sites)
end

function parse_constituent(text::AbstractString)
    strs = split(text)
    @assert match(r"CONSTITUENT"i, strs[1]) ≠ ""
    phas = split(strs[2], ':')[1] # "Liquid", "Bcc", "Fcc", etc.
    cons_text = join(strs[3:end]) # ":Cu,Zn", ":Cu,Zn:Cu,Zn:", etc.
    cons_text = format_constitution(cons_text)
    cons_strs = split(cons_text, ':')
    cons_strs = filter(e -> e ≠ "", cons_strs)
    cons = map(text -> split(text, ','), cons_strs)
    return phas, cons
end

function parse_parameter(text::AbstractString)
    # Parameters
    strs = split(text)
    @assert match(r"PAR"i, strs[1]) ≠ ""

    text = join(strs[2:end], ' ') # ex. G(BCC2,...

    k = findfirst('(', text)
    l = findfirst(')', text)
    name = replace(text[1:l], " " => "") # G(BCC_B2,Cu:Cu,Zn;0)
    name = format_constitution(name)
    funcname = getfuncname(name) # GBCC_B2CuCuZn0
    symbol = split(name, '(')[1] # G

    # Skip MatCalc specific parameters
    symbol == "HMVA" && return nothing, nothing
    symbol == "SE" && return nothing, nothing

    symbol = symbol[1]

    strs = split(text[l+1:end]) # 298.15...

    text = text[k+1:l-1] # BCC_B2,Cu,Zn:Cu,Zn;0
    text = format_constitution(text)

    k = findfirst(',', text)
    phas = replace(text[1:k-1], " " => "") # BCC_B2

    text = text[k+1:end] # Cu:Cu,Zn;0
    comb_text, order_text = split(text, ';') # Cu:Cu,Zn, 0
    comb = parse_combination(comb_text)
    order = parse(Int, order_text)

    temp, funcstr = parse_function(funcname, strs)

    return phas, Parameter(name, symbol, comb, order, temp, funcstr)
end

function parse_combination(str::AbstractString)
    comb_strs = split(str, ':') # ["Cu", "Cu,Zn"]
    comb = map(text -> split(text, ','), comb_strs) # [[Cu], [Cu,Zn]]
end

function getfuncname(str::AbstractString)
    filter(c -> !in(c, ['(', ')', ',', ':', ';']), str)
end

function localname(param::Parameter)
    str = ""
    S = length(param.comb)
    for s in 1:S
    latt = param.comb[s]
        for cons in latt
            str *= subscript(cons)
            if cons ≠ latt[end]
                str *= ','
            end
        end
        if s ≠ S
            str *= ':'
        end
    end
    str *= ';' * string(param.order)
end

function format_constitution(str::AbstractString)
    # ex. str = "G(FCC,Cu:Cu,Zn,P;0)"
    replace(str, r"(?<=,|;|:).+?(?=,|;|:)" => format_specie)
end

function format_specie(str::AbstractString)
    # ex. str = "Cu2S1", "CU2O"
    str = replace(str, "%" => "")
    str = replace(str, r"(?<=\D)1(?=\D|$)" => "")
    replace(str, r"\D{1,2}(?=\d{1,}|$)" => format_element)
end

function format_element(str::AbstractString)
    N = length(str)
    N > 2 && error("Invalid element name: $str")
    N == 2 && return uppercase(str[1]) * lowercase(str[2])
    N == 1 && return uppercase(str)
end

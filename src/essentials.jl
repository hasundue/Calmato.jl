function subscript(c::Char)
    i = tryparse(Int, string(c))
    if isnothing(i)
        return c
    else
        return Char(0x2080 + parse(Int, c))
    end
end

function subscript(s::AbstractString)
    return map(subscript, s)
end

function subscript(i::Int)
    str = string(i)
    return map(subscript, str)
end

function subscript(x::Float64)
    str = string(x)
    sub = ""
    for c in str
        if c == '.'
            sub *= c
        else
            sub *= subscript(c)
        end
    end
    return sub
end

function composedof(spec::AbstractString)
    rms = collect(eachmatch(r"\D{1,2}(?=\d{1,}|$)", spec))
    isempty(rms) && @error "Unrecognized specie name"
    return map(rm -> rm.match, rms)
end

function iscomposedof(spec::AbstractString, elems::Vector{<:AbstractString})
    comps = composedof(spec)
    return all(elem -> elem in elems, comps)
end

function stoichiometry(spec::AbstractString)
    rms = collect(eachmatch(r"(?=\D{1,2}(\d{1,}|$))\D{1,2}", spec))
    isempty(rms) && @error "Unrecognized specie name"
    elems = map(rm -> rm.match, rms)
    stois = map(rm -> rm.captures, rms)
    dict = Dict{AbstractString,Int}()
    N = length(elems)
    @assert length(stois) == N
    for i in 1:N
        @assert length(stois[i]) == 1
        stoi = tryparse(Int, stois[i][1])
        stoi = isnothing(stoi) ? 1 : stoi
        push!(dict, elems[i] => stoi)
    end
    return dict
end

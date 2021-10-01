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

const Sublattice = Vector{AbstractString}
const Constitution = Vector{Sublattice}

mutable struct Parameter <: AbstractGFunction
    name::AbstractString
    symbol::Char
    comb::Constitution
    order::Int
    temp::Tuple{<:Real,<:Real}
    funcstr::AbstractString
end

mutable struct Phase
    name::AbstractString
    state::Char
    type::AbstractString
    sites::Vector{<:Real}
    cons::Constitution
    params::Vector{Parameter}

    function Phase(name::AbstractString,
                   state::Char,
                   type::AbstractString,
                   sites::Vector{<:Real},
                   cons::Constitution,
                   params::Vector{Parameter})
        return new(name, state, type, sites, cons, params)
    end

    function Phase(name::AbstractString,
                   state::Char,
                   type::AbstractString,
                   sites::Vector{<:Real})
        phas = new()
        phas.name = name
        phas.state = state
        phas.type = type
        phas.sites = sites
        phas.params = Parameter[]
        return phas
    end
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
        println(constitutionstring(phas))
        for param in phas.params
            println("\t\t\t" * param.symbol * '(' * localname(param) * ')')
        end
    end
end

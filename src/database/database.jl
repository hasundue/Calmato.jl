struct Element
    name::AbstractString
    refstate::AbstractString # reference state
    mass::Float64
    H298::Float64
    S298::Float64
end

abstract type AbstractGFunction end

mutable struct GFunction <: AbstractGFunction
    name::AbstractString
    temp::Tuple{<:Real,<:Real}
    temps::Vector{<:Real}
    exprs::Vector{AbstractString}
    funcname::AbstractString
    funcstr::AbstractString
end

function GFunction(name::AbstractString,
                   temp::Tuple{<:Real,<:Real},
                   funcstr::AbstractString)
    return GFunction(name, temp, Real[], AbstractString[], getfuncname(name), funcstr)
end

function GFunction(name::AbstractString,
                   temps::Vector{<:Real},
                   exprs::Vector{<:AbstractString})
    temp = (minimum(temps), maximum(temps))
    return GFunction(name, temp, temps, exprs, getfuncname(name), "")
end

const Sublattice = Vector{AbstractString}
const Constitution = Vector{Sublattice}

mutable struct Parameter <: AbstractGFunction
    name::AbstractString
    symbol::Char
    comb::Constitution
    order::Int
    temp::Tuple{<:Real,<:Real}
    temps::Vector{<:Real}
    exprs::Vector{AbstractString}
    funcname::AbstractString
    funcstr::AbstractString
end

function Parameter(name::AbstractString,
                   symbol::Char,
                   comb::Vector{<:Vector},
                   order::Int,
                   temps::Vector{<:Real},
                   exprs::Vector{<:AbstractString},
                   funcname::AbstractString)
    temp = (minimum(temps), maximum(temps))
    return Parameter(name, symbol, comb, order, temp, temps, exprs, funcname, "")
end

function Parameter(name::AbstractString,
                   symbol::Char,
                   comb::Vector{<:Vector},
                   order::Int,
                   temp::Tuple{<:Real,<:Real},
                   funcstr::AbstractString)
    return Parameter(name, symbol, comb, order, temp,
                     Real[], AbstractString[], getfuncname(name), funcstr)
end

mutable struct Phase
    name::AbstractString
    state::Char
    type::AbstractString
    sites::Vector{<:Real}
    cons::Constitution
    params::Vector{Parameter}
end

function Phase(name::AbstractString,
               state::Char,
               sites::Vector{<:Real})
    return Phase(name, state, "%", sites, Sublattice[], Parameter[])
end

function Phase(name::AbstractString,
               state::Char,
               sites::Vector{<:Real},
               cons::Vector{<:Vector},
               params::Vector{Parameter})
    return Phase(name, state, "%", sites, cons, params)
end

function Phase(name::AbstractString,
               state::Char,
               type::AbstractString,
               sites::Vector{<:Real})
    phas = Phase(name, state, type, sites, Sublattice[], Parameter[])
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

function Database(elems::Vector{Element},
                  funcs::Vector{GFunction},
                  phass::Vector{Phase})
    return Database(elems, funcs, phass, TypeDefinition[])
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

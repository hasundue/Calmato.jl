module Calmato

using Printf

using JuMP
using McCormick; xlogx
using RecipesBase
import EAGO
import ForwardDiff
import Ipopt
import GLPK
import MathOptInterface
import HTTP
import JSON
import PeriodicTable
import ZipFile

MOI = MathOptInterface

export read_tdb, search_db, select
export init_system
export equilib
export solidify

const R = 8.3144598 # gas constant
const P = 1.0 # pressure

function McCormick.xlogx(d::T) where T <: ForwardDiff.Dual
    return d * log(d)
end

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

include("database.jl")
include("tdbdb.jl")
include("system.jl")
include("equilib.jl")
include("solidify.jl")

end

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

# function Base.Float64(x::McCormick.MC{1, McCormick.NS})
#     error(string(x))
# end

include("database.jl")
include("tdbdb.jl")
include("system.jl")
include("equilib.jl")
include("solidify.jl")

end

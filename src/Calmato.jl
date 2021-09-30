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

export read_tdb, search_db, select, merge, merge!
export init_system
export equilib
export solidify

const R = 8.3144598 # gas constant
const P = 1.0 # pressure

function McCormick.xlogx(d::T) where T <: ForwardDiff.Dual
    return d * log(d)
end

include("essentials.jl")
include(joinpath("database", "database.jl"))
include(joinpath("database", "utils.jl"))
include(joinpath("database", "tdb.jl"))
include(joinpath("database", "operations.jl"))
include(joinpath("database", "tdbdb.jl"))
include("system.jl")
include("equilib.jl")
include("solidify.jl")

end

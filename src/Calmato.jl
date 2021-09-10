module Calmato

using Printf

using JuMP
using RecipesBase
import EAGO
import McCormick
import ForwardDiff
import Ipopt
import GLPK
import MathOptInterface as MOI
import HTTP
import JSON
import PeriodicTable
import ZipFile

export read_tdb, search_db
export init_system
export equilib
export solidify

const R = 8.3144598 # gas constant
const P = 1.0 # pressure

function McCormick.xlogx(d::T) where T <: ForwardDiff.Dual
    return d * log(d)
end

include("database.jl")
include("tdbdb.jl")
include("system.jl")
include("equilib.jl")
include("solidify.jl")

end

module Calmato

using Printf

using JuMP
import EAGO
import McCormick
import ForwardDiff
import Ipopt
import GLPK
import MathOptInterface as MOI

export read_tdb
export init_system
export equilib

const R = 8.3144598 # gas constant
const P = 1.0 # pressure

function McCormick.xlogx(d::T) where T <: ForwardDiff.Dual
    return d * log(d)
end

include("database.jl")
include("system.jl")
include("equilib.jl")

end

module Calmato

using Printf

using JuMP
import EAGO
import GLPK
import MathOptInterface as MOI

export read_tdb
export init_system
export equilib

const R = 8.3144598 # gas constant
const P = 1.0 # pressure

include("database.jl")
include("system.jl")
include("equilib.jl")

end

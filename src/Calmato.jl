module Calmato

using JuMP
import EAGO

export read_tdb
export init_system
export equilib

const R = 8.314 # gas constant

include("database.jl")
include("system.jl")
include("equilib.jl")

end

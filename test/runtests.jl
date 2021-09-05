using Calmato
using Test

@testset "Calmato.jl" begin
    db = read_tdb(joinpath("tdb", "cuzn_liang.tdb"))
    sys = init_system(db)
    res = equilib(sys)
end

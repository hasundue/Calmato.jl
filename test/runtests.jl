using Calmato
using Test

@testset "Calmato.jl" begin
    @test_nowarn db = read_tdb(joinpath("test", "tdb", "cuzn_liang.tdb"))
    @test_nowarn print(db)
end

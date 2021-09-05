using Calmato
using Test

@testset "Calmato.jl" begin
    @testset "Database" begin
        db = read_tdb(joinpath("tdb", "cuzn_liang.tdb"))
        print(db)
    end
    @testset "Equilib" begin
        sys = init_system(db)
        res = equilib(sys)
    end
end

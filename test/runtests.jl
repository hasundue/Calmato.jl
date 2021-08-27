using Calmato
using Test

@testset "Calmato.jl" begin
    db = read_tdb(joinpath("tdb", "cuzn_liang.tdb"))
    print(db)
    print_func(db)
end

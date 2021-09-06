using Calmato
using Test

@testset "Calmato.jl" begin
    db = read_tdb(joinpath("tdb", "cuzn_liang.tdb"))
    print(db)
    sys = init_system(db)
    print(sys)
    res = solidify(sys)
    print(res)
end

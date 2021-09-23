using Calmato
using Test

@testset "Calmato.jl" begin
    db = search_db("Cu Zn Liang Hsiao 2015")
    print(db)
    sys = init_system(db)
    print(sys)
    res = solidify(sys)
    print(res)

    db = search_db("Cu H O S P Magnusson")
    print(db)
    db = select(db, "Cu S")
    print(db)
    sys = init_system(db)
    print(sys)
    res = equilib(sys, [0.99, 0.01], 900)
    print(res)
end

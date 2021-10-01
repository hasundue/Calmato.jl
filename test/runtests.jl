using Calmato
using Test

@testset "Calmato.jl" begin
    db1 = search_db("Cu Zn Liang Hsiao 2015")
    print(db1)
    sys = init_system(db1)
    print(sys)
    res = solidify(sys)
    print(res)

    db2 = search_db("Cu H O S P Magnusson")
    db2 = select(db2, "Cu S", *, -1)
    print(db2)

    db3 = merge(db1, db2, 1 => 1, 2 => 4, *)
    print(db3)
    sys = init_system(db3)
    print(sys)
    X = [0.7, 0.3, 0.01]
    Ts = 600:100:1300
    res = solidify(sys, X, Ts)
    print(res)
end

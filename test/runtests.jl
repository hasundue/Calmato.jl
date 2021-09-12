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
    db = select(db, "Cu S P")
    print(db)
end

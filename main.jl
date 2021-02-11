include("posets.jl")
include("iposets.jl")

#ps = genPosets(6)
#savePosets(ps, "posets/6nodes")

#g = SimpleDiGraph(4)
#add_edge!(g, 2, 1)
#add_edge!(g, 3, 1)
#add_edge!(g, 3, 4)
#p = Iposet((2,3), (4,), g)
#
#s = SimpleDiGraph(4)
#add_edge!(s, 1, 2)
#add_edge!(s, 2, 3)
#add_edge!(s, 2, 4)
#q = Iposet((1,), (4,3), s)
#
#draw(PNG("/tmp/Iposet1.png", 16cm, 16cm), igplot(p))
#draw(PNG("/tmp/Iposet2.png", 16cm, 16cm), igplot(q))
#draw(PNG("/tmp/Iposet3.png", 16cm, 16cm), igplot(parallel(p, q)))

println(length(genGpIposets(3)))
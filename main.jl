include("posets.jl")
include("iposets.jl")
using StatProfilerHTML

#ps = genPosets(6)
#savePosets(ps, "posets/6nodes")

#g = SimpleDiGraph(5)
#add_edge!(g, 1, 2)
#add_edge!(g, 1, 3)
#add_edge!(g, 3, 4)
#add_edge!(g, 4, 5)
#p = Iposet((2,), (1,), g)
#
#s = SimpleDiGraph(4)
#add_edge!(s, 1, 2)
#add_edge!(s, 1, 3)
#add_edge!(s, 3, 4)
#add_edge!(s, 2, 4)
#q = Iposet((1,), (4,), s)

#draw(PNG("/tmp/Iposet1.png", 16cm, 16cm), igplot(p))
@time genPosets(8)
#println(length(genGpIposets(6)))

#isIsoIposet(p, q) ? println(true) : println(false)
#genPosets(2)
#@profilehtml genGpIposets(4)
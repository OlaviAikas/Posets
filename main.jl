include("posets.jl")
include("iposets.jl")
using StatProfilerHTML
using Traceur

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
#add_edge!(s, 3, 4)
#q = Iposet((1,), (4,), s)
#println(almostConnected(q))

@time genGpIposets(2)
@time println(length(genGpIposets(6)))

#draw(PNG("/tmp/Iposet1.png", 16cm, 16cm), igplot(p))
#res_gp = zeros(Int, 6, 7, 7)
#res_weak = zeros(Int, 6, 7, 7)
#for n in 1:6
#    a = genGpIposets(n)
#    for iposet in a
#        res_gp[n, length(iposet.s) + 1, length(iposet.t) + 1] += 1
#        if almostConnected(iposet) && ne(hasseDiagram(iposet.poset)) > 0
#            res_weak[n, length(iposet.s) + 1, length(iposet.t) + 1] += 1
#        end
#    end
#end
#
#open("gpiposets.txt", "w") do file
#    for n in 1:6
#        for k in 0:n
#            for j in 0:n
#                r = res_gp[n, k+1, j+1]
#                write(file, "$n, $k, $j -> $r\n")
#            end
#        end
#    end
#end
#
#open("almost_connected_no_discrete.txt", "w") do file
#    for n in 1:6
#        for k in 0:n
#            for j in 0:n
#                r = res_weak[n, k+1, j+1]
#                write(file, "$n, $k, $j -> $r\n")
#            end
#        end
#    end
#end
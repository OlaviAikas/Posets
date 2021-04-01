include("posets.jl")
include("iposets.jl")
include("uli.jl")
using StatProfilerHTML

#ps = genPosets(6)
#savePosets(ps, "posets/6nodes")

#g = SimpleDiGraph(12)
#add_edge!(g, 1, 2)
#add_edge!(g, 1, 3)
#add_edge!(g, 3, 4)
#add_edge!(g, 10, 11)
#p = Iposet((1,2,3), (4,5,12), g)
#println(toString(p))
#
#printPoset(p.poset, "/tmp/p.png")
#printPoset(inversion(p).poset, "/tmp/p_inv.png")
#
#s = SimpleDiGraph(4)
#add_edge!(s, 1, 2)
#add_edge!(s, 3, 4)
#q = Iposet((1,), (4,), s)
#println(glue(p, q))

@time gpPosets(2)
@time println(length(gpPosets(7)))

#@time gpiPosets(2)
#@time println(length(gpiPosets(7)))
#alliposets = Array{Array{Tuple{Iposet, Array{Tuple{Int, Int}}}}}(undef, (6*(6-1))÷2 + 1, 6 + 1, 6 + 1, 6)
#for i in eachindex(alliposets)
#    alliposets[i] = []
#end
#for s in [(), (1,)]
#    for t in [(), (1,)]
#        push!(alliposets[1, length(s) + 1, length(t) + 1, 1], (Iposet(s, t, SimpleDiGraph(1)), [(0, 0)]))
#    end
#end
#filled = zeros(Bool, 7, 7, 6)
#filled[1:2, 1:2, 1] = [true true; true true]
#a = gpiPosets(6, 4, 0, alliposets, filled)
#g = SimpleDiGraph(6)
#add_edge!(g, 1, 3)
#add_edge!(g, 2, 3)
#add_edge!(g, 4, 6)
#add_edge!(g, 5, 6)
#ip = Iposet((1, 4, 2, 5), (), g)
#itransitiveclosure!(ip, false)
#seen = false
#for iposet in a
#    if isIsoIposet(iposet, ip)
#        global seen = true
#        break
#    end
#end
#println(false)

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

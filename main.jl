include("posets.jl")
include("iposets.jl")
using StatProfilerHTML

#ps = genPosets(6)
#savePosets(ps, "posets/6nodes")

#g = SimpleDiGraph(4)
#add_edge!(g, 1, 2)
#add_edge!(g, 1, 3)
#add_edge!(g, 3, 4)
#p = Iposet((2,), (1,), g)
#
#s = SimpleDiGraph(4)
#add_edge!(s, 1, 2)
#add_edge!(s, 2, 3)
#add_edge!(s, 3, 4)
#q = Iposet((1,), (4,), s)

#println(minFiltration(g))
#println(maxFiltration(g))

#draw(PNG("/tmp/Iposet1.png", 16cm, 16cm), igplot(p))
#draw(PNG("/tmp/Iposet2.png", 16cm, 16cm), igplot(q))
#draw(PNG("/tmp/Iposet3.png", 16cm, 16cm), igplot(parallel(p, q)))

#isIsoIposet(p, q) ? println(true) : println(false)
genGpIposets(3)
#@profilehtml genGpIposets(5)
@time(genGpIposets(5))
#74 gp iposets to check
#uli_posets = Array{Iposet}(undef, 74)
#
#for i in 1:74
#    open("../UliPosets/3gpiposets/iposet" * string(i) * ".ipo", "r") do file
#        ipo = readlines(file)
#        s = Array{Int}(undef, 0)
#        t = Array{Int}(undef, 0)
#        poset = SimpleDiGraph(3)
#        stage = 1
#        for str in ipo
#            if str == "s"
#                stage = 1
#                continue
#            elseif str == "t"
#                stage = 2
#                continue
#            elseif str == "edges"
#                stage = 3
#                continue
#            end
#            if stage == 1
#                push!(s, parse(Int, str))
#            elseif stage == 2
#                push!(t, parse(Int, str))
#            elseif stage == 3
#                add_edge!(poset, parse(Int, str[1]), parse(Int, str[4]))
#            end
#        end
#        uli_posets[i] = Iposet(Tuple(s), Tuple(t), poset)
#    end
#end
#
#posets = genGpIposets(3)
#outliers = Array{Iposet}(undef, 0)
#
#for iposet in uli_posets
#    seen = false
#    for iq in posets
#        if isIso(iposet, iq)
#            println("Was iso")
#            seen = true
#            break
#        end
#    end
#    if !seen
#        push!(outliers, iposet)
#    end
#end
#
#
#comp1 = Iposet((1,), (1,), transitiveclosure(SimpleDiGraph(1), true))
#comp2 = Iposet((1,), (1,), transitiveclosure(SimpleDiGraph(1), true))
#comp3 = Iposet((1,), (), transitiveclosure(SimpleDiGraph(1), true))
#comp4 = Iposet((), (1,), transitiveclosure(SimpleDiGraph(1), true))
#lp = parallel(comp1, comp3)
#rp = parallel(comp4, comp2)
#pcomp = glue(lp, rp)
#
#pcomp = parallel(comp2, comp1)
#
#println(outliers[1].s)
#println(outliers[1].t)
#println(pcomp.s)
#println(pcomp.t)
#
#println(collect(edges(outliers[1].poset)))
#println(collect(edges(pcomp.poset)))
#
#
#draw(PNG("/tmp/Pcomp.png", 16cm, 16cm), igplot(iHasse(pcomp)))
#
#isIso(outliers[1], pcomp) ? println(true) : println(false)


#for i in 1:74
#    draw(PNG("/tmp/UliIposet" * string(i) * ".png", 16cm, 16cm), igplot(iHasse(uli_posets[i])))
#    if i < 74
#        draw(PNG("/tmp/Iposet" * string(i) * ".png", 16cm, 16cm), igplot(iHasse(posets[i])))
#    end
#end


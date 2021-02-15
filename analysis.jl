include("posets.jl")

posets7 = loadPosets("posets/7nodes")
gpposets7 = loadPosets("gpposets/7nodes")

println(length(posets7))
println(length(gpposets7))

exclusion = Array{SimpleDiGraph}(undef, 0)
for g in posets7
    seen = false
    for s in gpposets7
        if isIsoV3(transitiveclosure(g, true), transitiveclosure(s, true))
            seen = true
            break
        end
    end
    if !seen
        push!(exclusion, g)
    end
end
println(length(exclusion))

graph_nn = SimpleDiGraph(6)
graph_m = SimpleDiGraph(6)
graph_w = SimpleDiGraph(6)
graph_3c = SimpleDiGraph(6)
graph_ln = SimpleDiGraph(6)

add_edge!(graph_nn, 1, 4)
add_edge!(graph_nn, 2, 5)
add_edge!(graph_nn, 3, 6)
add_edge!(graph_nn, 2, 4)
add_edge!(graph_nn, 3, 5)

add_edge!(graph_m, 1, 3)
add_edge!(graph_m, 1, 4)
add_edge!(graph_m, 2, 4)
add_edge!(graph_m, 2, 5)
add_edge!(graph_m, 3, 6)

add_edge!(graph_w, 1, 4)
add_edge!(graph_w, 2, 5)
add_edge!(graph_w, 3, 5)
add_edge!(graph_w, 3, 6)
add_edge!(graph_w, 4, 6)

add_edge!(graph_3c, 1, 4)
add_edge!(graph_3c, 2, 5)
add_edge!(graph_3c, 3, 6)
add_edge!(graph_3c, 2, 4)
add_edge!(graph_3c, 3, 5)
add_edge!(graph_3c, 1, 6)

add_edge!(graph_ln, 1, 3)
add_edge!(graph_ln, 3, 5)
add_edge!(graph_ln, 2, 4)
add_edge!(graph_ln, 4, 6)
add_edge!(graph_ln, 2, 5)

forbiddens = [graph_nn, graph_w, graph_m, graph_3c, graph_ln]

novelties = Array{SimpleDiGraph}(undef, 0)
for g in exclusion
    seen = false
    for f in forbiddens
        if hasSubgraph(g, f)
            seen = true
            break
        end
    end
    if !seen
        push!(novelties, g)
    end
end
println(length(novelties))

for i in 1:length(novelties)
    draw(PNG("/tmp/novelty" * string(i) * ".png", 16cm, 16cm), gplot(hasseDiagram(novelties[i])))
end
include("posets.jl")

posets7 = loadPosets("posets/7nodes")
gpposets7 = loadPosets("gpposets/7nodes")

println(length(posets7))
println(length(gpposets7))

exclusion = Array{SimpleDiGraph}(undef, 0)
for g in posets7
    transitiveclosure!(g, true)
end
for g in gpposets7
    transitiveclosure!(g, true)
end
for g in posets7
    seen = false
    for s in gpposets7
        if isIso(g, s)
            seen = true
            break
        end
    end
    if !seen
        push!(exclusion, g)
    end
end

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
transitiveclosure!(graph_nn, true)

add_edge!(graph_m, 1, 3)
add_edge!(graph_m, 1, 4)
add_edge!(graph_m, 2, 4)
add_edge!(graph_m, 2, 5)
add_edge!(graph_m, 3, 6)
transitiveclosure!(graph_m, true)


add_edge!(graph_w, 1, 4)
add_edge!(graph_w, 2, 5)
add_edge!(graph_w, 3, 5)
add_edge!(graph_w, 3, 6)
add_edge!(graph_w, 4, 6)
transitiveclosure!(graph_w, true)


add_edge!(graph_3c, 1, 4)
add_edge!(graph_3c, 2, 5)
add_edge!(graph_3c, 3, 6)
add_edge!(graph_3c, 2, 4)
add_edge!(graph_3c, 3, 5)
add_edge!(graph_3c, 1, 6)
transitiveclosure!(graph_3c, true)


add_edge!(graph_ln, 1, 3)
add_edge!(graph_ln, 3, 5)
add_edge!(graph_ln, 2, 4)
add_edge!(graph_ln, 4, 6)
add_edge!(graph_ln, 2, 5)
transitiveclosure!(graph_ln, true)

forbiddens = [graph_nn, graph_w, graph_m, graph_3c, graph_ln]
flabels = ["NN", "W", "M", "3C", "LN"]
sgraphs = Array{Tuple{Int, String}}(undef, length(exclusion))

for i in 1:length(exclusion)
    for f in 1:length(forbiddens)
        sg = hasSubgraph(exclusion[i], forbiddens[f])
        if sg[1]
            rnode = 0
            for v in 1:nv(exclusion[i])
                if !(v in sg[2])
                    rnode = v
                    break
                end
            end
            sgraphs[i] = (rnode, flabels[f])
            break
        end
    end
end

println(sgraphs)

open("NonGPgraphs.gs", "w") do file
    for i in 1:length(exclusion)
        write(file, "Graph " * string(i) * "\n")
        edge_strings = []
        for e in edges(exclusion[i])
            push!(edge_strings, string(src(e)) * ", " * string(dst(e)))
        end
        for s in edge_strings
            write(file, s * "\n")
        end
        write(file, "rem " * string(sgraphs[i][1]) * " for " * sgraphs[i][2] * "\n")
        write(file, "end\n\n")
    end
end

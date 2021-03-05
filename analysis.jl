include("posets.jl")

function readPosetList()
    txt::String = ""
    sn = 0
    dn = 0
    res = Array{SimpleDiGraph}(undef, 0)
    open("8posetlist.txt", "r") do file
        txt = read(file, String)
    end
    for p in txt[8:end-2]
        if p == '{'
            push!(res, SimpleDiGraph(8))
        elseif  p == ']'
            add_edge!(res[end], sn, dn)
            sn = 0
            dn = 0
        elseif isdigit(p)
            if sn == 0
                sn = parse(Int, p)
            else
                dn = parse(Int, p)
            end
        end
    end
    return res
end


"""Read a file of posets in McKay's format and return them in an array"""
function readPosetsMcKay(posetlist)
    txt::Array{String} = []
    open(posetlist, "r") do file
        txt = readlines(file)
    end
    res = Array{SimpleDiGraph}(undef, 0)
    for line in txt
        nums = split(line)
        push!(res, SimpleDiGraph(parse(Int, nums[1])))
        for edges_num in nums[3:end]
            add_edge!(res[end], parse(Int, edges_num[1]) + 1, parse(Int, edges_num[2]) + 1)
        end
    end
    return res
end

posets9 = readPosetsMcKay("hasse9.txt")
gpposets9 = loadPosets("gpposets/9nodes")

println(length(posets9))
println(length(gpposets9))

exclusion = Array{SimpleDiGraph}(undef, 0)
for g in posets9
    transitiveclosure!(g, true)
end
for g in gpposets9
    transitiveclosure!(g, true)
end
for g in posets9
    seen = false
    for s in gpposets9
        if isIso(g, s)
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
graph_bf = SimpleDiGraph(8)

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

add_edge!(graph_bf, 1, 2)
add_edge!(graph_bf, 2, 3)
add_edge!(graph_bf, 4, 5)
add_edge!(graph_bf, 5, 6)
add_edge!(graph_bf, 4, 7)
add_edge!(graph_bf, 7, 3)
add_edge!(graph_bf, 1, 8)
add_edge!(graph_bf, 8, 6)
transitiveclosure!(graph_bf, true)

forbiddens = [graph_nn, graph_w, graph_m, graph_3c, graph_ln, graph_bf]
flabels = ["NN", "W", "M", "3C", "LN", "BF"]
new_forbiddens = Array{SimpleDiGraph}(undef, 0)

for i in 1:length(exclusion)
    seen = false
    for f in 1:length(forbiddens)
        sg = hasSubgraph(exclusion[i], forbiddens[f])
        if sg[1]
            seen = true
            break
        end
    end
    if !seen
        push!(new_forbiddens, exclusion[i])
    end
end

println(length(new_forbiddens))
if length(new_forbiddens) > 0
    for i in 1:length(new_forbiddens)
        printPoset(new_forbiddens[i], "NewForbiddens/newForbidden$(i)_9nodes.png")
    end
end

#open("NonGPgraphs.gs", "w") do file
#    for i in 1:length(exclusion)
#        write(file, "Graph " * string(i) * "\n")
#        edge_strings = []
#        for e in edges(exclusion[i])
#            push!(edge_strings, string(src(e)) * ", " * string(dst(e)))
#        end
#        for s in edge_strings
#            write(file, s * "\n")
#        end
#        write(file, "rem " * string(sgraphs[i][1]) * " for " * sgraphs[i][2] * "\n")
#        write(file, "end\n\n")
#    end
#end

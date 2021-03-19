import Cairo, Fontconfig
using GraphPlot, Compose
using LightGraphs
using Colors
using Printf
using Combinatorics
import Base.==

"""Check if the graph g is antisymmetric, assuming that it is already acyclic."""
function isAntiSymmetric(g::SimpleDiGraph)
    for e in edges(g)
        if e != reverse(e) && has_edge(g, reverse(e))
            return false
        end
    end
    return true
end

"""Compute the minimum filtration (e.g starting from minimal elements) of a
poset. Returns an array where each index corresponding to a node in the graph
is the label in the filtration"""
function minFiltration(g::SimpleDiGraph)
    g = transitiveclosure(g, true)
    res = -1*ones(Int, nv(g))
    deleted = zeros(Bool, nv(g))
    level = 0
    while -1 in res
        just_deleted = zeros(Bool, nv(g))
        for v in vertices(g)
            if deleted[v]
                continue
            end
            minimal = true
            for neighbour in inneighbors(g, v)
                if neighbour != v && (!deleted[neighbour] || just_deleted[neighbour])
                    minimal = false
                end
            end
            if minimal
                res[v] = level
                deleted[v] = true
                just_deleted[v] = true
            end
        end
        level += 1
    end
    return res
end

"""Compute the maximum filtration (e.g starting from maximal elements) of a
poset. Returns an array where each index corresponding to a node in the graph
is the label in the filtration"""
function maxFiltration(g::SimpleDiGraph)
    g = transitiveclosure(g, true)
    res = -1*ones(Int, nv(g))
    deleted = zeros(Bool, nv(g))
    level = 0
    while -1 in res
        just_deleted = zeros(Bool, nv(g))
        for v in vertices(g)
            if deleted[v]
                continue
            end
            maximal = true
            for neighbour in outneighbors(g, v)
                if neighbour != v && (!deleted[neighbour] || just_deleted[neighbour])
                    maximal = false
                end
            end
            if maximal
                res[v] = level
                deleted[v] = true
                just_deleted[v] = true
            end
        end
        level += 1
    end
    return res
end

"""Produce a new SimpleDiGraph object, that represents the hasse diagram of g.
This is the minimal graph gh such that transitiveclosure!(gh) = g. Used to draw
nicer diagrams -- Attention! This is buggy and will be soon reworked
with the help of the minFiltration function"""
function hasseDiagram(g::SimpleDiGraph)
    s::SimpleDiGraph = deepcopy(g)
    for i in 1:nv(s)
        rem_edge!(s, i, i)
    end
    redundant_e = Array{LightGraphs.SimpleGraphs.SimpleEdge}(undef, 0)
    for e in edges(s)
        for pos_intermediary in outneighbors(s, src(e))
            if pos_intermediary == dst(e)
                continue
            end
            if dst(e) in outneighbors(s, pos_intermediary)
                push!(redundant_e, e)
                break
            end
        end
    end
    for e in redundant_e
        rem_edge!(s, e)
    end
    return s
end

"""Check if the array a is an isomorphism for the graphs g and s"""
function isIso(g::SimpleDiGraph, s::SimpleDiGraph, a::Array{Int})
    #Easy stuff first
    if nv(g) != nv(s) || length(a) != nv(g)
        return false
    end
    for i in 1:length(a) #Confirm that "a" is an isomorphism
        if !(i in a)
            return false
        end
    end
    for vertex in 1:length(a)
        neighbours_in_g = inneighbors(g, vertex)
        neighbours_in_s = inneighbors(s, a[vertex])
        for neighbour in neighbours_in_g
            if !(a[neighbour] in neighbours_in_s)
                return false
            end
        end
        neighbours_out_g = outneighbors(g, vertex)
        neighbours_out_s = outneighbors(s, a[vertex])
        for neighbour in neighbours_out_g
            if !(a[neighbour] in neighbours_out_s)
                return false
            end
        end
    end
    return true
end

"""Check if the graphs g and s are isomorphic to each other"""
function isIso(g::SimpleDiGraph, s::SimpleDiGraph)
    #Start with the easy stuff
    n = nv(g)
    if n != nv(s) || ne(g) != ne(s)
        return false
    end
    #Make arrays of "graph invariants" e.g integer tuples that describe the nodes
    g_v_profiles = Array{Tuple{Int, Int}}(undef, n)
    s_v_profiles = Array{Tuple{Int, Int}}(undef, n)
    for v in 1:n
        g_v_profiles[v] = (inHash(g, v), outHash(g, v))
        s_v_profiles[v] = (inHash(s, v), outHash(s, v))
    end
    # Check that the vertex profiles match in both graphs
    used = zeros(Bool, n)
    for i in 1:n
        found = false
        for j in 1:n
            if g_v_profiles[i] == s_v_profiles[j] && !used[j]
                found = true
                break
            end
        end
        if !found
            return false
        end
    end
    # Construct classes of vertices that can be isomorphic
    v_classes = Array{Tuple{Array{Int}, Array{Int}}}(undef, 0)
    indexed_g = zeros(Bool, n)
    indexed_s = zeros(Bool, n)
    for i in 1:n
        if indexed_g[i]
            continue
        end
        class = (Array{Int}(undef, 1), Array{Int}(undef, 0))
        class[1][1] = i
        for v in 1:n
            if !indexed_s[v] && s_v_profiles[v] == g_v_profiles[i]
                push!(class[2], v)
                indexed_s[v] = true
            end
            if !indexed_g[v] && g_v_profiles[v] == g_v_profiles[i] && v != i
                push!(class[1], v)
                indexed_g[v] = true
            end
        end
        push!(v_classes, class)
        indexed_g[i] = true
    end
    isom_parts = Array{Array{Array{Tuple{Int, Int}}}}(undef, length(v_classes))
    for i in 1:length(v_classes)
        isom_parts[i] = []
        lh = v_classes[i][1] # Nodes in g
        rh = v_classes[i][2] # Nodes in s
        for perm in permutations(rh)
            push!(isom_parts[i], collect(zip(lh, perm)))
        end
    end
    for perm in Iterators.product(isom_parts...)
        pos_isom = Array{Int}(undef, n)
        mapped = zeros(Bool, n)
        # Now we use the vertex profiles to see if we can fill in pos_isom
        #println(perm)
        for info_array in perm
            if !(false in mapped)
                break
            end
            for mapping in info_array
                pos_isom[mapping[1]] = mapping[2]
            end
        end
        #println(perm)
        #println(pos_isom)
        if isIso(g, s, pos_isom)
            return true
        end
    end
    return false
end

"""Check if the tuple a is an isomorphism for the graphs g and s without checking
for number of vertices or that a is bijective"""
function isIsoX(g::SimpleDiGraph, s::SimpleDiGraph, a::Tuple{Vararg{Int}})
    for vertex in 1:length(a)
        neighbours_in_g = inneighbors(g, vertex)
        neighbours_in_s = inneighbors(s, a[vertex])
        for neighbour in neighbours_in_g
            if !(a[neighbour] in neighbours_in_s)
                return false
            end
        end
        neighbours_out_g = outneighbors(g, vertex)
        neighbours_out_s = outneighbors(s, a[vertex])
        for neighbour in neighbours_out_g
            if !(a[neighbour] in neighbours_out_s)
                return false
            end
        end
    end
    return true
end

"""Same as isIsoIposetX but for normal posets, use with caution"""
function isIsoX(g::SimpleDiGraph, g_v_profiles::Array{Tuple{Int, Int}}, 
                s::SimpleDiGraph, s_v_profiles::Array{Tuple{Int, Int}})
    #Start with the easy stuff
    n = nv(g)
    # Check that the vertex profiles match in both graphs
    used = zeros(Bool, n)
    @inbounds for i in 1:n
        found = false
        for j in 1:n
            if g_v_profiles[i] == s_v_profiles[j] && !used[j]
                found = true
                break
            end
        end
        if !found
            return false
        end
    end
    # Construct an array/mapping node -> list of possible nodes it can map to
    targets = Array{Array{Int}}(undef, n)
    @inbounds for i in 1:n
        targets[i] = []
    end
    @inbounds for v in 1:n
        for i in 1:n
            if g_v_profiles[v] == s_v_profiles[i]
                push!(targets[v], i)
            end
        end
    end
    # Now we see which combinations of those targets can give us bijective
    # mappings, and hope they're graph isomorphisms
    for pos_isom in Iterators.product(targets...)
        #Check bijectivity
        bij = true
        @inbounds for i in 1:n
            if !(i in pos_isom)
                bij = false
                break
            end
        end
        if !bij
            continue
        end
        if isIsoX(g, s, pos_isom)
            return true
        end
    end
    return false
end

"""Generate all of the posets of size n, return them as an array"""
function genPosets(n::Int)
    #Quickly handle special cases 0 and 1 because they're easy
    if n == 0
        return [SimpleDiGraph(0)]
    elseif n == 1
        return [SimpleDiGraph(1)]
    end
    #Movin on...
    max_edges = (n*(n-1))÷2 #÷ is the integer division character
    all_edges = Array{Tuple{Int, Int}}(undef, max_edges)
    en::UInt = 1
    for i in 1:(n-1)
        for j in (i+1):n
            all_edges[en] = (i, j)
            en += 1
        end
    end
    println(all_edges)
    res = Array{Array{Tuple{SimpleDiGraph, Array{Tuple{Int, Int}}}}}(undef, max_edges + 1)
    for i in eachindex(res)
        res[i] = []
    end
    for numedges in 0:max_edges
        #Get all the ways to distribute the edges
        for distribution in combinations(all_edges, numedges)
            g = SimpleDiGraph(n) #New graph with n nodes
            #Now we put the edges in
            for e in distribution
                add_edge!(g, e[1], e[2])
            end
            transitiveclosure!(g, false) #Take the transitive closure
            pne = ne(g)
            vprof = Array{Tuple{Int, Int}}(undef, n)
            @inbounds for v in 1:n
                vprof[v] = (inHash(g, v), outHash(g, v))
            end
            if isAntiSymmetric(g)
                is_iso = false
                for s in res[pne + 1]
                    if isIsoX(g, vprof, s[1], s[2])
                        is_iso = true
                        break
                    end
                end
                if is_iso
                    continue
                end
                push!(res[pne + 1], (g, vprof))
            end
        end
    end
    vc = vcat(res...)
    res = Array{SimpleDiGraph}(undef, length(vc))
    @inbounds for i in 1:length(vc)
        res[i] = vc[i][1]
    end
    return res
end

"""Print a poset into the given filename in a nice way such that the nodes are
organised from left to right according to the minFiltration"""
function printPoset(g::SimpleDiGraph, filename::String)
    filtration = minFiltration(g)
    columns = Array{Array{Int}}(undef, max(filtration...) + 1)
    for i in 1:length(columns)
        columns[i] = []
    end
    for i in 1:length(filtration)
        push!(columns[filtration[i] + 1], i)
    end
    y_vals = Array{Float64}(undef, length(filtration))
    for col in columns
        if length(col) == 1
            y_vals[col[1]] = 0.5
        else
            spacing = 1/(length(col) - 1)
            for i in 1:length(col)
                y_vals[col[i]] = (i - 1)*spacing
            end
        end
    end
    draw(PNG(filename, 16cm, 16cm), gplot(g, float(filtration), y_vals))
end

"""Print an array of Posets into path for inspection. Enter a path with a prefix
of the desired filename, for example /tmp/poset will result in files like
/tmp/poset1.png, /tmp/poset2.png, etc. Prints the Hasse diagram of posets"""
function printPosets(a::Array{SimpleDiGraph}, path::String)
    i = 0
    for g in a
        gh = hasseDiagram(hasseDiagram(g))
        i += 1
        printPoset(gh, path * string(i) * ".png")
    end
end

"""Save an array of posets to disk, 1 file per poset. Takes a list of posets
and a directory to save them to as arguments"""
function savePosets(a::Array{SimpleDiGraph}, dir::String)
    counter = 1
    for g in a
        savegraph(dir * "/poset" * string(counter) * ".lgz", hasseDiagram(g))
        counter += 1
    end
end

"""From a directory of posets on disk, return them all in an array of
SimpleDiGraphs"""
function loadPosets(dir::String)
    res = Array{SimpleDiGraph}(undef, 0)
    counter = 1
    still_files = true
    while still_files
        try
            push!(res, loadgraph(dir * "/poset" * string(counter) * ".lgz"))
            counter += 1
        catch SystemError
            still_files = false
        end
    end
    if length(res) == 0
        println("No graphs loaded, make sure you typed the directory right")
    end
    return res
end

"""Find minimal nodes in a poset"""
function minNodes(p::SimpleDiGraph)
    res = []
    for v in vertices(p)
        n = inneighbors(p, v)
        if length(n) == 0 || n == [v]
            push!(res, v)
        end
    end
    return res
end

"""Find maximal nodes in an iposet"""
function maxNodes(p::SimpleDiGraph)
    res = []
    for v in vertices(p)
        n = outneighbors(p, v)
        if length(n) == 0 || n == [v]
            push!(res, v)
        end
    end
    return res
end

"""Compute an "in-hash" of a vertex in a graph g, which is basically a hashing
on the in-edges of the vertex and its neighbouring vertices, that can be used
as an isomorphism invariant"""
function inHash(g::SimpleDiGraph, v::Int)
    n = inneighbors(g, v)
    if length(n) == 0 || n == [v]
        return 0
    else
        nhashes = 0
        for neighbour in n
            if neighbour != v
                nhashes += inHash(g, neighbour)
            end
        end
        return length(n) + nhashes << 4
    end
end

"""Same as inHash but for outgoing endges"""
function outHash(g::SimpleDiGraph, v::Int)
    n = outneighbors(g, v)
    if length(n) == 0 || n == [v]
        return 0
    else
        nhashes = 0
        for neighbour in n
            if neighbour != v
                nhashes += outHash(g, neighbour)
            end
        end
        return length(n) + 16*nhashes
    end
end

"""Check if some desired graph s is a subgraph (or isomorphic to some subgraph)
of the graph g. Return a tuple (Bool, Array{Int}) to signify success and the
subset of nodes that give the induced subgraph."""
function hasSubgraph(g::SimpleDiGraph, s::SimpleDiGraph)
    subnodes = combinations(1:nv(g), nv(s))
    for nodes in subnodes
        subgraph = induced_subgraph(g, nodes)[1]
        if isIso(subgraph, s)
            return (true, nodes)
        end
    end
    return (false, [])
end

"""Check if some desired graph s is a subgraph (or isomorphic to some subgraph)
of the graph g, and print the graph with the nodes of the subgraph highlighted,
and the subgraph with the filenames filenamef (full) and filenames (sub)"""
function printSubgraph(g::SimpleDiGraph, s::SimpleDiGraph, filenamef::String, filenames::String)
    subnodes = combinations(1:nv(g), nv(s))
    for nodes in subnodes
        subgraph = induced_subgraph(g, nodes)[1]
        if isIso(subgraph, s)
            println("Found subgraph, printing...")
            filtration = minFiltration(g)
            columns = Array{Array{Int}}(undef, max(filtration...) + 1)
            for i in 1:length(columns)
                columns[i] = []
            end
            for i in 1:length(filtration)
                push!(columns[filtration[i] + 1], i)
            end
            y_vals = Array{Float64}(undef, length(filtration))
            for col in columns
                if length(col) == 1
                    y_vals[col[1]] = 0.5
                else
                    spacing = 1/(length(col) - 1)
                    for i in 1:length(col)
                        y_vals[col[i]] = (i - 1)*spacing
                    end
                end
            end
            colors = Array{ColorTypes.AbstractRGB}(undef, nv(g))
            for i in 1:length(colors)
                if i in nodes
                    colors[i] = colorant"red"
                else
                    colors[i] = colorant"blue"
                end
            end
            draw(PNG(filenamef, 16cm, 16cm), gplot(g, float(filtration), y_vals, nodefillc=colors))
            printPoset(s, filenames)
            return true
        end
    end
    println("Subgraph not found")
    return false
end

"""Overload equality for Array{SimpleDiGraph} so you can easily check if two
arrays have the same posets up to isomorphism"""
function ==(a::Array{SimpleDiGraph}, b::Array{SimpleDiGraph})
    for ps1 in a
        seen = false
        for ps2 in b
            if isIso(ps1, ps2)
                seen = true
                break
            end
        end
        if !seen
            return false
        end
    end
    return true
end
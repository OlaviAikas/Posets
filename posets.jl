import Cairo, Fontconfig
using GraphPlot, Compose
using LightGraphs
using Printf
using Combinatorics

"""Check if the graph g is antisymmetric, assuming that it is already acyclic."""
function isAntiSymmetric(g::SimpleDiGraph)
    for e in edges(g)
        if e != reverse(e) && has_edge(g, reverse(e))
            return false
        end
    end
    return true
end

"""Produce a new SimpleDiGraph object, that represents the hasse diagram of g.
This is the minimal graph gh such that transitiveclosure!(gh) = g. Used to draw
nicer diagrams -- Attention! This is buggy and will be soon reworked
with the help of the minFiltration function"""
function hasseDiagram(g::SimpleDiGraph)
    for i in 1:nv(g)
        rem_edge!(g, i, i)
    end
    s::SimpleDiGraph = g
    for v in vertices(g)
        for u in outneighbors(g, v)
            for w in outneighbors(g, v)
                if u != w && has_path(g, u, w)
                    rem_edge!(s, v, w)
                end
            end
        end
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
    #Take the transitive closures in case they weren't good before
    #transitiveclosure!(g, true)
    #transitiveclosure!(s, true)
    #Make arrays of "graph invariants" e.g integer tuples that describe the nodes
    g_v_profiles = Array{Tuple{Int, Int}}(undef, n)
    s_v_profiles = Array{Tuple{Int, Int}}(undef, n)
    for v in 1:n
        g_v_profiles[v] = (length(inneighbors(g, v)), length(outneighbors(g, v)))
        s_v_profiles[v] = (length(inneighbors(s, v)), length(outneighbors(s, v)))
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
    indexed = zeros(Bool, n)
    for i in 1:n
        if indexed[i]
            continue
        end
        class = (Array{Int}(undef, 1), Array{Int}(undef, 0))
        class[1][1] = i
        for v in 1:n
            if s_v_profiles[v] == g_v_profiles[i]
                push!(class[2], v)
            end
            if g_v_profiles[v] == g_v_profiles[i] && v != i
                push!(class[1], v)
                indexed[v] = true
            end
        end
        push!(v_classes, class)
        indexed[i] = true
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
        # Now we use the vertex profiles to see if we can fill in pos_isom
        #println(perm)
        for info_array in perm
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

"""Generate all of the posets of size n, return them as a set"""
function genPosets(n::Int)
    #Quickly handle special cases 0 and 1 because they're easy
    if n == 0
        return SimpleDiGraph(0)
    elseif n == 1
        return SimpleDiGraph(1)
    end
    #Movin on...
    graph_counter = 0
    max_edges = (n*(n-1))÷2 #÷ is the integer division character
    representatives = Array{SimpleDiGraph}(undef, 0)
    all_edges = Array{Tuple{Int, Int}}(undef, max_edges)
    en::UInt = 1
    for i in 1:(n-1)
        for j in (i+1):n
            all_edges[en] = (i, j)
            en += 1
        end
    end
    println(all_edges)
    for numedges in 0:max_edges
        #@printf("Considering %d edges\n", numedges)
        #Get all the ways to distribute the edges
        distributions = collect(combinations(all_edges, numedges))
        ld = length(distributions)
        for i in 1:ld
            #@printf("Checking distribution %d/%d\n", i, ld)
            g = SimpleDiGraph(n) #New graph with n nodes
            #Now we put the edges in
            for e in distributions[i]
                add_edge!(g, e[1], e[2])
            end
            transitiveclosure!(g, true) #Take the transitive closure
            if isAntiSymmetric(g)
                is_iso = false
                for s in representatives
                    if isIso(g, s)
                        is_iso = true
                        break
                    end
                end
                if is_iso
                    continue
                end
                graph_counter += 1
                println("Found a representative, " * string(graph_counter) * " so far")
                push!(representatives, g)
            end
        end
    end
    return representatives
end

"""Print an array of Posets into /tmp/ for inspection"""
function printPosets(a::Array{SimpleDiGraph})
    i::Uint = 0
    for g in a
        gh = hasseDiagram(hasseDiagram(g))
        i += 1
        draw(PNG("/tmp/Hasse" * string(i) * ".png", 16cm, 16cm), gplot(gh))
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

"""Compute the minimum filtration (e.g starting from minimal elements) of a
poset. Returns an array where each index corresponding to a node in the graph
is the label in the filtration"""
function minFiltration(g::SimpleDiGraph)
    transitiveclosure!(g, true)
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
    transitiveclosure!(g, true)
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

"""Compute an "in-hash" of a vertex in a graph g, which is basically a hashing
on the in-edges of the vertex and its neighbouring vertices, that can be used
as an isomorphism invariant"""
function inHash(g::SimpleDiGraph, v::Int)
    res = 0
    for n in inneighbors(g, v)
        res += length(inneighbors(g, n))
    end
    return res*length(inneighbors(g, v))
end

"""Same as inHash but for outgoing endges"""
function outHash(g::SimpleDiGraph, v::Int)
    res = 0
    for n in outneighbors(g, v)
        res += length(outneighbors(g, n))
    end
    return res*length(outneighbors(g, v))
end

"""Check if some desired graph s is a subgraph (or isomorphic to some subgraph)
of the graph g"""
function hasSubgraph(g::SimpleDiGraph, s::SimpleDiGraph)
    subnodes = combinations(1:nv(g), nv(s))
    for nodes in subnodes
        subgraph = induced_subgraph(g, nodes)[1]
        if isIso(subgraph, s)
            return true
        end
    end
    return false
end
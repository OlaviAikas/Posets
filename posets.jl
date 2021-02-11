import Cairo, Fontconfig
using GraphPlot, Compose
using LightGraphs
using Printf
using Combinatorics

"""Check if the graph g is antisymmetric. Assumes that it is already acyclic.
Returns true or false"""
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
nicer diagrams"""
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

"""Check if the mapping d::Dict is an isomorphism from g to s"""
function isIso(g::SimpleDiGraph, s::SimpleDiGraph, d::Dict{Int, Int})
    for i in 1:(length(nv(g))) #Check that the values make a permutation of the nodes (bijection)
        if !(i in values(d)) || !(i in keys(d))
            return false
        end
    end
    if nv(g) != length(keys(d))
        return false
    end
    for vertex in keys(d)
        neighbours_in_g = inneighbors(g, vertex)
        neighbours_in_s = inneighbors(s, d[vertex])
        for neighbour in neighbours_in_g
            if !(d[neighbour] in neighbours_in_s)
                return false
            end
        end
        neighbours_out_g = outneighbors(g, vertex)
        neighbours_out_s = outneighbors(s, d[vertex])
        for neighbour in neighbours_out_g
            if !(d[neighbour] in neighbours_out_s)
                return false
            end
        end
    end
    return true
end

"""Check if graphs g and s are isomorphic to each other. Returns a boolean"""
function isIso(g::SimpleDiGraph, s::SimpleDiGraph)
    #Start with the easy stuff
    n = nv(g)
    if n != nv(s) || ne(g) != ne(s)
        return false
    end
    for permutation in collect(permutations(1:n))
        pos_isom = Dict(zip(1:n, permutation))
        if isIso(g, s, pos_isom)
            return true
        end
    end
    return false
end

"""Better verion of isIso which uses a more complicated, but more efficient
method to iterate over possible isomorphism."""
function isIsoV2(g::SimpleDiGraph, s::SimpleDiGraph)
    #Start with the easy stuff
    n = nv(g)
    if n != nv(s) || ne(g) != ne(s)
        return false
    end
    g_nb_in_per_node = Array{Array{Int,1},1}(undef, n + 1)
    s_nb_in_per_node = Array{Array{Int,1},1}(undef, n + 1)
    g_nb_out_per_node = Array{Array{Int,1},1}(undef, n + 1)
    s_nb_out_per_node = Array{Array{Int,1},1}(undef, n + 1)
    for i in 1:n #Construct arrays where the ith index is an array of all the
                 #nodes that have i in/out neighbours in g and s
                 #Reserve index 1 for 0 neighbours and push everything else to
                 #the right
        gnib = length(inneighbors(g, i))
        if isassigned(g_nb_in_per_node, gnib + 1)
            push!(g_nb_in_per_node[gnib + 1], i)
        else 
            g_nb_in_per_node[gnib + 1] = [i]
        end
        snib = length(inneighbors(s, i))
        if isassigned(s_nb_in_per_node, snib + 1)
            push!(s_nb_in_per_node[snib + 1], i)
        else 
            s_nb_in_per_node[snib + 1] = [i]
        end
        gnob = length(outneighbors(g, i))
        if isassigned(g_nb_out_per_node, gnob + 1)
            push!(g_nb_out_per_node[gnob + 1], i)
        else 
            g_nb_out_per_node[gnob + 1] = [i]
        end
        snob = length(outneighbors(s, i))
        if isassigned(s_nb_out_per_node, snob + 1)
            push!(s_nb_out_per_node[snob + 1], i)
        else 
            s_nb_out_per_node[snob + 1] = [i]
        end
    end
    for i in 1:(n+1)
        if isassigned(g_nb_in_per_node, i) && isassigned(s_nb_in_per_node, i)
            if length(g_nb_in_per_node[i]) != length(s_nb_in_per_node[i])
                return false
            end
        elseif !isassigned(g_nb_in_per_node, i) && !isassigned(s_nb_in_per_node, i)
            #Do nothing
        else
            return false
        end
        if isassigned(g_nb_out_per_node, i) && isassigned(s_nb_out_per_node, i)
            if length(g_nb_out_per_node[i]) == length(s_nb_out_per_node[i])
                continue
            else
                return false
            end
        elseif !isassigned(g_nb_out_per_node, i) && !isassigned(s_nb_out_per_node, i)
            continue
        else
            return false
        end
    end
    in_isom_parts = Array{Array{Dict{Int, Int},1},1}(undef, n + 1)
    out_isom_parts = Array{Array{Dict{Int, Int},1},1}(undef, n + 1)
    for i in 1:(n+1)
        if isassigned(g_nb_in_per_node, i)
            for permutation in collect(permutations(s_nb_in_per_node[i]))
                if isassigned(in_isom_parts, i)
                    push!(in_isom_parts[i], Dict(zip(g_nb_in_per_node[i], permutation)))
                else
                    in_isom_parts[i] = [Dict(zip(g_nb_in_per_node[i], permutation))]
                end
            end
        else
            in_isom_parts[i] = [Dict()]
        end
        if isassigned(g_nb_out_per_node, i)
            for permutation in collect(permutations(s_nb_out_per_node[i]))
                if isassigned(out_isom_parts, i)
                    push!(out_isom_parts[i], Dict(zip(g_nb_out_per_node[i], permutation)))
                else
                    out_isom_parts[i] = [Dict(zip(g_nb_out_per_node[i], permutation))]
                end
            end
        else
            out_isom_parts[i] = [Dict()]
        end
    end
    for permutation in Iterators.product(in_isom_parts..., out_isom_parts...)
        pos_isom = Dict{Int, Int}()
        for isom_part in permutation
            pos_isom = merge(pos_isom, isom_part)
        end
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
                    if isIsoV2(g, s)
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

"""Find minimal nodes in an iposet"""
function minNodes(p::SimpleDiGraph)
    res = []
    for v in vertices(p)
        if length(inneighbors(p, v)) <= 1
            push!(res, v)
        end
    end
    return res
end

"""Find maximal nodes in an iposet"""
function maxNodes(p::SimpleDiGraph)
    res = []
    for v in vertices(p)
        if length(outneighbors(p, v)) <= 1
            push!(res, v)
        end
    end
    return res
end
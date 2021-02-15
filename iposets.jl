using LightGraphs
using Colors
using GraphPlot, Compose
using Combinatorics
import Cairo, Fontconfig
include("posets.jl")

struct GluingInterfacesDontMatchError <: Exception end

"""Check that the nodes in the tuple s are minimal in the poset p"""
function is_min(p::SimpleDiGraph, s::Tuple{Vararg{Int}})
    for v in s
        if length(inneighbors(p, v)) > 1
            return false
        end
    end
    return true
end

"""Check that the nodes in the tuple t are maximal in the poset p"""
function is_max(p::SimpleDiGraph, t::Tuple{Vararg{Int}})
    for v in t
        if length(outneighbors(p, v)) > 1
            return false
        end
    end
    return true
end

"""Define Iposet with tuples for the interfaces and an associated poset"""
struct Iposet
    s::Tuple{Vararg{Int}}
    t::Tuple{Vararg{Int}}
    poset::SimpleDiGraph
    Iposet(s, t, p) = is_min(p, s) && is_max(p, t) ? new(s, t, p) : error("Interfaces not minimal/maximal")
end

"""Hasse diagram of an Iposet"""
function iHasse(g::Iposet)
    return Iposet(g.s, g.t, hasseDiagram(g.poset))
end

"""Get the gplot object for an Iposet with the interfaces coloured
The key is:
green - node in source
red - node in target
yellow - node in both source and target
blue - normal node (not in source or target)"""
function igplot(p::Iposet)
    colors = Array{ColorTypes.AbstractRGB}(undef, nv(p.poset))
    for i in 1:length(colors)
        if i in p.s && i in p.t
            colors[i] = colorant"#C7C795"
        elseif i in p.s
            colors[i] = colorant"#95C7AE"
        elseif i in p.t
            colors[i] = colorant"#C795AE"
        else
            colors[i] = colorant"#95AEC7"
        end
    end
    return gplot(p.poset, nodefillc=colors)
end

"""Return gluing product of two Iposets"""
function glue(g::Iposet, p::Iposet)
    if length(g.t) != length(p.s)
        throw(GluingInterfacesDontMatchError)
    end
    rposet = SimpleDiGraph(nv(g.poset) + nv(p.poset) - length(g.t))
    count::UInt, oldnew1, oldnew2 = 0, Dict{Int, Int}(), Dict{Int, Int}()
    for v in vertices(g.poset)
        count += 1
        oldnew1[v] = count
    end
    for v in vertices(p.poset)
        if v in p.s
            pos = findfirst(x -> x == v ? true : false, p.s)
            corrtarget = g.t[pos]
            oldnew2[v] = oldnew1[corrtarget]
        else
            count += 1
            oldnew2[v] = count
        end
    end
    for e in edges(g.poset)
        !add_edge!(rposet, oldnew1[src(e)], oldnew1[dst(e)])
    end
    for e in edges(p.poset)
        !add_edge!(rposet, oldnew2[src(e)], oldnew2[dst(e)])
    end
    for x in vertices(g.poset)
        if !(x in g.t)
            for y in vertices(p.poset)
                if !(y in p.s)
                    add_edge!(rposet, oldnew1[x], oldnew2[y])
                end
            end
        end
    end
    rs, rt = Array{Int}(undef, length(g.s)), Array{Int}(undef, length(p.t))
    for i in 1:length(g.s)
        rs[i] = oldnew1[g.s[i]]
    end
    for i in 1:length(p.t)
        rt[i] = oldnew2[p.t[i]]
    end
    return Iposet(Tuple(rs), Tuple(rt), rposet)
end

"""Parallel composition of two Iposets"""
function parallel(g::Iposet, s::Iposet)
    gsize = nv(g.poset)
    rposet = SimpleDiGraph(gsize + nv(s.poset))
    for e in edges(g.poset)
        add_edge!(rposet, e)
    end
    for e in edges(s.poset)
        add_edge!(rposet, src(e) + gsize, dst(e) + gsize)
    end
    ns = Array{Int}(undef, length(g.s) + length(s.s))
    nt = Array{Int}(undef, length(g.t) + length(s.t))
    for i in 1:(length(g.s) + length(s.s))
        if i <= length(g.s)
            ns[i] = g.s[i]
        else
            ns[i] = s.s[i - length(g.s)] + gsize
        end
    end
    for i in 1:(length(g.t) + length(s.t))
        if i <= length(g.t)
            nt[i] = g.t[i]
        else
            nt[i] = s.t[i - length(g.t)] + gsize
        end
    end
    return Iposet(Tuple(ns), Tuple(nt), rposet)
end

"""Check if the array a is an isomorphism for the iposets g and s"""
function isIsoIposet(g::Iposet, s::Iposet, a::Array{Int})
    #Easy stuff first
    if nv(g.poset) != nv(s.poset) || length(a) != nv(g.poset)
        return false
    end
    for i in 1:length(a) #Confirm that "a" is an isomorphism
        if !(i in a)
            return false
        end
    end
    # Check that the interfaces are mapped nicely
    for i in 1:length(g.s)
        if a[g.s[i]] != s.s[i]
            return false
        end
    end
    for i in 1:length(g.t)
        if a[g.t[i]] != s.t[i]
            return false
        end
    end
    for vertex in 1:length(a)
        neighbours_in_g = inneighbors(g.poset, vertex)
        neighbours_in_s = inneighbors(s.poset, a[vertex])
        for neighbour in neighbours_in_g
            if !(a[neighbour] in neighbours_in_s)
                return false
            end
        end
        neighbours_out_g = outneighbors(g.poset, vertex)
        neighbours_out_s = outneighbors(s.poset, a[vertex])
        for neighbour in neighbours_out_g
            if !(a[neighbour] in neighbours_out_s)
                return false
            end
        end
    end
    return true
end

"""Check if the iposets g and s are isomorphic to each other"""
function isIsoIposet(g::Iposet, s::Iposet)
    #Start with the easy stuff
    n = nv(g.poset)
    if n != nv(s.poset) || ne(g.poset) != ne(s.poset)
        return false
    end
    if length(g.s) != length(s.s) || length(g.t) != length(s.t)
        return false
    end
    #Make arrays of "graph invariants" e.g integer tuples that describe the nodes
    g_v_profiles = Array{Tuple{Int, Int}}(undef, n)
    s_v_profiles = Array{Tuple{Int, Int}}(undef, n)
    for v in 1:n
        g_v_profiles[v] = (length(inneighbors(g.poset, v)), length(outneighbors(g.poset, v)))
        s_v_profiles[v] = (length(inneighbors(s.poset, v)), length(outneighbors(s.poset, v)))
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
    # Each class is a tuple with 2 arrays such that the nodes in the first array
    # are vertices in g that can potentially be mapped to those vertices of s in
    # the other array
    indexed = Array{Bool}(undef, n)
    for i in 1:n
        if i in g.s || i in g.t
            indexed[i] = true # We already know where the interfaces are mapped
        else                  # to, so they don't need to be put in classes
            indexed[i] = false
        end
    end
    for i in 1:n
        if indexed[i]
            continue
        end
        class = (Array{Int}(undef, 1), Array{Int}(undef, 0))
        class[1][1] = i # Construct a new class with all of the nodes with the
                        # same invariant as "i"
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
    # Here we construct all the possible mappings of the nodes in each v-class
    for i in 1:length(v_classes)
        isom_parts[i] = []
        lh = v_classes[i][1] # Nodes in g
        rh = v_classes[i][2] # Nodes in s
        for perm in permutations(rh)
            push!(isom_parts[i], collect(zip(lh, perm)))
        end
    end
    # Now we look at select permutations of all the nodes by taking the cartessian
    # product of those mappings of the v-classes we just constructed 
    for perm in Iterators.product(isom_parts...)
        pos_isom = Array{Int}(undef, n)
        # The interfaces must map nicely so this is sure information
        for i in 1:length(g.s)
            pos_isom[g.s[i]] = s.s[i]
        end
        for i in 1:length(g.t)
            pos_isom[g.t[i]] = s.t[i]
        end
        # Now we use the vertex profiles to see if we can fill in the rest of
        # pos_isom
        for info_array in perm
            for mapping in info_array
                pos_isom[mapping[1]] = mapping[2]
            end
        end
        #println(perm)
        #println(pos_isom)
        if isIsoIposet(g, s, pos_isom)
            return true
        end
    end
    return false
end

"""Compute the size of a potential gluing of g and s. Return -1 if there's a
mismatch with the interfaces"""
function potGlue(g::Iposet, s::Iposet)
    if length(g.t) != length(s.s)
        return -1
    end
    return nv(g.poset) + nv(s.poset) - length(g.t)
end

"""Generate all iposets with n points, up to isomorphism"""
function genGpIposets(n::Int)
    if n == 0
        return [Iposet((), (), SimpleDiGraph(0))]
    end
    alliposets = Array{Array{Iposet, 1}, 1}(undef, n)
    for i in 1:length(alliposets)
        alliposets[i] = []
    end
    for s in [(), (1,)]
        for t in [(), (1,)]
            push!(alliposets[1], Iposet(s, t, transitiveclosure(SimpleDiGraph(1), true)))
        end
    end
    if n == 1
        return alliposets[1]
    end
    for numpoints in 2:n
        a = Array{Int}(undef, n)
        for i in 1:n
            a[i] = length(alliposets[i])
        end
        println(a)
        for n1 in 1:(numpoints-1)
            for n2 in 1:(numpoints-1)
                for ip1 in alliposets[n1]
                    for ip2 in alliposets[n2]
                        pg = potGlue(ip1, ip2)
                        if pg == numpoints
                            ip = glue(ip1, ip2)
                            seen = false
                            for iq in alliposets[numpoints]
                                if isIsoIposet(iq, ip)
                                    seen = true
                                    break
                                end
                            end
                            if !seen
                                push!(alliposets[numpoints], ip)
                            end
                        else
                            continue
                        end
                    end
                end
            end
        end
        for n1 in 1:(numpoints-1)
            n2 = numpoints - n1
            ips1 = alliposets[n1]
            ips2 = alliposets[n2]
            for ip1 in ips1
                for ip2 in ips2
                    ip = parallel(ip1, ip2)
                    seen = false
                    for iq in alliposets[numpoints]
                        if isIsoIposet(iq, ip)
                            seen = true
                            break
                        end
                    end
                    if !seen
                        push!(alliposets[numpoints], ip)
                    end
                end
            end
        end
    end                  
    return alliposets[n]
end

"""From an array of Iposets, get the different posets associated to them,
up to isomorphism"""
function removeInterfaces(a::Array{Iposet})
    b = Array{SimpleDiGraph}(undef, length(a))
    for i in 1:length(a)
        b[i] = a[i].poset
    end
    res = Array{SimpleDiGraph}(undef, 0)
    for poset in b
        seen = false
        for g in res
            if isIsoV2(poset, g)
                seen = true
            end
        end
        if !seen
            push!(res, poset)
        end
    end
    return res
end

"""Generate all of the possible iposets with n nodes, up to isomorphism"""
function genAllIposets(n::Int)
    println("Getting posets")
    posets = genPosets(n)
    iposets = Array{Iposet}(undef, 0)
    println("Found posets, computing possible interfaces")
    for poset in posets
        mins = minNodes(poset)
        maxes = maxNodes(poset)
        #println(mins)
        for s_length in 0:length(mins)
            for t_length in 0:length(maxes)
                pos_sources = collect(permutations(mins, s_length))
                pos_targets = collect(permutations(maxes, t_length))
                for combination in Iterators.product(pos_sources, pos_targets)
                    #println(combination)
                    push!(iposets, Iposet(Tuple(combination[1]), Tuple(combination[2]), poset))
                end
            end
        end
    end
    println("Got interfaces, checking for isomorphisms")
    representatives = Array{Iposet}(undef, 0)
    ipl = length(iposets)
    for i in 1:length(iposets)
        @printf("Checking iposet %d/%d\n", i, ipl)
        seen = false
        for iq in representatives
            if isIsoIposet(iposets[i], iq)
                seen = true
                break
            end
        end
        if !seen
            push!(representatives, iposets[i])
        end
    end
    println("isomorphisms done")
    return representatives
end
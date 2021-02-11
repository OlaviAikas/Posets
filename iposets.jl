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

"""This is copy-paste from posets.jl so I don't have to import all the functions"""
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

"""Hasse diagram of an Iposet"""
function iHasse(g::Iposet)
    return Iposet(g.s, g.t, hasseDiagram(hasseDiagram(g.poset)))
end

"""Get the gplot object for an Iposet with the interfaces coloured"""
function igplot(p::Iposet)
    colors = Array{ColorTypes.AbstractRGB}(undef, nv(p.poset))
    for i in 1:length(colors)
        if i in p.s
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
            pos = getindex(p.s, v)
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

"""Check if graphs g and s are isomorphic to each other by checking every
permutation"""
function isIsoPerm(g::Iposet, s::Iposet)
    #Start with the easy stuff
    n = nv(g.poset)
    if n != nv(s.poset) || ne(g.poset) != ne(s.poset)
        return false
    end
    if length(g.s) != length(s.s) || length(g.t) != length(s.t)
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

"""Check if the mapping d::Dict is an isomorphism from g to s"""
function isIso(g::Iposet, s::Iposet, d::Dict{Int, Int})
    for i in 1:(length(nv(g.poset))) #Check that the values make a permutation of the nodes (bijection)
        if !(i in values(d)) || !(i in keys(d))
            return false
        end
    end
    if nv(g.poset) != length(keys(d))
        return false
    end
    # Check that the interfaces are mapped nicely
    for i in 1:length(g.s)
        if d[g.s[i]] != s.s[i]
            return false
        end
    end
    for i in 1:length(g.t)
        if d[g.t[i]] != s.t[i]
            return false
        end
    end
    for vertex in keys(d)
        neighbours_in_g = inneighbors(g.poset, vertex)
        neighbours_in_s = inneighbors(s.poset, d[vertex])
        for neighbour in neighbours_in_g
            if !(d[neighbour] in neighbours_in_s)
                return false
            end
        end
        neighbours_out_g = outneighbors(g.poset, vertex)
        neighbours_out_s = outneighbors(s.poset, d[vertex])
        for neighbour in neighbours_out_g
            if !(d[neighbour] in neighbours_out_s)
                return false
            end
        end
    end
    return true
end

"""Check if an isomorphism exists between g and s"""
function isIso(g::Iposet, s::Iposet)
    #Start with the easy stuff
    n = nv(g.poset)
    if n != nv(s.poset) || ne(g.poset) != ne(s.poset)
        return false
    end
    if length(g.s) != length(s.s) || length(g.t) != length(s.t)
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
        gnib = length(inneighbors(g.poset, i))
        if isassigned(g_nb_in_per_node, gnib + 1)
            push!(g_nb_in_per_node[gnib + 1], i)
        else 
            g_nb_in_per_node[gnib + 1] = [i]
        end
        snib = length(inneighbors(s.poset, i))
        if isassigned(s_nb_in_per_node, snib + 1)
            push!(s_nb_in_per_node[snib + 1], i)
        else 
            s_nb_in_per_node[snib + 1] = [i]
        end
        gnob = length(outneighbors(g.poset, i))
        if isassigned(g_nb_out_per_node, gnob + 1)
            push!(g_nb_out_per_node[gnob + 1], i)
        else 
            g_nb_out_per_node[gnob + 1] = [i]
        end
        snob = length(outneighbors(s.poset, i))
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
            push!(alliposets[1], Iposet(s, t, SimpleDiGraph(1)))
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
                ips1 = alliposets[n1]
                ips2 = alliposets[n2]
                for ip1 in ips1
                    for ip2 in ips2
                        ip = Iposet((),(),SimpleDiGraph(0))
                        try
                            ip = glue(ip1, ip2)
                        catch GluingInterfacesDontMatchError
                            continue
                        end
                        ipsize = nv(ip.poset)
                        if ipsize <= n
                            seen = false
                            for iq in alliposets[ipsize]
                                if isIso(iq, ip)
                                    seen = true
                                    break
                                end
                            end
                            if !seen
                                push!(alliposets[ipsize], ip)
                            end
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
                        if isIso(iq, ip)
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

"""Generate all gp-iposets without memoization, debugging purposes only"""
function genGpIposetsNoMemo(n::Int)
    # Base cases
    if n == 0
        return Iposet((),(),SimpleDiGraph(0))
    end
    if n == 1
        alliposets = Array{Iposet}(undef, 0)
        for s in [(), (1,)]
            for t in [(), (1,)]
                push!(alliposets, Iposet(s, t, SimpleDiGraph(1)))
            end
        end
        return alliposets
    end
    # Inductive cases
    alliposets = Array{Iposet}(undef, 0)
    for n1 in 1:(n-1)
        for n2 in 1:(n-1)
            ips1 = genGpIposetsNoMemo(n1)
            ips2 = genGpIposetsNoMemo(n2)
            for ip1 in ips1
                for ip2 in ips2
                    ip = Iposet((),(),SimpleDiGraph(0))
                    try
                        ip = glue(ip1, ip2)
                    catch GluingInterfacesDontMatchError
                        continue
                    end
                    if nv(ip.poset) == n
                        seen = false
                        for iq in alliposets
                            if isIso(ip, iq)
                                seen = true
                                break
                            end
                        end
                        if !seen
                            push!(alliposets, ip)
                        end
                    end
                end
            end
        end
    end
    for n1 in 1:(n-1)
        n2 = n - n1
        ips1 = genGpIposetsNoMemo(n1)
        ips2 = genGpIposetsNoMemo(n2)
        for ip1 in ips1
            for ip2 in ips2
                ip = parallel(ip1, ip2)
                seen = false
                for iq in alliposets
                    if isIso(ip, iq)
                        seen = true
                        break
                    end
                end
                if !seen
                    push!(alliposets, ip)
                end
            end
        end
    end
    return alliposets
end

"""Different style of generating gp-iposets, in case this works. Takes an array
of Iposets, and an integer n, then returns the gp closure of that array with
all graphs up to n nodes, up to iso"""
function gpClosure(a::Array{Iposet}, n::Int)
    old_length::Int  = length(a)
    new_length::Int = 0
    while old_length != new_length
        println(old_length)
        println(new_length)
        for ip1 in a
            for ip2 in a
                gluing = Iposet((), (), SimpleDiGraph(0))
                glued::Bool = false
                try
                    gluing = glue(ip1, ip2)
                    glued = true
                catch GluingInterfacesDontMatchError
                    glued = false
                end
                if glued && nv(gluing.poset) <= n
                    seen::Bool = false
                    for iposet in a
                        if isIso(iposet, gluing)
                            seen = true
                            break
                        end
                    end
                    if !seen
                        push!(a, gluing)
                    end
                end
                if nv(ip1.poset) + nv(ip2.poset) <= n
                    pcomp = parallel(ip1, ip2)
                    seen = false
                    for iposet in a
                        if isIso(iposet, pcomp)
                            seen = true
                            break
                        end
                    end
                    if !seen
                        push!(a, pcomp)
                    end
                end
            end
        end
        old_length = new_length
        new_length = length(a)
    end
    return a
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
            if isIso(iposets[i], iq)
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
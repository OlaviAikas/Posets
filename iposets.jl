using LightGraphs
using Colors
using GraphPlot, Compose
using Combinatorics
import Cairo, Fontconfig
import Base.Threads.@spawn
import Base.==
include("posets.jl")

struct GluingInterfacesDontMatchError <: Exception end

"""Check that the nodes in the tuple s are minimal in the poset p"""
function is_min(p::SimpleDiGraph, s::Tuple{Vararg{Int}})
    for v in s
        if length(inneighbors(p, v)) > 1
            return false
        end
        if length(inneighbors(p, v)) == 1 && inneighbors(p, v) != [v]
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
        if length(outneighbors(p, v)) == 1 && outneighbors(p, v) != [v]
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

"""Transitive closure of an Iposet, doesn't change the original"""
function itransitiveclosure(g::Iposet, b::Bool)
    return Iposet(g.s, g.t, transitiveclosure(g.poset, b))
end

"""Transitive closure of an Iposet, changes the original"""
function itransitiveclosure!(g::Iposet, b::Bool)
    transitiveclosure!(g.poset, b)
end

"""Return the inversion an iposet, which is the same iposet but all of the
edges are reversed and what used to be the target is now the source and
vice versa."""
function inversion(g::Iposet)
    rposet = SimpleDiGraph(nv(g.poset))
    for e in edges(g.poset)
        add_edge!(rposet, dst(e), src(e))
    end
    return Iposet(g.t, g.s, rposet)
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
    il = length(g.t)
    if il != length(p.s)
        throw(GluingInterfacesDontMatchError)
    end
    np = nv(p.poset)
    ng = nv(g.poset)
    rposet = SimpleDiGraph(ng + np - il)
    for e in edges(g.poset)
        add_edge!(rposet, e)
    end
    pmap = zeros(Int, np) # pmap[vertex in p.poset] = vertex in rposet
    for i in 1:il
        pmap[p.s[i]] = g.t[i]
    end
    counter = 1
    for v in 1:np
        if pmap[v] == 0
            pmap[v] = ng + counter
            counter += 1
        end
    end
    for e in edges(p.poset)
        add_edge!(rposet, pmap[src(e)], pmap[dst(e)])
    end
    for v1 in vertices(g.poset)
        if !(v1 in g.t)
            for v2 in (ng + 1):(ng + np - il)
                add_edge!(rposet, v1, v2)
            end
        end
    end
    tl = length(p.t)
    rt = Array{Int}(undef, tl)
    for i in 1:tl
        rt[i] = pmap[p.t[i]]
    end
    #foreach(println, edges(rposet))
    return Iposet(g.s, Tuple(rt), rposet)
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
        g_v_profiles[v] = (inHash(g.poset, v), outHash(g.poset, v))
        s_v_profiles[v] = (inHash(s.poset, v), outHash(s.poset, v))
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
    indexed_g = Array{Bool}(undef, n)
    indexed_s = Array{Bool}(undef, n)
    for i in 1:n
        if i in g.s || i in g.t
            indexed_g[i] = true # We already know where the interfaces are mapped
        else                    # to, so they don't need to be put in classes
            indexed_g[i] = false
        end
        if i in s.s || i in s.t
            indexed_s[i] = true # We already know where the interfaces are mapped
        else                    # to, so they don't need to be put in classes
            indexed_s[i] = false
        end   
    end
    for i in 1:n
        if indexed_g[i]
            continue
        end
        class = (Array{Int}(undef, 1), Array{Int}(undef, 0))
        class[1][1] = i # Construct a new class with all of the nodes with the
                        # same invariant as "i"
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
        mapped = zeros(Bool, n)
        # The interfaces must map nicely so this is sure information
        for i in 1:length(g.s)
            pos_isom[g.s[i]] = s.s[i]
            mapped[g.s[i]] = true
        end
        for i in 1:length(g.t)
            pos_isom[g.t[i]] = s.t[i]
            mapped[g.t[i]] = true
        end
        # Now we use the vertex profiles to see if we can fill in the rest of
        # pos_isom
        for info_array in perm
            if !(false in mapped)
                break
            end
            for mapping in info_array
                pos_isom[mapping[1]] = mapping[2]
            end
        end
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

"""Generate all gp-iposets with n points, up to isomorphism"""
function genGpIposets(n::Int)
    if n == 0
        return [Iposet((), (), SimpleDiGraph(0))]
    end
    # The indexing format is [num_edges + 1, size_s + 1, size_t + 1, numpoints]
    alliposets = Array{Array{Tuple{Iposet, Array{Tuple{Int, Int}}}}}(undef, (n*(n-1))÷2 + 1, n + 1, n + 1, n)
    @inbounds for i in eachindex(alliposets)
        alliposets[i] = []
    end
    for s in [(), (1,)]
        for t in [(), (1,)]
            push!(alliposets[1, length(s) + 1, length(t) + 1, 1], (Iposet(s, t, SimpleDiGraph(1)), [(0, 0)]))
        end
    end
    if n == 1
        vc = vcat(alliposets[:,:,:,1]...)
        res = Array{Iposet}(undef, length(vc))
        @inbounds for i in 1:length(vc)
            res[i] = vc[i][1]
        end
        return res
    end
    for numpoints in 2:n
        for n1 in 1:(numpoints-1)
            for n2 in 1:(numpoints-1)
                u = n1 + n2 - numpoints
                if u < 0 || u >= n1 || u >= n2
                    continue
                end
                for ip1 in vcat(alliposets[:, :, u + 1, n1]...)
                    for ip2 in vcat(alliposets[:, u + 1, :, n2]...)
                        ip = itransitiveclosure(glue(ip1[1], ip2[1]), false)
                        ipe = ne(ip.poset)
                        ips = length(ip.s)
                        ipt = length(ip.t)
                        vprof = Array{Tuple{Int, Int}}(undef, numpoints)
                        @inbounds for v in 1:numpoints
                            vprof[v] = (inHash(ip.poset, v), outHash(ip.poset, v))
                        end
                        seen = false
                        for iq in alliposets[ipe + 1, ips + 1, ipt + 1, numpoints]
                            if isIsoIposetX(iq[1], iq[2], ip, vprof)
                                seen = true
                                break
                            end
                        end
                        if !seen
                            push!(alliposets[ipe + 1, ips + 1, ipt + 1, numpoints], (ip, vprof))
                        end
                    end
                end
            end
        end
        for n1 in 1:(numpoints-1)
            n2 = numpoints - n1
            ips1 = vcat(alliposets[:, :, :, n1]...)
            ips2 = vcat(alliposets[:, :, :, n2]...)
            for ip1 in ips1
                for ip2 in ips2
                    ip = itransitiveclosure(parallel(ip1[1], ip2[1]), false)
                    ipe = ne(ip.poset)
                    ips = length(ip.s)
                    ipt = length(ip.t)
                    vprof = Array{Tuple{Int, Int}}(undef, numpoints)
                    for v in 1:numpoints
                        vprof[v] = (inHash(ip.poset, v), outHash(ip.poset, v))
                    end
                    seen = false
                    for iq in alliposets[ipe + 1, ips + 1, ipt + 1, numpoints]
                        if isIsoIposetX(iq[1], iq[2], ip, vprof)
                            seen = true
                            break
                        end
                    end
                    if !seen
                        push!(alliposets[ipe + 1, ips + 1, ipt + 1, numpoints], (ip, vprof))
                    end
                end
            end
        end
    end
    vc = vcat(alliposets[:, :, :, n]...)
    res = Array{Iposet}(undef, length(vc))
    @inbounds for i in 1:length(vc)
        res[i] = vc[i][1]
    end
    return res
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
            if isIso(poset, g)
                seen = true
            end
        end
        if !seen
            push!(res, poset)
        end
    end
    return res
end

"""Recursively compute the set (array) of gp-iposets with n nodes."""
function gpiPosets(n)
    alliposets = Array{Array{Tuple{Iposet, Array{Tuple{Int, Int}}}}}(undef, (n*(n-1))÷2 + 1, n + 1, n + 1, n)
    @inbounds for i in eachindex(alliposets)
        alliposets[i] = []
    end
    for s in [(), (1,)]
        for t in [(), (1,)]
            push!(alliposets[1, length(s) + 1, length(t) + 1, 1], (Iposet(s, t, SimpleDiGraph(1)), [(0, 0)]))
        end
    end
    filled = zeros(Bool, n + 1, n + 1, n)
    filled[1:2, 1:2, 1] = [true true; true true]
    locks = Array{ReentrantLock}(undef, (n*(n-1))÷2 + 2, n + 1, n + 1, n)
    for i in eachindex(locks)
        locks[i] = ReentrantLock()
    end
    elems = Array{Array{Iposet}}(undef, n + 1, n + 1)
    for i in eachindex(elems)
        elems[i] = []
    end
    Threads.@threads for i in 1:(n+1)
        Threads.@threads for j in 1:(n+1)
            elems[i, j] = gpiPosets(n, i - 1, j - 1, alliposets, filled, locks)
        end
    end
    return vcat(elems...)
end

"""Recursively compute the set (array) of gp-posets with n nodes."""
function gpPosets(n)
    alliposets = Array{Array{Tuple{Iposet, Array{Tuple{Int, Int}}}}}(undef, (n*(n-1))÷2 + 1, n + 1, n + 1, n)
    @inbounds for i in eachindex(alliposets)
        alliposets[i] = []
    end
    for s in [(), (1,)]
        for t in [(), (1,)]
            push!(alliposets[1, length(s) + 1, length(t) + 1, 1], (Iposet(s, t, SimpleDiGraph(1)), [(0, 0)]))
        end
    end
    filled = zeros(Bool, n + 1, n + 1, n)
    filled[1:2, 1:2, 1] = [true true; true true]
    locks = Array{ReentrantLock}(undef, (n*(n-1))÷2 + 2, n + 1, n + 1, n)
    for i in eachindex(locks)
        locks[i] = ReentrantLock()
    end
    a = gpiPosets(n, 0, 0, alliposets, filled, locks)
    res = Array{SimpleDiGraph}(undef, length(a))
    @inbounds for i in 1:length(a)
        res[i] = a[i].poset
    end
    return res
end

"""Recursively compute the set of gp-iposets with k sources and l targets.
The arrays alliposets and filled are for memoisation"""
function gpiPosets(n, k, l, alliposets, filled, locks)
    #@printf("Fill ratio: %f\n", count(x -> x, filled)/length(filled))
    lock(locks[end, k + 1, l + 1, n])
    if filled[k + 1, l + 1, n]
        unlock(locks[end, k + 1, l + 1, n])
        vc = vcat(alliposets[:, k + 1, l + 1, n]...)
        res = Array{Iposet}(undef, length(vc))
        @inbounds for i in 1:length(vc)
            res[i] = vc[i][1]
        end
        return res
    end
    if filled[l + 1, k + 1, n]
        for ne in 1:size(alliposets, 1)
            nes = length(alliposets[ne, l + 1, k + 1, n])
            iposets = Array{Tuple{Iposet, Array{Tuple{Int, Int}}}}(undef, nes)
            for i in 1:nes
                nip = inversion(alliposets[ne, l + 1, k + 1, n][i][1])
                vproflen = length(alliposets[ne, l + 1, k + 1, n][i][2])
                nvprof = Array{Tuple{Int, Int}}(undef, vproflen)
                for j in 1:vproflen
                    nvprof[j] = (alliposets[ne, l + 1, k + 1, n][i][2][j][2], alliposets[ne, l + 1, k + 1, n][i][2][j][1])
                end
                iposets[i] = (nip, nvprof)
            end
            alliposets[ne, k + 1, l + 1, n] = iposets
        end
        filled[k + 1, l + 1, n] = true
        unlock(locks[end, k + 1, l + 1, n])
        vc = vcat(alliposets[:, k + 1, l + 1, n]...)
        res = Array{Iposet}(undef, length(vc))
        @inbounds for i in 1:length(vc)
            res[i] = vc[i][1]
        end
        return res
    end
    Threads.@threads for n1 in 1:(n-1)
        Threads.@threads for n2 in 1:(n-1)
            u = n1 + n2 - n
            if u < 0
                continue
            end
            for ip1 in gpiPosets(n1, k, u, alliposets, filled, locks)
                for ip2 in gpiPosets(n2, u, l, alliposets, filled, locks)
                    ip = glue(ip1, ip2)
                    ipe = ne(ip.poset)
                    vprof = Array{Tuple{Int, Int}}(undef, n)
                    @inbounds for v in 1:n
                        vprof[v] = (inHash(ip.poset, v), outHash(ip.poset, v))
                    end
                    lock(locks[ipe + 1, k + 1, l + 1, n])
                    seen = false
                    for iq in alliposets[ipe + 1, k + 1, l + 1, n]
                        if isIsoIposetX(iq[1], iq[2], ip, vprof)
                            seen = true
                            break
                        end
                    end
                    if !seen
                        push!(alliposets[ipe + 1, k + 1, l + 1, n], (ip, vprof))
                    end
                    unlock(locks[ipe + 1, k + 1, l + 1, n])
                end
            end
        end
    end
    Threads.@threads for n1 in 1:(n-1)
        n2 = n - n1
        for ip1s in 0:k
            ip2s = k - ip1s
            for ip1t in 0:l
                ip2t = l - ip1t
                for ip1 in gpiPosets(n1, ip1s, ip1t, alliposets, filled, locks)
                    for ip2 in gpiPosets(n2, ip2s, ip2t, alliposets, filled, locks)
                        ip = parallel(ip1, ip2)
                        ipe = ne(ip.poset)
                        vprof = Array{Tuple{Int, Int}}(undef, n)
                        for v in 1:n
                            vprof[v] = (inHash(ip.poset, v), outHash(ip.poset, v))
                        end
                        lock(locks[ipe + 1, k + 1, l + 1, n])
                        seen = false
                        for iq in alliposets[ipe + 1, k + 1, l + 1, n]
                            if isIsoIposetX(iq[1], iq[2], ip, vprof)
                                seen = true
                                break
                            end
                        end
                        if !seen
                            push!(alliposets[ipe + 1, k + 1, l + 1, n], (ip, vprof))
                        end
                        unlock(locks[ipe + 1, k + 1, l + 1, n])
                    end
                end
            end
        end
    end
    filled[k + 1, l + 1, n] = true
    unlock(locks[end, k + 1, l + 1, n])
    vc = vcat(alliposets[:, k + 1, l + 1, n]...)
    res = Array{Iposet}(undef, length(vc))
    @inbounds for i in 1:length(vc)
        res[i] = vc[i][1]
    end
    return res
end


"""Generate all of the possible iposets with n nodes, up to isomorphism"""
function genAllIposets(n::Int)
    println("Getting posets")
    posets = genPosets(n)
    println("Got posets")
    np = length(posets)
    iposets = Array{Array{Iposet}}(undef, n + 1, n + 1, np)
    for i in eachindex(iposets)
        iposets[i] = []
    end
    vprofs = Array{Array{Tuple{Int, Int}}}(undef, length(posets))
    for i in 1:np
        vprofs[i] = Array{Tuple{Int, Int}}(undef, n)
        for v in vertices(posets[i])
            vprofs[i][v] = (inHash(posets[i], v), outHash(posets[i], v))
        end
    end
    for i in 1:np
        println("Interfacing poset $i")
        Threads.@threads for k in 0:n
            Threads.@threads for l in 0:n
                mins = minNodes(posets[i])
                maxes = maxNodes(posets[i])
                pos_sources = collect(permutations(mins, k))
                pos_targets = collect(permutations(maxes, l))
                for combination in Iterators.product(pos_sources, pos_targets)
                    nip = Iposet(Tuple(combination[1]), Tuple(combination[2]), posets[i])
                    seen = false
                    for ip in iposets[k + 1, l + 1, i]
                        if isIsoIposetX(ip, vprofs[i], nip, vprofs[i])
                            seen = true
                            break
                        end
                    end
                    if !seen
                        push!(iposets[k + 1, l + 1, i], nip)
                    end
                end
            end
        end
    end
    return vcat(iposets...)
end

"""Generate all of the iposets with n nodes, k source interfaces and l target
interfaces"""
function genIposets(n::Int, k::Int, l::Int)
    posets = genPosets(n)
    iposets = Array{Array{Iposet}}(undef, length(posets))
    vprofs = Array{Array{Tuple{Int, Int}}}(undef, length(posets))
    for i in 1:length(iposets)
        iposets[i] = []
    end
    for i in 1:length(posets)
        vprofs[i] = Array{Tuple{Int, Int}}(undef, n)
        for v in vertices(posets[i])
            vprofs[i][v] = (inHash(posets[i], v), outHash(posets[i], v))
        end
    end
    for i in 1:length(posets)
        mins = minNodes(posets[i])
        maxes = maxNodes(posets[i])
        pos_sources = collect(permutations(mins, k))
        pos_targets = collect(permutations(maxes, l))
        for combination in Iterators.product(pos_sources, pos_targets)
            nip = Iposet(Tuple(combination[1]), Tuple(combination[2]), posets[i])
            seen = false
            for ip in iposets[i]
                if isIsoIposetX(ip, vprofs[i], nip, vprofs[i])
                    seen = true
                    break
                end
            end
            if !seen
                push!(iposets[i], nip)
            end
        end
    end
    return vcat(iposets...)
end

"""Special isomorphisms checking function that uses tuples instead of arrays
for the more streamlined isIsoIposetX checker"""
function isIsoIposetX(g::Iposet, s::Iposet, a::Tuple{Vararg{Int}})
    # Check that the interfaces are mapped nicely
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

"""Special isomorphism function that takes the vertex profiles in as arguments
so they're not recomputed every time. Also does not check for the number of
nodes or edges or source interfaces or target interfaces, since those are 
handled by the caller. ATTENTION: this function can crash or produce incorrect
results if you do not check that the number of nodes, edges and interfaces
match in the two graphs"""
function isIsoIposetX(g::Iposet, g_v_profiles::Array{Tuple{Int, Int}}, s::Iposet, s_v_profiles::Array{Tuple{Int, Int}})
    #Start with the easy stuff
    n = nv(g.poset)
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
    indexed_g = zeros(Bool, n)
    indexed_s = zeros(Bool, n)
    @inbounds for i in 1:length(g.s)
        push!(targets[g.s[i]], s.s[i])
        indexed_g[g.s[i]] = true
        indexed_s[s.s[i]] = true
    end
    @inbounds for i in 1:length(g.t)
        if !(s.t[i] in targets[g.t[i]])
            push!(targets[g.t[i]], s.t[i])
            indexed_g[g.t[i]] = true
            indexed_s[s.t[i]] = true
        end
    end
    for a in targets
        if length(a) > 1
            return false
        end
    end
    @inbounds for v in 1:n
        if indexed_g[v]
            continue
        end
        for i in 1:n
            if g_v_profiles[v] == s_v_profiles[i] && !indexed_s[i]
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
        if isIsoIposetX(g, s, pos_isom)
            return true
        end
    end
    return false
end

"""Check if an iposet is "almost weakly connected", e.g have one weakly connected
component and 0 or more individual points"""
function almostConnected(g::Iposet)
    comps = weakly_connected_components(g.poset)
    big_one = false
    for comp in comps
        if length(comp) > 1 && big_one == false
            big_one = true
        elseif length(comp) <= 1
            continue
        else
            return false
        end
    end
    return true
end

"""Return a string representation of an Iposet to be printed to file"""
function toString(g::Iposet)
    n_vertices = string(nv(g.poset), base=16)
    n_edges = string(ne(g.poset), base=16)

    res = "$n_vertices $n_edges "

    e_list = [string(src(e), base=16)*string(dst(e), base=16) for e in edges(g.poset)]
    for e in e_list
        res *= e * " "
    end

    if length(g.s) == 0 && length(g.t) == 0
        return res[1:end-1]
    end
    if length(g.s) == 0
        res *= "- "
        for t in g.t
            res *= string(t, base=16)
        end
        return res
    end
    for s in g.s
        res *= string(s, base=16)
    end
    if length(g.t) == 0
        return res
    end
    res *= " "
    for t in g.t
        res *= string(t, base=16)
    end
    return res
end

"""Save an array of iposets into the given filename"""
function saveIposets(a::Array{Iposet}, filename::String)
    open(filename, "w") do file
        for ips in a
            write(file, toString(ips)*"\n")
        end
    end
end

"""Read a file of iposets and return them in an array"""
function loadIposets(filename::String)
    res = Array{Iposet}(undef, 0)
    open(filename, "r") do file
        for line in eachline(file)
            nums = split(line)
            np = parse(Int, nums[1], base=16) # number of points
            ne = parse(Int, nums[2], base=16) # number of edges
            #println(nums, np, ne)
            dg = SimpleDiGraph(np)
            for edges_num in nums[3:ne+2]
                add_edge!(dg, parse(Int, edges_num[1], base=16), parse(Int, edges_num[2], base=16))
            end
            if length(nums) == ne+2 # no sources neither targets
                s, t = (), ()
            elseif length(nums) == ne+3 # sources, but no targets
                s = Tuple([parse(Int, x, base=16) for x in nums[ne+3]])
                t = ()
            elseif length(nums) == ne+4 && nums[ne+3] == "-" # targets, but no sources
                s = ()
                t = Tuple([parse(Int, x, base=16) for x in nums[ne+4]])
            else # sources and targets
                s = Tuple([parse(Int, x, base=16) for x in nums[ne+3]])
                t = Tuple([parse(Int, x, base=16) for x in nums[ne+4]])
            end
            push!(res, Iposet(s, t, dg))
    end
        return res
    end
end

"""Overload equality for Array{Iposet} so you can easily check if two arrays
have the same iposets up to isomorphism"""
function ==(a::Array{Iposet}, b::Array{Iposet})
    for ips1 in a
        seen = false
        for ips2 in b
            if isIsoIposet(ips1, ips2)
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
# Generate gp-iposets
n = 9
k = 0
l = 5

using LightGraphs

"""Iposet: source tuple; target tuple; poset.
Poset is transitively closed, and without reflexive edges(!)
Sources must be minimal, targets maximal: this is not checked!
"""
struct Iposet
    s::Tuple{Vararg{Int}}
    t::Tuple{Vararg{Int}}
    poset::SimpleDiGraph
#    Iposet(s, t, p) = is_min(p, s) && is_max(p, t) ? new(s, t, p) : error("Interfaces not minimal/maximal")
    Iposet(s, t, p) = new(s, t, p)
end

function itransitiveclosure(g::Iposet)
    return Iposet(g.s, g.t, transitiveclosure(g.poset, false))
end

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
    pmap = zeros(Int, np)
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
    return Iposet(g.s, Tuple(rt), rposet)
end

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

function iHasse(g::Iposet)
    return Iposet(g.s, g.t, hasseDiagram(g.poset))
end

function iposetToString(g::Iposet)
    g = iHasse(g)
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

function saveIposets(a::Array{Iposet}, filename::String)
    open(filename, "w") do file
        for ips in a
            write(file, iposetToString(ips)*"\n")
        end
    end
end

function iposetFromString(line::String)
    nums = split(line)
    np = parse(Int, nums[1], base=16)
    ne = parse(Int, nums[2], base=16)
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
    return itransitiveclosure(Iposet(s, t, dg))
end

function loadIposets(filename::String)
    res = Array{Iposet}(undef, 0)
    open(filename, "r") do file
        for line in eachline(file)
            push!(res, iposetFromString(line))
        end
    end
    return res
end

function isIsoIposetX(g::Iposet, s::Iposet, a::Tuple{Vararg{Int}})
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

function isIsoIposetX(g::Iposet, g_v_profiles::Array{Tuple{Int, Int}}, s::Iposet, s_v_profiles::Array{Tuple{Int, Int}})
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

function inversion(g::Iposet)
    rposet = SimpleDiGraph(nv(g.poset))
    for e in edges(g.poset)
        add_edge!(rposet, dst(e), src(e))
    end
    return Iposet(g.t, g.s, rposet)
end

function gpiPosets(n, k, l, alliposets, filled, locks)
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
        println("Generating gpi-", n, "-", k, "-", l, " by inversion")
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
        println("Generating gpi-", n, "-", k, "-", l, " by inversion: Done")
        return res
    end
    println("Generating gpi-", n, "-", k, "-", l, " recursively")
    Threads.@threads for n1 in 1:(n-1)
        Threads.@threads for n2 in 1:(n-1)
            u = n1 + n2 - n
            if u < 0 || u > n1 || u > n2 || k > n1 || l > n2
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
                if n1 < ip1s || n1 < ip1t || n2 < ip2s || n2 < ip2t
                    continue
                end
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

function gpiPosetsMemo(n, k, l)
    infname = string("gpiupto",n-1,".ips")
    alliposets = Array{Array{Tuple{Iposet, Array{Tuple{Int, Int}}}}}(undef, (n*(n-1))÷2 + 1, n + 1, n + 1, n)
    @inbounds for i in eachindex(alliposets)
        alliposets[i] = []
    end
    filled = zeros(Bool, n + 1, n + 1, n)
    println("Loading gpi cache")
    open(infname, "r") do file
        for line in eachline(file)
            ip = iposetFromString(line)
            np = nv(ip.poset)
            ipe = ne(ip.poset)
            kp = length(ip.s)
            lp = length(ip.t)
            vprof = Array{Tuple{Int, Int}}(undef, np)
            @inbounds for v in 1:np
                vprof[v] = (inHash(ip.poset, v), outHash(ip.poset, v))
            end
            push!(alliposets[ipe + 1, kp + 1, lp + 1, np], (ip, vprof))
            # if we see one, we see all
            filled[kp + 1, lp + 1, np] = true
        end
    end
    locks = Array{ReentrantLock}(undef, (n*(n-1))÷2 + 2, n + 1, n + 1, n)
    for i in eachindex(locks)
        locks[i] = ReentrantLock()
    end
    println("Start generating")
    return gpiPosets(n, k, l, alliposets, filled, locks)
end

function ggpi(n, k, l)
    outfname = string("gpi",n,"-",k,"-",l,".ips")
    println("Generating ", outfname)
    ips = gpiPosetsMemo(n, k, l)
    println("Done: ", length(ips))
    saveIposets(ips, outfname)
end

#if n > 7 && n < 10
ggpi(n, k, l)
#end

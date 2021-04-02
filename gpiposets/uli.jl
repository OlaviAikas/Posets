include("ggpi.jl")

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

function all_min(p::SimpleDiGraph)
    res = Array{Int}(undef, 0)
    for v in vertices(p)
        if length(inneighbors(p, v)) == 0
            push!(res, v)
        end
    end
    return res
end

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

function all_max(p::SimpleDiGraph)
    res = Array{Int}(undef, 0)
    for v in vertices(p)
        if length(outneighbors(p, v)) == 0
            push!(res, v)
        end
    end
    return res
end

function is_left_Winkowski(g::Iposet)
    for v in all_min(g.poset)
        if !( v in g.s)
            return false
        end
    end
    return true
end

function is_right_Winkowski(g::Iposet)
    for v in all_max(g.poset)
        if !( v in g.t)
            return false
        end
    end
    return true
end

function is_Winkowski(g::Iposet)
    return is_left_Winkowski(g) && is_right_Winkowski(g)
end

function is_starter(g::Iposet)
    return (ne(g.poset) == 0) && is_right_Winkowski(g)
end    

function is_terminator(g::Iposet)
    return (ne(g.poset) == 0) && is_left_Winkowski(g)
end    

# extract Winskowskis
function extractW(n)
    infile = string("gpi",n,".ips")
    outfile = string("gpwi",n,".ips")
    ips = loadIposets(infile)
    res = Array{Iposet}(undef, 0)
    for ip in ips
        if is_Winkowski(ip)
            push!(res, ip)
        end
    end
    saveIposets(res, outfile)
end

function extractW8(k, l)
    infile = string("gpi8-",k,"-",l,".ips")
    outfile = string("gpwi8-",k,"-",l,".ips")
    ips = loadIposets(infile)
    res = Array{Iposet}(undef, 0)
    for ip in ips
        if is_Winkowski(ip)
            push!(res, ip)
        end
    end
    saveIposets(res, outfile)
end

# extract starters
function extractS(n)
    infile = string("gpi",n,".ips")
    outfile = string("gpsi",n,".ips")
    ips = loadIposets(infile)
    res = Array{Iposet}(undef, 0)
    for ip in ips
        if is_starter(ip)
            push!(res, ip)
        end
    end
    saveIposets(res, outfile)
end

# extract terminators
function extractT(n)
    infile = string("gpi",n,".ips")
    outfile = string("gpti",n,".ips")
    ips = loadIposets(infile)
    res = Array{Iposet}(undef, 0)
    for ip in ips
        if is_terminator(ip)
            push!(res, ip)
        end
    end
    saveIposets(res, outfile)
end

function gpwstiPosets(n, k, l, allwiposets, allsiposets, alltiposets, filled, locks)
    lock(locks[end, k + 1, l + 1, n])
    if filled[k + 1, l + 1, n]
        unlock(locks[end, k + 1, l + 1, n])
        vcw = vcat(allwiposets[:, k + 1, l + 1, n]...)
        vcs = vcat(allsiposets[:, k + 1, l + 1, n]...)
        vct = vcat(alltiposets[:, k + 1, l + 1, n]...)
        resw = Array{Iposet}(undef, length(vcw))
        ress = Array{Iposet}(undef, length(vcs))
        rest = Array{Iposet}(undef, length(vct))
        @inbounds for i in 1:length(vcw)
            resw[i] = vcw[i][1]
        end
        @inbounds for i in 1:length(vcs)
            ress[i] = vcs[i][1]
        end
        @inbounds for i in 1:length(vct)
            rest[i] = vct[i][1]
        end
        return resw, ress, rest
    end
    if filled[l + 1, k + 1, n]
        for ne in 1:size(allwiposets, 1)
            nes = length(allwiposets[ne, l + 1, k + 1, n])
            iposets = Array{Tuple{Iposet, Array{Tuple{Int, Int}}}}(undef, nes)
            for i in 1:nes
                nip = inversion(allwiposets[ne, l + 1, k + 1, n][i][1])
                vproflen = length(allwiposets[ne, l + 1, k + 1, n][i][2])
                nvprof = Array{Tuple{Int, Int}}(undef, vproflen)
                for j in 1:vproflen
                    nvprof[j] = (allwiposets[ne, l + 1, k + 1, n][i][2][j][2], allwiposets[ne, l + 1, k + 1, n][i][2][j][1])
                end
                iposets[i] = (nip, nvprof)
            end
            allwiposets[ne, k + 1, l + 1, n] = iposets
        end
        for ne in 1:size(allsiposets, 1)
            nes = length(allsiposets[ne, l + 1, k + 1, n])
            iposets = Array{Tuple{Iposet, Array{Tuple{Int, Int}}}}(undef, nes)
            for i in 1:nes
                nip = inversion(allsiposets[ne, l + 1, k + 1, n][i][1])
                vproflen = length(allsiposets[ne, l + 1, k + 1, n][i][2])
                nvprof = Array{Tuple{Int, Int}}(undef, vproflen)
                for j in 1:vproflen
                    nvprof[j] = (allsiposets[ne, l + 1, k + 1, n][i][2][j][2], allsiposets[ne, l + 1, k + 1, n][i][2][j][1])
                end
                iposets[i] = (nip, nvprof)
            end
            alltiposets[ne, k + 1, l + 1, n] = iposets
        end
        for ne in 1:size(alltiposets, 1)
            nes = length(alltiposets[ne, l + 1, k + 1, n])
            iposets = Array{Tuple{Iposet, Array{Tuple{Int, Int}}}}(undef, nes)
            for i in 1:nes
                nip = inversion(alltiposets[ne, l + 1, k + 1, n][i][1])
                vproflen = length(alltiposets[ne, l + 1, k + 1, n][i][2])
                nvprof = Array{Tuple{Int, Int}}(undef, vproflen)
                for j in 1:vproflen
                    nvprof[j] = (alltiposets[ne, l + 1, k + 1, n][i][2][j][2], alltiposets[ne, l + 1, k + 1, n][i][2][j][1])
                end
                iposets[i] = (nip, nvprof)
            end
            allsiposets[ne, k + 1, l + 1, n] = iposets
        end
        filled[k + 1, l + 1, n] = true
        unlock(locks[end, k + 1, l + 1, n])
        vcw = vcat(allwiposets[:, k + 1, l + 1, n]...)
        vcs = vcat(allsiposets[:, k + 1, l + 1, n]...)
        vct = vcat(alltiposets[:, k + 1, l + 1, n]...)
        resw = Array{Iposet}(undef, length(vcw))
        ress = Array{Iposet}(undef, length(vcs))
        rest = Array{Iposet}(undef, length(vct))
        @inbounds for i in 1:length(vcw)
            resw[i] = vcw[i][1]
        end
        @inbounds for i in 1:length(vcs)
            ress[i] = vcs[i][1]
        end
        @inbounds for i in 1:length(vct)
            rest[i] = vct[i][1]
        end
        return resw, ress, rest
    end
    # generate Winkowskis
    Threads.@threads for n1 in 1:(n-1)
        Threads.@threads for n2 in 1:(n-1)
            u = n1 + n2 - n
            if u < 0 || u >= n1 || u >= n2
                continue
            end
            # ... from Winkowski gluings
            for ip1 in gpwstiPosets(n1, k, u, allwiposets, allsiposets, alltiposets, filled, locks)[1]
                for ip2 in gpwstiPosets(n2, u, l, allwiposets, allsiposets, alltiposets, filled, locks)[1]
                    ip = glue(ip1, ip2)
                    ipe = ne(ip.poset)
                    vprof = Array{Tuple{Int, Int}}(undef, n)
                    @inbounds for v in 1:n
                        vprof[v] = (inHash(ip.poset, v), outHash(ip.poset, v))
                    end
                    lock(locks[ipe + 1, k + 1, l + 1, n])
                    seen = false
                    for iq in allwiposets[ipe + 1, k + 1, l + 1, n]
                        if isIsoIposetX(iq[1], iq[2], ip, vprof)
                            seen = true
                            break
                        end
                    end
                    if !seen
                        push!(allwiposets[ipe + 1, k + 1, l + 1, n], (ip, vprof))
                    end
                    unlock(locks[ipe + 1, k + 1, l + 1, n])
                end
            end
            # ... from T*S gluings
            for ip1 in gpwstiPosets(n1, k, u, allwiposets, allsiposets, alltiposets, filled, locks)[3]
                for ip2 in gpwstiPosets(n2, u, l, allwiposets, allsiposets, alltiposets, filled, locks)[2]
                    ip = glue(ip1, ip2)
                    ipe = ne(ip.poset)
                    vprof = Array{Tuple{Int, Int}}(undef, n)
                    @inbounds for v in 1:n
                        vprof[v] = (inHash(ip.poset, v), outHash(ip.poset, v))
                    end
                    lock(locks[ipe + 1, k + 1, l + 1, n])
                    seen = false
                    for iq in allwiposets[ipe + 1, k + 1, l + 1, n]
                        if isIsoIposetX(iq[1], iq[2], ip, vprof)
                            seen = true
                            break
                        end
                    end
                    if !seen
                        push!(allwiposets[ipe + 1, k + 1, l + 1, n], (ip, vprof))
                    end
                    unlock(locks[ipe + 1, k + 1, l + 1, n])
                end
            end
        end
    end
    # ... from Winkowski parallels
    Threads.@threads for n1 in 1:(n-1)
        n2 = n - n1
        for ip1s in 0:k
            ip2s = k - ip1s
            for ip1t in 0:l
                ip2t = l - ip1t
                for ip1 in gpwstiPosets(n1, ip1s, ip1t, allwiposets, allsiposets, alltiposets, filled, locks)[1]
                    for ip2 in gpwstiPosets(n2, ip2s, ip2t, allwiposets, allsiposets, alltiposets, filled, locks)[1]
                        ip = parallel(ip1, ip2)
                        ipe = ne(ip.poset)
                        vprof = Array{Tuple{Int, Int}}(undef, n)
                        for v in 1:n
                            vprof[v] = (inHash(ip.poset, v), outHash(ip.poset, v))
                        end
                        lock(locks[ipe + 1, k + 1, l + 1, n])
                        seen = false
                        for iq in allwiposets[ipe + 1, k + 1, l + 1, n]
                            if isIsoIposetX(iq[1], iq[2], ip, vprof)
                                seen = true
                                break
                            end
                        end
                        if !seen
                            push!(allwiposets[ipe + 1, k + 1, l + 1, n], (ip, vprof))
                        end
                        unlock(locks[ipe + 1, k + 1, l + 1, n])
                    end
                end
            end
        end
    end
    # generate starters
    Threads.@threads for n1 in 1:(n-1)
        n2 = n - n1
        for ip1s in 0:k
            ip2s = k - ip1s
            for ip1t in 0:l
                ip2t = l - ip1t
                for ip1 in gpwstiPosets(n1, ip1s, ip1t, allwiposets, allsiposets, alltiposets, filled, locks)[2]
                    for ip2 in gpwstiPosets(n2, ip2s, ip2t, allwiposets, allsiposets, alltiposets, filled, locks)[2]
                        ip = parallel(ip1, ip2)
                        ipe = ne(ip.poset)
                        vprof = Array{Tuple{Int, Int}}(undef, n)
                        for v in 1:n
                            vprof[v] = (inHash(ip.poset, v), outHash(ip.poset, v))
                        end
                        lock(locks[ipe + 1, k + 1, l + 1, n])
                        seen = false
                        for iq in allsiposets[ipe + 1, k + 1, l + 1, n]
                            if isIsoIposetX(iq[1], iq[2], ip, vprof)
                                seen = true
                                break
                            end
                        end
                        if !seen
                            push!(allsiposets[ipe + 1, k + 1, l + 1, n], (ip, vprof))
                        end
                        unlock(locks[ipe + 1, k + 1, l + 1, n])
                    end
                end
            end
        end
    end
    # generate terminators
    Threads.@threads for n1 in 1:(n-1)
        n2 = n - n1
        for ip1s in 0:k
            ip2s = k - ip1s
            for ip1t in 0:l
                ip2t = l - ip1t
                for ip1 in gpwstiPosets(n1, ip1s, ip1t, allwiposets, allsiposets, alltiposets, filled, locks)[3]
                    for ip2 in gpwstiPosets(n2, ip2s, ip2t, allwiposets, allsiposets, alltiposets, filled, locks)[3]
                        ip = parallel(ip1, ip2)
                        ipe = ne(ip.poset)
                        vprof = Array{Tuple{Int, Int}}(undef, n)
                        for v in 1:n
                            vprof[v] = (inHash(ip.poset, v), outHash(ip.poset, v))
                        end
                        lock(locks[ipe + 1, k + 1, l + 1, n])
                        seen = false
                        for iq in alltiposets[ipe + 1, k + 1, l + 1, n]
                            if isIsoIposetX(iq[1], iq[2], ip, vprof)
                                seen = true
                                break
                            end
                        end
                        if !seen
                            push!(alltiposets[ipe + 1, k + 1, l + 1, n], (ip, vprof))
                        end
                        unlock(locks[ipe + 1, k + 1, l + 1, n])
                    end
                end
            end
        end
    end
    filled[k + 1, l + 1, n] = true
    unlock(locks[end, k + 1, l + 1, n])
    vcw = vcat(allwiposets[:, k + 1, l + 1, n]...)
    vcs = vcat(allsiposets[:, k + 1, l + 1, n]...)
    vct = vcat(alltiposets[:, k + 1, l + 1, n]...)
    resw = Array{Iposet}(undef, length(vcw))
    ress = Array{Iposet}(undef, length(vcs))
    rest = Array{Iposet}(undef, length(vct))
    @inbounds for i in 1:length(vcw)
        resw[i] = vcw[i][1]
    end
    @inbounds for i in 1:length(vcs)
        ress[i] = vcs[i][1]
    end
    @inbounds for i in 1:length(vct)
        rest[i] = vct[i][1]
    end
    return resw, ress, rest
end

function gpiWinkowskiMemo8(k, l)
    n = 8
    allwiposets = loadIposetsSplit("gpwiupto7.ips", n)
    allsiposets = loadIposetsSplit("gpsiupto7.ips", n)
    alltiposets = loadIposetsSplit("gptiupto7.ips", n)
    filled = zeros(Bool, n + 1, n + 1, n)
    for np in 1:7
        for kp in 0:np
            for lp in 0:np
                filled[kp + 1, lp + 1, np] = true
            end
        end
    end
    locks = Array{ReentrantLock}(undef, (n*(n-1))÷2 + 2, n + 1, n + 1, n)
    for i in eachindex(locks)
        locks[i] = ReentrantLock()
    end
    return gpwstiPosets(n, k, l, allwiposets, allsiposets, alltiposets, filled, locks)
end

function ggpwi8(k, l)
    wfname = string("gpwi8-",k,"-",l,".ips")
    sfname = string("gpsi8-",k,"-",l,".ips")
    tfname = string("gpti8-",k,"-",l,".ips")
    println(wfname)
    gpwi, gpsi, gpti = gpiWinkowskiMemo8(k, l)
    println("Done: W: ", length(gpwi), ", S: ", length(gpsi), ", T: ", length(gpti))
    saveIposets(gpwi, wfname)
    saveIposets(gpsi, sfname)
    saveIposets(gpti, tfname)
end

# 8-4-6 is wrong?!
# 8-1-6 is unclear
n = 8
k = 1
l = 6

ggpwi8(k, l)

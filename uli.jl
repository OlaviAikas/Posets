include("posets.jl")
include("iposets.jl")
#using StatProfilerHTML

"""Generate all gp-iposets with n nodes, k sources, l targets,
using files in gpiposets/ for memo initialization"""
function gpiPosetsMemo(n, k, l)
    alliposets = Array{Array{Tuple{Iposet, Array{Tuple{Int, Int}}}}}(undef, (n*(n-1))÷2 + 1, n + 1, n + 1, n)
    @inbounds for i in eachindex(alliposets)
        alliposets[i] = []
    end
    filled = zeros(Bool, n + 1, n + 1, n)
    locks = Array{ReentrantLock}(undef, (n*(n-1))÷2 + 2, n + 1, n + 1, n)
    for i in eachindex(locks)
        locks[i] = ReentrantLock()
    end
    for np in 1:7
        fname = string("gpiposets/gpi",np,".ips")
        ips = loadIposets(fname)
        for ip in ips
            ipe = ne(ip.poset)
            kp = length(ip.s)
            lp = length(ip.t)
            vprof = Array{Tuple{Int, Int}}(undef, np)
            @inbounds for v in 1:np
                vprof[v] = (inHash(ip.poset, v), outHash(ip.poset, v))
            end
            push!(alliposets[ipe + 1, kp + 1, lp + 1, np], (ip, vprof))
        end
        for kp in 0:np
            for lp in 0:np
                filled[kp + 1, lp + 1, np] = true
            end
        end
    end
    return gpiPosets(n, k, l, alliposets, filled, locks)
end


#saveIposets(gpiPosets(7), "gpi7.ips")
#ips = loadIposets("gpiposets/gpi0.ips")
#println(ips)

n = 8
k = 0
l = 1
fname = string("gpi",n,"-",k,"-",l,".ips")
println(fname)
ips = gpiPosetsMemo(n, k, l)
println(length(ips))
saveIposets(ips, fname)

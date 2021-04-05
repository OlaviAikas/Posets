include("posets.jl")
include("iposets.jl")
#using StatProfilerHTML

function count()
    alliposets = Array{Array{Iposet}}(undef, 8, 8, 7)
    @inbounds for i in eachindex(alliposets)
        alliposets[i] = []
    end
    ips = loadIposets("gpiposets/gpiupto7.ips")
    for ip in ips
        np = nv(ip.poset)
        kp = length(ip.s)
        lp = length(ip.t)
        push!(alliposets[kp + 1, lp + 1, np], ip)
    end
    for n in 1:7
        for k in 0:n
            for l in k:n
                println(n, k, l, " : ", length(alliposets[k+1, l+1, n]))
            end
        end
    end
end

count()





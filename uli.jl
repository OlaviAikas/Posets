include("posets.jl")
include("iposets.jl")
using StatProfilerHTML

"""Read a file of iposets and return them in an array"""
function readIposets(filename)
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

#@time genAllIposets(2)
#@time println(length(genAllIposets(4)))
#ip = readIposets("test.txt")
#println(ip[end])

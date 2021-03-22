include("posets.jl")
include("iposets.jl")
using StatProfilerHTML

"""Read a file of iposets and return them in an array"""
function readIposets(filename)
    res = Array{Iposet}(undef, 0)
    open(filename, "r") do file
        for line in eachline(file)
            nums = split(line)
            np = parse(Int, nums[1]) # number of points
            ne = parse(Int, nums[2]) # number of edges
            #println(nums, np, ne)
            dg = SimpleDiGraph(np)
            for edges_num in nums[3:ne+2]
                add_edge!(dg, parse(Int, edges_num[1]) + 1, parse(Int, edges_num[2]) + 1)
            end
            if length(nums) == ne+2 # no sources neither targets
                s, t = (), ()
            elseif length(nums) == ne+3 # sources, but no targets
                s = Tuple([parse(Int, x)+1 for x in nums[ne+3]])
                t = ()
            elseif length(nums) == ne+4 && nums[ne+3] == "-" # targets, but no sources
                s = ()
                t = Tuple([parse(Int, x)+1 for x in nums[ne+4]])
            else # sources and targets
                s = Tuple([parse(Int, x)+1 for x in nums[ne+3]])
                t = Tuple([parse(Int, x)+1 for x in nums[ne+4]])
            end
            push!(res, Iposet(s, t, dg))
    end
        return res
    end
end

"""Turn iposet into string representation a la McKay"""
function iposetToString(ip)
    np = nv(ip.poset)
    ne = length(edges(ip.poset))
    println(np, ne)
end


"""Write a file of iposets"""
function writeIposets(ips, filename)
    open(filename, "w") do file
        write(file, "Hello\n\nHell\n")
    end
end

ips = readIposets("test.txt")
println(length(ips), " ", ips[end])

iposetToString(ips[end])
#writeIposets(ip, "out.txt")





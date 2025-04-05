using Plots

include("kdTrees.jl")
using .kdTrees


# x = [8:-1:1...]
x = Vector{Float64}[]
num_pts = 128
for i = 1:num_pts 
    push!(x, rand(2))
end


tree = kdTree(x, num_leaf_pts = 4);


p = plotPartitions(tree)
p = plotPartitions(tree, watertight=false)

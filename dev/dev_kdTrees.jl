using Plots

using NeighborGraphs.kdTrees


# x = [8:-1:1...]
x = Vector{Float64}[]
num_pts = 1280
for i = 1:num_pts 
    push!(x, randn(2))
end

X = spatial2Data(x)


tree = kdTree(X, num_leaf_pts = 10);


p = treePlot(tree, plot_text=false)

p = partitionPlot(tree)
display(p)
# p = partitionPlot(tree, watertight=false)



 
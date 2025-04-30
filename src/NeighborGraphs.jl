# Graphs to support:
# - Epsilon-ball 
# - KNN
# - NUN/kNUN (double-cone)
# - Tau rule 
# - Beta-skeletons, Gabriel
# - Theta, Yao, Semi-Yao
# - Background obstruction field

# Multiclass/
# - SVM
# - Triangulations, Edge-tracing
# - KNN
module NeighborGraphs

    using SparseArrays

    import Base.getindex

    const DEFAULT_INDEX_SEARCH = true
    const DEFAULT_RADIUS_TOL = 1e-12
    const DEFAULT_DIST_CACHE_SPARSITY_THRESHOLD = 1000
    const DEFAULT_USE_DIST_CACHE = true


    # Submodules
    export GeometricPrimitives
    include("GeometricPrimitives.jl")
    using .GeometricPrimitives: SearchableGeometry, BoundingVolume, Ball, intersects, getIntersection, isContained

    export kdTrees
    include("kdTrees.jl")
    using .kdTrees: DataPoint, IndexRange, kdNode, kdTree, getContainingNode, getLeaves, search


    # Connectivity types:
    export Edge, UnweightedEdge, WeightedEdge
    export ConnectivityType, Edgelist, AdjacencyList, DirectedAdjacencyList, AdjacencyMatrix, WeightMatrix, DistanceMatrix

    # NeighborRule types:
    export NeighborRule, EpsilonBall, NNG, KNN, TauRule, NUN, kNUN, Gabriel, BetaSkeleton, Theta, Yao, SemiYao, SVM
    include("graphs.jl")

end
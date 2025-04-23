# Graphs to support:
# - Epsilon-ball 
# - kNN
# - NUN/kNUN (double-cone)
# - Tau rule 
# - Beta-skeletons, Gabriel
# - Theta, Yao, Semi-Yao
# - Background obstruction field

# Multiclass/
# - SVM
# - Triangulations, Edge-tracing
# - kNN
module NeighborGraphs

    # Data types
    # export EpsilonBall, NNG, kNN, TauRule, NUN, kNUN, Gabriel, BetaSkeleton, Theta, Yao, SemiYao, SVM
    #export Delaunay, Fiber, MulticlassSVM
    
    # Functions

    # Submodules
    export GeometricPrimitives, kdTrees

    # Submodules
    include("GeometricPrimitives.jl")
    using .GeometricPrimitives

    include("kdTrees.jl")
    using .kdTrees

    abstract type NeighborRule end

    # TODO: add help dialogs for each kind of graph that explain their parameters
    struct NNG <: NeighborRule end

    struct kNN <:NeighborRule 
        k::Integer
    end 
    struct EpsilonBall <:NeighborRule end
    struct TauRule <:NeighborRule end

    struct NUN <:NeighborRule  end
    struct kNUN <:NeighborRule  end

    struct Gabriel <:NeighborRule end
    struct BetaSkeleton <:NeighborRule  end

    struct Theta <:NeighborRule end
    struct Yao <: NeighborRule end
    struct SemiYao <:NeighborRule end


    # TODO: functionality for multiclass rules
    function getLocalNeighbors(pt::DP, tree::kdTree{T, VDP}, search_type::NeighborRule; directed_edges=false::Bool, index_search=false::Bool) where {T, DP<:DataPoint{T}, VDP<:Vector{DataPoint{T}}}
        
    end

    function getLocalNeighbors!(nbr_list::VDP, pt::DP, tree::kdTree{T, VDP}, search_type) where {T, DP<:DataPoint{T}, VDP<:Vector{DataPoint{T}}}
        push!(nbr_list, nbr)
    end


end
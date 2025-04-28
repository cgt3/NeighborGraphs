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
    using StaticArrays

    import Base.getindex

    const DEFAULT_INDEX_SEARCH = true
    const DEFAULT_RADIUS_TOL = 1e-12

    # Data types
    # export EpsilonBall, NNG, KNN, TauRule, NUN, kNUN, Gabriel, BetaSkeleton, Theta, Yao, SemiYao, SVM
    #export Delaunay, Fiber, MulticlassSVM
    
    # Functions

    # Submodules
    export GeometricPrimitives, kdTrees

    # Submodules
    include("GeometricPrimitives.jl")
    using .GeometricPrimitives

    include("kdTrees.jl")
    using .kdTrees

    abstract type ConnectivityType end # ===============================================
    abstract type Edge end

    # TODO: many graphs will have already calculated the distance/weight between nodes during
    #       their construction routines. Can these values be saved to avoid recomputing?
    struct UnweightedEdge{T} <: Edge
        src::T 
        dst::T
    end    
    
    struct UnweightedEdgeList{T} <:ConnectivityType
        edges::Vector{UnweightedEdge{T}}
    end

    struct WeightedEdge{T_v, T_weight} <: Edge
        src::T_v
        dst::T_v
        weight::T_weight
    end

    struct WeightedEdgeList{T_v, T_wt} <:ConnectivityType
        edges::Vector{WeightedEdge{T_v, T_wt}}
        weightFunc::Function
    end

    struct AdjacencyList{T} <:ConnectivityType
        nbrs::Vector{Vector{T}}
    end

    struct DirectedAdjacencyList{T} <: ConnectivityType
        outgoing::Vector{Vector{T}}
        incoming::Vector{Vector{T}}
    end

    struct AdjacencyMatrix <: ConnectivityType
        is_edge::Union{Matrix{Bool}, SparseArraysCSC}
        num_pts::Integer
        is_sparse::Bool
        is_directed::Bool
    end

    struct WeightMatrix <: ConnectivityType
        edge_weight::Union{Matrix, SparseArraysCSC}
        num_pts::Integer
        is_sparse::Bool
        is_directed::Bool
        weightFunc::Function
    end

    struct DistanceMatrix <: ConnectivityType
        edge_length::Union{Matrix, SparseArraysCSC}
        num_pts::Integer
        is_sparse::Bool
        is_directed::Bool
        p::Real
    end


    abstract type NeighborRule end # ===============================================

    # TODO: add help dialogs for each kind of graph that explain their parameters

    struct NNG <: NeighborRule 
        p::Real
        include_duplicates::Bool
    end
    NNG(; p = 2::Integer, include_duplicates=false::Bool) = NNG(p, include_duplicates)
    NNG(; n::Integer, include_duplicates=false::Bool) = NNG(floor(Int64, log10(n)), include_duplicates)

    struct KNN <:NeighborRule 
        k::Function
        p::Real
        include_duplicates::Bool
    end 
    KNN(k::Integer; p=2::Integer, include_duplicates=false::Bool) = KNN((x->k), p, include_duplicates)
    KNN(k::Integer; n::Integer, include_duplicates=false::Bool)   = KNN((x->k), floor(Int64, log10(n)), include_duplicates)

    struct EpsilonBall <:NeighborRule 
        epsilon::Function
        p::Real
        include_boundary::Bool
    end
    EpsilonBall(epsilon; p=2::Real, include_boundary=true::Bool)  = EpsilonBall((x->epsilon), p, include_boundary)
    EpsilonBall(epsilon; n::Integer, include_boundary=true::Bool) = EpsilonBall((x->epsilon), floor(Int64, log10(n)), include_boundary)

    struct TauRule <:NeighborRule end

    struct NUN <:NeighborRule  end
    struct kNUN <:NeighborRule  end

    struct Gabriel <:NeighborRule end
    struct BetaSkeleton <:NeighborRule  end

    struct Theta <:NeighborRule end
    struct Yao <: NeighborRule end
    struct SemiYao <:NeighborRule end

    # TODO: Multiclass Neighbor Rules: =====================


    # For building graphs
    function addNeighbor!(L::UnweightedEdgeList, current_pt::Integer, new_nbr::Integer, data::VecDP) where {T, VecDP<:Vector{DataPoint{T}}}
        push!(L.edges, UnweightedEdge(current_pt, new_nbr))
    end

    function addNeighbor!(L::WeightedEdgeList, current_pt::Integer, new_nbr::Integer, data::VecDP) where {T, VecDP<:Vector{DataPoint{T}}}
        push!(L.edges, WeightedEdge(current_pt, new_nbr, L.weightFunc(current_pt, new_nbr, data)))
    end

    function addNeighbor!(L::AdjacencyList, current_pt::Integer, new_nbr::Integer, data::VecDP) where {T, VecDP<:Vector{DataPoint{T}}}
        push!(L[current_pt], new_nbr)
    end

    function addNeighbor!(L::DirectedAdjacencyList, current_pt::Integer, new_nbr::Integer, data::VecDP) where {T, VecDP<:Vector{DataPoint{T}}}
        push!(L[current_pt].outgoing, new_nbr)
        push!(L[new_nbr].incoming, current_pt)
    end

    function addNeighbor!(A::AdjacencyMatrix, current_pt::Integer, new_nbr::Integer, data::VecDP) where {T, VecDP<:Vector{DataPoint{T}}}
        if A.is_directed
            A.is_edge[current_pt, new_nbr] = true
        else
            if new_nbr.y < current_pt
                A.is_edge[current_pt, new_nbr] = true
            else
                A.is_edge[new_nbr, current_pt] = true
            end
        end
    end

    function addNeighbor!(W::WeightMatrix, current_pt::Integer, new_nbr::Integer, data::VecDP) where {T, VecDP<:Vector{DataPoint{T}}}
        edge_weight = W.weightFunc(data[current_pt].x, data[new_nbr].x)
        if A.is_directed
            A.edge_length[current_pt, new_nbr.y] = edge_weight
        else
            if new_nbr.y < current_pt
                A.edge_length[current_pt, new_nbr.y] = edge_weight
            else
                A.edge_length[new_nbr.y, current_pt] = edge_weight
            end
        end
    end

    function addNeighbor!(D::DistanceMatrix, current_pt::Integer, new_nbr::Integer, data::VecDP) where {T, VecDP<:Vector{DataPoint{T}}}
        edge_length = norm(data[current_pt].x .- data[new_nbr].x, D.p)
        if A.is_directed
            A.edge_length[current_pt, new_nbr.y] = edge_length
        else
            if new_nbr.y < current_pt
                A.edge_length[current_pt, new_nbr.y] = edge_length
            else
                A.edge_length[new_nbr.y, current_pt] = edge_length
            end
        end
    end


    # Finding local neighbors: =============================================================
    # TODO: can T_int be restricted to being an integer type?
    function findLocalNeighbors(nbr_rule::NeighborRule, pt::DP, tree::kdTree{T, VecDP}; params, 
        index_search = DEFAULT_INDEX_SEARCH::Bool,
               T_int = Int64::Type ) where {T, DP<:DataPoint{T}, VecDP<:Vector{DataPoint{T}}}

        nbr_list = index_search ? T_int[] : DataPoint{typeof(pt.y)}[]
        findLocalNeighbors!(nbr_rule, nbr_list, pt, tree, params=params, index_search=index_search)
        return nbr_list
    end

    # Epsilon ball:
    function findLocalNeighbors!(epsilonBall::EpsilonBall, nbr_list::Union{VecInt, VecDP}, pt::DP, tree::kdTree{T, VDP}; params, index_search=DEFAULT_INDEX_SEARCH::Bool) where {T, DP<:DataPoint{T}, VecDP<:Vector{DataPoint{T}}, T_int::Integer, VecInt<:Vector{T_int}}
        ball = Ball(pt.x, epsilonBall.epsilon(pt.x), p=epsilonBall.p)
        neighbor_nodes = search(tree, ball, include_boundary=epsilonBall.include_boundary, index_search=true)
        for node_index_range in neighbor_nodes
            for j in node_index_range.first:node_index_range.last
                if isContained(ball, data[j].x, include_boundary=epsilonBall.include_boundary)
                    push!(nbr_list, index_search ? j : data[j])
                end
            end
        end
    end

    # NNG
    function findLocalNeighbors!(nng::NNG, nbr_list::Union{VecInt, VecDP}, pt::DP, tree::kdTree{T, VDP}; params, index_search = DEFAULT_INDEX_SEARCH::Bool ) where {T, DP<:DataPoint{T}, VecDP<:Vector{DataPoint{T}}, T_int::Integer, VecInt<:Vector{T_int}}
        containing_node = !hasproperty(params, :containing_node) ? getContainingNode(tree, pt, pt_tol=pt_tol) : params.containing_node
        if isnothing(containing_node)
            return
        end
        query_index = !hasproperty(params, :query_index) ? 0 : params.query_index
        rad_tol = !hasproperty(params, :rad_tol) ? DEFAULT_RADIUS_TOL : params.rad_tol


        # Find the distance to the closest point in this node
        min_dist = Inf
        if query_index != 1
            closest_pt = 1
        else
            closest_pt = 2
        end

        for j in containing_node.index_range.first:containing_node.index_range.last
            if j != query_index
                dist = norm(pt.x .- data[j].x, nng.p)
                if dist < min_dist
                    min_dist = dist 
                    closest_pt = j
                end
            end
        end

        # Search the tree for all points within a ball of radius min_dist
        ball = Ball(pt.x, min_dist, p=nng.p)
        neighboring_nodes = search(tree, ball, index_search=true, include_boundary=nng.include_duplicates)

        # Find the closest points
        candidate_pts = [closest_pt]
        for index_range in neighboring_nodes
            for j in index_range.first:index_range.last
                if j != pt_index 
                    dist = norm(pt.x .- data[j].x, nng.p)
                    if nng.include_duplicates && abs(dist - min_dist) < rad_tol
                        push!(candidate_pts, j)
                    elseif dist < min_dist
                        candidate_pts = [j]
                        min_dist = dist
                    end
                end
            end
        end

        # Return the closest point(s) in nbr_list
        for j in candidate_pts
            push!(nbr_list, index_search ? j : data[j])
        end
    end


    # KNN
    function findLocalNeighbors!(knn:KNN, nbr_list::Union{VecInt, VecDP}, pt::DP, tree::kdTree{T, VDP}; params, index_search = DEFAULT_INDEX_SEARCH::Bool ) where {T, DP<:DataPoint{T}, VecDP<:Vector{DataPoint{T}}, T_int::Integer, VecInt<:Vector{T_int}}
        containing_node = !hasproperty(params, :containing_node) ? getContainingNode(tree, pt, pt_tol=pt_tol) : params.containing_node
        if isnothing(containing_node)
            return
        end
        query_index = !hasproperty(params, :query_index) ? 0 : params.query_index
        rad_tol = !hasproperty(params, :rad_tol) ? DEFAULT_RADIUS_TOL : params.rad_tol

        # Find the distance to all points inside this node and sort
        distances = zeros(Float64, query_index ? node.index_range.n : node.index_range.n - 1)
        for j in node.index_range.first:node.index_range.last
            if j != query_index
                distances[i] = norm(pt.x .- data[j].x, knn.p)
            end
        end
        J_sorted = sortperm(distances)

        # If there are more than k points, set R to the distance to the kth nearest point; If
        # there are less than k points, set R to the distance to the furthest point in the BV.
        R = node.index_range.n > knn.k(pt.x) ? distances[J_sorted[knn.k(pt.x)]] : norm(pt .- getFurthestPoint(node.bv_watertight, pt), knn.p)

        # Search the tree for points within a ball of radius R 
        ball = Ball(pt, R, p=knn.p)
        neighboring_nodes = search(tree, ball, include_boundary=knn.include_duplicates, index_search=true)

        # For points that are within R, sort and take the first k
        candidate_pts = Int64[]
        distances = Float64[]
        for index_range in neighboring_nodes
            for j in index_range.first:index_range.last
                dist = norm(pt .- data[j].x, knn.p)
                if dist < ball.radius || ( knn.include_duplicates && dist < ball.radius + rad_tol )
                    push!(candidate_pts, j) 
                    push!(distances, dsit)
                end
            end
        end

        J_sorted = sortperm(distances)
        if knn.include_duplicates
            j0 = J_sorted[knn.k(pt.x)]
            j_end = j0 + 1
            while j_end < length(J_sorted) && abs(distances[J_sorted[j0]] - distances[J_sorted[j_end]]) < rad_tol
                j_end += 1
            end
        else
            j_end = knn.k(pt.x)
        end

        for j = 1:j_end 
            push!(nbr_list, index_search ? J_sorted[j] : data[J_sorted[j]])
        end
    end


    # Building graphs ========================================================================================
    struct Graph{T_int<:Integer, T_data, T_rule<:NeighborRule, T_connectivity<:ConnectivityType}
        dim::T_int
        tree::kdTree{T_data, Vector{DataPoint{T_data}}}
        nbr_rule::T_rule 
        connectivity::T_connectivity

        function Graph(data::VecDP, nbr_rule::NeighborRule, connectivity::ConnectivityType;
            num_leaf_pts = DEFAULT_NUM_LEAF_PTS::Integer ) where {T, VecDP<:Vector{DataPoint{T}}}
  
            dim = length(data[1].x)
          
            # Build the k-d tree for fast(er) search operations:
            tree = kdTree(data, num_leaf_pts=num_leaf_pts)
  
            # Build the graph
            buildGraph!(connectivity, tree, nbr_rule)

            return new(dim, tree, nbr_rule, connectivity)
        end
    end 


    function buildGraph!(connectivity::ConnectivityType, epsilonBall::EpsilonBall, tree::kdTree) where {T, VecDP<:Vector{DataPoint{T}}}
        # For each data point, find its local/outgoing neighbors
        for i in eachindex(tree.data)
            nbrs = findLocalNeighbors(epsilonBall, tree.data[i], tree, index_search=true)
            for j in nbrs
                addNeighbor!(connectivity, i, j, tree.data)
            end
        end
    end

    # NNG and kNNs
    function buildGraph!(connectivity::ConnectivityType, nbr_rule::Union{NNG, KNN}, tree::kdTree; params) where {T, VecDP<:Vector{DataPoint{T}}}
        rad_tol = !hasproperty(params, :rad_tol) ? DEFAULT_RADIUS_TOL : params.rad_tol

        # NOTE: points are ordered by node, so process them by node to avoid unnecessary searches in the tree
        leaf_nodes = getLeaves(tree)
        for node in leaf_nodes
            for i in node.index_range.first:node.index_range.last
                params = (; query_index=i, containing_node=node, rad_tol)
                nbrs = findLocalNeighbors(nng, tree.data[i], tree, params=params, index_search=true)
                for j in nbrs 
                    addNeighbor!(connectivity, i, j, tree.data)
                end
            end # points in each node
        end # leaf nodes
    end

end
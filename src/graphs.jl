
abstract type ConnectivityType end # ===============================================
abstract type AbstractEdge end

# TODO: many graphs will have already calculated the distance/weight between nodes during
#       their construction routines. Can these values be saved to avoid recomputing?
struct UnweightedEdge{T} <: AbstractEdge
    src::T 
    dst::T
end

struct WeightedEdge{T_pt, T_wt} <: AbstractEdge 
    src::T_pt
    dst::T_pt 
    weight::T_wt
end

struct EdgeList{T} <:ConnectivityType where T<:AbstractEdge
    edges::Vector{T}
    is_weighted::Bool
    weightFunc::Union{Function, Nothing}
end

struct AdjacencyList{T} <:ConnectivityType
    nbrs::Vector{Vector{T}}
end

struct DirectedAdjacencyList{T} <: ConnectivityType
    outgoing::Vector{Vector{T}}
    incoming::Vector{Vector{T}}
end

struct AdjacencyMatrix <: ConnectivityType
    is_edge::Union{Matrix{Bool}, SparseMatrixCSC}
    num_pts::Integer
    is_sparse::Bool
    is_directed::Bool
end

struct WeightMatrix <: ConnectivityType
    edge_weight::Union{Matrix, SparseMatrixCSC}
    num_pts::Integer
    is_sparse::Bool
    is_directed::Bool
    weightFunc::Function
end

struct DistanceMatrix <: ConnectivityType
    edge_length::Union{Matrix, SparseMatrixCSC}
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
# NNG(; n::Integer, include_duplicates=false::Bool) = NNG(floor(Int64, log10(n)), include_duplicates)

struct KNN <:NeighborRule 
    k::Function
    p::Real
    include_duplicates::Bool
end 
KNN(k::Integer; p=2::Integer, include_duplicates=false::Bool) = KNN((x->k), p, include_duplicates)
# KNN(k::Integer; n::Integer, include_duplicates=false::Bool)   = KNN((x->k), floor(Int64, log10(n)), include_duplicates)

struct EpsilonBall <:NeighborRule 
    epsilon::Function
    p::Real
    include_boundary::Bool
end
EpsilonBall(epsilon; p=2::Real, include_boundary=true::Bool)  = EpsilonBall((x->epsilon), p, include_boundary)
# EpsilonBall(epsilon; n::Integer, include_boundary=true::Bool) = EpsilonBall((x->epsilon), floor(Int64, log10(n)), include_boundary)

struct TauRule <:NeighborRule end

struct NUN <:NeighborRule  end
struct kNUN <:NeighborRule  end

struct Gabriel <:NeighborRule end
struct BetaSkeleton <:NeighborRule  end

struct Theta <:NeighborRule end
struct Yao <: NeighborRule end
struct SemiYao <:NeighborRule end

struct SVM <: NeighborRule end

# TODO: Multiclass Neighbor Rules: =====================


# For building graphs
function addNeighbor!(L::EdgeList, current_pt::Integer, new_nbr::Integer, data::VecDP) where {T, VecDP<:Vector{DataPoint{T}}}
    if L.is_weighted
        push!(L.edges, WeightedEdge(current_pt, new_nbr, L.weightFunc(data[current_pt], data[new_nbr])) )
    else
        push!(L.edges, UnweightedEdge(current_pt, new_nbr))
    end
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

# Functions for efficiency in high-dim spaces and other helper functions: -----------------------------------------------------------
function distanceCache(n::Integer, cache_distances::Bool; sparse_cache=(n > DEFAULT_DIST_CACHE_SPARSITY_THRESHOLD)::Bool)
    if !cache_distances
        return nothing 
    end

    if sparse_cache
        return spzeros(Float64, n, n)
    else
        return zeros(Float64, n, n)
    end
end

function getDistance!(i1::Vector, i2::Vector, p1::Vector, p2::Vector, distances::Union{SparseMatrixCSC, Matrix, Nothing}; p=2::Real)
    if i1 == 0 || isnothing(distances)
        return norm(p1 .- p2, p)
    end

    if i1 < i2
        swap = i2
        i2 = i1
        i1 = swap 
    end

    if i1 == i2 || distances[i1, i2]
        return 0
    elseif distances[i1, i2] != 0
        return distances[i1, i2]
    else
        distances[i1, i2] = norm(data[i1] .= data[i2], p)
        return distances[i1, i2]
    end
end

function getQueryIndex(params)
    return !hasproperty(params, :query_index) ? 0 : params.query_index
end

function getRadTol(params)
    return !hasproperty(params, :rad_tol) ? DEFAULT_RADIUS_TOL : params.rad_tol
end



# Finding local neighbors: ==========================================================================================
# TODO: can T_int be restricted to being an integer type?
function findLocalNeighbors(nbr_rule::NeighborRule, pt::DP, tree::kdTree{T, VecDP}; params, 
       index_search = DEFAULT_INDEX_SEARCH::Bool,
    known_distances = nothing::Union{SparseMatrixCSC, Matrix, Nothing},
           T_int = Int64::Type ) where {T, DP<:DataPoint{T}, VecDP<:Vector{DataPoint{T}}}

    nbr_list = index_search ? T_int[] : DataPoint{typeof(pt.y)}[]
    findLocalNeighbors!(nbr_rule, nbr_list, pt, tree, params=params, index_search=index_search, known_distances=known_distances)
    return nbr_list
end

# Epsilon ball:
function findLocalNeighbors!(epsilonBall::EpsilonBall, nbr_list::Union{VecInt, VecDP}, pt::T_DP, tree::kdTree{T, VecDP}; params, 
       index_search = DEFAULT_INDEX_SEARCH::Bool,
    known_distances = nothing::Union{SparseMatrixCSC, Matrix, Nothing} ) where {T, T_DP<:DataPoint{T}, VecDP<:Vector{DataPoint{T}}, T_int<:Integer, VecInt<:Vector{T_int}}
    
    query_index = getQueryIndex(params)

    ball = Ball(pt.x, epsilonBall.epsilon(pt.x), p=epsilonBall.p)
    neighbor_nodes = search(tree, ball, include_boundary=epsilonBall.include_boundary, index_search=true)
    for node_index_range in neighbor_nodes
        for j in node_index_range.first:node_index_range.last
            if j != params.query_index 
                R = getDistance!(query_index, j, ball.center, data[j].x, known_distances, p=ball.p)
                if R < ball.radius || (epsilon_ball.include_boundary && R == ball.radius)
                    push!(nbr_list, index_search ? j : data[j])
                end
            end
        end
    end
end

# NNG
function findLocalNeighbors!(nng::NNG, nbr_list::Union{VecInt, VecDP}, pt::DP, tree::kdTree{T, VecDP};  params,
       index_search = DEFAULT_INDEX_SEARCH::Bool, 
    known_distances = nothing::Union{SparseMatrixCSC, Matrix, Nothing} ) where {T, DP<:DataPoint{T}, VecDP<:Vector{DataPoint{T}}, T_int<:Integer, VecInt<:Vector{T_int}}
    
    containing_node = !hasproperty(params, :containing_node) ? getContainingNode(tree, pt, pt_tol=pt_tol) : params.containing_node
    if isnothing(containing_node)
        return
    end
    query_index = getQueryIndex(params)
    rad_tol = getRadTol(params)

    # Find the distance to the closest point in this node
    min_dist = Inf
    closest_pt = query_index != 1 ? 1 : 2
    for j in containing_node.index_range.first:containing_node.index_range.last
        if j != query_index
            dist = getDistance!(query_index, j, pt.x, data[j].x, known_distances, p=nng.p)
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
                dist = getDistance!(query_index, j, pt.x, data[j].x, known_distances, p=nng.p)
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
function findLocalNeighbors!(knn::KNN, nbr_list::Union{VecInt, VecDP}, pt::DP, tree::kdTree{T, VecDP}; params, 
       index_search = DEFAULT_INDEX_SEARCH::Bool, 
    known_distances = nothing::Union{SparseMatrixCSC, Matrix, Nothing} ) where {T, DP<:DataPoint{T}, VecDP<:Vector{DataPoint{T}}, T_int<:Integer, VecInt<:Vector{T_int}}

    containing_node = !hasproperty(params, :containing_node) ? getContainingNode(tree, pt, pt_tol=pt_tol) : params.containing_node
    if isnothing(containing_node)
        return
    end
    query_index = getQueryIndex(params)
    rad_tol = getRadTol(params)

    # Find the distance to all points inside this node and sort
    distances = zeros(Float64, query_index ? node.index_range.n : node.index_range.n - 1)
    for j in node.index_range.first:node.index_range.last
        if j != query_index
            distances[i] = getDistance!(query_index, j, pt.x, data[j].x, known_distanes, p=knn.p)
        end
    end
    J_sorted = sortperm(distances)

    # If there are more than k points, set R to the distance to the kth nearest point; If
    # there are less than k points, set R to the distance to the furthest point of the BV.
    R = node.index_range.n > knn.k(pt.x) ? distances[J_sorted[knn.k(pt.x)]] : norm(pt .- getFurthestPoint(node.bv_watertight, pt), knn.p)

    # Search the tree for points within a ball of radius R 
    ball = Ball(pt, R, p=knn.p)
    neighboring_nodes = search(tree, ball, include_boundary=knn.include_duplicates, index_search=true)

    # For points that are within R, sort and take the first k
    candidate_pts = Int64[]
    distances = Float64[]
    for index_range in neighboring_nodes
        for j in index_range.first:index_range.last
            dist = getDistance!(query_index, j, pt.x, data[j].x, known_distances, p=knn.p)
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
struct Graph{T_data, T_rule<:NeighborRule, T_connectivity<:ConnectivityType}
    dim::Integer
    tree::kdTree{T_data, Vector{DataPoint{T_data}}}
    nbr_rule::T_rule 
    connectivity::T_connectivity
    cache_distances::Bool
    known_distances::Union{SparseMatrixCSC, Matrix, Nothing}

    function Graph(data::VecDP, nbr_rule::NeighborRule, connectivity::ConnectivityType; params,
           num_leaf_pts = DEFAULT_NUM_LEAF_PTS::Integer,
        cache_distances = DEFAULT_USE_DIST_CACHE::Bool,
           sparse_cache = (length(data[1].x) > DEFAULT_DIST_CACHE_SPARSITY_THRESHOLD)::Bool ) where {T, VecDP<:Vector{DataPoint{T}}}

        dim = length(data[1].x)
        
        # Build the k-d tree for fast(er) search operations:
        tree = kdTree(data, num_leaf_pts=num_leaf_pts)

        # If caching is to be used, set up the matrix of distance values
        known_distances = distanceCache(dim, cache_distances, sparse_cache=sparse_cache)

        # Build the graph
        buildGraph!(connectivity, tree, nbr_rule, params=params, known_distances=known_distances)

        return new{typeof(data[1].y), typeof(nbr_rule), typeof(connectivity) }(dim, tree, nbr_rule, connectivity, cache_distances, known_distances)
    end
end 


function buildGraph!(connectivity::ConnectivityType, epsilonBall::EpsilonBall, tree::kdTree; params,
    known_distances = nothing::Union{SparseMatrixCSC, Matrix, Nothing} ) where {T, VecDP<:Vector{DataPoint{T}}}

    # For each data point, find its local/outgoing neighbors
    for i in eachindex(tree.data)
        nbrs, known_distances = findLocalNeighbors(epsilonBall, tree.data[i], tree, index_search=true, known_distances=known_distances)
        for j in nbrs
            addNeighbor!(connectivity, i, j, tree.data)
        end
    end
end

# NNG and kNNs
function buildGraph!(connectivity::ConnectivityType, nbr_rule::Union{NNG, KNN}, tree::kdTree; params,
    known_distances = nothing::Union{SparseMatrixCSC, Matrix, Nothing} )

    rad_tol = getRadTol(params)

    # NOTE: points are ordered by node, so process them by node to avoid unnecessary searches in the tree
    leaf_nodes = getLeaves(tree)
    for node in leaf_nodes
        for i in node.index_range.first:node.index_range.last
            params = (; query_index=i, containing_node=node, rad_tol)
            nbrs = findLocalNeighbors(nng, tree.data[i], tree, params=params, index_search=true, known_distances=known_distances)
            for j in nbrs 
                addNeighbor!(connectivity, i, j, tree.data)
            end
        end # points in each node
    end # leaf nodes
end
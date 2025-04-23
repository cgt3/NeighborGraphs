module kdTrees

    using ..GeometricPrimitives: Ball, BoundingVolume, intersects, isContained, getIntersection

    using Plots
    using Printf

    const DEFAULT_NUM_LEAF_PTS = 40
    const DEFAULT_PT_TOL = 1e-12

    const DEFAULT_PLOT_WIDTH_PER_NODE = 100
    const DEFAULT_PLOT_HEIGHT_PER_NODE = 200
    const DEFAULT_LEAF_FONTSIZE = 8
    const DEFAULT_MAX_TREE_PLOT_FONTSIZE = 24

    const MIN_POINT_PLOT_SIZE = 2
    const MAX_POINT_PLOT_SIZE = 16
    const MIN_SPLIT_LINE_WIDTH = 3
    const MAX_SPLIT_LINE_WIDTH = 10
    const TREE_MARKER_SIZE = 20


    const DEFAULT_PARTITION_PALETTE = Plots.palette([:darkred,       :orangered3,  :darkorange2, :goldenrod1,
                                               :chartreuse3,   :forestgreen, :darkgreen,      #:olivedrab2
                                               :navy,          :blue,        :deepskyblue,    #:blue
                                               :mediumpurple1, :purple3,     :purple4  ])     #:purple1



    # Data types
    export IndexRange, DataPoint, kdNode, kdTree

    # Functions
    export spatial2Data, getSplit, contains, getContainingNode, getLeaves, partitionPlot, treePlot
   
   
    # TODO: there is most definitely something in Julia that replicates the below, but if
    #     there is true need to also have n and not just i0,in, than it may be ok to keep.
    struct IndexRange
        first::Integer
        last::Integer
        n::Integer

        IndexRange(first; n) = new(first, first + n - 1, n)
        IndexRange(; first, last) = new(first, last, last - first + 1)
    end

    struct DataPoint{T}
        x::AbstractVector
        y::T
    end

    function spatial2Data(data_num::Vector{Vector{T}}) where T <: Real
        return DataPoint.(data_num, nothing)
    end


    function Base.getindex(p::DataPoint{T}, i::Integer) where T
        return p.x[i]
    end

    function BoundingVolume(data::VDP, index_range::IndexRange) where {T, VDP<:Vector{DataPoint{T}}}
        ub = copy(data[index_range.first].x)
        lb = copy(data[index_range.first].x)
        for i = index_range.first:index_range.last
            d_ub = ub .< data[i].x
            ub[d_ub] = data[i].x[d_ub]

            d_lb = lb .> data[i].x
            lb[d_lb] = data[i].x[d_lb]
        end
        return BoundingVolume(lb, ub)
    end


    # Convention: use <= goes to the left, > goes to the right
    mutable struct kdNode{T, VDP<:Vector{DataPoint{T}}}
        data::VDP
        original_indices::Vector{Int64}
        is_leaf::Bool
        index_range::IndexRange
        split_dim::Integer
        split_val::Real
        bv::BoundingVolume
        bv_watertight::BoundingVolume
        parent::Union{kdNode, Nothing}
        l_child::Union{kdNode, Nothing}
        r_child::Union{kdNode, Nothing}
        depth::Integer

        # Default constructor
        kdNode(data::VDP, original_indices, is_leaf, index_range, split_dim, split_val, bv, bv_wt, parent, l_child, r_child, depth) where {T, VDP<:Vector{DataPoint{T}}} = 
            new{typeof(data[index_range.first].y), typeof(data)}(data, original_indices, is_leaf, index_range, split_dim, split_val, bv, bv_wt, parent, l_child, r_child, depth) 

        # The recursive, top-level constructor
        function kdNode(data::VDP; original_indices=[1:length(data)...]::Vector{Int64}, 
                  parent = nothing::Union{kdNode, Nothing}, 
             index_range = IndexRange(1, n=length(data))::IndexRange, 
                   depth = 0::Integer, 
            num_leaf_pts = DEFAULT_NUM_LEAF_PTS::Integer) where {T, VDP<:Vector{DataPoint{T}}}

            if index_range.n == 0
                return nothing
            elseif index_range.n <= num_leaf_pts 
                bv = BoundingVolume(data, index_range)
                bv_wt = getWatertightBoundingVolume(parent, data, index_range)
                return kdNode(data, original_indices, true, index_range, 0, NaN, bv, bv_wt, parent, nothing, nothing, depth)
            end

            node = new{typeof(data[index_range.first].y), typeof(data)}()
            node.data = data
            node.original_indices = original_indices
            node.is_leaf = false
            node.index_range = index_range
            node.depth=depth

            node.bv = BoundingVolume(data, index_range)
            node.bv_watertight = getWatertightBoundingVolume(parent, data, index_range)

            # Find the splitting dimension and value
            node.split_dim, node.split_val = getSplit(node.bv, node.bv_watertight, data, index_range)

            # Partition the data
            i_split = partitionData!(data, original_indices, index_range, node.split_dim, node.split_val)

            # For the children nodes
            i0_left  = node.index_range.first
            i0_right = i_split + 1
            n_left = i_split - i0_left + 1
            n_right = node.index_range.n - n_left 

            # Attach family nodes
            node.parent = parent
            node.l_child = kdNode(data, original_indices=original_indices, parent=node, index_range=IndexRange(i0_left,  n=n_left),  depth=depth+1, num_leaf_pts=num_leaf_pts)
            node.r_child = kdNode(data, original_indices=original_indices, parent=node, index_range=IndexRange(i0_right, n=n_right), depth=depth+1, num_leaf_pts=num_leaf_pts)
            return node
        end
    end # kdNode struct declaration


    function getCentroid(bv::BoundingVolume, data::VDP, index_range::IndexRange) where {T, VDP<:Vector{DataPoint{T}}}
        x0 = data[1].x
        centroid = zeros(eltype(x0), size(x0))
        for i in index_range.first:index_range.last
            centroid += data[i].x
        end
        centroid *= 1/index_range.n

        return centroid
    end

    function getWatertightBoundingVolume(parent::Union{kdNode, Nothing}, data, index_range::IndexRange)
        if isnothing(parent)
            return BoundingVolume(data, index_range)
        else
            parent_bv = parent.bv_watertight
            if parent.index_range.first == index_range.first # we are constructing the left child
                bv_ub = copy(parent_bv.ub)
                bv_ub[parent.split_dim] = parent.split_val
                return BoundingVolume(copy(parent_bv.lb), bv_ub)
            else
                bv_lb = copy(parent_bv.lb)
                bv_lb[parent.split_dim] = parent.split_val
                return BoundingVolume(bv_lb, copy(parent_bv.ub))
            end
        end
    end

    # TODO: change this function to allow different schemes for choosing the split value
    function getSplit(bv::BoundingVolume, bv_wt::BoundingVolume, data, index_range::IndexRange)
        centroid = getCentroid(bv, data, index_range)
        # mid = 0.5 * (bv.ub - bv.lb) + bv.lb
        # dist = centroid - mid

        # # Choose the dimension where the centroid is closest to the midpoint
        # split_dim = argmin(dist)
        # split_val = centroid[split_dim]

        # Choose the dimension that splits the watertight BV into the most square pieces
        bv_length = bv_wt.ub .- bv_wt.lb
        split_dim = argmax(bv_length)
        split_val = centroid[split_dim]

        return split_dim, split_val
    end


    function partitionData!(data::VDP, original_indices, index_range::IndexRange, split_dim, split_val) where {T, VDP<:Vector{DataPoint{T}}}
        i_start = index_range.first
        i_end = index_range.last

        while i_start < i_end 
            is_valid_L = data[i_start][split_dim] <= split_val
            is_valid_R = data[i_end][split_dim] > split_val

            if !is_valid_L && !is_valid_R
                swap          = data[i_end]
                data[i_end]   = data[i_start]
                data[i_start] = swap

                swap = original_indices[i_end]
                original_indices[i_end]   = original_indices[i_start]
                original_indices[i_start] = swap
            end

            if is_valid_L && !is_valid_R
                i_start += 1
            elseif !is_valid_L && is_valid_R
                i_end -= 1
            else
                i_start += 1
                i_end -=1
            end
        end

        if data[i_start][split_dim] > split_val
            i_start -=1
        end

        return i_start
    end

    struct kdTree{T, VDP<:Vector{DataPoint{T}}}
        data::VDP
        original_indices::Vector{Int64}
        root::kdNode
        max_depth::Integer
        num_leaves::Integer
        num_pts::Integer

        function kdTree(data::VDP; num_leaf_pts::Integer=DEFAULT_NUM_LEAF_PTS,
            original_indices = [1:length(data)...]::Vector{Int64},
            root = kdNode(data, original_indices=original_indices, num_leaf_pts=num_leaf_pts) ) where {T, VDP<:Vector{DataPoint{T}}}

            max_depth, num_leaves, num_pts = getTreeStats(root)
            return new{typeof(data[root.index_range.first].y), typeof(data)}(data, original_indices, root, max_depth, num_leaves, num_pts)
        end
    end

    function getTreeStats(node::Union{kdNode, Nothing})
        if isnothing(node)
            return 0, 0, 0
        elseif node.is_leaf
            return node.depth, 1, node.index_range.n
        else
            depth_L, num_leaves_L, num_pts_L = getTreeStats(node.l_child)
            depth_R, num_leaves_R, num_pts_R = getTreeStats(node.r_child)

            return max(depth_L, depth_R), num_leaves_L + num_leaves_R, num_pts_L + num_pts_R
        end
    end

    # User-exposed functions -------------------------------------------------------------
    function getContainingNode(tree::kdTree, query_pt; pt_tol=DEFAULT_NUM_LEAF_PTS)
        return getContainingNode(tree.root, query_pt, pt_tol=pt_tol)
    end

    function getContainingNode(node::Union{kdNode, Nothing}, query_pt; pt_tol=DEFAULT_PT_TOL)
        if isnothing(kd_node)
            return nothing
        elseif node.depth == 0 # Check if the point in question is even in the tree's BV
            if any(query_pt .- node.bv_watertight.ub .> pt_tol) || any( query_pt .- node.bv_watertight.lb .< pt_tol)
                return nothing
            end
        end

        if node.is_leaf == true
            return node
        elseif query_pt[kd_node.split_dim] <= node.split_val
            return getContainingNode(node.l_child, query_pt, pt_tol=pt_tol)
        else
            return getContainingNode(node.r_child, query_pt, pt_tol=pt_tol)
        end
    end

    function contains(tree::kdTree, query_pt; pt_tol=DEFAULT_PT_TOL)
        node = getContainingNode(query_pt, tree.root, pt_tol=pt_tol)
        return contains(node, query_pt, pt_tol=pt_tol)
    end

    function contains(node::Union{kdNode, Nothing}, query_pt; pt_tol=DEFAULT_PT_TOL)
        if isnothing(node)
            return false
        end

        for i = node.index_range.first:node.index_range.last
            if all(abs.(query_pt .- node.data[i].x)) .< pt_tol
                return true
            end
        end

        return false
    end

    function treePlot(tree::kdTree;
        x_spacing=DEFAULT_PLOT_WIDTH_PER_NODE,
        y_spacing=DEFAULT_PLOT_HEIGHT_PER_NODE,
        plot_text=true,
        leaf_fontsize=DEFAULT_LEAF_FONTSIZE,
        max_fontsize=DEFAULT_MAX_TREE_PLOT_FONTSIZE)

        p = plot(
            size=( x_spacing*(2^tree.max_depth + 1), y_spacing*(tree.max_depth + 3) ), 
            ylims=(-1, tree.max_depth + 1), 
            xlims=(0, 2^tree.max_depth + 1),
            axis=([], false),
            legend=false,
            rmargin=0Plots.px,
            lmargin=0Plots.px
            )
        yflip!(true)
        treePlot!(p, 1, tree.max_depth, tree.root, plot_text=plot_text, leaf_fontsize=leaf_fontsize, max_fontsize=max_fontsize)
        return p
    end

    function treePlot!(p, i_level::Integer, max_depth::Integer, node::Union{kdNode, Nothing}; 
        plot_text=true,
        leaf_fontsize=DEFAULT_TREE_PLOT_FONTSIZE,
        max_fontsize=DEFAULT_MAX_TREE_PLOT_FONTSIZE )

        x_pos = 0.5 + (i_level - 0.5)*( 2^(max_depth - node.depth) )
        y_pos = node.depth

        fontsize = min(leaf_fontsize + 2*(max_depth - node.depth), max_fontsize)

        if isnothing(node)
            if plot_text
                annotate!(p, x_pos, y_pos, text("(empty leaf)", :center, fontsize))
            else
                scatter!(p, [x_pos], [y_pos], markercolor=:gray80, markerstrokewidth=0, markersize=TREE_MARKER_SIZE)
            end
        elseif !node.is_leaf
            # Plot lines to child nodes
            x_spacing = 2^(max_depth - node.depth)
            x_pos_L, x_pos_R = x_pos - x_spacing/4 , x_pos + x_spacing / 4
            y_pos_child = y_pos + 1
            plot!(p, [x_pos, x_pos_L], [y_pos, y_pos_child], linecolor=:black, linewidth=3)
            plot!(p, [x_pos, x_pos_R], [y_pos, y_pos_child], linecolor=:black, linewidth=3)

            # Add text for this node
            if plot_text
                annotate!(p, x_pos, y_pos, text(@sprintf("Split dim: %d\nVal: %.4e", node.split_dim, node.split_val), :center, fontsize))
            else
                scatter!(p, [x_pos], [y_pos], markercolor=:black, markerstrokewidth=0, markersize=TREE_MARKER_SIZE)
            end

            treePlot!(p, 2*i_level - 1, max_depth, node.l_child, plot_text=plot_text, leaf_fontsize=leaf_fontsize, max_fontsize=max_fontsize)
            treePlot!(p, 2*i_level,     max_depth, node.r_child, plot_text=plot_text, leaf_fontsize=leaf_fontsize, max_fontsize=max_fontsize)
        else
            if plot_text
                annotate!(p, x_pos, y_pos, text("Leaf: i0=$(node.index_range.first),\n n=$(node.index_range.n)", :center, fontsize))
            else
                scatter!(p, [x_pos], [y_pos], markercolor=:red, markerstrokewidth=0, markersize=TREE_MARKER_SIZE)
            end
        end
    end

    function partitionPlot(tree::kdTree; 
        size=(1500,1000), 
        fontsize=20,
        index=(1,2), 
        watertight=true,
        color_palette=DEFAULT_PARTITION_PALETTE,
        linewidth=4)

        p = plot(size=size,
            xtickfont=font(fontsize), 
            ytickfont=font(fontsize),
            leg=:false
        )
        partitionPlot!(p, tree.root, watertight=watertight, index=index, color_palette=color_palette, linewidth=linewidth)
        return p
    end

    # TODO: adapt to plot 1D points as well
    function partitionPlot!(p, node::Union{kdNode, Nothing}; 
            index=(1,2),
            watertight=true,
            color_palette=DEFAULT_PARTITION_PALETTE,
            linewidth=4)

        if isnothing(node)
            return
        elseif node.is_leaf
            I = node.index_range.first : node.index_range.last
            markersize = max(MAX_POINT_PLOT_SIZE - 1.5*node.depth, MIN_POINT_PLOT_SIZE)
            scatter!(p, getindex.(node.data[I], index[1]), getindex.(node.data[I], index[2]), 
                markercolor=:darkgrey, markerstrokewidth=0, markersize=markersize, markeralpha=0.75)
        else
            i_color = node.depth % length(color_palette) + 1 # Note on indexing: node depth is already 0-indexed

            if watertight 
                split_bv = node.bv_watertight
            else
                split_bv = node.bv 
            end

            partitionPlot!(p, node.l_child, watertight=watertight, index=index, color_palette=color_palette, linewidth=linewidth)
            partitionPlot!(p, node.r_child, watertight=watertight, index=index, color_palette=color_palette, linewidth=linewidth)

            if !watertight || node.depth == 0
                # Plot the parent BV:
                plot!(p, node.bv.lb[index[1]]*[1, 1], [node.bv.lb[index[2]], node.bv.ub[index[2]]], linecolor=color_palette[i_color], linestyle=:dash, linewidth=linewidth/2) # left
                plot!(p, node.bv.ub[index[1]]*[1, 1], [node.bv.lb[index[2]], node.bv.ub[index[2]]], linecolor=color_palette[i_color], linestyle=:dash, linewidth=linewidth/2) # right
                plot!(p, [node.bv.lb[index[1]], node.bv.ub[index[1]]], node.bv.lb[index[2]]*[1, 1], linecolor=color_palette[i_color], linestyle=:dash, linewidth=linewidth/2) # bottom
                plot!(p, [node.bv.lb[index[1]], node.bv.ub[index[1]]], node.bv.ub[index[2]]*[1, 1], linecolor=color_palette[i_color], linestyle=:dash, linewidth=linewidth/2) # top
            end

            # Plot the splitting plane
            if node.split_dim == index[1]
                plot!(p, node.split_val*[1, 1], [split_bv.lb[index[2]], split_bv.ub[index[2]]], linecolor=color_palette[i_color], linewidth=max(MIN_SPLIT_LINE_WIDTH, MAX_SPLIT_LINE_WIDTH - node.depth))
            elseif node.split_dim == index[2]
                plot!(p, [split_bv.lb[index[1]], split_bv.ub[index[1]]], node.split_val*[1, 1], linecolor=color_palette[i_color], linewidth=max(MIN_SPLIT_LINE_WIDTH, MAX_SPLIT_LINE_WIDTH - node.depth))
            end
        end
    end


    function getLeaves(tree::kdTree; index_search=false::Bool)
        if index_search
            leaves = IndexRange[]
        else
            leaves = kdNode[]
        end
        getLeaves!(leaves, tree.root,index_search=index_search)

        return leaves
    end

    function getLeaves!(leaves, node::Union{kdNode, Nothing}; index_search=false::Bool)
        if isnothing(node)
            return
        elseif node.is_leaf
            if index_search
                push!(leaves, node.index_range)
            else
                push!(leaves, node)
            end
        else
            getLeaves!(leaves, node.l_child, index_search=index_search)
            getLeaves!(leaves, node.r_child, index_search=index_search)
        end
    end

    # Non-point search related functions
    mutable struct SearchNode
        search::Union{Bool, Nothing}
        parent::Union{SearchNode, Nothing}
        l_child::Union{SearchNode, Nothing}
        r_child::Union{SearchNode, Nothing}

        SearchNode() = new(nothing, nothing, nothing, nothing)

        SearchNode(search::Bool, parent::Union{SearchNode, Nothing}, l_child::Union{SearchNode, Nothing}, r_child::Union{SearchNode, Nothing}) = new(search, parent, l_child, r_child)
    end

    function initializeNodeList(index_search::Bool)
        return index_search == true ? IndexRange[] : kdNode[]
    end

    function addNode!(nodes::Vector{kdNode{T, VecDP}}, node::kdNode) where {T, VecDP<:Vector{DataPoint{T}}}
        push!(nodes, node)
    end

    function addNode!(nodes::Vector{IndexRange}, node::kdNode)
        push!(nodes, node.index_range)
    end

    # Top level search-dispatch function
    # TODO: add functionality for search-and-mark to this function? Need: 
    function search(tree::kdTree, query<:GeometricPrimitive; 
        include_boundary = true::Bool, 
              watertight = false::Bool,
         fully_contained = false::Bool, 
            index_search = false::Bool,
             lazy_search = false::Bool,
             search_tree = false::Bool )

        nodes = initializeNodeList(index_search)
        if search_tree
            search_tree_root = SearchNode()
        else
            search_tree_root = nothing
        end
        search!(nodes, tree.root, query, include_boundary=include_boundary, watertight=watertight, fully_contained=fully_contained, lazy_search=lazy_search, search_node=search_tree_root)
        
        return nodes, search_tree_root
    end

    # Bounding volume search
    function search!(nodes, node::Union{kdNode, Nothing}, query_bv::BoundingVolume; 
        include_boundary = true::Bool,
              watertight = false::Bool,
         fully_contained = false::Bool,
            index_search = false::Bool,
             lazy_search = false::Bool,
             search_node = nothing::Union{SearchNode, Nothing} )

        if isnothing(node)
            return
        else 
            node_bv = watertight ? node.bv_watertight : node.bv

            # TODO: watertight searches: For efficiency, only check the BV against the necessary dimensions to save computation
            #       non-watertight: do the reduced check on the watertight BV first, then the real BV if necessary?
            cropped_query_bv = getIntersection(node_bv, query_bv)

            if !cropped_query_bv.is_empty
                if node.is_leaf
                    if ( fully_contained && isContained(query_bv, node_bv, include_boundary=include_boundary) ) || !fully_contained
                        addNode(nodes, node)
                    end
                else
                    if lazy_search && ( ( fully_contained && isContained(query_bv, node_bv, include_boundary=include_boundary) ) || !fully_contained )
                        addNode(nodes, node)
                    else
                        if (include_boundary && cropped_query_bv.lb[node.split_dim] <= node.split_val ) || cropped_query_bv.lb[node.split_dim] < node.split_val
                            search!(nodes, node.l_child, cropped_query_bv, include_boundary=include_boundary, watertight=watertight, fully_contained=fully_contained, index_search=index_search)
                        end

                        if (include_boundary && cropped_query_bv.ub[node.split_dim] >= node.split_val ) || cropped_query_bv.ub[node.split_dim] > node.split_val
                            search!(nodes, node.r_child, cropped_query_bv, include_boundary=include_boundary, watertight=watertight, fully_contained=fully_contained, index_search=index_search)
                        end
                    end
                end
            end
        end
    end

    # Ball search
    function search!(nodes, node::Union{kdNode, Nothing}, query_ball::Ball;
        include_boundary = true::Bool,
              watertight = false::Bool,
         fully_contained = false::Bool,
            index_search = false::Bool,
             lazy_search = false::Bool,
             search_node = nothing::Union{SearchNode, Nothing} )

        if isnothing(node)
            return
        else
            node_bv = watertight ? node.bv_watertight : node.bv

            # TODO: watertight searches: For efficiency, only check the BV against the necessary dimensions to save computation
            #       non-watertight: do the reduced check on the watertight BV first, then the real BV if necessary?
            cropped_query_bv = getIntersection(node_bv, query_ball)

            if !cropped_query_bv.is_empty
                if node.is_leaf
                    if ( fully_contained && isContained(query_ball, node_bv, include_boundary=include_boundary) ) || !fully_contained
                        addNode(nodes, node)
                    end
                else
                    if lazy_search && ( ( fully_contained && isContained(query_ball, node_bv, include_boundary=include_boundary) ) || !fully_contained )
                        addNode(nodes, node)
                    else
                        if (include_boundary && cropped_query_bv.lb[node.split_dim] <= node.split_val ) || cropped_query_bv.lb[node.split_dim] < node.split_val
                            search!(nodes, node.l_child, ball, include_boundary=include_boundary, watertight=watertight, fully_contained=fully_contained, index_search=index_search)
                        end

                        if (include_boundary && cropped_query_bv.ub[node.split_dim] >= node.split_val ) || cropped_query_bv.ub[node.split_dim] > node.split_val
                            search!(nodes, node.r_child, ball, include_boundary=include_boundary, watertight=watertight, fully_contained=fully_contained, index_search=index_search)
                        end
                    end
                end
            end
        end
    end

end # kdTree module

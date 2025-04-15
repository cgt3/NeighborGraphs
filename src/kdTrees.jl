module kdTrees

    using LinearAlgebra
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


    const DEFAULT_PARTITION_PALETTE = palette([:darkred,       :orangered3,  :darkorange2, :goldenrod1,
                                               :chartreuse3,   :forestgreen, :darkgreen,      #:olivedrab2
                                               :navy,          :blue,        :deepskyblue,    #:blue
                                               :mediumpurple1, :purple3,     :purple4  ])     #:purple1

                            
    # Data types
    export Ball, Cone, BoundingVolume, kdNode,  kdTree, DataPoint, IndexRange

    # Functions
    export spatial2Data, getSplit, contains, getContainingNode, getLeaves, partitionPlot, treePlot


    struct DataPoint
        x::AbstractVector
        y
    end

    function spatial2Data(data_num::Vector{Vector{T}}) where T <: Real
        return DataPoint.(data_num, nothing)
    end


    function Base.getindex(p::DataPoint, i::Integer)
        return p.x[i]
    end


    # TODO: there is most definitely something in Julia that replicates the below, but if
    #     there is true need to also have n and not just i0,in, than it may be ok to keep.
    struct IndexRange
        first::Integer
        last::Integer
        n::Integer

        IndexRange(first; last) = new(first, last, last - first + 1)
        IndexRange(first; n) = new(first, first + n - 1, n)
    end

    # Geometric primitives
    struct Ball
        center
        radius::Real
        p::Real
    end

    struct BoundingVolume
        min::AbstractArray
        max::AbstractArray
        is_empty::Bool

        function BoundingVolume()
            return new(Inf, -Inf, true)
        end

        function BoundingVolume(min, max)
            if any(min .> max)
                @error "kdTree::BoundingVolume: Cannot construct bounding volume with min (=$min) > max (=$max)"
            elseif length(min) != length(max)
                @error "kdTree::BoundingVolume: length(min) = $(length(min)) != length(max) = $(length(max))"
            end

            return new(min, max, false)
        end
    end

    function BoundingVolume(ball::Ball)
        return BoundingVolume(ball.center .- ball.radius, ball.center .+ radius)
    end

    function BoundingVolume(data::Vector{DataPoint}, index_range::IndexRange)
        max = copy(data[index_range.first].x)
        min = copy(data[index_range.first].x)
        for i = index_range.first:index_range.last
            d_max = max .< data[i].x
            max[d_max] = data[i].x[d_max]

            d_min = min .> data[i].x
            min[d_min] = data[i].x[d_min]
        end
        return BoundingVolume(min, max)
    end


    # Convention: use <= goes to the left, > goes to the right
    mutable struct kdNode
        data::Vector{DataPoint}
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
        kdNode(data, original_indices, is_leaf, index_range, split_dim, split_val, bv, bv_wt, parent, l_child, r_child, depth) = 
            new(data, original_indices, is_leaf, index_range, split_dim, split_val, bv, bv_wt, parent, l_child, r_child, depth)

        # The recursive, top-level constructor
        function kdNode(data::Vector{DataPoint}; original_indices=[1:length(data)...]::Vector{Int64}, 
                  parent = nothing::Union{kdNode, Nothing}, 
             index_range = IndexRange(1, n=length(data))::IndexRange, 
                   depth = 0::Integer, 
            num_leaf_pts = DEFAULT_NUM_LEAF_PTS::Integer)

            if index_range.n == 0
                return nothing
            elseif index_range.n <= num_leaf_pts 
                bv = BoundingVolume(data, index_range)
                bv_wt = getWatertightBoundingVolume(parent, data, index_range)
                return kdNode(data, original_indices, true, index_range, 0, NaN, bv, bv_wt, parent, nothing, nothing, depth)
            end

            node = new()
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


    function getCentroid(bv::BoundingVolume, data::Vector{DataPoint}, index_range::IndexRange)
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
                bv_max = copy(parent_bv.max)
                bv_max[parent.split_dim] = parent.split_val
                return  BoundingVolume(copy(parent_bv.min), bv_max)
            else
                bv_min = copy(parent_bv.min)
                bv_min[parent.split_dim] = parent.split_val
                return BoundingVolume(bv_min, copy(parent_bv.max))
            end
        end
    end

    # TODO: change this function to allow different schemes for choosing the split value
    function getSplit(bv::BoundingVolume, bv_wt::BoundingVolume, data, index_range::IndexRange)
        centroid = getCentroid(bv, data, index_range)
        # mid = 0.5 * (bv.max - bv.min) + bv.min
        # dist = centroid - mid

        # # Choose the dimension where the centroid is closest to the midpoint
        # split_dim = argmin(dist)
        # split_val = centroid[split_dim]

        # Choose the dimension that splits the watertight BV into the most square pieces
        bv_length = bv_wt.max .- bv_wt.min
        split_dim = argmax(bv_length)
        split_val = centroid[split_dim]

        return split_dim, split_val
    end


    function partitionData!(data::Vector{DataPoint}, original_indices, index_range::IndexRange, split_dim, split_val)
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

    struct kdTree
        data::Vector{DataPoint}
        original_indices::Vector{Int64}
        root::kdNode
        max_depth::Integer
        num_leaves::Integer
        num_pts::Integer

        function kdTree(data::Vector{DataPoint}; num_leaf_pts::Integer=DEFAULT_NUM_LEAF_PTS)
            original_indices = [1:length(data)...]
            root = kdNode(data, original_indices=original_indices, num_leaf_pts=num_leaf_pts)
            max_depth, num_leaves, num_pts = getTreeStats(root)

            return new(data, original_indices, root, max_depth, num_leaves, num_pts)
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
            if any(query_pt .- node.bv_watertight.max .> pt_tol) || any( query_pt .- node.bv_watertight.min .< pt_tol)
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
            i_color = node.depth % length(color_palette) + 1

            if watertight 
                split_bv = node.bv_watertight
            else
                split_bv = node.bv 
            end

            partitionPlot!(p, node.l_child, watertight=watertight, index=index, color_palette=color_palette, linewidth=linewidth)
            partitionPlot!(p, node.r_child, watertight=watertight, index=index, color_palette=color_palette, linewidth=linewidth)

            if !watertight || node.depth == 0
                # Plot the parent BV:
                plot!(p, node.bv.min[index[1]]*[1, 1], [node.bv.min[index[2]], node.bv.max[index[2]]], linecolor=color_palette[i_color], linestyle=:dash, linewidth=linewidth/2) # left
                plot!(p, node.bv.max[index[1]]*[1, 1], [node.bv.min[index[2]], node.bv.max[index[2]]], linecolor=color_palette[i_color], linestyle=:dash, linewidth=linewidth/2) # right
                plot!(p, [node.bv.min[index[1]], node.bv.max[index[1]]], node.bv.min[index[2]]*[1, 1], linecolor=color_palette[i_color], linestyle=:dash, linewidth=linewidth/2) # bottom
                plot!(p, [node.bv.min[index[1]], node.bv.max[index[1]]], node.bv.max[index[2]]*[1, 1], linecolor=color_palette[i_color], linestyle=:dash, linewidth=linewidth/2) # top
            end

            # Plot the splitting plane
            if node.split_dim == index[1]
                plot!(p, node.split_val*[1, 1], [split_bv.min[index[2]], split_bv.max[index[2]]], linecolor=color_palette[i_color], linewidth=max(MIN_SPLIT_LINE_WIDTH, MAX_SPLIT_LINE_WIDTH - node.depth))
            elseif node.split_dim == index[2]
                plot!(p, [split_bv.min[index[1]], split_bv.max[index[1]]], node.split_val*[1, 1], linecolor=color_palette[i_color], linewidth=max(MIN_SPLIT_LINE_WIDTH, MAX_SPLIT_LINE_WIDTH - node.depth))
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

    function intersects(bv1::BoundingVolume, bv2::BoundingVolume; include_boundary=true::Bool)
        if (  include_boundary && ( any(bv1.min .> bv2.max)  || any(bv1.max .<  bv2.min) ) ) ||
           ( !include_boundary && ( any(bv1.min .>= bv2.max) || any(bv1.max .<= bv2.min) ) )
            return false
        else
            return true
        end
    end

    function isContained(bv::BoundingVolume, query_bv::BoundingVolume; include_boundary=true::Bool)
        if ( !include_boundary && ( all(query_bv.max .<  bv.max) && all(query_bv.min .>  bv.min) ) ) || 
           (  include_boundary && ( all(query_bv.max .<= bv.max) && all(query_bv.min .>= bv.min) ) )
            return true
        else
            return false
        end
    end

    function getIntersection(bv1::BoundingVolume, bv2::BoundingVolume)
        if bv1.is_empty || bv2.is_empty
            return BoundingVolume()
        end

        new_min = max.(bv1.min, bv2.min)
        new_max = min.(bv1.max, bv2.max)
        if any(new_min .> new_max)
            return BoundingVolume()
        end
        return BoundingVolume(new_min, new_max)
    end

    function search(tree::kdTree, query_bv::BoundingVolume; 
        include_boundary = true::Bool, 
              watertight = false::Bool,
         fully_contained = false::Bool, 
            index_search = false::Bool )

        if index_search
            nodes = IndexRange[]
        else
            nodes = kdNode[]
        end
        search!(nodes, tree.root, query_bv, include_boundary=include_boundary, watertight=watertight, fully_contained=fully_contained)
        return nodes
    end

    # TODO: Add a go_to_leaf/lazy_search parameter/flag to allow returning early
    function search!(nodes, node::Union{kdNode, Nothing}, query_bv::BoundingVolume; 
        include_boundary = true::Bool,
              watertight = false::Bool,
         fully_contained = false::Bool,
            index_search = false::Bool )

        if isnothing(node)
            return
        else 
            node_bv = watertight ? node.bv_watertight : node.bv
            cropped_query_bv = getIntersection(node_bv, query_bv)

            if !cropped_query_bv.is_empty
                if node.is_leaf
                    if ( fully_contained && isContained(query_bv, node_bv, include_boundary=include_boundary) ) || !fully_contained
                        if index_search
                            push!(nodes, node.index_range)
                        else 
                            push!(nodes, node)
                        end
                    end
                else
                    search!(nodes, node.l_child, cropped_query_bv, include_boundary=include_boundary, watertight=watertight, fully_contained=fully_contained, index_search=index_search)
                    search!(nodes, node.r_child, cropped_query_bv, include_boundary=include_boundary, watertight=watertight, fully_contained=fully_contained, index_search=index_search)
                end
            end
        end
    end

end # kdTree module

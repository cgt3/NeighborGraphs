
module kdTrees

    using Plots
    using Printf

    # Data types
    export BoundingVolume, kdNode,  kdTree

    # Constucting methods
    export getSplit

    # Processing methods
    export contains, getContainingNode, getLeaves, partitionPlot, treePlot

    const DEFAULT_NUM_LEAF_PTS = 40
    const DEFAULT_PT_TOL = 1e-12

    const DEFAULT_PLOT_WIDTH_PER_NODE = 100
    const DEFAULT_PLOT_HEIGHT_PER_NODE = 200
    const DEFAULT_LEAF_FONTSIZE = 8
    const DEFAULT_MAX_TREE_PLOT_FONTSIZE = 20


    const DEFAULT_PARTITION_PALETTE = palette(:RdYlBu_10) # palette([])

    struct BoundingVolume{T}
        min::T
        max::T
    end

    function BoundingVolume(data::Vector{T}, i0, n) where T <:Number
        max = copy(data[i0])
        min = copy(data[i0])
        for i = i0:i0+n-1
            if max .< data[i]
                max = data[i]
            end

            if min .> data[i]
                min = data[i]
            end
        end
        return BoundingVolume(min, max)
    end

    function BoundingVolume(data, i0, n)
        max = copy(data[i0])
        min = copy(data[i0])
        for i = i0:i0+n-1
            d_max = max .< data[i]
            max[d_max] = data[i][d_max]

            d_min = min .> data[i]
            min[d_min] = data[i][d_min]
        end
        return BoundingVolume(min, max)
    end

    function getCentroid(bv::BoundingVolume, data, i0, n)
        # TODO: will not work for 1D data
        centroid = zeros(eltype(data[1]), size(data[1]))
        for i in i0:i0+n-1
            centroid += data[i]
        end
        centroid *= 1/n

        return centroid
    end

    # Convention: use <= goes to the left, > goes to the right
    mutable struct kdNode
        data::AbstractArray
        is_leaf::Bool
        i0::Integer
        n::Integer
        split_dim::Integer
        split_val::Real
        bv::BoundingVolume
        bv_watertight::BoundingVolume
        parent::Union{kdNode, Nothing}
        l_child::Union{kdNode, Nothing}
        r_child::Union{kdNode, Nothing}
        depth::Integer

        # Default constructor
        kdNode(data, is_leaf, i0, n, split_dim, split_val, bv, bv_wt, parent, l_child, r_child, depth) = new(data, is_leaf, i0, n, split_dim, split_val, bv, bv_wt, parent, l_child, r_child, depth)

        # The recursive, top-level constructor
        function kdNode(data::AbstractArray; parent=nothing::Union{kdNode, Nothing}, i0=1::Integer, n=length(data)::Integer, depth=0::Integer, num_leaf_pts=DEFAULT_NUM_LEAF_PTS::Integer)
            if n == 0
                return nothing
            elseif n <= num_leaf_pts 
                bv = BoundingVolume(data, i0, n)
                bv_wt = getWatertightBoundingVolume(parent, data, i0, n)
                return kdNode(data, true, i0, n, 0, NaN, bv, bv_wt, parent, nothing, nothing, depth)
            end

            node = new()
            node.data = data
            node.is_leaf = false
            node.i0 = i0
            node.n = n
            node.depth=depth

            node.bv = BoundingVolume(data, i0, n)
            node.bv_watertight = getWatertightBoundingVolume(parent, data, i0, n)

            # Find the splitting dimension and value
            node.split_dim, node.split_val = getSplit(node.bv, node.bv_watertight, data, i0, n)

            # Partition the data
            i_split = partitionData!(data, i0, n, node.split_dim, node.split_val)

            # For the children nodes
            i0_left  = node.i0
            i0_right = i_split + 1
            n_left = i_split - node.i0 + 1
            n_right = node.n - n_left 

            # Attach family nodes
            node.parent = parent
            node.l_child = kdNode(data, parent=node, i0=i0_left,  n=n_left,  depth=depth+1, num_leaf_pts=num_leaf_pts)
            node.r_child = kdNode(data, parent=node, i0=i0_right, n=n_right, depth=depth+1, num_leaf_pts=num_leaf_pts)
            return node
        end
    end # kdNode struct declaration


    function getWatertightBoundingVolume(parent::Union{kdNode, Nothing}, data, i0::Integer, n::Integer)
        if isnothing(parent)
            return BoundingVolume(data, i0, n)
        else
            parent_bv = parent.bv_watertight
            if parent.i0 == i0 # we are constructing the left child
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
    function getSplit(bv::BoundingVolume, bv_wt::BoundingVolume, data, i0, n)
        centroid = getCentroid(bv, data, i0, n)
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


    function partitionData!(data, i0, n, split_dim, split_val)
        i_start = i0
        i_end = i0 + n - 1

        while i_start < i_end 
            is_valid_L = data[i_start][split_dim] <= split_val
            is_valid_R = data[i_end][split_dim] > split_val

            if !is_valid_L && !is_valid_R
                swap          = data[i_end]
                data[i_end]   = data[i_start]
                data[i_start] = swap
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
        root::kdNode
        max_depth::Integer
        num_leaves::Integer
        num_pts::Integer

        function kdTree(data; num_leaf_pts::Integer=DEFAULT_NUM_LEAF_PTS)
            root = kdNode(data, num_leaf_pts=num_leaf_pts)
            max_depth, num_leaves, num_pts = getTreeStats(root)

            return new(root, max_depth, num_leaves, num_pts)
        end
    end

    function getTreeStats(node::Union{kdNode, Nothing})
        if isnothing(node)
            return 0, 0, 0
        elseif node.is_leaf
            return node.depth, 1, node.n
        else
            depth_L, num_leaves_L, num_pts_L = getTreeStats(node.l_child)
            depth_R, num_leaves_R, num_pts_R = getTreeStats(node.r_child)

            return max(depth_L, depth_R), num_leaves_L + num_leaves_R, num_pts_L + num_pts_R
        end
    end

    # User-exposed functions -------------------------------------------------------------
    function getContainingNode(query_pt, tree::kdTree; pt_tol=DEFAULT_NUM_LEAF_PTS)
        return getContainingNode(query_pt, tree.root, pt_tol=pt_tol)
    end

    function getContainingNode(query_pt, node::Union{kdNode, Nothing}; pt_tol=DEFAULT_PT_TOL)
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
            return getContainingNode(query_pt, node.l_child, pt_tol=pt_tol)
        else
            return getContainingNode(query_pt, node.r_child, pt_tol=pt_tol)
        end
    end

    function contains(query_pt, tree::kdTree; pt_tol=DEFAULT_PT_TOL)
        node = getContainingNode(query_pt, tree.root, pt_tol=pt_tol)
        return contains(query_pt, node, pt_tol=pt_tol)
    end

    function contains(query_pt, node::Union{kdNode, Nothing}; pt_tol=DEFAULT_PT_TOL)
        if isnothing(node)
            return false
        end

        i_end = node.i0 + node.n - 1
        for i = i0:i_end 
            if all(abs.(query_pt .- node.data[i])) .< pt_tol
                return true
            end
        end

        return false
    end

    function treePlot(tree::kdTree;
        x_spacing=DEFAULT_PLOT_WIDTH_PER_NODE,
        y_spacing=DEFAULT_PLOT_HEIGHT_PER_NODE,
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
        treePlot!(p, 1, tree.max_depth, tree.root, leaf_fontsize=leaf_fontsize, max_fontsize=max_fontsize)
        return p
    end

    function treePlot!(p, i_level::Integer, max_depth::Integer, node::Union{kdNode, Nothing}; 
        leaf_fontsize=DEFAULT_TREE_PLOT_FONTSIZE,
        max_fontsize=DEFAULT_MAX_TREE_PLOT_FONTSIZE )

        x_pos = 0.5 + (i_level - 0.5)*( 2^(max_depth - node.depth) )
        y_pos = node.depth
        
        fontsize = min(leaf_fontsize + 2*(max_depth - node.depth), max_fontsize)

        if isnothing(node)
            annotate!(p, x_pos, y_pos, text("(empty leaf)", :center, fontsize))
        elseif !node.is_leaf
            # Plot lines to child nodes
            x_spacing = 2^(max_depth - node.depth)
            x_pos_L, x_pos_R = x_pos - x_spacing/4 , x_pos + x_spacing / 4
            y_pos_child = y_pos + 1
            plot!(p, [x_pos, x_pos_L], [y_pos, y_pos_child], linecolor=:black, linewidth=3)
            plot!(p, [x_pos, x_pos_R], [y_pos, y_pos_child], linecolor=:black, linewidth=3)

            # Add text for this node
            annotate!(p, x_pos, y_pos, text(@sprintf("Split dim: %d\nVal: %.4e", node.split_dim, node.split_val), :center, fontsize))


            treePlot!(p, 2*i_level - 1, max_depth, node.l_child, leaf_fontsize=leaf_fontsize, max_fontsize=max_fontsize)
            treePlot!(p, 2*i_level,     max_depth, node.r_child, leaf_fontsize=leaf_fontsize, max_fontsize=max_fontsize)
        else
            annotate!(p, x_pos, y_pos, text("Leaf: i0=$(node.i0),\n n=$(node.n)", :center, fontsize))
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
            i_end = node.i0 + node.n - 1
            scatter!(p, getindex.(node.data[node.i0:i_end], index[1]), getindex.(node.data[node.i0:i_end], index[2]), markercolor=color_palette[1], ms=5)
        else
            i_color = node.depth % length(color_palette) + 1
            if !watertight || node.depth == 0
                # Plot the parent BV:
                plot!(p, node.bv.min[index[1]]*[1, 1], [node.bv.min[index[2]], node.bv.max[index[2]]], linecolor=color_palette[i_color], linestyle=:dash, linewidth=linewidth/2) # left
                plot!(p, node.bv.max[index[1]]*[1, 1], [node.bv.min[index[2]], node.bv.max[index[2]]], linecolor=color_palette[i_color], linestyle=:dash, linewidth=linewidth/2) # right
                plot!(p, [node.bv.min[index[1]], node.bv.max[index[1]]], node.bv.min[index[2]]*[1, 1], linecolor=color_palette[i_color], linestyle=:dash, linewidth=linewidth/2) # bottom
                plot!(p, [node.bv.min[index[1]], node.bv.max[index[1]]], node.bv.max[index[2]]*[1, 1], linecolor=color_palette[i_color], linestyle=:dash, linewidth=linewidth/2) # top
            end

            if watertight 
                split_bv = node.bv_watertight
            else
                split_bv = node.bv 
            end

            # Plot the splitting plane
            if node.split_dim == index[1]
                plot!(p, node.split_val*[1, 1], [split_bv.min[index[2]], split_bv.max[index[2]]], linecolor=color_palette[i_color], linewidth=linewidth)
            elseif node.split_dim == index[2]
                plot!(p, [split_bv.min[index[1]], split_bv.max[index[1]]], node.split_val*[1, 1], linecolor=color_palette[i_color], linewidth=linewidth)
            end

            partitionPlot!(p, node.l_child, watertight=watertight, index=index, color_palette=color_palette, linewidth=linewidth)
            partitionPlot!(p, node.r_child, watertight=watertight, index=index, color_palette=color_palette, linewidth=linewidth)
        end
    end


    function getLeaves(tree::kdTree)
        leaves = kdNode[]
        getLeaves!(leaves, tree.root)

        return leaves
    end

    function getLeaves!(leaves, node::Union{kdNode, Nothing})
        if isnothing(node)
            return
        elseif node.is_leaf
            push!(leaves, node)
        else
            getLeaves!(leaves, node.l_child)
            getLeaves!(leaves, node.r_child)
        end
    end
end # kdTree module
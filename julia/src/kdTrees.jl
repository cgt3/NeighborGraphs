
module kdTrees

    using Plots

    export BoundingVolume
    export kdTree, kdNode, search, plotPartitions

    DEFAULT_NUM_LEAF_PTS=40

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
        function kdNode(data::AbstractArray; parent=nothing, i0=1, n=length(data), depth=0, num_leaf_pts=DEFAULT_NUM_LEAF_PTS)
            if n <= num_leaf_pts 
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


    # User-exposed functions -------------------------------------------------------------
    function kdTree(data; num_leaf_pts::Integer=DEFAULT_NUM_LEAF_PTS)
        return kdNode(data, num_leaf_pts=num_leaf_pts)
    end


    function search(query_pt, kd_node::kdNode)
        if kd_node.is_leaf == true
            return kd_node
        elseif query_pt[kd_node.split_dim] <= kd_node.split_val
            return search(query_pt, kd_node.l_child)
        else
            return search(query_pt, kd_node.r_child)
        end
    end

    function plotTree(tree::kdNode)
    end

    function plotTree!(node::kdNode)
    end

    function plotPartitions(tree::kdNode; 
        watertight=true,
        index=(1,2), 
        size=(1500,1000), 
        fontsize=20,
        color_palette=palette(:RdYlBu_10), # TODO: create a custom palette for plotting partitions
        linewidth=4)

        p = plot(size=size,
            xtickfont=font(fontsize), 
            ytickfont=font(fontsize),
            leg=:false
        )
        plotPartitions!(p, tree, watertight=watertight, index=index, color_palette=color_palette, linewidth=linewidth)
        return p
    end

    # TODO: adapt to plot 1D points as well
    function plotPartitions!(p, node::kdNode; index=(1,2),
            watertight,
            color_palette,
            linewidth)

        if node.is_leaf
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

            plotPartitions!(p, node.l_child, watertight=watertight, index=index, color_palette=color_palette, linewidth=linewidth)
            plotPartitions!(p, node.r_child, watertight=watertight, index=index, color_palette=color_palette, linewidth=linewidth)
        end
    end


    function getLeaves(tree::kdNode)

    end

    function getLeaves!(node::kdNode)

    end
end # kdTree module
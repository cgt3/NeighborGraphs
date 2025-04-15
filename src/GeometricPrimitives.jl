module GeometricPrimitives

    # Data types
    export Ball, Cone, BoundingVolume#, Cone, Hyperplane, Simplex

    # Functions
    export intersects, isContained, getIntersection

    struct Ball
        center::Vector{Real}
        radius::Real
        p::Real
    end

    # TODO: Finish these classes
    struct Cone end 
    struct Hyperplane{Integer} end
    struct Simplex{Integer} end

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

    # Other primitives to BVs:
    function BoundingVolume(ball::Ball)
        return BoundingVolume(ball.center .- ball.radius, ball.center .+ radius)
    end
end
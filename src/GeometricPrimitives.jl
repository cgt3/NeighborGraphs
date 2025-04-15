module GeometricPrimitives

    export Ball, Cone, BoundingVolume#, Cone, Hyperplane, Simplex

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

    struct Ball
        center::Vector{Real}
        radius::Real
        p::Real
    end

    function BoundingVolume(ball::Ball)
        return BoundingVolume(ball.center .- ball.radius, ball.center .+ radius)
    end

    # TODO: Finish these classes
    struct Cone end 
    struct Hyperplane{Integer} end
    struct Simplex{Integer} end
end
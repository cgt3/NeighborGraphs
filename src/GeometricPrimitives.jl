module GeometricPrimitives

    using LinearAlgebra

    # Data types
    export Ball, Cone, BoundingVolume#, Cone, Hyperplane, Simplex

    # Functions
    export intersects, isContained, getIntersection

    const DEFAULT_BV_POINT_TOL = 1e-15

    struct Ball
        center::Vector
        radius::Real
        p::Real
        dim::Integer
        active_dim::Vector
        inactive_dim::Vector
        is_active::Vector{Bool}
        embedding_dim::Integer

        function Ball(center::VecReal, radius::Real; p=2::Real, active_dim=true::Bool, indices=(active_dim ? [eachindex(center)...] : VecInt[])::VecInt) where {T_real<:Real, T_int<:Integer, VecReal<:Vector{T_real}, VecInt<:Vector{T_int}}
            unique_indices=unique(indices)
            if active_dim
                is_active = zeros(Bool, length(center))
                is_active[unique_indices] = true
    
                all_dim = [eachindex(center)...]
                inactive_dim = all_dim[!is_active]
                dim = length(unique_indices)
                return new(center, radius, p, dim, unique_indices, inactive_dim, is_active, length(center))
            else
                is_active = ones(Bool, length(center))
                is_active[unique_indices] = false

                all_dim = [eachindex(center)...]
                active_dim = all_dim[is_active]
                dim = length(unique_indices)
                return new(center, radius, p, dim, active_dim, unique_indices, is_active, length(center))
            end
        end
    end # struct

    # TODO: Finish these classes
    struct Cone end 
    struct Hyperplane end
    struct Simplex end

    struct BoundingVolume{}
        lb::Vector
        ub::Vector
        is_empty::Bool
        dim::Integer
        active_dim::Vector
        inactive_dim::Vector
        is_active::Vector{Bool}

        # For invalid/emtpy BV's (note an empty BV differs from a non-empty BV with dimension 0 (a point))
        function BoundingVolume()
            return new(Inf, -Inf, true, 0, [], [], [])
        end

        function BoundingVolume(lb::Vector{T1}, ub::Vector{T2}; tol=DEFAULT_BV_POINT_TOL::Real) where {T1<:Real, T2<:Real}
            if length(lb) != length(ub)
                @error "kdTree::BoundingVolume: lb (length=$(length(lb))) and ub (length=$(length(lb))) points have different dimensions"
            elseif any(lb .> ub)
                @error "kdTree::BoundingVolume: Cannot construct bounding volume with lb (=$lb) > ub (=$ub)"
            elseif length(lb) != length(ub)
                @error "kdTree::BoundingVolume: length(lb) = $(length(lb)) != length(ub) = $(length(ub))"
            end

            dim = length(lb)
            is_active = ones(Bool, length(lb))
            for d in eachindex(lb)
                if abs(lb[d] - ub[d]) < tol
                    dim -= 1
                    is_active[d] = false
                end
            end
            all_dim = [eachindex(lb)...]
            active_dim   = all_dim[is_active]
            inactive_dim = all_dim[ is_active .!= true ] # Cannot do !x it seems

            return new(lb, ub, false, dim, active_dim, inactive_dim, is_active)
        end
    end

    function getClosestPoint(bv::BoundingVolume, query_pt::Array)
        closest_pt = copy(query_pt)

        I_lb = query_pt .< bv.lb
        I_ub = query_pt .> bv.ub

        closest_pt[I_lb] = bv.lb[I_lb]
        closest_pt[I_ub] = bv.ub[I_ub]

        return closest_pt
    end

    function getFurthestPoint(bv::BoundingVolume, query_pt::Array)
        furthest_pt = similar(query_pt)        

        ub_is_closer = bv.ub - query_pt .<= query_pt - bv.lb
        I_lb = query_pt .> bv.ub .|| (bv.lb .<= query_pt &&  ub_is_closer)
        I_ub = query_pt .< bv.lb .|| (query_pt .<= bv.ub && !ub_is_closer)

        furthest_pt[I_lb] = bv.lb[I_lb]
        furthest_pt[I_ub] = bv.ub[I_ub]

        return furthest_pt
    end

    function isContained(bv::BoundingVolume, query_pt::Array; include_boundary=true::Bool)
        if ( include_boundary  && all(bv.lb .<= query_pt .<= bv.ub) ) ||
           ( !include_boundary && all(bv.lb .<  query_pt .<  bv.ub) )
           return true
        else
            return false
        end
    end

    function isContained(bv::BoundingVolume, query_bv::BoundingVolume; include_boundary=true::Bool)
        if ( !include_boundary && ( all(query_bv.ub .<  bv.ub) && all(query_bv.lb .>  bv.lb) ) ) || 
           (  include_boundary && ( all(query_bv.ub .<= bv.ub) && all(query_bv.lb .>= bv.lb) ) )
            return true
        else
            return false
        end
    end

    function intersects(bv1::BoundingVolume, bv2::BoundingVolume; include_boundary=true::Bool)
        if (  include_boundary && ( any(bv1.lb .> bv2.ub)  || any(bv1.ub .<  bv2.lb) ) ) ||
           ( !include_boundary && ( any(bv1.lb .>= bv2.ub) || any(bv1.ub .<= bv2.lb) ) )
            return false
        else
            return true
        end
    end

    function getIntersection(bv1::BoundingVolume, bv2::BoundingVolume; tol=DEFAULT_BV_POINT_TOL::Real)
        if bv1.is_empty || bv2.is_empty
            return BoundingVolume()
        end

        new_lb = max.(bv1.lb, bv2.lb)
        new_ub = min.(bv1.ub, bv2.ub)
        if any(new_lb .> new_ub)
            return BoundingVolume()
        end
        return BoundingVolume(new_lb, new_ub, tol=tol)
    end

    function faceIndex2SpatialIndex(face_index::Integer, num_dim::Integer)
        return face_index <= num_dim ? face_index : face_index - num_dim
    end

    function getFaceBoundingVolume(face_index::Integer, bv::BoundingVolume; tol=DEFAULT_BV_POINT_TOL)
        face_lb, face_ub = copy(bv.lb), copy(bv.ub)

        if face_index <= length(bv.lb) # Lower bound face
            face_ub[face_index] = bv.lb[face_index]
        else # Upper bound face
            d = face_index - length(bv.lb)
            face_lb[d] = bv.ub[d]
        end

        return BoundingVolume(face_lb, face_ub, tol=tol)
    end


    # Balls ----------------------------------------------------------------------
    function BoundingVolume(ball::Ball; tol=DEFAULT_BV_POINT_TOL::Real)
        return BoundingVolume(ball.center .- ball.radius, ball.center .+ radius, tol=tol)
    end

    function isContained(ball::Ball, query_pt::Pt; include_boundary=true::Bool, tol=DEFAULT_BV_POINT_TOL) where {T<:Real, Pt<:Vector{T}}
        if ball.dim < length(ball.center) # The ball does not have full dimension
            for d_fixed in ball.inactive_dim
                if abs(query_pt[d_fixed] - ball.center[d_fixed]) > tol
                    return false
                end
            end
            R_query = norm(query_pt[ball.active_dim] - ball.center[ball.active_dim], ball.p)
        else # The ball has full-dimension
            R_query = norm(query_pt - ball.center, ball.p)
        end

        return include_boundary ? R_query <= ball.radius : R_query < ball.radius
    end

    function isContained(bv::BoundingVolume, query_ball::Ball; include_boundary=true::Bool, tol=DEFAULT_BV_POINT_TOL::Real)
        return isContained(bv, BoundingVolume(query_ball, tol=tol), include_boundary=include_boundary)
    end

    function isContained(ball::Ball, query_bv::BoundingVolume; include_boundary=true::Bool)
        furthest_pt = getFurthestPoint(query_bv, ball.center)
        R_max = norm(furthest_pt, ball.p)
        if (  include_boundary && R_max <= ball.radius ) ||
           ( !include_boundary && R_max <  ball.radius )
            return true
        end

        return false
    end

    function intersects(bv::BoundingVolume, ball::Ball; include_boundary=true::Bool, tol=DEFAULT_BV_POINT_TOL::Real)
        bv_ball = BoundingVolume(ball, tol=tol)
        cropped_bv = getIntersection(bv, bv_ball, tol=tol)

        # First, do the easy checks against the ball's BV:
        if cropped_bv.is_empty 
            # The two are completely disjoint
            return false
        elseif isContained(bv, ball, include_boundary=include_boundary) 
            # The ball is completely contained by the BV
            # Note this case also covers the special case where the ball's radius is zero
            return true
        end

        # Second, check if the center is contained in the bv
        if isContained(bv, ball.center, include_boundary=include_boundary)
            return true
        end
        
        # If the easy checks fail, check against the ball itself:
        #  Find the point in the  BV closest to the ball's center
        #  and check if it is within the ball.
        closest_pt = getClosestPoint(bv, ball.center)
        R_closest_pt = norm(closest_pt .- ball.center, ball.p)
        if (  include_boundary && R_closest_pt <= ball.radius ) ||
           ( !include_boundary && R_closest_pt <  ball.radius )
            return true
        end

        return false
    end

    function getReducedDimBall(removal_dim::Integer, x_d::Real, ball::Ball)
        if x_d < ball.center[removal_dim] - ball.radius || ball.center[removal_dim] + ball.radius < x_d
            @error "GeometricPrimitives:getReducedDimBall: coordinate plane defined by x_$removal_dim = $x_d does not intersect the ball (center=$(ball.center), radius=$(ball.radius)"
        end
        new_center = copy(ball.center)
        new_center[removal_dim] = x_d

        new_radius = (ball.radius^ball.p - abs(x_d - ball.center)^ball.p )^(1/p)
        inactive_dim = [ball.inactive_dim..., removal_dim]
        return Ball(new_center, new_radius, p=ball.p, active_dim=false, indices=inactive_dim)
    end


    function tightenBVBounds!(bv::BoundingVolume, ball::Ball; tol=DEFAULT_BV_POINT_TOL)
        if ball.dim == 1
            d = ball.active_dim[1]
            lb_ball = ball.center[d] - ball.radius
            ub_ball = ball.center[d] + ball.radius

            if bv.lb[d] < lb_ball < bv.ub[d]
                bv.lb[d] = lb_ball
                altered_lb_indices = [d]
            else
                altered_lb_indices = []
            end

            if bv.lb[d] < ub_ball < bv.ub[d]
                bv.ub[d] = ub_ball
                altered_ub_indices = [d]
            else
                altered_ub_indices = []
            end
            return altered_lb_indices, altered_ub_indices
        end

        # For every face with no intersection with the ball, recurse
        altered_lb_indices = []
        altered_ub_indices = []
        num_dim = length(ball.center)
        for f_target in 1:2*num_dim # for each face
            face_bv = getFaceBoundingVolume(f_target, bv, tol=tol)

            if !intersects(face_bv, ball, include_boundary=true, tol=tol)
                adjacent_faces = [1:f_target-1..., f_target+1:2*num_dim...]

                for f_adjacent in adjacent_faces
                    adjacent_face_bv = getFaceBoundingVolume(f_adjacent, bv, tol=tol)

                    if intersects(adjacent_face_bv, ball, include_boundary=true)
                        d_fixed = faceIndex2SpatialIndex(f_adjacent, length(bv.lb))
                        reduced_ball = getReducedDimBall(d_fixed, adjacent_face_bv.lb[d_fixed], ball)

                         # This will modify face_adjacent's bounds 
                        altered_lb_indices_new, altered_ub_indices_new = tightenBounds!(adjacent_face_bv, reduced_ball, tol=tol)

                        # Update the higher-dim BV with the new bounds on face_adjacent
                        bv.lb[altered_lb_indices_new] .= adjacent_face_bv.lb[altered_lb_indices_new]
                        bv.ub[altered_ub_indices_new] .= adjacent_face_bv.ub[altered_ub_indices_new]

                        altered_lb_indices = vcat(altered_lb_indices, altered_lb_indices_new)
                        altered_ub_indices = vcat(altered_ub_indices, altered_ub_indices_new)
                    end
                end
            end
        end

        return altered_lb_indices, altered_ub_indices
    end

    function getIntersection(bv::BoundingVolume, ball::Ball; tol=DEFAULT_BV_POINT_TOL::Real)
        if !intersects(bv, ball, include_boundary=true)
            return BoundingVolume()
        end

        bv_ball = BoundingVolume(ball, tol=tol)
        cropped_bv = getIntersection(bv, bv_ball, tol=tol)

        # Check if the ball's center is in the BV or the BV is completely contained in the ball
        if isContained(bv, ball.center, include_boundary=true) || isContained(ball, bv, include_boundary=true)
            return cropped_bv
        end

        # The ball's center is not contained in the BV, so it
        # may be possible to crop the BV further
        tightenBVBounds!(cropped_bv, ball, tol=tol)
        return cropped_bv
    end
end
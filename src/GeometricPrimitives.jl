module GeometricPrimitives

    using LinearAlgebra

    # Data types
    export SearchableGeometry, BoundingVolume, Ball, Cone, Hyperplane, Simplex

    # General Functions:
    export intersects, isContained, getIntersection

    # BV only functions:
    export getClosestPoint, getFurthestPoint, faceIndex2SpatialIndex, getFaceBoundingVolume

    # Ball functions:
    export tightenBVBounds!, getReducedDimBall


    const DEFAULT_BV_POINT_TOL = 1e-15

    # TODO: allow sparse vectors where vectors are allowed

    struct Line 
        dir::Vector
        source::Vector

        function Line(dir, source) 
            if length(dir) != length(source)
                throw("GeometricPrimitives.Line: Direction and point vector have different dimensions.")
            end

            return new(dir / norm(dir), source)
        end
    end

    function (::Line)(s::Real)
        return source + dir*s
    end

    abstract type SearchableGeometry end

    struct BoundingVolume <: SearchableGeometry
        lb::Vector
        ub::Vector
        is_empty::Bool
        dim::Integer
        active_dim::Vector
        inactive_dim::Vector
        is_active::Vector{Bool}

        # For invalid/emtpy BV's (note an empty BV differs from a non-empty BV with dimension 0 (a point))
        function BoundingVolume()
            return new([Inf], [-Inf], true, 0, [], [], Bool[])
        end

        function BoundingVolume(lb::Vector{T1}, ub::Vector{T2}; tol=DEFAULT_BV_POINT_TOL::Real) where {T1<:Real, T2<:Real}
            if length(lb) != length(ub)
                throw("GeometricPrimitives.BoundingVolume: lb (length=$(length(lb))) and ub (length=$(length(lb))) points have different dimensions")
            elseif any(lb .> ub)
                throw("GeometricPrimitives.BoundingVolume: Cannot construct bounding volume with lb (=$lb) > ub (=$ub)")
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


    struct Ball <: SearchableGeometry
        center::Vector
        radius::Real
        p::Real
        dim::Integer
        active_dim::Vector
        inactive_dim::Vector
        is_active::Vector{Bool}
        embedding_dim::Integer

        function Ball(center::VecReal, radius::Real; p=2::Real, active_indices=true::Bool, indices=[eachindex(center)...]::Vector{Int64}) where {T_real<:Real, VecReal<:Vector{T_real}}
            if radius < 0
                throw("GeometricPrimitives.Ball: Cannot construct ball with negative radius.")
            elseif radius == 0
                return new(center, radius, p, 0, [], [eachindex(center)...], zeros(Bool, length(center)), length(center))
            end

            unique_indices=unique(indices)
            if active_indices
                is_active = zeros(Bool, length(center))
                is_active[unique_indices] .= true
    
                all_dim = [eachindex(center)...]
                inactive_dim = all_dim[is_active .== false]
                dim = sum(is_active)
                return new(center, radius, p, dim, unique_indices, inactive_dim, is_active, length(center))
            else
                is_active = ones(Bool, length(center))
                is_active[unique_indices] .= false

                all_dim = [eachindex(center)...]
                active_dim = all_dim[is_active]
                dim = sum(is_active)
                return new(center, radius, p, dim, active_dim, unique_indices, is_active, length(center))
            end
        end
    end # struct

    struct Hyperplane <: SearchableGeometry 
        point::Vector
        n::Vector
        dim::Integer
        embedding_dim::Integer
        active_dim::Vector
        inactive_dim::Vector
        is_active::Vector{Bool}

        function Hyperplane(point::Vector, n::Vector)
            if length(point) != length(n)
                throw("GeometricPrimitives.Hyperplane: Point and normal vector do not have matching dimensions.")
            end

            all_indices = [eachindex(point)]
            is_active = n .!= 0
            return new(point, n ./ norm(n), sum(n .!= 0), length(point), all_indices[is_active], all_indices[n .== 0], is_active) 
        end
    end

    struct Cone <: SearchableGeometry
        vertex::Vector 
        axis::Vector
        slope::Real

        function Cone(vertex::Vector, axis::Vector, slope::Real)
            if length(vertex) != length(axis)
                throw("GeometricPrimitives.Cone: Vertex and normal vector do not have the same dimensions (dim(vertex)=$(length(vertex)), dim(axis)=$(length(n)))")
            elseif slope < 0 
                throw("GeometricPrimitives.Cone: Cannot construct cone with negative slope.")
            end
            return new(vertex, axis ./ norm(axis), slope)
        end
    end

    struct Simplex end


    # BoundingVolume functions ====================================================================
    import Base.==
    function Base.:(==)(bv1::BoundingVolume, bv2::BoundingVolume)
        return           all(bv1.lb .== bv2.lb) &&
                         all(bv1.ub .== bv2.ub) &&
                       bv1.is_empty .== bv2.is_empty &&
                            bv1.dim .== bv2.dim &&
                 all(bv1.active_dim .== bv2.active_dim) &&
               all(bv1.inactive_dim .== bv2.inactive_dim) &&
                  all(bv1.is_active .== bv2.is_active)    
    end

    # TODO: do furthest/closest generalize to all lp-norms with p > 1?
    function getClosestPoint(bv::BoundingVolume, query_pt::Vector)
        closest_pt = copy(query_pt)

        I_lb = query_pt .< bv.lb
        I_ub = query_pt .> bv.ub

        closest_pt[I_lb] = bv.lb[I_lb]
        closest_pt[I_ub] = bv.ub[I_ub]

        return closest_pt
    end

    function getFurthestPoint(bv::BoundingVolume, query_pt::Vector)
        furthest_pt = similar(query_pt)        

        ub_is_closer = 0.5*(bv.ub + bv.lb) .<= query_pt
        lb_is_closer = ub_is_closer .== false
        furthest_pt[ub_is_closer] = bv.lb[ub_is_closer]
        furthest_pt[lb_is_closer] = bv.ub[lb_is_closer]

        return furthest_pt
    end

    function isContained(bv::BoundingVolume, query_pt::Vector; include_boundary=true::Bool)
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
        if length(query_pt) != ball.embedding_dim
            throw("GeometricPrimitive.isContained(Ball,pt): Point dim(=$(length(query_pt))) does not match ball embedding dim(=$(ball.embedding_dim))")
        end
        
        if ball.dim < ball.embedding_dim # The ball does not have full dimension
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
        return isContained(ball, furthest_pt, include_boundary=include_boundary)
    end

    function intersects(bv::BoundingVolume, ball::Ball; include_boundary=true::Bool, tol=DEFAULT_BV_POINT_TOL::Real)
        # First, do the easy checks against the ball's BV:
        if !intersects(bv, BoundingVolume(ball, tol=tol), include_boundary=include_boundary)
            # The two are completely disjoint
            return false
        elseif isContained(bv, ball.center, include_boundary=include_boundary)
            return true
        end
        
        # Note: the below check could catch all of the above cases but its use of the
        #       lp-norm makes it more expensive. The above checks are to avoid having
        #       to compute norms.

        # If the easy checks fail, check against the ball itself
        closest_pt = getClosestPoint(bv, ball.center)
        return isContained(ball, closest_pt, include_boundary=include_boundary)
    end

    function getReducedDimBall(removal_dim::Integer, x_d::Real, ball::Ball)
        if x_d < ball.center[removal_dim] - ball.radius || ball.center[removal_dim] + ball.radius < x_d
            throw("GeometricPrimitives.getReducedDimBall: coordinate plane defined by x_$removal_dim = $x_d does not intersect the ball (center=$(ball.center), radius=$(ball.radius)")
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

        # For non-simple intersections
        ub_pt_projected = ball.center
        lb_pt_projected = ball.center

        # For every face with no intersection with the ball, recurse
        altered_lb_indices = []
        altered_ub_indices = []
        num_dim = length(ball.center)
        for f_target in 1:2*num_dim # for each face
            face_bv = getFaceBoundingVolume(f_target, bv, tol=tol)

            non_simple = false
            if !intersects(face_bv, ball, include_boundary=true, tol=tol)
                adjacent_faces = [1:f_target-1..., f_target+1:2*num_dim...]
                d_target = faceIndex2SpatialIndex(f_target, num_dim)

                # Check if this face needs to be updated using a non-simple intersection
                if f_target < num_dim # f_target is a lb face
                    lb_pt_projected[d_target] = face_bv.lb[d_target]
                    if isContained(fave_bv, lb_pt_projected, include_boundary=true)
                        push!(altered_lb_indices, d_target)
                        bv.lb[d_target] = ball.center[d_target] - ball.radius
                        non_simple = true
                    end
                    lb_pt_projected[d_target] = ball.center[d_target]
                else # f_target is an ub face
                    ub_pt_projected[d_target] = face_bv.ub[d_target]
                    if isContained(fave_bv, ub_pt_projected, include_boundary=true)
                        push!(altered_ub_indices, d_target)
                        bv.ub[d_target] = ball.center[d_target] + ball.radius
                        non_simple = true
                    end
                    ub_pt_projected[d_target] = ball.center[d_target]
                end

                # For simple intersections
                if !non_simple
                    for f_adjacent in adjacent_faces
                        adjacent_face_bv = getFaceBoundingVolume(f_adjacent, bv, tol=tol)

                        if intersects(adjacent_face_bv, ball, include_boundary=true)
                            d_fixed = faceIndex2SpatialIndex(f_adjacent, num_dim)
                            reduced_ball = getReducedDimBall(d_fixed, adjacent_face_bv.lb[d_fixed], ball)

                            # This will modify face_adjacent's bounds 
                            altered_lb_indices_new, altered_ub_indices_new = tightenBounds!(adjacent_face_bv, reduced_ball, tol=tol)

                            # Update the higher-dim BV with the new bounds on face_adjacent
                            bv.lb[altered_lb_indices_new] .= adjacent_face_bv.lb[altered_lb_indices_new]
                            bv.ub[altered_ub_indices_new] .= adjacent_face_bv.ub[altered_ub_indices_new]

                            altered_lb_indices = vcat(altered_lb_indices, altered_lb_indices_new)
                            altered_ub_indices = vcat(altered_ub_indices, altered_ub_indices_new)
                        end
                    end # for
                end # if simple
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

    # Hyperplane Functions: ==================================================================================
    # TODO: this functions is not guaranteed to have a unique output; provide a means to break ties
    function getClosestPoint(bv::BoundingVolume, query_plane::Hyperplane)
    end

    # TODO: this functions is not guaranteed to have a unique output; provide a means to break ties
    function getFurthestPoint(bv::BoundingVolume, query_plane::Hyperplane)
    end

    function isContained(plane::Hyperplane, query_pt::Vector)
    end

    # Note: hyperplanes are unbounded geometries, so isContained(bv, query_plane) is ill-defined
    function isContained(plane::Hyperplane, query_bv::BoundingVolume)
    end


    function intersects(bv::BoundingVolume, query_plane::Hyperplane; include_boundary=true::Bool)
    end

    function getIntersection(bv::BoundngVolume, query_plane::Hyperplane; tol=DEFAULT_BV_POINT_TOL::Real)
    end

    # Cone Functions: ========================================================================================
    # TODO: this functions is not guaranteed to have a unique output; provide a means to break ties
    function getClosestPoint(bv::BoundingVolume, query_line::Line)
        # TODO: Finish
    end

    # TODO: this functions is not guaranteed to have a unique output; provide a means to break ties
    function getFurthestPoint(bv::BoundingVolume, query_line::Line)
        # TODO: Finish
    end

    function getMajorRadius(cone::Cone, p::Vector) 
        return (cone.vertex .- p)' * cone.axis
    end

    function getRadii(cone::Cone, p::Vector) 
        dist2Vertex = cone.vertex .- p
        R = dist2Vertex' * cone.axis
        r = norm( dist2Vertex .- R*cone.axis )

        return R, r
    end

    function getConeBounds(cone::Cone)
        if cone.slope == 0
            return Line(cone.vertex, cone.axis)
        end

        ub_scalings = similar(cone.vertex)
        lb_scalings = similar(cone.vertex)
        e_i = spzeros(length(cone.vertex))
        for i in eachindex(cone.vertex)
            e_i[dim] = 1
            if abs(e_i' * cone.axis) == 1
                lb_scalings[i] = 0.0
                ub_scalings[i] = 0.0
            else
                u0 = cone.vertex + e_i
                R0, r = getRadii(cone, u0)

                R_bound = r / cone.slope
                u_ub_i = u0[i] + (R_bound - R0)*cone.axis[i]
                u_lb_i = (cone.vertex - e_i)[i] + (R_bound + R0)*cone.axis[i]

                lb_scalings[i] = ( u_lb_i - cone.vertex[i] ) / R_bound
                ub_scalings[i] = ( u_ub_i - cone.vertex[i] ) / R_bound
            end
        end
        
        return Line(cone.vertex, lb_scalings), Line(cone.vertex, ub_scalings)
    end

    function getBoundingRadii(cone::Cone, query_bv::BoundingVolume)
        closest_R_pt  = getClosestPoint(query_bv, Hyperplane(cone.vertex, cone.axis))
        furthest_R_pt = getFurthestPoint(bv, Hyperplane(cone.vertex, cone.axis))
        return getMajorRadius(cone, closest_R_pt), getMajorRadius(cone, furthest_R_pt)
    end

    function BoundingVolume(ub_func, lb_func, s1::Real, s2::Real; tol=DEFAULT_BV_POINT_TOL::Real)
        lb1 = lb_func(s1)
        ub1 = ub_func(s1)

        lb2 = lb_func(s2)
        ub2 = ub_func(s2)

        return BoundingVolume(min.(lb1, lb2), max.(ub1, ub2), tol=tol)
    end

    function BoundingVolume(cone::Cone, R_max::Real, R_min=0::Real; tol=DEFAULT_BV_POINT_TOL::Real)
        if R_max < R_min 
            throw("GeometricPrimitives.BoundingVolume(Cone): R_max cannot be < R_min")
        end

        if R_min < 0
            throw("GeometricPrimitives.BoundingVolume(Cone): R_min and R_max do no satisfy 0 <= R_min <= R_max")
        end
        return BoundingVolume(cone_ub, cone_lb, R_max, R_min, tol=tol)
    end

    # Note: Cones are unbounded geometries, so isContained(bv, cone) is ill-defined
    function isContained(cone::Cone, query_pt::Vector; include_boundary=true::Bool )
        R, r = getRadii(cone, query_pt)
        
        if include_boundary
            return R >= 0 ? r <= cone.slope*R : false
        else
            return R > 0 ? r < cone.slope*R : false
        end
    end


    # TODO: this function is incomplete
    function isContained(cone::Cone, query_bv::BoundingVolume; include_boundary=true::Bool )

        # Cyclindrical BV condition: Check the cylindrical BV defined by R_min and r_max:
        # TODO: what is there are two or more points tied to be closest point? In the cases
        # of two points we need the negative one here.
        smallest_R_pt = getClosestPoint(bv, Hyperplane(cone.vertex, cone.axis))
        R_min = getMajorRadius(cone, smallest_R_pt)
        if R_min < 0
            return false
        end

        # If R_min = 0, the only way a BV can fit in the cone is if it:
        #   1) Is a line that passes throught the cone's vertex
        #   2) Has a corner at the vertex and the cone has a slope
        #      corresponding to 45deg or more

        furthest_pt = getFurthestPoint(query_bv, Line(cone.vertex, cone.axis))
        _, r_max = getRadii(cone, furthest_pt)
        

        # Note: the cylindrical BV condition is a strong condition than induced by the
        #       cone itself: it will reject some BV's that are in fact inside the cone. 
        #       However, if it states a BV is in the cone, that BV must be in the cone.
        if (include_boundary && r_max <= cone.slope*R_min) || (r_max < cone.slope*R_min)
            return true
        end
        return isContained(cone, bounding_pt, include_boundary=include_boundary)
    end

    # TODO: this function is not fully defined since getClosestPoint can return duplicates
    function intersects(bv::BoundingVolume, cone::Cone; include_boundary=true::Bool) 
        closest_pt = getClosestPoint(query_bv, Line(cone.vertex, cone.axis))
        return isContained(cone, closest_pt, include_boundary=include_boundary)
    end

    function getIntersection(bv::BoundingVolume, cone::Cone; tol=DEFAULT_BV_POINT_TOL::Real) 
        intersection = bv
        cone_lbs, cone_ubs = getConeBounds(cone)
        while !intersection.is_empty
            R_min, R_max = getBoundingRadii(cone, intersection)
            bv_cone = BoundingVolume(cone_ub, cone_lb, R_max, R_min, tol=tol)

            old_intersection = intersection 
            intersection = getIntersection(bv_cone, bv_cone, tol=tol)

            # TODO: Should have a tolerance instead of strict '=='?
            if old_intersection == intersection 
                return intersection 
            end
        end

        # If we reach here intersection is empty
        return intersection
    end
end # submodule


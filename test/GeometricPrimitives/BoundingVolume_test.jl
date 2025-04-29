
using NeighborGraphs.GeometricPrimitives

# Constructing BVs: --------------------------------------------------------------
@testset "    Constructing BVs: Invalid bounds" begin
    lb = [1, 2]
    ub = [-1, -2]
    @test_throws "Cannot construct bounding volume" BoundingVolume(lb, ub)

    ub = [3]
    @test_throws "points have different dimensions" BoundingVolume(lb, ub)
end

@testset "    Empty BV" begin
    bv = BoundingVolume()

    @test bv.is_empty
    @test bv.dim == 0
    @test bv.lb[1] == Inf
    @test bv.ub[1] == -Inf
    @test length(bv.active_dim) == 0
    @test length(bv.inactive_dim) == 0
    @test length(bv.is_active) == 0
end

@testset "    Full-dimension BV" begin
    bv = BoundingVolume([1,2,3], [4,5,6])

    @test bv.dim == 3
    @test bv.is_empty == false
    @test bv.active_dim == [1,2,3]
    @test length(bv.inactive_dim) == 0
    @test bv.is_active == ones(Bool, 3)
end

@testset "    Low-dimension BV" begin
    bv = BoundingVolume([1, 2, 3], [4, 2, 5])

    @test bv.is_empty == false
    @test bv.dim == 2
    @test bv.active_dim == [1, 3]
    @test bv.inactive_dim == [2]
    @test bv.is_active == [true, false, true]
end

@testset "    Point BV" begin
    bv = BoundingVolume([1, 2, 3], [1, 2, 3])

    @test bv.is_empty == false
    @test bv.dim == 0
    @test bv.active_dim == [ ]
    @test bv.inactive_dim == [1, 2, 3]
    @test bv.is_active == [false, false, false]
end


# `intersects`: ------------------------------------------------------------------
@testset "    intersects(BV, BV): Empty Intersection" begin
    bv1 = BoundingVolume([1, 2], [3, 4])
    bv2 = BoundingVolume([-3, -4], [-1, -2])

    @test intersects(bv1, bv2) == false
    @test intersects(bv1, bv2, include_boundary=true) == false
end

@testset "    intersects(BV, BV): Low-Dim Intersection" begin
    bv1 = BoundingVolume([1, 2], [3, 4])
    bv2 = BoundingVolume([0, 0], [1, 2])

    @test intersects(bv1, bv2, include_boundary=true) == true
    @test intersects(bv1, bv2, include_boundary=false) == false
end

@testset "    intersects(BV, BV): Full-Dim Intersection" begin
    bv1 = BoundingVolume([1, 2], [3, 4])
    bv2 = BoundingVolume([0, 0], [2, 3])

    @test intersects(bv1, bv2) == true
    @test intersects(bv1, bv2, include_boundary=true) == true
end


# `isContained`: -----------------------------------------------------------------
@testset "    isContained(BV, pt)" begin
    bv = BoundingVolume([0,0], [1,1])
    interior_pt = [0.5, 0.5]
    boundary_pt = [1, 0]
    exterior_pt = [2, 2]

    @test isContained(bv, interior_pt, include_boundary=true) == true
    @test isContained(bv, interior_pt, include_boundary=false) == true

    @test isContained(bv, boundary_pt, include_boundary=true) == true
    @test isContained(bv, boundary_pt, include_boundary=false) == false

    @test isContained(bv, exterior_pt, include_boundary=true) == false
    @test isContained(bv, exterior_pt, include_boundary=false) == false
end

@testset "    isContained(BV, BV): Empty Intersection" begin
    bv = BoundingVolume([0, 0], [1, 1])
    bv_query = BoundingVolume([-2, -2], [-1, -1])

    @test isContained(bv, bv_query, include_boundary=true) == false
    @test isContained(bv, bv_query, include_boundary=false) == false
end

@testset "    isContained(BV, BV): Partial Intersections" begin
    bv = BoundingVolume([0, 0], [1, 1])
    bv_query = BoundingVolume([-1, -1], [0.5, 0.5])

    # Full-dim intersection
    @test isContained(bv, bv_query, include_boundary=true) == false
    @test isContained(bv, bv_query, include_boundary=false) == false

    # Low-dimension intersection
    bv_query = BoundingVolume([-1, 0], [0, 1])
    @test isContained(bv, bv_query, include_boundary=true) == false
    @test isContained(bv, bv_query, include_boundary=false) == false
end

@testset "    isContained(BV, BV): Contained" begin
    bv = BoundingVolume([0, 0], [1, 1])
    bv_query = BoundingVolume([0.25, 0.25], [0.75, 0.75])

    @test isContained(bv, bv_query, include_boundary=true) == true
    @test isContained(bv, bv_query, include_boundary=false) == true
end

@testset "    isContained(BV, BV): Strictly/fully Contained" begin
    bv = BoundingVolume([0, 0], [1, 1])
    bv_query = BoundingVolume([0.25, 0.25], [1, 1])

    @test isContained(bv, bv_query, include_boundary=true) == true
    @test isContained(bv, bv_query, include_boundary=false) == false
end

# `getIntersection`: -------------------------------------------------------------
# `getClosestPoint`: -------------------------------------------------------------
# `getFurthestPoint`: ------------------------------------------------------------
# `getFaceBoundingVolume`: -------------------------------------------------------
# `faceIndex2SpatialIndex`: ------------------------------------------------------



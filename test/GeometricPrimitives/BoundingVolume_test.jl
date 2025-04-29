
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

# `getIntersection`: -------------------------------------------------------------
# `getClosestPoint`: -------------------------------------------------------------
# `getFurthestPoint`: ------------------------------------------------------------
# `getFaceBoundingVolume`: -------------------------------------------------------
# `faceIndex2SpatialIndex`: ------------------------------------------------------



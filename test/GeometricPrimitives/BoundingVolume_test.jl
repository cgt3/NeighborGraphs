
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

# @testset "Low-dimension BVs"


# `intersects`: ------------------------------------------------------------------
# `isContained`: -----------------------------------------------------------------
# `getIntersection`: -------------------------------------------------------------
# `getClosestPoint`: -------------------------------------------------------------
# `getFurthestPoint`: ------------------------------------------------------------
# `getFaceBoundingVolume`: -------------------------------------------------------
# `faceIndex2SpatialIndex`: ------------------------------------------------------



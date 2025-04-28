
using NeighborGraphs.GeometricPrimitives

# Constructing BVs: --------------------------------------------------------------
@testset "Constructing BVs: Invalid bounds" begin
    lb = [1, 2]
    ub = [-1, -2]
    @test_throws "Cannot construct bounding volume" BoundingVolume(lb, ub)

    ub = [3]
    @test_throws "points have different dimensions" BoundingVolume(lb, ub)
end

# Invalid bounds
# Empty BVs
# Low-dimension BVs


# `intersects`: ------------------------------------------------------------------
# `isContained`: -----------------------------------------------------------------
# `getIntersection`: -------------------------------------------------------------
# `getClosestPoint`: -------------------------------------------------------------
# `getFurthestPoint`: ------------------------------------------------------------
# `getFaceBoundingVolume`: -------------------------------------------------------
# `faceIndex2SpatialIndex`: ------------------------------------------------------



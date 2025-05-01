using NeighborGraphs.GeometricPrimitives

# Ball Constructors -----------------------------------------------------
@testset "    Constructing Balls: Invalid radius" begin
    @test_throws "Cannot construct ball with negative radius" Ball([0,0], -1)
end

# @testset "    Constructing Balls: Zero Radius" begin 
# end

# @testset "    Constructing Balls: Low-Dimension Ball" begin 
# end


# @testset "    Constructing Balls: Fill-Dim, Euclidean Distance" begin 
# end

# @testset "    Constructing Balls: Full-Dim, abritrary p" begin  
# end


# Ball -> BV ------------------------------------------------------------
# `isContained(Ball, pt)` -----------------------------------------------
# `isContained(BV, Ball)` -----------------------------------------------
# `isContained(Ball, BV)` -----------------------------------------------
# `intersects(BV, Ball)` ------------------------------------------------
# `getReducedDimBall` ---------------------------------------------------
# `tightenBVBounds!` ----------------------------------------------------
# `getIntersection(BV, Ball)` -------------------------------------------
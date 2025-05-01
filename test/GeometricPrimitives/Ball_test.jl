using NeighborGraphs.GeometricPrimitives

# Ball Constructors -----------------------------------------------------
@testset "    Constructing Balls: Invalid radius" begin
    @test_throws "Cannot construct ball with negative radius" Ball([0,0], -1)
end

@testset "    Constructing Balls: Zero Radius" begin 
    ball = Ball([0,0], 0)
    @test ball.radius == 0
    @test ball.dim == 0
    @test ball.embedding_dim == 2
    @test all(ball.active_dim .== [])
    @test all(ball.inactive_dim .== [1,2])
    @test all(ball.is_active .== [false, false])
end

@testset "    Constructing Balls: Low-Dimension Ball" begin 
    ball = Ball([0,0,0], 1, p=2, active_indices=true, indices=[1, 3])
    @test ball.radius == 1
    @test ball.p == 2
    @test ball.dim == 2
    @test ball.embedding_dim == 3
    @test all(ball.active_dim .== [1,3])
    @test all(ball.inactive_dim .== [2])
    @test all(ball.is_active .== [true, false, true])

    ball = Ball([0,0,0], 1, p=2, active_indices=false, indices=[2])
    @test ball.radius == 1
    @test ball.p == 2
    @test ball.dim == 2
    @test ball.embedding_dim == 3
    @test all(ball.active_dim .== [1,3])
    @test all(ball.inactive_dim .== [2])
    @test all(ball.is_active .== [true, false, true])
end


@testset "    Constructing Balls: Fill-Dim, Euclidean Distance" begin 
    ball = Ball([0,0,0], 1)

    @test ball.radius == 1
    @test ball.p == 2
    @test ball.dim == 3
    @test ball.embedding_dim == 3
    @test all(ball.active_dim .== [1,2,3])
    @test all(ball.inactive_dim .== [])
    @test all(ball.is_active .== [true, true, true])
end

# @testset "    Constructing Balls: Full-Dim, abritrary p" begin  
    ball = Ball([0,0,0], 1, p=5)

    @test ball.radius == 1
    @test ball.p == 5
    @test ball.dim == 3
    @test ball.embedding_dim == 3
    @test all(ball.active_dim .== [1,2,3])
    @test all(ball.inactive_dim .== [])
    @test all(ball.is_active .== [true, true, true])
# end


# Ball -> BV ------------------------------------------------------------
# `isContained(Ball, pt)` -----------------------------------------------
# `isContained(BV, Ball)` -----------------------------------------------
# `isContained(Ball, BV)` -----------------------------------------------
# `intersects(BV, Ball)` ------------------------------------------------
# `getReducedDimBall` ---------------------------------------------------
# `tightenBVBounds!` ----------------------------------------------------
# `getIntersection(BV, Ball)` -------------------------------------------
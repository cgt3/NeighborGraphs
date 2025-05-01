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

@testset "    Constructing Balls: Full-Dim, abritrary p" begin  
    ball = Ball([0,0,0], 1, p=Inf)

    @test ball.radius == 1
    @test ball.p == Inf
    @test ball.dim == 3
    @test ball.embedding_dim == 3
    @test all(ball.active_dim .== [1,2,3])
    @test all(ball.inactive_dim .== [])
    @test all(ball.is_active .== [true, true, true])
end


# Ball -> BV ------------------------------------------------------------


# `isContained(Ball, pt)` -----------------------------------------------
@testset "    isContained(Ball, pt): Point was wrong dim" begin
    ball = Ball([0,0,0], 1)
    pt = [0,0]
    
    @test_throws "Point dimension(=2) does not match ball embedding dimension(=3)" isContained(ball, pt)
end

@testset "    isContained(Ball, pt): Zero-Dim Ball" begin
    ball_p1   = Ball([0,0], 0, p=1)
    ball_p2   = Ball([0,0], 0, p=2)
    ball_pInf = Ball([0,0], 0, p=Inf)

    pt = [0.0, 0.0]
    @test isContained(ball_p1,   pt, include_boundary=true)
    @test isContained(ball_p2,   pt, include_boundary=true)
    @test isContained(ball_pInf, pt, include_boundary=true)

    @test !isContained(ball_p1,   pt, include_boundary=false)
    @test !isContained(ball_p2,   pt, include_boundary=false)
    @test !isContained(ball_pInf, pt, include_boundary=false)
end

@testset "    isContained(Ball, pt): Low-Dim Ball - Interior point" begin
    ball_p1   = Ball([0,0,0], 1, p=1, active_indices=true, indices=[1,2])
    ball_p2   = Ball([0,0,0], 1, p=2, active_indices=true, indices=[1,2])
    ball_pInf = Ball([0,0,0], 1, p=Inf, active_indices=true, indices=[1,2])

    pt = [0.25, 0.25, 0]

    @test isContained(ball_p1,   pt, include_boundary=true)
    @test isContained(ball_p2,   pt, include_boundary=true)
    @test isContained(ball_pInf, pt, include_boundary=true)

    @test isContained(ball_p1,   pt, include_boundary=false)
    @test isContained(ball_p2,   pt, include_boundary=false)
    @test isContained(ball_pInf, pt, include_boundary=false)
end

@testset "    isContained(Ball, pt): Low-Dim Ball - Boundary point" begin
    ball_p1   = Ball([0,0,0], 1, p=1, active_indices=true, indices=[1,2])
    ball_p2   = Ball([0,0,0], 1, p=2, active_indices=true, indices=[1,2])
    ball_pInf = Ball([0,0,0], 1, p=Inf, active_indices=true, indices=[1,2])

    pt = [0, 1, 0]

    @test isContained(ball_p1,   pt, include_boundary=true)
    @test isContained(ball_p2,   pt, include_boundary=true)
    @test isContained(ball_pInf, pt, include_boundary=true)

    @test !isContained(ball_p1,   pt, include_boundary=false)
    @test !isContained(ball_p2,   pt, include_boundary=false)
    @test !isContained(ball_pInf, pt, include_boundary=false)
end


@testset "    isContained(Ball, pt): Low-Dim Ball - Exterior point" begin
    ball_p1   = Ball([0,0,0], 1, p=1, active_indices=true, indices=[1,2])
    ball_p2   = Ball([0,0,0], 1, p=2, active_indices=true, indices=[1,2])
    ball_pInf = Ball([0,0,0], 1, p=Inf, active_indices=true, indices=[1,2])

    # Out-of-plane point
    pt = [0, 0, 1]
    @test !isContained(ball_p1,   pt, include_boundary=true)
    @test !isContained(ball_p2,   pt, include_boundary=true)
    @test !isContained(ball_pInf, pt, include_boundary=true)

    @test !isContained(ball_p1,   pt, include_boundary=false)
    @test !isContained(ball_p2,   pt, include_boundary=false)
    @test !isContained(ball_pInf, pt, include_boundary=false)


    # In-plane exterior point
    pt = [2, 2, 0]
    @test !isContained(ball_p1,   pt, include_boundary=true)
    @test !isContained(ball_p2,   pt, include_boundary=true)
    @test !isContained(ball_pInf, pt, include_boundary=true)

    @test !isContained(ball_p1,   pt, include_boundary=false)
    @test !isContained(ball_p2,   pt, include_boundary=false)
    @test !isContained(ball_pInf, pt, include_boundary=false)
end

@testset "    isContained(Ball, pt): Full-Dim Ball - Interior point" begin
    ball_p1   = Ball([0,0], 1, p=1)
    ball_p2   = Ball([0,0], 1, p=2)
    ball_pInf = Ball([0,0], 1, p=Inf)

    pt = [0.25, 0.25]

    @test isContained(ball_p1,   pt, include_boundary=true)
    @test isContained(ball_p2,   pt, include_boundary=true)
    @test isContained(ball_pInf, pt, include_boundary=true)

    @test isContained(ball_p1,   pt, include_boundary=false)
    @test isContained(ball_p2,   pt, include_boundary=false)
    @test isContained(ball_pInf, pt, include_boundary=false)
end

@testset "    isContained(Ball, pt): Full-Dim Ball - Boundary point" begin
    ball_p1   = Ball([0,0], 1, p=1)
    ball_p2   = Ball([0,0], 1, p=2)
    ball_pInf = Ball([0,0], 1, p=Inf)

    pt = [0, 1]

    @test isContained(ball_p1,   pt, include_boundary=true)
    @test isContained(ball_p2,   pt, include_boundary=true)
    @test isContained(ball_pInf, pt, include_boundary=true)

    @test !isContained(ball_p1,   pt, include_boundary=false)
    @test !isContained(ball_p2,   pt, include_boundary=false)
    @test !isContained(ball_pInf, pt, include_boundary=false)
end


@testset "    isContained(Ball, pt): Full-Dim Ball - Exterior point" begin
    ball_p1   = Ball([0,0], 1, p=1)
    ball_p2   = Ball([0,0], 1, p=2)
    ball_pInf = Ball([0,0], 1, p=Inf)

    pt = [2,2]

    @test !isContained(ball_p1,   pt, include_boundary=true)
    @test !isContained(ball_p2,   pt, include_boundary=true)
    @test !isContained(ball_pInf, pt, include_boundary=true)

    @test !isContained(ball_p1,   pt, include_boundary=false)
    @test !isContained(ball_p2,   pt, include_boundary=false)
    @test !isContained(ball_pInf, pt, include_boundary=false)
end

@testset "    isContained(Ball, pt): Differing p" begin
    ball_p1   = Ball([0,0], 1, p=1)
    ball_p2   = Ball([0,0], 1, p=2)
    ball_pInf = Ball([0,0], 1, p=Inf)

    pt = [0.75, 0.75]
    @test !isContained(ball_p1,  pt, include_boundary=true)
    @test !isContained(ball_p2,  pt, include_boundary=true)
    @test isContained(ball_pInf, pt, include_boundary=true)

    pt = [0.51, 0.51]
    @test !isContained(ball_p1, pt, include_boundary=true)
    @test isContained(ball_p2, pt, include_boundary=true)
end

# `isContained(BV, Ball)` -----------------------------------------------
# `isContained(Ball, BV)` -----------------------------------------------
# `intersects(BV, Ball)` ------------------------------------------------
# `getReducedDimBall` ---------------------------------------------------
# `tightenBVBounds!` ----------------------------------------------------
# `getIntersection(BV, Ball)` -------------------------------------------
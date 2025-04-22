using SafeTestsets

println("Testing submodule GeometricPrimitives:")
@safetestset "  BoundingVolumes" begin include("GeometricPrimitives/BoundingVolume_test.jl") end
@safetestset "  Balls" begin include("GeometricPrimitives/Ball_test.jl") end

println("Testing submodule kdTrees:")
@safetestset "  Base Functionality" begin include("kdTrees/kdTrees_test.jl") end
@safetestset "  Searches" begin include("kdTrees/search_test.jl") end
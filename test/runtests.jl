using SafeTestsets

println("Testing submodule GeometricPrimitives:")
@safetestset "  GeometricPrimitive: BoundingVolumes:" begin include("GeometricPrimitives/BoundingVolume_test.jl") end
@safetestset "  GeometricPrimitive: Balls" begin include("GeometricPrimitives/Ball_test.jl") end

println("\nTesting submodule kdTrees:")
@safetestset "  kdTrees: Base Functionality" begin include("kdTrees/kdTrees_test.jl") end
@safetestset "  kdTrees: Searches" begin include("kdTrees/search_test.jl") end
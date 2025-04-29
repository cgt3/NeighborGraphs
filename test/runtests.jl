using SafeTestsets

@safetestset "GeometricPrimitives:" begin 
    @safetestset "  BoundingVolumes:" begin include("GeometricPrimitives/BoundingVolume_test.jl") end
    @safetestset "  Balls:" begin include("GeometricPrimitives/Ball_test.jl") end
end

@safetestset "kdTrees:" begin 
    @safetestset "  Base functionality: " begin include("kdTrees/kdTrees_test.jl") end
    @safetestset "  Searches" begin include("kdTrees/search_test.jl") end
end

using BubbleBath
using Distributions: Uniform
using LinearAlgebra: norm
using Test

@testset "BubbleBath.jl" begin
    @testset "Spheres" begin
        radius = 1
        pos = ntuple(_ -> 5.0, 3)
        # sphere dimensionality must match pos length
        @test_throws MethodError Sphere{2}(pos, radius)
        # if sphere dimensionality is not specified, inherit from pos
        sphere = Sphere(pos, radius)
        @test sphere isa Sphere{3}
        # check fields are assigned correctly
        @test sphere.pos == pos
        @test sphere.radius == radius
        # should work identically when pos is NTuple{D,Float64} or NTuple{D,Int}
        pos = ntuple(_ -> 5, 3)
        sphere2 = Sphere(pos, radius)
        @test sphere2.pos == sphere.pos
    end
    @testset "BubbleBath algorithm" begin
        L = 10
        extent = ntuple(_ -> L, 3)
        r = 4.0
        bath = bubblebath([r], extent)
        # should be a vector with only one sphere
        @test bath isa Vector{Sphere{3}}
        @test length(bath) == 1
        # sphere radius should be r
        @test bath[1].radius == r
        # sphere position should not overlap with domain boundaries
        for i in 1:3
            @test r ≤ bath[1].pos[i] ≤ L-r
        end
        # packing fraction should be volume of the sphere over domain volume
        ϕ = packing_fraction(bath, extent)
        @test ϕ ≈ (4π*r^3/3)/prod(extent)

        L = 10
        extent = ntuple(_ -> L, 3)
        r = 8.0
        # with a large sphere radius (r > L/2),
        # if through_boundaries = false (default), sphere can't be inserted
        bath = bubblebath([r], extent)
        @test isempty(bath)
        # if through_boundaries = true, it will surely cross all domain boundaries
        bath = bubblebath([r], extent; through_boundaries=true)
        for i in 1:3
            @test !(r ≤ bath[1].pos[i] ≤ L-r)
        end

        extent = (10, 15, 12)
        radius_pdf = Uniform(1, 2)
        ϕ_max = 0.4
        bath = bubblebath(radius_pdf, ϕ_max, extent)
        # all the spheres should have radius between 1 and 2
        @test all(map(s -> 1 ≤ s.radius ≤ 2, bath))
        # no sphere should cross domain boundaries
        @test all(map(s -> BubbleBath.is_inside_boundaries(s.pos, s.radius, extent), bath))
        # packing fraction is below ϕ_max
        ϕ = packing_fraction(bath, extent)
        @test ϕ ≤ ϕ_max

        extent = ntuple(_ -> 10, 3)
        radius_pdf = [2.0]
        ϕ_max = 0.4
        min_distance = 0.5
        bath = bubblebath(radius_pdf, ϕ_max, extent; min_distance)
        # all spheres should be at distance > min_distance (between their surfaces)
        surface_distances = vec([
            norm(bath[i].pos .- bath[j].pos) - (bath[i].radius + bath[j].radius)
            for i in eachindex(bath), j in eachindex(bath) if j>i
        ])
        @test all(surface_distances .≥ min_distance)
    end
end

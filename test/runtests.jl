using BubbleBath
using Distributions: Uniform
using LinearAlgebra: norm
using Test

function capture_stderr(f, args, kwargs)
    original_stderr = stderr
    out_read, out_write = redirect_stderr()
    f(args...; kwargs...)
    # without this the program hangs if f does not write to stderr
    @info "dummy text"
    close(out_write)
    data = readavailable(out_read)
    close(out_read)
    s = String(copy(data))
    redirect_stderr(original_stderr)
    return s
end

@testset "BubbleBath.jl" begin
    @testset "Spheres" begin
        # sphere dimensionality is always inherited from pos
        radius = 1
        pos = ntuple(_ -> 5.0, 2)
        sphere = Sphere(pos, radius)
        @test sphere isa Sphere{2}
        pos = ntuple(_ -> 5.0, 3)
        sphere = Sphere(pos, radius)
        @test sphere isa Sphere{3}
        # check fields are assigned correctly
        @test sphere.pos == pos
        @test sphere.radius == radius
        # should work identically when pos is NTuple{D,Float64} or NTuple{D,Int}
        pos = ntuple(_ -> 5, 3)
        sphere2 = Sphere(pos, radius)
        @test sphere2.pos == sphere.pos
        # if pos is an AbstractVector it should be converted to a tuple
        pos = rand(3)
        sphere = Sphere(pos, 1)
        @test sphere.pos == Tuple(pos)
        pos = 1:5
        sphere = Sphere(pos, 1)
        @test sphere.pos == Tuple(pos)
        # non-positive radius not allowed
        @test_throws ArgumentError Sphere((5,5), 0)
        @test_throws ArgumentError Sphere((5,5), -1)
    end

    @testset "Bubblebath" begin
        L = 10
        extent = ntuple(_ -> L, 3)
        r = 4.0
        # negative radii not allowed
        @test_throws ArgumentError bubblebath([-r], extent)
        # packing fraction must be ϕ∈(0,1]
        @test_throws ArgumentError bubblebath([r], -0.1, extent)
        @test_throws ArgumentError bubblebath([r], 1.1, extent)
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
        # if max fails are reached, should print info message
        msg = capture_stderr(bubblebath,
            ([r,r], extent),
            (max_tries=0, max_fails=0)
        )
        @test contains(msg, "Reached max. number of tries")
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

        # test the verbose keyword
        L = 10
        extent = (L,L)
        radius_pdf = Uniform(1,2)
        ϕ_max = 0.3
        msg = capture_stderr(bubblebath,
            (radius_pdf, ϕ_max, extent),
            (verbose=true,) # default
        )
        @test contains(msg, r"Generated \d+ spheres")
        @test contains(msg, r"\d+/\d+ new spheres inserted")
        msg = capture_stderr(bubblebath,
            (radius_pdf, ϕ_max, extent),
            (verbose=false,)
        )
        @test ~contains(msg, r"Generated \d+ spheres")
        @test ~contains(msg, r"\d+/\d+ new spheres inserted")
    end

    @testset "In-place Bubbleath" begin
        L = 50
        extent = (L,L)
        # initialize with one sphere of radius 3
        spheres = [Sphere((L/2,L/2), 3.0)]
        # add 10 more spheres of radius 1
        r = 1.0
        radii = repeat([r], 10)
        bubblebath!(spheres, radii, extent)
        # should be a total of 11 spheres
        @test length(spheres) == 11
        @test count(map(s -> s.radius==3, spheres)) == 1
        @test count(map(s -> s.radius==1, spheres)) == 10
        # the new spheres should not overlap with the original one
        overlaps = [
            norm(spheres[1].pos .- s.pos) ≤ spheres[1].radius + s.radius
            for s in spheres[2:end]
        ]
        @test !(any(overlaps))
        
        L = 10
        extent = (L,L)
        r = L/2
        spheres_old = [Sphere((L/2,L/2), r)] # largest sphere to fit extent
        spheres_new = copy(spheres_old)
        # should add nothing since r is too large to fit another sphere
        bubblebath!(spheres_new, [r], extent)
        @test spheres_new == spheres_old
        # if max fails are reached, should print info message
        max_tries = 1
        max_fails = 1
        msg = capture_stderr(bubblebath!,
            (spheres_new, [r,r], extent),
            (max_tries=0, max_fails=0)
        )
        @test contains(msg, "Reached max. number of tries")
    end

    @testset "Walkmap" begin
        extent = (10, 10)
        pos = (5,5)
        r = 3
        spheres = [Sphere(pos, r)]
        res = 0.1
        wm = walkmap(spheres, extent, res)
        # walkmap has same dimensions as extent
        @test ndims(wm) == length(extent)

        extent = (5.07, 12.0, 8.15)
        pos = (3, 3, 3)
        r = 1
        spheres = [Sphere(pos, r)]
        res = 0.1
        wm = walkmap(spheres, extent, res)
        # res defines the size of wm
        n_nodes = @. floor(Int, extent / res)
        @test size(wm) == n_nodes

        extent = (10, 10)
        pos = (5, 5)
        r = 3
        spheres = [Sphere(pos, r)]
        res = 0.1
        wm = walkmap(spheres, extent, res)
        xs = range(res/2, extent[1]-res/2; step=res)
        ys = range(res/2, extent[2]-res/2; step=res)
        grid = (xs, ys)
        function getidx(pos, grid, res)
            ntuple(i ->
                findfirst(j -> grid[i][j]-res/2 ≤ pos[i] ≤ grid[i][j]+res/2,
                    eachindex(grid[i])
                ),
                length(pos)
            )
        end
        # wm should be 0 in positions occupied by spheres
        p₁ = pos
        p₂ = pos .+ (r-res/2, 0)
        p₃ = pos .+ (-r/6, r/6)
        I₁, I₂, I₃ = map(p -> CartesianIndex(getidx(p,grid,res)), (p₁,p₂,p₃))
        @test ~wm[I₁]
        @test ~wm[I₂]
        @test ~wm[I₃]
        # 1 in unoccupied positions
        p₄ = pos .+ (r+res/2, 0)
        p₅ = pos .+ (0, r+1)
        I₄, I₅ = map(p -> CartesianIndex(getidx(p,grid,res)), (p₄,p₅))
        @test wm[I₄]
        @test wm[I₅]
        # if probe_radius>0 the occupied region increases
        wm₂ = walkmap(spheres, extent, res, 1.0)
        I = CartesianIndex(getidx(pos.+(r+0.5,0), grid, res))
        @test wm[I] && ~wm₂[I]

        # put a sphere at the boundary
        extent = (10, 10)
        pos = (0, 5)
        r = 4
        spheres = [Sphere(pos, r)]
        res = 0.1
        xs = range(res/2, extent[1]-res/2; step=res)
        ys = range(res/2, extent[2]-res/2; step=res)
        grid = (xs, ys)
        wm₁ = walkmap(spheres, extent, res) # boundaries = :cut
        wm₂ = walkmap(spheres, extent, res; boundaries=:wrap)
        # if boundaries=:cut the opposite edge is free
        I = CartesianIndex(getidx((10-res,5), grid, res))
        @test wm₁[I]
        # if boundaries=:wrap the oppposite edge is occupied
        @test ~wm₂[I]
    end

    @testset "Packing fraction" begin
        # packing fraction should match theoretical values
        L = 10
        extent = (L,L)
        r = 4
        bath = bubblebath([r], extent)
        @test packing_fraction(bath, extent) ≈ π*r^2 / L^2 
        extent = (L,L,L)
        bath = bubblebath([r], extent)
        @test packing_fraction(bath, extent) ≈ 4π*r^3/3 / L^3
        # test same behavior on radii vector instead of bath
        @test packing_fraction([r], (L,L)) ≈ π*r^2 / L^2
        @test packing_fraction([r], (L,L,L)) ≈ 4π*r^3/3 / L^3
        # an empty collection should give 0 packing fraction
        @test packing_fraction(Sphere{2}[], (L,L)) == 0
        @test packing_fraction(Float64[], (L,L)) == 0

        # packing fractions produced by bubblebath should never be > ϕ_max
        extent = (8, 10)
        radius_pdf = Uniform(2,5)
        ϕ_max = 0.2
        bath = bubblebath(radius_pdf, ϕ_max, extent)
        @test packing_fraction(bath, extent) ≤ ϕ_max
        extent = (8, 10, 12)
        radius_pdf = Uniform(2,5)
        ϕ_max = 0.2
        bath = bubblebath(radius_pdf, ϕ_max, extent)
        @test packing_fraction(bath, extent) ≤ ϕ_max

        # pre-initialize a bath
        extent = (15, 15)
        r = 3
        bath = [Sphere((7.5,7.5), r)]
        ϕ₀ = packing_fraction(bath, extent) # ≈ 0.126
        # fill with more spheres
        radius_pdf = [0.1]
        ϕ_max = 0.3
        bubblebath!(bath, radius_pdf, ϕ_max, extent)
        # final ϕ will be ϕ_max+ϕ₀ ≈ 0.426
        @test packing_fraction(bath, extent) ≈ ϕ_max+ϕ₀ atol=0.02
        bath = [Sphere((7.5,7.5), r)]
        bubblebath!(bath, radius_pdf, ϕ_max-ϕ₀, extent)
        # now final ϕ will be ϕ_max
        @test packing_fraction(bath, extent) ≈ ϕ_max atol=0.02

        # packing fraction of a walkmap
        extent = (10, 10)
        r = 3
        spheres = [Sphere((5,5), r)]
        ϕ₁ = packing_fraction(spheres, extent)
        res = 0.1
        wm = walkmap(spheres, extent, res)
        ϕ₂ = packing_fraction(wm)
        @test ϕ₂ ≈ ϕ₁ atol=res^2
        # if half of a sphere goes through a boundary
        # packing fraction of the walkmap is correct
        spheres = [Sphere((10,5), r)]
        ϕ = (π*r^2 / prod(extent)) / 2 # exact packing fraction
        ϕ₀ = packing_fraction(spheres, extent)
        wm₁ = walkmap(spheres, extent, res)
        wm₂ = walkmap(spheres, extent, res; boundaries=:wrap)
        ϕ₁ = packing_fraction(wm₁)
        ϕ₂ = packing_fraction(wm₂)
        @test ϕ₀ ≈ 2ϕ
        @test ϕ₁ ≈ ϕ atol=res^2
        @test ϕ₂ ≈ 2ϕ atol=res^2
    end
end

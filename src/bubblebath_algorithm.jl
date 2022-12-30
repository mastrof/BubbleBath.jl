export bubblebath, bubblebath!

"""
    bubblebath(
        radius_pdf, ϕ_max::Real, extent::NTuple{D,Real};
        through_boundaries = false,
        max_tries = 10000, max_fails = 100,
        verbose = true
    ) where D
Generate a bath of spheres in the domain `extent`,
extracting radii from `radius_pdf` trying to reach a target
packing fraction `ϕ_max`.
The domain is filled with spheres in order of decreasing radius.

## Keywords
* `min_distance = 0.0`: Minimum allowed distance between spheres.
* `through_boundaries = false`: Whether spheres can cross through the domain boundaries.
* `max_tries = 10000`: Maximum number of insertion tries for each sphere.
    When `max_tries` is reached, the sphere is discarded and the algorithm moves on
    to the next one.
* `max_fails = 100`: Maximum number of failures (i.e. discarded spheres) allowed.
    Once `max_fails` is reached, the program halts.
* `verbose = true`: Whether info logs should be printed.
"""
function bubblebath(
    radius_pdf,
    ϕ_max::Real,
    extent::NTuple{D,Real};
    min_distance::Real = 0.0,
    through_boundaries = false,
    max_tries = 10000,
    max_fails = 100,
    verbose = true
)::Vector{Sphere{D}} where D
    radii = generate_radii(radius_pdf, ϕ_max, extent; max_tries, verbose)
    bubblebath(radii, extent;
        min_distance, through_boundaries, max_tries, max_fails, verbose
    )
end

"""
    generate_radii(
        radius_pdf, ϕ_max::Real, extent::NTuple{D,Real};
        max_tries = 10000, verbose = true
    ) where D
Generate a vector of radii from the `radius_pdf` distribution,
with a limit packing fraction `ϕ_max` in the domain `extent`.
"""
function generate_radii(
    radius_pdf,
    ϕ_max::Real,
    extent::NTuple{D,Real};
    max_tries = 10000,
    verbose = true
)::Vector{Float64} where D
    radii = Float64[]
    V₀ = prod(extent)
    tries = 0
    while true
        r = rand(radius_pdf)
        V = volume(r, D)
        ϕ = (sum(volume.(radii,D))+V)/V₀
        if ϕ ≤ ϕ_max
            push!(radii, r)
        else
            tries += 1
            tries > max_tries && break
        end
    end
    if verbose
        @info "Generated $(length(radii)) spheres."
    end
    return radii
end

"""
    bubblebath(
        radii::Vector{<:Real}, extent::NTuple{D,Real};
        min_distance = 0.0,
        through_boundaries = false,
        max_tries = 10000, max_fails = 100,
        verbose = true
    ) where D
Generate a bath of spheres with radii `radii` in the domain `extent`.
"""
function bubblebath(
    radii::Vector{<:Real},
    extent::NTuple{D,Real};
    min_distance::Real = 0.0,
    through_boundaries = false,
    max_tries = 10000,
    max_fails = 100,
    verbose = true
)::Vector{Sphere{D}} where D
    spheres = Sphere{D}[]
    bubblebath!(spheres, radii, extent;
        min_distance, through_boundaries, max_tries, max_fails, verbose
    )
    return spheres
end

"""
    is_overlapping(p::NTuple{D,Real}, r::Real, spheres::Vector{Sphere{D}}) where D
Test if a sphere of radius `r` centered at `p` overlaps with any sphere in `spheres`.
Surface contact is not counted as overlap.
"""
function is_overlapping(
    p::NTuple{D,Real}, r::Real,
    spheres::Vector{Sphere{D}}
)::Bool where D
    for sphere in spheres
        if is_overlapping(p, r, sphere)
            return true
        end
    end
    return false
end

"""
    is_overlapping(p₁::NTuple{D,Real}, r₁::Real, s₂::Sphere{D}) where D
Test if a sphere of radius `r₁` centered at `p₁` overlaps with sphere `s₂`.
Surface contact is not counted as overlap.
"""
@inline function is_overlapping(
    p₁::NTuple{D,Real}, r₁::Real,
    s₂::Sphere{D}
)::Bool where D
    is_overlapping(p₁, r₁, s₂.pos, s₂.radius)
end

"""
    is_overlapping(p₁::NTuple{D,Real}, r₁::Real, p₂::NTuple{D,Real}, r₂::Real) where D
Test if two spheres with radii `r₁` and `r₂`, centered at `p₁` and `p₂` respectively,
are overlapping. Surface contact is not counted as overlap.
"""
@inline function is_overlapping(
    p₁::NTuple{D,Real}, r₁::Real,
    p₂::NTuple{D,Real}, r₂::Real
)::Bool where D
    norm(p₁ .- p₂) < r₁ + r₂
end

"""
    is_inside_boundaries(pos::NTuple{D,Real}, radius::Real, extent::NTuple{D,Real}) where D
Check if a sphere of radius `radius` centered at `pos` is within domain `extent`.
"""
function is_inside_boundaries(
    pos::NTuple{D,Real}, radius::Real, extent::NTuple{D,Real}
)::Bool where D
    for i in 1:D
        if !(radius ≤ pos[i] ≤ extent[i]-radius)
            return false
        end
    end
    return true
end


"""
    bubblebath!(
        spheres::Vector{Sphere{D}},
        radius_pdf, ϕ_max::Real, extent::NTuple{D,Real};
        min_distance::Real = 0.0, through_boundaries = false,
        max_tries = 10000, max_fails = 100,
        verbose = true
    ) where D
In-place version of `bubblebath`, adds new spheres to the `spheres`
vector, which can be already populated.

Here, `ϕ_max` does **not** account for spheres that might
already be present in the `spheres` vector.
E.g. if `packing_fraction(spheres, extent)` is 0.2 and `ϕ_max=0.3`,
then the algorithm generates new spheres for a packing fraction of 0.3,
which upon insertion will (try to) add up to a total packing fraction of 0.5.
To account for the pre-initialized spheres, decrease `ϕ_max` accordingly
(`ϕ_max = 0.3-packing_fraction(spheres,extent)` giving 0.1 for this example).
"""
function bubblebath!(
    spheres::Vector{Sphere{D}},
    radius_pdf,
    ϕ_max::Real,
    extent::NTuple{D,Real};
    min_distance::Real = 0.0,
    through_boundaries = false,
    max_tries = 10000,
    max_fails = 100,
    verbose = true
)::Nothing where D
    radii = generate_radii(radius_pdf, ϕ_max, extent; max_tries, verbose)
    bubblebath!(spheres, radii, extent;
        min_distance, through_boundaries, max_tries, max_fails, verbose
    )
end

"""
    bubblebath!(
        spheres::Vector{Sphere{D}},
        radii::Vector{<:Real}, extent::NTuple{D,Real};
        min_distance::Real = 0.0, through_boundaries = false,
        max_tries = 10000, max_fails = 100,
        verbose = true
    ) where D
In-place version of `bubblebath`, adds new spheres to the
`spheres` vector (which can be already populated).
"""
function bubblebath!(
    spheres::Vector{Sphere{D}},
    radii::Vector{<:Real},
    extent::NTuple{D,Real};
    min_distance::Real = 0.0,
    through_boundaries = false,
    max_tries = 10000,
    max_fails = 100,
    verbose = true
)::Nothing where D
    n₀ = length(spheres)
    sizehint!(spheres, length(spheres)+n₀)
    fails = 0
    for radius in sort(radii, rev=true)
        tries = 0
        Δ = through_boundaries ? 0.0 : radius
        while true
            if tries > max_tries
                if fails > max_fails
                    @warn "Reached max. number of tries. Interrupting."
                    @goto packing_complete
                else
                    fails += 1
                    break
                end
            end
            pos = Δ .+ Tuple(rand(D)) .* (extent .- 2Δ)
            isvalid_pos = (
                !is_overlapping(pos, radius+min_distance, spheres) &&
                (through_boundaries || is_inside_boundaries(pos, radius, extent))
            )
            if isvalid_pos
                push!(spheres, Sphere(pos, radius))
                break
            else
                tries += 1
            end
        end
    end
    @label packing_complete
    if verbose
        @info "$(length(spheres)-n₀)/$(length(radii)) new spheres inserted."
    end
    return nothing
end
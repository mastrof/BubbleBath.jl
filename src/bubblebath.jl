export bubblebath

"""
    bubblebath(
        radius_pdf, ϕ_max::Real, extent::NTuple{D,Real};
        through_boundaries = false,
        max_tries = 10000, max_fails = 100
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
"""
function bubblebath(
    radius_pdf,
    ϕ_max::Real,
    extent::NTuple{D,Real};
    min_distance::Real = 0.0,
    through_boundaries = false,
    max_tries = 10000,
    max_fails = 100
)::Vector{Sphere{D}} where D
    radii = generate_radii(radius_pdf, ϕ_max, extent)
    bubblebath(radii, extent;
        min_distance, through_boundaries, max_tries, max_fails
    )
end

"""
    generate_radii(
        radius_pdf, ϕ_max::Real, extent::NTuple{D,Real}
    ) where D
Generate a vector of radii from the `radius_pdf` distribution,
with a limit packing fraction `ϕ_max` in the domain `extent`.
"""
function generate_radii(
    radius_pdf,
    ϕ_max::Real,
    extent::NTuple{D,Real}
)::Vector{Float64} where D
    radii = Float64[]
    V₀ = prod(extent)
    while true
        r = rand(radius_pdf)
        V = volume(r, D)
        ϕ = isempty(radii) ? 0.0 : (sum(volume.(radii,D))+V)/V₀
        if ϕ ≤ ϕ_max
            push!(radii, r)
        else
            break
        end
    end
    return radii
end

"""
    bubblebath(
        radii::Vector{<:Real}, extent::NTuple{D,Real};
        min_distance = 0.0,
        through_boundaries = false,
        max_tries = 10000, max_fails = 100
    ) where D
Generate a bath of spheres with radii `radii` in the domain `extent`.
"""
function bubblebath(
    radii::Vector{<:Real},
    extent::NTuple{D,Real};
    min_distance::Real = 0.0,
    through_boundaries = false,
    max_tries = 10000,
    max_fails = 100
)::Vector{Sphere{D}} where D
    spheres = Sphere{D}[]
    sizehint!(spheres, length(radii))
    fails = 0
    for radius in sort(radii, rev=true)
        tries = 0
        Δ = through_boundaries ? 0.0 : radius
        while true
            if tries > max_tries
                if fails > max_fails
                    @info "Reached max. number of tries. Interrupting."
                    @goto packing_complete
                else
                    fails += 1
                    break
                end
            end
            pos = Δ .+ Tuple(rand(D)) .* (extent .- 2Δ)
            if is_overlapping(pos, radius+min_distance, spheres)
                tries += 1
            else
                push!(spheres, Sphere(pos, radius))
                break
            end
        end
    end
    @label packing_complete
    @info "$(length(spheres))/$(length(radii)) spheres inserted."
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
export walkmap

"""
    walkmap(
        spheres::AbstractVector{Sphere{D}}, extent::NTuple{D,<:Real},
        resolution::Real, probe_radius::Real = 0;
        boundaries::Symbol = :cut
    ) where D
Generate a walkmap for the given configuration of `spheres` in the domain
`extent`, with the desired `resolution.`
A positive `probe_radius` restricts the walkable space assuming that the
"probe" has a finite size. 

`boundaries` can be set to `:cut` or `:wrap` to define how to deal with
spheres that cross through the domain boundaries (in case there is any).
"""
function walkmap(
    spheres::AbstractVector{Sphere{D}},
    extent::NTuple{D,<:Real},
    resolution::Real,
    probe_radius::Real = 0;
    boundaries::Symbol = :cut
)::BitArray{D} where D
    grid = ntuple(i -> (resolution/2):resolution:(extent[i]-resolution/2), D)
    itr = Iterators.product(grid...)
    BitArray{D}([is_walkable(pos, probe_radius, spheres, extent, boundaries) for pos in itr])
end

"""
    is_walkable(
        pos::NTuple{D,<:Real}, r::Real, spheres::AbstractVector{Sphere{D}},
        extent::NTuple{D,<:Real}, boundaries::Symbol
    ) where D
Determines whether an object of size `r` can occupy position `pos`
in a domain `extent` filled by `spheres`.
"""
function is_walkable(
    pos::NTuple{D,<:Real}, r::Real,
    spheres::AbstractVector{Sphere{D}},
    extent::NTuple{D,<:Real}, boundaries::Symbol
)::Bool where D
    if boundaries == :cut
        return is_walkable(pos, r, spheres)
    elseif boundaries == :wrap
        return is_walkable_periodic(pos, r, spheres, extent)
    else
        throw(ArgumentError(
            "Mode $(boundaries) unrecognized. Choose between `:cut` and `:wrap`."
        ))
    end
end

function is_walkable(
    pos::NTuple{D,<:Real}, r::Real, spheres::AbstractVector{Sphere{D}}
)::Bool where D
    for sphere in spheres
        if ~is_walkable(pos, r, sphere)
            return false
        end
    end
    return true
end
function is_walkable(pos::NTuple{D,<:Real}, r::Real, sphere::Sphere{D})::Bool where D
    return norm(pos .- sphere.pos) ≥ r + sphere.radius
end

function is_walkable_periodic(
    pos::NTuple{D,<:Real}, r::Real,
    spheres::AbstractVector{Sphere{D}}, extent::NTuple{D,<:Real}
)::Bool where D
    for sphere in spheres
        if ~is_walkable_periodic(pos, r, sphere, extent)
            return false
        end
    end
    return true
end
function is_walkable_periodic(
    pos::NTuple{D,<:Real}, r::Real,
    sphere::Sphere{D}, extent::NTuple{D,<:Real}
)::Bool where D
    a = @. (pos - sphere.pos) / extent
    d = @. (a - round(a)) * extent # minimum-image distance
    return norm(d) ≥ r + sphere.radius
end
export packing_fraction

"""
    packing_fraction(spheres::Vector{Sphere{D}}, extent::NTuple{D,Float64}) where D
Evaluate the packing fraction of `spheres` in domain `extent`.
"""
function packing_fraction(
    spheres::Vector{Sphere{D}},
    extent::NTuple{D,Float64}
)::Float64 where D
    V₀ = prod(extent)
    V = 0.0
    for sphere in spheres
        V += volume(sphere)
    end
    return V/V₀
end

"""
    volume(sphere::Sphere{2})
Evaluate volume of two-dimensional sphere, i.e. the area of a circle (πr²).
"""
volume(sphere::Sphere{2})::Float64 = π * sphere.radius^2
"""
    volume(sphere::Sphere{3})
Evaluate volume of a three-dimensional sphere (4πr³/3)
"""
volume(sphere::Sphere{3})::Float64 = 4π/3 * sphere.radius^3

"""
    packing_fraction(radii::Vector{Real}, extent::NTuple{D,Float64}) where D
Evaluate the packing fraction of a collection of spheres with radii `radii`
in domain `extent`.
"""
function packing_fraction(
    radii::Vector{Real},
    extent::NTuple{D, Float64}
)::Float64 where D
    V₀ = prod(extent)
    V = 0.0
    for radius in radii
        V += volume(radius, D)
    end
    return V/V₀
end

"""
    volume(r::Real, D::Int)
Evaluate volume of a sphere of radius `r` in `D` dimensions.
Only D=2 and D=3 currently supported.
"""
function volume(r::Real, D::Int)::Float64
    if D == 2
        return π*r^2
    elseif D == 3
        return 4π/3*r^3
    end
end
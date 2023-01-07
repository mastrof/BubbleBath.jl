module BubbleBath

using LinearAlgebra: norm

export Sphere

"""
    Sphere(pos::NTuple{D,Real}, radius::Real) where D
Create a `Sphere{D}` object centered at `pos` with radius `radius`.
"""
struct Sphere{D}
    pos::NTuple{D,Float64}
    radius::Float64
    function Sphere(pos::NTuple{D,Real}, radius::Real) where D
        if radius ≤ 0
            throw(ArgumentError("Sphere radius must be a non-negative real number."))
        end
        new{D}(Float64.(pos), Float64(radius))
    end
end
Sphere(pos::AbstractVector{<:Real}, radius::Real) = Sphere(Tuple(pos), radius)

include("bubblebath_algorithm.jl")
include("walkmaps.jl")
include("packing_fraction.jl")

end

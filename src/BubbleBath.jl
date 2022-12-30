module BubbleBath

using LinearAlgebra: norm

export Sphere

struct Sphere{D}
    pos::NTuple{D,Float64}
    radius::Float64
end
"""
    Sphere(pos::NTuple{D,Real}, radius::Real) where D
Create a `Sphere{D}` object centered at `pos` with radius `radius`.
"""
function Sphere(pos::NTuple{D,Real}, radius::Real) where D
    if radius ≤ 0
        throw(ArgumentError("Sphere radius must be a non-negative real number."))
    end
    Sphere{D}(Float64.(pos),Float64(radius))
end

include("bubblebath_algorithm.jl")
include("packing_fraction.jl")

end

module BubbleBath

using LinearAlgebra: norm

export Sphere

struct Sphere{D}
    pos::NTuple{D,Float64}
    radius::Float64
end
"""
    Sphere(pos::NTuple{D,Float64}, radius::Real) where D
Create a `Sphere{D}` object centered at `pos` with radius `radius`.
"""
Sphere(pos::NTuple{D,Float64}, radius::Real) where D = Sphere{D}(pos,Float64(radius))

include("bubblebath.jl")
include("packing_fraction.jl")

end

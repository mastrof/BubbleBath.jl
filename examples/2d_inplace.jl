using BubbleBath
using Distributions: Exponential
using Plots

function circle(x,y,r,n=100)
    θ = LinRange(0, 2π, n)
    x .+ r.*sin.(θ), y .+ r.*cos.(θ)
end # function

Lx = 400
Ly = 400
extent = (Lx,Ly)
R = 50
D = 60
spheres = [
    Sphere((Lx/2-D,Ly/2-D), R),
    Sphere((Lx/2+D,Ly/2-D), R),
    Sphere((Lx/2,Ly/2+3D/4), R)
]

radius_pdf = Exponential(2.0)
ϕ_max = 0.25 - packing_fraction(spheres, extent)
min_distance = 2.0
bubblebath!(spheres, radius_pdf, ϕ_max, extent; min_distance)

plot(
    xlims=(0,extent[1]), ylims=(0,extent[2]),
    ratio=1, legend=false,
    grid=false, axis=false,
    bgcolor=:transparent,
    size=(Lx,Ly)
)
for i in eachindex(spheres)
    s = spheres[i]
    plot!(circle(s.pos..., s.radius),
        seriestype = :shape, lw = 0,
        fillcolor = i ≤ 3 ? i : :gray,
        alpha = i ≤ 3 ? 1.0 : 0.7,
    )
end
plot!()
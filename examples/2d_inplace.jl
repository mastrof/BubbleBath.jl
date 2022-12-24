using BubbleBath
using Distributions: Exponential
using Plots

function circle(x,y,r,n=500)
    θ = LinRange(0, 2π, n)
    x .+ r.*sin.(θ), y .+ r.*cos.(θ)
end # function

radius_pdf = Exponential(0.5)
L = 100
extent = (L,L)
R = 14
spheres = [
    Sphere((L/3,L/3), R),
    Sphere((2L/3,L/3), R),
    Sphere((L/2,2L/3), R)
]
ϕ_max = 0.3 - packing_fraction(spheres, extent)
min_distance = 1.0
bubblebath!(spheres, radius_pdf, ϕ_max, extent; min_distance)


plot(
    xlims=(0,extent[1]), ylims=(0,extent[2]),
    ratio=1, legend=false,
    grid=false, axis=false,
    bgcolor=:transparent,
    size=(400,400)
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
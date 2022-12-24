using BubbleBath
using Distributions: Uniform
using Plots

function circle(x,y,r,n=500)
    θ = LinRange(0, 2π, n)
    x .+ r.*sin.(θ), y .+ r.*cos.(θ)
end # function

radius_pdf = Uniform(1,5)
extent = (100, 50)
ϕ_max = 0.4
spheres = bubblebath(radius_pdf, ϕ_max, extent)

ϕ = round(packing_fraction(spheres, extent), digits=3)

plot(
    xlims=(0,extent[1]), ylims=(0,extent[2]),
    ratio=1, legend=false, grid=false,
    title="target ϕ=$(ϕ_max); actual ϕ=$(ϕ)"
)
for s in spheres
    plot!(circle(s.pos..., s.radius), seriestype=:shape)
end
plot!()
using BubbleBath
using Distributions: Uniform
using Plots

function circle(x₀,y₀,z₀,r,n=20)
    θ = LinRange(-π, π, n)
    φ = LinRange(0, π, n)
    x = x₀ .+ r.*cos.(θ).*sin.(φ)'
    y = y₀ .+ r.*sin.(θ).*sin.(φ)'
    z = z₀ .+ r.*ones(n).*cos.(φ)'
    return x, y, z
end # function

radius_pdf = Uniform(10, 25)
extent = (100, 100, 100)
ϕ_max = 0.3
min_distance = 10.0
spheres = bubblebath(radius_pdf, ϕ_max, extent; min_distance)

ϕ = round(packing_fraction(spheres, extent), digits=3)

plot(
    xlims=(0,extent[1]), ylims=(0,extent[2]), zlims=(0,extent[3]),
    size=(400,400), aspect_ratio=:equal,
    legend=false, colorbar=false,
    ticks=0:50:100,
    camera=(20,20)
)
for s in spheres
    surface!(circle(s.pos..., s.radius), alpha=0.7)
end
plot!()
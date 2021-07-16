begin #imports

include("Orrery.jl")

using GLMakie
using GeometryBasics
using LinearAlgebra
using Unitful, Unitful.DefaultSymbols

end #imports

begin #theme

marker_theme = Theme(
	Scatter=(
		markersize=20,
		marker=:+,
		color=:white,
	)
)

update_theme!(theme_dark())
update_theme!(marker_theme)

end #theme

begin #bodies

earth = Orrery.Body(
	μ=3.986004418e+14 * u"m^3/s^2",
	orbit=Orrery.FixedOrbit(),
)

moon = Orrery.Body(
	μ=4.9048695e+12 * u"m^3/s^2",
	orbit=Orrery.ClosedOrbit(
		a=384_748u"km",
		ẽ=Vec2(0.0549006, 0.),
		around=earth,
	),
)

satellite = Orrery.ClosedOrbit(
	a = 50_000u"km",
	ẽ = Vec2(0.25, 0.),
	around=moon
)

end #bodies

begin #observables

t = Node(0.0s)

moon_state = lift(t) do t
	θ = Orrery.propagate(moon.orbit, t)
	Orrery.GeometricState(moon.orbit, θ)
end

moon_position = @lift(Orrery.position($moon_state))

satellite_state = lift(t) do t
	θ = Orrery.propagate(satellite, t)
	Orrery.GeometricState(satellite, θ)
end

satellite_position = @lift(Orrery.position($satellite_state) + $moon_position)

end #observables

scene = Scene(scale_plot=true, show_axis=false)
earth_marker = scatter!(scene, Vec2(0), color=:grey)
moon_marker = scatter!(scene, @lift(ustrip.($moon_position/6e5)), color=:grey)
satellite_marker = scatter!(scene, @lift(ustrip.($satellite_position/6e5)), color=:grey)

T = 2 * 27.322 * 24 * 60 * 60s
Δt = 2 * 60 * 60s

t[] = 0.0s

while t[] < T
	t[] += Δt
	sleep(1//30)
end

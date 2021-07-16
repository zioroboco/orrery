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

text_theme = Theme(
	Text=(
		textsize=24,
		color=:white,
		align=(:left, :center),
	)
)

update_theme!(theme_dark())
update_theme!(marker_theme)
update_theme!(text_theme)

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

moon_velocity = lift(moon_position) do r̃
	a = ustrip(moon.orbit.a)
	μ = ustrip(moon.orbit.around.μ)
	r = ustrip(Orrery.magnitude(r̃))
  sqrt(μ * (2/r - 1/a)) * u"km/s"
end

satellite_state = lift(t) do t
	θ = Orrery.propagate(satellite, t)
	Orrery.GeometricState(satellite, θ)
end

satellite_position = @lift(Orrery.position($satellite_state) + $moon_position)

end #observables

scene = Scene(scale_plot=true, show_axis=false, resolution=(800, 800))

earth_marker = scatter!(scene, Vec2(0), color=:grey)
moon_marker = scatter!(scene, @lift(ustrip.($moon_position/6e5)), color=:white)
satellite_marker = scatter!(scene, @lift(ustrip.($satellite_position/6e5)), color=:grey)

moon_velocity_stats = text!(scene, @lift(string(round(u"km/s", $moon_velocity, digits=3))), position=Vec2(-0.9, 0.9))
elapsed_time_stats = text!(scene, @lift(string("$(ustrip(round(u"d", $t, digits=1))) days")), position=Vec2(-0.9, 0.8))

T = 2 * 27.322 * 24 * 60 * 60s
Δt = 2 * 60 * 60s

t[] = 0.0s

while t[] < T
	t[] += Δt
	sleep(1//30)
end

begin # Imports
	include("Orrery.jl")
	import .Orrery
	using GLMakie
	using GeometryBasics
	using LinearAlgebra
	using Unitful
end

begin # Theme
	theme = Theme(
		Scatter=(
			markersize=20,
			marker=:+,
			color=:grey,
		),
		Text=(
			textsize=24,
			color=:grey,
		),
	)
	update_theme!(theme_dark())
	update_theme!(theme)
end

begin # Entities

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
		a = 30_000u"km",
		ẽ = Vec2(0.20, 0.),
		around=moon
	)

	return
end

begin # Observables

	t = Node(0.0u"s")

	moon_state = lift(t) do t
		θ = Orrery.propagate(moon.orbit, t)
		Orrery.StateElements(moon.orbit, θ)
	end

	satellite_state = lift(t) do t
		θ = Orrery.propagate(satellite, t)
		Orrery.StateElements(satellite, θ)
	end

	particle_position = lift(t) do t
		particle_state[].r̃
	end

	moon_position = lift(moon_state) do state
		Orrery.position(state)
	end

	satellite_position = lift(satellite_state, moon_state) do s, m
		Orrery.position(s) + Orrery.position(m)
	end

	moon_velocity = lift(moon_position) do r̃
		a = ustrip(moon.orbit.a)
		μ = ustrip(moon.orbit.around.μ)
		r = ustrip(Orrery.magnitude(r̃))
		sqrt(μ * (2/r - 1/a)) * u"km/s"
	end

	return
end

ZOOM = 5e5
ZOOM_RECT = FRect(-1*ZOOM, -1*ZOOM, 2*ZOOM, 2*ZOOM)

scene = Scene(
	scale_plot=true,
	show_axis=false,
	resolution=(800, 800),
	limits=ZOOM_RECT,
)

## ---

satellite_marker = scatter!(scene, @lift(ustrip.($satellite_position)))
moon_marker = scatter!(scene, @lift(ustrip.($moon_position)))
earth_marker = scatter!(scene, Vec2(0))

update_cam!(scene, ZOOM_RECT)

begin # Run!

	T = 55 * 24 * 60 * 60u"s"
	Δt = 4 * 60 * 60u"s"

	t[] = 0.0u"s"

	while t[] < T
		t[] += Δt
		sleep(1//30)
	end

end

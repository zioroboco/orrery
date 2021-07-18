begin # Imports

	include("Orrery.jl")
	import .Orrery

	using GLMakie
	using GeometryBasics
	using LinearAlgebra
	using Unitful

end

begin # Scene

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

	t = Node(0.0u"s")

	scene = Scene(scale_plot=true, show_axis=false, resolution=(800, 800))

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
		a = 50_000u"km",
		ẽ = Vec2(0.25, 0.),
		around=moon
	)

	particle_initial_state = Orrery.StateVectors(
		Vec2(100_000u"km", 0.0u"km"),
		Vec2(0.0u"km/s", 2.0u"km/s"),
	)

	return
end

begin # Observables

	particle_state = lift(t) do t
		if @isdefined(particle_state)
			return Orrery.update(particle_state[], Δt, around=earth)
		else
			return particle_initial_state
		end
	end

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

begin # Markers

	particle_marker = scatter!(scene, @lift(ustrip.($particle_position/6e5)), color=:orangered)

	earth_marker = scatter!(scene, Vec2(0), color=:grey)
	moon_marker = scatter!(scene, @lift(ustrip.($moon_position/6e5)), color=:white)
	satellite_marker = scatter!(scene, @lift(ustrip.($satellite_position/6e5)), color=:grey)

	moon_velocity_stats = text!(scene, @lift(string(round(u"km/s", $moon_velocity, digits=3))), position=Vec2(-0.9, 0.9))
	elapsed_time_stats = text!(scene, @lift(string("$(ustrip(round(u"d", $t, digits=1))) days")), position=Vec2(-0.9, 0.8))

	return
end

begin # Run!

	T = 2 * 27.322 * 24 * 60 * 60u"s"
	Δt = 2 * 60 * 60u"s"

	t[] = 0.0u"s"

	while t[] < T
		t[] += Δt
		sleep(1//30)
	end

end

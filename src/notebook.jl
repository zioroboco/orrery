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

	Rot(θ) = [
		cos(θ) -sin(θ);
		sin(θ)  cos(θ)
	]

	function draw_orbit!(scene::Scene, orbit::Orrery.Orbit)
		θs = range(0*u"rad", stop=2π*u"rad", length=100)
		r̃s = map(θᵢ -> Orrery.position(Orrery.StateElements(orbit, θᵢ)), θs)
		θₑ = angle(orbit.ẽ[1] + orbit.ẽ[2]*im)
		r̃ₑs = Vec2.(Rot.(θₑ + π) .* r̃s)
		lines!(scene, ustrip.(r̃ₑs), color=:orangered)
	end

	moon_position_rotated = lift(moon_state, moon_position) do state, r̃
		θₑ = angle(state.orbit.ẽ[1] + state.orbit.ẽ[2]*im)
		Vec2(Rot(θₑ + π) * r̃)
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

satellite_marker = scatter!(scene, @lift([ustrip.($satellite_position)]), color=:grey)
moon_marker = scatter!(scene, @lift([ustrip.($moon_position_rotated)]))
earth_marker = scatter!(scene, Vec2(0))

moon_plot = Node([ustrip.(moon_position[])])
moon_line = lines!(scene, @lift(ustrip.($plot)))

draw_orbit!(scene, moon.orbit)

update_cam!(scene, ZOOM_RECT)

begin # Run!

	T = 55 * 24 * 60 * 60u"s"
	Δt = 4 * 60 * 60u"s"

	t[] = 0.0u"s"

	while t[] < T
		t[] += Δt
		moon_plot[] = push!(moon_plot[], ustrip.(moon_position[]))
		sleep(1//30)
	end

end

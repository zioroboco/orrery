module Orrery

begin # Imports
	using Base: @kwdef
	using GeometryBasics: Vec2
	using LinearAlgebra: norm, normalize
	using Memoize: @memoize
	using Unitful
end

begin # Units

	const Dimensionless = Float64
	const Anomaly = typeof(1.0u"rad")
	const Distance = typeof(1.0u"km")
	const Velocity = typeof(1.0u"km/s")
	const Gravitation = typeof(1.0u"km^3/s^2")

end

begin # Orbit types

	"""
	Supertype of all concrete orbit types.
	"""
	abstract type Orbit end

	"""
	Gravitationally significant massive body.
	"""
	@kwdef struct Body
		# TODO: Model SOIs!
		orbit::Orbit # around parent body (or fixed)
		μ::Gravitation
	end

	"""
	Sentinel type indicating the root of the tree of bodies.

	Alternatively: stationary in this reference frame (with obvious caveats).
	"""
	struct FixedOrbit <: Orbit end

	"""
	Keplerian hyperbolic trajectory.
	"""
	@kwdef struct OpenOrbit <: Orbit
		a::Distance
		ẽ::Vec2{Dimensionless}
		around::Body
	end

	"""
	Keplerian elliptical trajectory.
	"""
	@kwdef struct ClosedOrbit <: Orbit
		a::Distance
		ẽ::Vec2{Dimensionless}
		around::Body
	end

	function Orbit(a::Distance, ẽ::Vec2{Dimensionless}; around::Body)::Orbit
		if magnitude(ẽ) > 1.0
			OpenOrbit(a * -sign(a), ẽ, around)
		else
			ClosedOrbit(a, ẽ, around)
		end
	end
end

begin # State types

	"""
	Supertype of all concrete state types.
	"""
	abstract type State end

	"""
	State of a particle in a Newtonian regime.
	"""
	@kwdef struct StateVectors
		r̃::Vec2{Distance}
		ṽ::Vec2{Velocity}
	end

	"""
	State of a particle in a Keplerian regime.
	"""
	@kwdef struct StateElements
		orbit::Orbit
		θ::Anomaly
	end

end

begin # Anomalies

	"""
	Mean anomaly (elliptic).
	"""
	function M(orbit::ClosedOrbit, t)
		-sqrt(orbit.around.μ / orbit.a^3) * t
	end

	"""
	Mean anomaly (hyperbolic).
	"""
	function M(orbit::OpenOrbit, t)
		-sqrt(orbit.around.μ / (-orbit.a)^3) * t
	end

	"""
	Eccentric anomaly.
	"""
	function E(orbit::ClosedOrbit, M; Eₖ=M, ϵ=eps(Float32))
		e = magnitude(orbit.ẽ)
		Eₖ₊₁ = Eₖ - (M - Eₖ + e*sin(Eₖ)) / (e*cos(Eₖ) - 1)
		if abs(M - Eₖ₊₁ + e*sin(Eₖ₊₁)) <= ϵ
			Eₖ₊₁
		else
			E(orbit, M, Eₖ=Eₖ₊₁, ϵ=ϵ)
		end
	end

	"""
	Hyperbolic anomaly.

	Uses the same name as eccentric anomaly to enable dispatching on the orbit.
	"""
	function E(orbit::OpenOrbit, M; Hₖ=M, ϵ=eps(Float32))
		e = magnitude(orbit.ẽ)
		Hₖ₊₁ = Hₖ + (M - e*sinh(Hₖ) + Hₖ) / (e*cosh(Hₖ) - 1)
		if abs(M - e*sinh(Hₖ₊₁) + Hₖ₊₁) <= ϵ
			Hₖ₊₁
		else
			E(orbit, M, Hₖ=Hₖ₊₁, ϵ=ϵ)
		end
	end

	"""
	True anomaly (elliptic).
	"""
	function f(orbit::ClosedOrbit, E)
		e = magnitude(orbit.ẽ)
		2 * atan(sqrt((1 + e) / (1 - e)) * tan(E / 2))
	end

	"""
	True anomaly (hyperbolic).
	"""
	function f(orbit::OpenOrbit, H)
		e = magnitude(orbit.ẽ)
		2 * atan(sqrt((e + 1) / (e - 1)) * tanh(H / 2))
	end

end

begin # Propagation

	function propagate(orbit::Orbit, t; ϵ=eps(Float32))
		f(orbit, E(orbit, M(orbit, t), ϵ=ϵ)) * u"rad"
	end

	function position(state::StateVectors)
		state.r̃
	end

	function position(state::StateElements)::Vec2
		e = Orrery.magnitude(state.orbit.ẽ)
		r = -state.orbit.a * (1 - e^2) / (1 + e*cos(state.θ))
		Vec2(r*cos(state.θ), r*sin(state.θ))
	end

	function update(state::Orrery.StateVectors, Δt)::Orrery.StateVectors
		ṽ = state.ṽ
		r̃ = state.r̃
		r̂ = Orrery.direction(r̃)
		r = Orrery.magnitude(r̃)
		# TODO Fix hard coded parent body
		ṽ′ = ṽ - r̂ * earth.μ / r^2 * Δt
		r̃′ = r̃ + ṽ′ * Δt
		Orrery.StateVectors(r̃′, ṽ′)
	end

end

begin # Utilities

	# Memoise using objectid, which checks hashed identity of objects (not values).
	@memoize magnitude(ã) = norm(ã)
	@memoize direction(ã) = normalize(ã)

end

export
	Anomaly,
	Body,
	ClosedOrbit,
	Dimensionless,
	Distance,
	Elements,
	FixedOrbit,
	Gravitation,
	OpenOrbit,
	Orbit,
	Orbit,
	State,
	StateElements,
	StateVectors,
	Velocity,
	magnitude,
	position,
	propagate,
	update

end

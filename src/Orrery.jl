module Orrery

begin #imports

using Base: @kwdef
using GeometryBasics: Vec2
using LinearAlgebra: norm, normalize
using Memoize: @memoize
using Unitful

end #imports

begin #units

const Dimensionless = Float64
const Anomaly = typeof(1.0u"rad")
const Distance = typeof(1.0u"km")
const Velocity = typeof(1.0u"km/s")
const Gravitation = typeof(1.0u"km^3/s^2")

end #units

begin #orbits

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

end #orbits

begin #state

"""
Supertype of all concrete state types.
"""
abstract type State end

"""
State of a particle in a Newtonian regime.
"""
@kwdef struct CartesianState
	r̃::Vec2{Distance}
	ṽ::Vec2{Velocity}
end

"""
State of a particle in a Keplerian regime.
"""
@kwdef struct GeometricState
	orbit::Orbit
	θ::Anomaly
end

end #state

begin #anomalies

"""
Mean anomaly (elliptic).
"""
function M(orbit::ClosedOrbit, Δt)
	-sqrt(orbit.around.μ / orbit.a^3) * Δt
end

"""
Mean anomaly (hyperbolic).
"""
function M(orbit::OpenOrbit, Δt)
	-sqrt(orbit.around.μ / (-orbit.a)^3) * Δt
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

end #anomalies

begin #propagation

function propagate(orbit::Orbit, Δt; ϵ=eps(Float32))
	f(orbit, E(orbit, M(orbit, Δt), ϵ=ϵ)) * u"rad"
end

function position(state::CartesianState)
	state.r̃
end

function position(state::GeometricState)::Vec2
	e = Orrery.magnitude(state.orbit.ẽ)
	r = -state.orbit.a * (1 - e^2) / (1 + e*cos(state.θ)) * u"km"
	Vec2(r*cos(state.θ), r*sin(state.θ))
end

end #propagation

begin #utilities

# Memoise using objectid, which checks hashed identity of objects (not values).
@memoize magnitude(ã) = norm(ã)
@memoize direction(ã) = normalize(ã)

end #utilities

begin #exports

export
	Anomaly,
	Body,
	CartesianState,
	ClosedOrbit,
	Dimensionless,
	Distance,
	Elements,
	FixedOrbit,
	GeometricState,
	Gravitation,
	OpenOrbit,
	Orbit,
	State,
	Velocity,
	magnitude,
	position,
	propagate

end #exports

end #module

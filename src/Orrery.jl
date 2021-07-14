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
function M(orbit::ClosedOrbit, Δt, μ)
	-sqrt(μ / orbit.a^3) * Δt
end

"""
Mean anomaly (hyperbolic).
"""
function M(orbit::OpenOrbit, Δt, μ)
	-sqrt(μ / (-orbit.a)^3) * Δt
end

"""
Eccentric anomaly.
"""
function EH(orbit::ClosedOrbit, M; Eₖ=M, ϵ=eps(Float32))
	e = magnitude(orbit.ẽ)
	Eₖ₊₁ = Eₖ - (M - Eₖ + e*sin(Eₖ)) / (e*cos(Eₖ) - 1)
	if abs(M - Eₖ₊₁ + e*sin(Eₖ₊₁)) <= ϵ
		Eₖ₊₁
	else
		EH(orbit, Mₕ, Hₖ=Hₖ₊₁, ϵ=ϵ)
	end
end

"""
Hyperbolic anomaly.
"""
function EH(orbit::OpenOrbit, Mₕ; Hₖ=Mₕ, ϵ=eps(Float32))
	e = magnitude(orbit.ẽ)
	Hₖ₊₁ = Hₖ + (Mₕ - e*sinh(Hₖ) + Hₖ) / (e*cosh(Hₖ) - 1)
	if abs(Mₕ - e*sinh(Hₖ₊₁) + Hₖ₊₁) <= ϵ
		Hₖ₊₁
	else
		EH(orbit, Mₕ, Hₖ=Hₖ₊₁, ϵ=ϵ)
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

begin #utilities

@memoize Dict magnitude(ã) = norm(ã)
@memoize Dict direction(ã) = normalize(ã)

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
	unit

end #exports

end #module

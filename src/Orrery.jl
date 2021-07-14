module Orrery

begin #imports

using Base: @kwdef
using GeometryBasics: Vec2
using LinearAlgebra: norm, normalize
using Memoize: @memoize
using Unitful: Quantity, FixedUnits, NoDims, 𝐋, 𝐓, km, s, rad

end #imports

begin #units

const Dimensionless = Float64
const Anomaly = Quantity{Float64, NoDims, FixedUnits{(rad,), NoDims, nothing}}
const Distance = Quantity{Float64, 𝐋, FixedUnits{(km,), 𝐋, nothing}}
const Velocity = Quantity{Float64, 𝐋/𝐓, FixedUnits{(km/s), 𝐋/𝐓, nothing}}
const Gravitation = Quantity{Float64, 𝐋^3/𝐓^2, FixedUnits{(km^3, s^-2), 𝐋^3/𝐓^2, nothing}}

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

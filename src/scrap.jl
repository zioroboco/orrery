using Base: @kwdef
using Chain: @chain
using Dictionaries
using GLMakie
using GeometryBasics
using LinearAlgebra

const Position = Point2{Float64}
const Velocity = Vec2{Float64}

abstract type AbstractSituation end

@kwdef struct StateVectors <: AbstractSituation
	position::Position
	velocity::Velocity
end

@kwdef struct Stationary <: AbstractSituation
	position::Position
end

@enum Body begin
	Earth = 1
	Moon
end

ledger = dictionary([
	:situation => dictionary([
		Earth => Stationary(position=[0.0, 0.0]),
		Moon => StateVectors(position=[0.5, 0.0], velocity=[0.0, 1.0]),
	]),
	:gravitation => dictionary([
		Earth => 1.0,
		Moon => 0.01,
	]),
	:parent => dictionary([
		Moon => Earth,
	]),
])

function query(symbols)::Vector{Body}
	@chain begin
		collect(symbols)
		map(
			body -> map(symbol -> all(body in keys(ledger[symbol])), _),
			instances(Body),
		)
		reduce(hcat, _)
		reduce(.&, eachrow(_))
		instances(Body)[_]
		collect(_)
	end
end

scene = Scene(
	scale_plot=true,
	show_axis=false,
	limits=FRect(-1, -1, 2, 2)
)

function get_position(situation::AbstractSituation)::Position
	return situation.position
end

positions = Node(zeros(Position, length(instances(Body))))
blips = scatter!(scene, positions)

positions[] = map(get_position, collect(values(ledger[:situation])))

const Δt = 0.01

function update(body::Body, parent::Body)
	r̃ = ledger[:situation][body].position
	ṽ = ledger[:situation][body].velocity

	r̃ₚ = r̃ - ledger[:situation][parent].position
	μₚ = ledger[:gravitation][parent]

	rₚ = norm(r̃ₚ)
	r̂ₚ = normalize(r̃ₚ)

	ṽ′ = ṽ - r̂ₚ * μₚ / rₚ^2 * Δt
	r̃′ = r̃ + ṽ′ * Δt

	return StateVectors(position=r̃′, velocity=ṽ′)
end

for t in 0:Δt:10.0
	for body in query([:situation, :parent])
		parent = ledger[:parent][body]
		ledger[:situation][body] = update(body, parent)
	end
	positions[] = map(s -> s.position, collect(values(ledger[:situation])))
	sleep(1/30)
end

using Base: @kwdef
using Chain: @chain
using Dictionaries
using GLMakie
using GeometryBasics

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
		Moon => StateVectors(position=[0.5, 0.0], velocity=[0.0, 0.5]),
	]),
	:parent => dictionary([
		Moon => Earth,
	]),
])

function query(symbols)::Set{Body}
	@chain begin
		collect(symbols)
		map(
			body -> map(symbol -> all(body in keys(ledger[symbol])), _),
			instances(Body),
		)
		reduce(hcat, _)
		reduce(.&, eachrow(_))
		instances(Body)[_]
		Set(_)
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

positions = Node(map(get_position, collect(values(ledger[:situation]))))
blips = scatter!(scene, positions)

const Δt = 0.1

function update(stationary::Stationary)
	return stationary
end

function update(vectors::StateVectors)
	return StateVectors(
		vectors.position + vectors.velocity * Δt,
		vectors.velocity,
	)
end

for t in 0:Δt:100.0
	for (body::Body, situation::AbstractSituation) in pairs(ledger[:situation])
		ledger[:situation][body] = update(situation)
	end
	positions[] = map(s -> s.position, collect(values(ledger[:situation])))
	sleep(1/30)
end

using Base: @kwdef
using Dictionaries
using GLMakie
using GeometryBasics

@enum Body begin
	Earth = 1
	Moon
end

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

scene = Scene(
	scale_plot=true,
	show_axis=false,
	limits=FRect(-1, -1, 2, 2)
)

positions = Node{Array{Position}}(zeros(Vec2, length(instances(Body))))
blips = scatter!(scene, positions)

world = Dictionary{Body, AbstractSituation}()
insert!(world, Earth, Stationary(position=[0.0, 0.0]))
insert!(world, Moon, StateVectors(position=[0.5, 0.0], velocity=[0.0, 0.5]))

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
	for (body::Body, situation::AbstractSituation) in pairs(world)
		world[body] = update(situation)
	end
	positions[] = map(s -> s.position, collect(values(world)))
	sleep(1/30)
end

using Base: @kwdef
using Chain: @chain
using Dictionaries
using GLMakie
using GeometryBasics
using LinearAlgebra

@enum Body begin
	Earth = 1
	Moon
end

@kwdef struct SOI
	parent::Body
	μ::Float64
end

ledger = dictionary([
	:position => dictionary([
		Earth => Point2([0.0, 0.0]),
		Moon => Point2([0.5, 0.0])
	]),
	:velocity => dictionary([
		Moon => Point2([0.0, 1.0]),
	]),
	:soi => dictionary([
		Moon => SOI(parent=Earth, μ=1.0),
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

positions = Node(zeros(Point2{Float64}, length(instances(Body))))
blips = scatter!(scene, positions)

positions[] = collect(values(ledger[:position]))

const Δt = 0.01

function update(body::Body)
	r̃ = ledger[:position][body]
	ṽ = ledger[:velocity][body]

	r̃ₚ = r̃ - ledger[:position][ledger[:soi][body].parent]
	μₚ = ledger[:soi][body].μ

	rₚ = norm(r̃ₚ)
	r̂ₚ = normalize(r̃ₚ)

	ṽ′ = ṽ - r̂ₚ * μₚ / rₚ^2 * Δt
	r̃′ = r̃ + ṽ′ * Δt

	ledger[:position][body] = r̃′
	ledger[:velocity][body] = ṽ′
end

for t in 0.0:Δt:10.0
	foreach(body -> update(body), query([:position, :velocity, :soi]))
	positions[] = collect(values(ledger[:position]))
	sleep(1/30)
end

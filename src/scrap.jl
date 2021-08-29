using Base: @kwdef
using Chain: @chain
using Dictionaries
using GLMakie
using GeometryBasics
using LinearAlgebra

magnitude = norm
direction = normalize

@enum Body begin
	Earth = 1
	Moon
	Satellite
	Comet
end

@kwdef struct SOI
	parent::Body
	μ::Float64
end

@kwdef struct Elements
	a::Float64
	ẽ::Vec2{Float64}
	f::Float64
end

ledger = dictionary([
	:position => dictionary([
		Earth => Point2([0.0, 0.0]),
		Moon => Point2([0.0, 0.75]),
		Satellite => Point2([-0.5, 0.0]),
		Comet => Point2([0.25, 0.25]),
	]),
	:velocity => dictionary([
		Satellite => Point2([0.0, -1.0]),
	]),
	:elements => dictionary([
		Moon => Elements(a=1.0, ẽ=Vec2(0.1, 0.0), f=0.0),
		Comet => Elements(a=0.05, ẽ=Vec2(2.5, 2.5), f=0.0),
	]),
	:soi => dictionary([
		Moon => SOI(parent=Earth, μ=1.0),
		Satellite => SOI(parent=Earth, μ=1.0),
		Comet => SOI(parent=Earth, μ=1.0),
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

function newtonian_update!(body::Body)
	r̃ = ledger[:position][body]
	ṽ = ledger[:velocity][body]

	r̃ₚ = r̃ - ledger[:position][ledger[:soi][body].parent]
	μₚ = ledger[:soi][body].μ

	rₚ = magnitude(r̃ₚ)
	r̂ₚ = direction(r̃ₚ)

	ṽ′ = ṽ - r̂ₚ * μₚ / rₚ^2 * Δt
	r̃′ = r̃ + ṽ′ * Δt

	ledger[:position][body] = r̃′
	ledger[:velocity][body] = ṽ′
end

function keplerian_update!(body::Body, t)
	function Mₑ(elements::Elements, soi::SOI, t)
		-sqrt(soi.μ / elements.a^3) * t
	end

	function Mₕ(elements::Elements, soi::SOI, t)
		-sqrt(soi.μ / (elements.a)^3) * t
	end

	function E(elements::Elements, M; Eₖ=M, ϵ=eps(Float32))
		e = magnitude(elements.ẽ)
		Eₖ₊₁ = Eₖ - (M - Eₖ + e*sin(Eₖ)) / (e*cos(Eₖ) - 1)
		if abs(M - Eₖ₊₁ + e*sin(Eₖ₊₁)) <= ϵ
			Eₖ₊₁
		else
			E(elements, M, Eₖ=Eₖ₊₁, ϵ=ϵ)
		end
	end

	function H(elements::Elements, Mₕ; Hₖ=Mₕ, ϵ=eps(Float32))
		e = magnitude(elements.ẽ)
		Hₖ₊₁ = Hₖ + (Mₕ - e*sinh(Hₖ) + Hₖ) / (e*cosh(Hₖ) - 1)
		if abs(Mₕ - e*sinh(Hₖ₊₁) + Hₖ₊₁) <= ϵ
			Hₖ₊₁
		else
			H(elements, Mₕ, Hₖ=Hₖ₊₁, ϵ=ϵ)
		end
	end

	function fₑ(elements::Elements, E)
		e = magnitude(elements.ẽ)
		2 * atan(sqrt((1 + e) / (1 - e)) * tan(E / 2))
	end

	function fₕ(elements::Elements, H)
		e = magnitude(elements.ẽ)
		2 * atan(sqrt((e + 1) / (e - 1)) * tanh(H / 2))
	end

	function propagate(elements::Elements, soi::SOI, t; ϵ=eps(Float32))
		if magnitude(elements.ẽ) < 1.0
			fₑ(elements, E(elements, Mₑ(elements, soi, t), ϵ=ϵ))
		else
			fₕ(elements, H(elements, Mₕ(elements, soi, t), ϵ=ϵ))
		end
	end

	function to_position(elements::Elements)::Vec2{Float64}
		e = magnitude(elements.ẽ)
		r = -elements.a * (1 - e^2) / (1 + e*cos(elements.f))
		Vec2(r*cos(elements.f), r*sin(elements.f))
	end

	elements = ledger[:elements][body]
	soi = ledger[:soi][body]

	ledger[:elements][body] = Elements(
		elements.a,
		elements.ẽ,
		propagate(elements, soi, t)
	)

	ledger[:position][body] = to_position(ledger[:elements][body])
end

for t in -5.0:Δt:5.0
	foreach(body -> newtonian_update!(body), query([:position, :velocity, :soi]))
	foreach(body -> keplerian_update!(body, t), query([:position, :elements, :soi]))
	positions[] = collect(values(ledger[:position]))
	sleep(1/30)
end

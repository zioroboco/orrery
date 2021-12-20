using Base: @kwdef
using Chain: @chain
using Dictionaries
using GLMakie
using GeometryBasics
using LinearAlgebra

const magnitude = norm
const direction = normalize

set_theme!(theme_dark())

const SCALE = 1e-9

scene = Scene(
	scale_plot=true,
	show_axis=false,
	limits=FRect(-1, -1, 2, 2),
)

@enum Body begin
	Earth = 1
	Moon
  Satellite
end

positions = Node(zeros(Point2{Float64}, length(instances(Body))))
blips = scatter!(scene, positions, marker=:+, markersize=20, color=:ivory)

stats_content = Node("...")
stats = text!(scene, stats_content, position=Point2([0.0, 0.8]), color=:ivory, align=(:center, :center))

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
		Moon => Point2([0.0, 0.0]),
		Satellite => Point2([2.0e8, 0.0]),
	]),
	:velocity => dictionary([
    Satellite => Vec2([0.0, 1.0e3]),
	]),
	:elements => dictionary([
		Moon => Elements(a=3.84748e+8, ẽ=Vec2(0.0549006, 0.0), f=0),
	]),
	:soi => dictionary([
		Moon => SOI(parent=Earth, μ=3.986004418e+14),
		Satellite => SOI(parent=Earth, μ=3.986004418e+14),
	]),
])

function query(ledger, symbols)::Vector{Body}
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

query(f::Function, symbols) = map(f, query(ledger, symbols))

positions[] = collect(values(ledger[:position]))

function newtonian_update!(ledger, body::Body, Δt)
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

Rot(θ) = [cos(θ) -sin(θ); sin(θ)  cos(θ)]
radial(a, e, f) = a * (1 - e^2) / (1 + e*cos(f))

function to_position(elements::Elements)::Vec2{Float64}
	r = radial(elements.a, magnitude(elements.ẽ), elements.f)
	θₑ = angle(elements.ẽ[1] + elements.ẽ[2]*im)
	Rot(θₑ) * Vec2(r * cos(elements.f), r * sin(elements.f))
end

function keplerian_update!(ledger, body::Body, t)
	function Mₑ(elements::Elements, soi::SOI, t)
		-sqrt(soi.μ / elements.a^3) * t # FIXME t=0 => M=0
	end

	function Mₕ(elements::Elements, soi::SOI, t)
		# TODO Why not (-a)³ in Mₕ?
		-sqrt(soi.μ / elements.a^3) * t
	end

	function E(elements::Elements, Mₑ; Eₖ=Mₑ, ϵ=eps(Float32))
		e = magnitude(elements.ẽ)
		Eₖ₊₁ = Eₖ - (Mₑ - Eₖ + e*sin(Eₖ)) / (e*cos(Eₖ) - 1)
		if abs(Mₑ - Eₖ₊₁ + e*sin(Eₖ₊₁)) <= ϵ
			Eₖ₊₁
		else
			E(elements, Mₑ, Eₖ=Eₖ₊₁, ϵ=ϵ)
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

	elements = ledger[:elements][body]
	soi = ledger[:soi][body]

	ledger[:elements][body] = Elements(
		elements.a,
		elements.ẽ,
		propagate(elements, soi, t)
	)

	ledger[:position][body] = to_position(ledger[:elements][body])
end

function draw_orbit!(scene::Scene, elements::Elements)
	θs = range(0, stop=2π, length=100)
	r̃s = map(θᵢ -> to_position(Elements(elements.a, elements.ẽ, θᵢ)), θs)
	lines!(scene, r̃s .* SCALE, color=:grey)
end

function to_elements(r̃, ṽ, μ)::Elements
	r = magnitude(r̃)
	v = magnitude(ṽ)
	ẽ = (1/μ) * ((v^2 - (μ/r)) * r̃ - dot(r̃, ṽ) * ṽ)
	Energy = v^2/2 - μ/r
	a = -μ / (2 * Energy)
	f = acos(dot(r̃, ẽ) / (r * magnitude(ẽ)))
	return Elements(a, ẽ, f)
end

function v²(elements::Elements, soi::SOI)
	r = radial(elements.a, magnitude(elements.ẽ), elements.f)
	return soi.μ * (2/r - 1/elements.a)
end

function v_display(ledger, body::Body)::String
	v = (
		body in keys(ledger[:elements])
			? sqrt(v²(ledger[:elements][body], ledger[:soi][body]))
			:	magnitude(ledger[:velocity][body])
	)
	return "$(round(v)) m/s"
end

function main(ledger)
	Δt = 60*60*2
	for t in 0.0:Δt:60*60*24*27.322
		stats_content[] = join([
			"t = $(round(t/60/60/24, digits=1)) days",
			"\n",
			"v_moon = $(v_display(ledger, Moon))",
			"v_satellite = $(v_display(ledger, Satellite))",
			"\n",
			"Regime: $(Satellite in keys(ledger[:elements]) ? "Keplerian" : "Newtonian")"
		], "\n")

		query([:position, :velocity]) do body
			newtonian_update!(ledger, body, Δt)
		end

		query([:position, :elements]) do body
			keplerian_update!(ledger, body, t)
		end

		positions[] = collect(values(ledger[:position])) * SCALE
		sleep(1/30)
	end
end

## ---

# Run the simulation
main(deepcopy(ledger))

# Draw the orbit
draw_orbit!(scene, ledger[:elements][Moon])

# Clear last drawn orbit
pop!(scene.plots); scene

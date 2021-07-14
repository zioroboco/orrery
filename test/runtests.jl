using GeometryBasics
using Test
using Unitful

include("../src/Orrery.jl")
@testset "Orrery.jl" begin

local earth::Orrery.Body = Orrery.Body(
	μ=3.986004418e11u"km^3/s^2",
	orbit=Orrery.FixedOrbit(),
)

local moon:: Orrery.Body = Orrery.Body(
	μ=4.9048695e9u"km^3/s^2",
	orbit=Orrery.ClosedOrbit(
		a=1.0u"km",
		ẽ=Vec2(0.0549, 0.),
		around=earth,
	),
)

@test typeof(earth.orbit) == Orrery.FixedOrbit
@test Orrery.magnitude(moon.orbit.ẽ) ≈ 0.0549 atol=0.0001
@test Orrery.direction(moon.orbit.ẽ) == [1.0, 0.0]

end;

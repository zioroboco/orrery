begin #imports

include("../src/Orrery.jl")
import .Orrery

using GeometryBasics
using Test
using Unitful, Unitful.DefaultSymbols

end #imports
begin #fixtures

earth = Orrery.Body(
	μ=3.986004418e+14 * u"m^3/s^2",
	orbit=Orrery.FixedOrbit(),
)

moon = Orrery.Body(
	μ=4.9048695e+12 * u"m^3/s^2",
	orbit=Orrery.ClosedOrbit(
		a=384_748u"km",
		ẽ=Vec2(0.0549006, 0.),
		around=earth,
	),
)

end #imports

@testset "Orrery.jl" begin

@test typeof(earth.orbit) == Orrery.FixedOrbit
@test Orrery.magnitude(moon.orbit.ẽ) ≈ 0.0549 atol=0.0001
@test Orrery.direction(moon.orbit.ẽ) == [1.0, 0.0]

@test typeof(Orrery.propagate(moon.orbit, 60s)) == Orrery.Anomaly

# TODO Investigate slight propagation error
@test Orrery.propagate(moon.orbit, 27.322u"d", ϵ=1.0e-15) ≈ 0rad atol=0.05

end #testset

module Orrery

using InlineTest

backwards(s) = reverse(s)

@testset "backwards" begin
	@test backwards("straw") == "warts"
end

end # module

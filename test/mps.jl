using IsoPEPS
using Test
@testset "mps" begin
    mps=generate_mps(3,10)
    result= IsoPEPS.LinearAlgebra.dot(mps,mps)
    @test  isless(0.0,result)
end

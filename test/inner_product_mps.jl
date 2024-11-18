using IsoPEPS
using Test

@testset "inner_product_mps" begin
    tensor1 = [1 1; 1 1]
    tensor2 = [1 1; 1 1]
    mps1 = IsoPEPS.MPS([tensor1, tensor2])
    mps2 = IsoPEPS.MPS([tensor1, tensor2])
    result = IsoPEPS.inner_product(mps1, mps2)
    expected_value = 16.0
    @test result ≈ expected_value
end
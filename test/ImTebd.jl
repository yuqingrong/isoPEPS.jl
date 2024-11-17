using Test, isoPEPS

@testset "ImTebd" begin
    nsites = 10
    J=1.0
    h=0.2
    t=20
    H = isoPEPS.transverse_ising(nsites;J, h) |> isoPEPS.cache
    p = 1.0

    reg = isoPEPS.rand_state(isoPEPS.nqubits(H))
    reg |> itime_groundstate!(H; t,p)
    EG = isoPEPS.real(isoPEPS.expect(H, reg)/nsites)
    @test isoPEPS.ishermitian(H)
    # using imaginary time Evolution
    @test isapprox(EG, -0.9120354170186685, atol=1e-1)
end

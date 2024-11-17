using Test, isoPEPS,Arpack

@testset "LanczosAlgorithm" begin
    nsites = 10
    J=1.0
    h=0.2
    H_L = isoPEPS.hamiltonian(nsites,J,h)
    ground_state_energy, ground_state_vector = isoPEPS.eigs(H_L, nev=1, which=:SR)
    E_L = ground_state_energy[1]/nsites
    @test isapprox(E_L, -0.9120354170186685, atol=1e-4)
end

using Test, isoPEPS
import Optimisers
@testset "KrylovkitYao" begin
    nsites=10
    J=1.0
    h=0.2
    hami = ising_hamiltonian(nsites,J,h)
    E, V = ed_groundstate(hami)
    @test isapprox(E/10, -0.9120354170186685, atol=1e-4)
end

@testset "VariationalCircuit" begin
    nsites=10
    J=1.0
    h=0.2
    hami = ising_hamiltonian(nsites,J,h)
    @test isoPEPS.ishermitian(hami)
    c = isoPEPS.variational_circuit(nsites)
    isoPEPS.dispatch!(c, :random)
    params = isoPEPS.parameters(c)
    optimizer = Optimisers.setup(Optimisers.Adam(0.01), params)
    niter = 1000
   
    for i = 1:niter
        grad_input, grad_params = isoPEPS.expect'(hami, isoPEPS.zero_state(nsites) => c)
        Optimisers.update!(optimizer, params, grad_params)
        isoPEPS.dispatch!(c, params)
          
    end
    EG = isoPEPS.expect(hami, isoPEPS.zero_state(nsites) |> c)/nsites
    @test isapprox(EG, -0.9120354170186685, atol=1e-2) # Q:EG used by VariationalCircuit can't reach 1e-4
end



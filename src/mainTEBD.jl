using Printf
using Yao
using Yao.EasyBuild
using LinearAlgebra
using CairoMakie
using Test

include("ImTebd.jl")
include("Lanczos.jl")

nsites = 15
J=1.0
h=0.2
t=50
H = transverse_ising(nsites;J, h) |> cache
span = LinRange(1.0,2.0, 2)
y = []

#exact result from Lanczos
H_L = ising_hamiltonian(nsites, J, h)
ground_state_energy, ground_state_vector = eigs(H_L, nev=1, which=:SR)
E_L = ground_state_energy[1]/nsites
println("Ground state energy: ", E_L)

for p in span
    reg = rand_state(nqubits(H))
    reg |> itime_groundstate!(H; t,p)
    EG = expect(H, reg)/nsites
    deltaE= sqrt((EG-E_L)^2)
    push!(y,deltaE)
    println(EG)
    println(y)
end

fig = Figure()
ax1 = Axis(fig[1, 1], title = "imaginary time evolution", xlabel = "dτ", ylabel = "δE/L")
lines!(ax1, span, y, color = :red, label = "sin")
axislegend(ax1)
display(fig)
save("line_plot.png", fig)





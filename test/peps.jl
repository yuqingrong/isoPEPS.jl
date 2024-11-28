using IsoPEPS
using Test


@testset "truncate" begin
    mps=generate_mps(8,10; d=2)
    mps0=deepcopy(mps)
    mps0.tensors[3] = cat(mps0.tensors[3], zeros(8, 2, 4), dims=3)
    mps0.tensors[4] = cat(mps0.tensors[4], zeros(4, 2, 8), dims=1)
    @show size(mps0.tensors[3]),size(mps0.tensors[4])
    mps1=IsoPEPS.truncated_bmps(mps0.tensors,8)
    @show size(mps1[3]),size(mps1[4]),size(mps1[5])
    @test size(mps.tensors[3],3)==size(mps1[3],3)
end

@testset "peps" begin
    peps=generate_peps(4,3,3)
    bra_ket=IsoPEPS.contract_2peps(peps,peps)
    @show size(bra_ket[1][1])
    @show size(bra_ket[2][1])
    dmax=4
    #bmps=IsoPEPS.truncated_bmps(bra_ket[1],dmax)
    
    result=IsoPEPS.overlap_peps(bra_ket,dmax)
    @test  isless(0.0,result)
end
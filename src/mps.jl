mutable struct MPS{T,AT<:AbstractArray{T,3}}
    tensors::Vector{AT}
    canonical_center::Int
    function MPS(tensors::Vector{AT},canonical_center::Int=0) where {T, AT<:AbstractArray{T,3}}
        n=length(tensors)
        physical_dim=size(tensors[1],2)
        @assert canonical_center > 0 && canonical_center <= n || canonical_center == 0
        @assert all(i->size(tensors[i],2)==physical_dim && size(tensors[i],3)==size(tensors[mod1(i+1,n)],1),1:n)
        new{T,AT}(tensors,canonical_center)
    end
end

tensors(mps::MPS) = mps.tensors
function generate_mps(::Type{T},bond_dim::Int,nsites::Int;d::Int=2) where T
    tensors=[randn(T,1,d,bond_dim)]
    for i in 2:(nsites-1)
        push!(tensors,randn(T,bond_dim,d,bond_dim))
    end
    push!(tensors,randn(T,bond_dim,d,1))
    return MPS(tensors)
end

generate_mps(bond_dim::Int,nsites::Int; d::Int=2) = generate_mps(Float64, bond_dim,nsites;d)


function code_dot(bra::MPS,ket::MPS;optimizer=GreedyMethod())
    store=IndexStore()
    index_bra=Vector{Int}[]
    index_ket=Vector{Int}[]
    firstidx_bra=newindex!(store)
    previdx_bra=firstidx_bra
    firstidx_ket=newindex!(store)
    previdx_ket=firstidx_ket
    nsites=length(bra.tensors)
    for k = 1:nsites
        physidx=newindex!(store)
        nextidx_bra=k==nsites ? firstidx_bra : newindex!(store)
        nextidx_ket=k==nsites ? firstidx_ket : newindex!(store)
        push!(index_bra,[previdx_bra,physidx,nextidx_bra])
        push!(index_ket,[previdx_ket,physidx,nextidx_ket])
        previdx_bra=nextidx_bra
        previdx_ket=nextidx_ket
    end
    ixs=[index_bra...,index_ket...]
    size_dict=OMEinsum.get_size_dict(ixs,[bra.tensors...,ket.tensors...])
    return optimize_code(DynamicEinCode(ixs,Int[]),size_dict,optimizer)
end

function LinearAlgebra.dot(mps1::MPS, mps2::MPS)
    code = code_dot(mps1, mps2)
    return code(conj.(mps1.tensors)..., mps2.tensors...)[]
end



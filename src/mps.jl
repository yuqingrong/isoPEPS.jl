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
    code=optimize_code(DynamicEinCode(ixs,Int[]),size_dict,optimizer)
    return code,code(conj.(bra.tensors)..., ket.tensors...)[]
end


function vec2mps(v::AbstractVector; d=2, Dmax=typemax(Int), atol=1e-10)
    state = reshape(v, 1, length(v))
    tensors = typeof(reshape(state, 1, length(v), 1))[]
    nsite = round(Int, log2(length(v)) ÷ log2(d))
    @assert d^nsite == length(v) 
    for _ = 1:nsite-1
        state = reshape(state, (d * size(state, 1), size(state, 2) ÷ d))
        u, s, v, err = truncated_svd(state, Dmax, atol)
        push!(tensors, reshape(u, size(u, 1) ÷ d, d, size(u, 2)))
        state = s .* v
    end
    push!(tensors, reshape(state, size(state, 1), d, 1))
    return MPS(tensors)
end

function code_mps2vec(mps; optimizer=GreedyMethod())
    store = IndexStore()
    ixs = Vector{Int}[]
    iy = Int[]
    firstidx = newindex!(store)
    previdx = firstidx
    for k = 1:length(mps)
        physical = newindex!(store)
        nextidx = k == length(mps) ? firstidx : newindex!(store)
        push!(ixs, [previdx, physical, nextidx])
        push!(iy, physical)
        previdx = nextidx
    end
    size_dict = OMEinsum.get_size_dict(ixs, mps)
    code=optimize_code(DynamicEinCode(ixs, iy), size_dict, optimizer)
    return vec(code(mps...))
end


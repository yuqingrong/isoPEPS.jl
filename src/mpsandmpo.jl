struct IndexStore{LT}
    indices::Vector{LT}
end
IndexStore() = IndexStore(Int[])
Base.length(store::IndexStore) = length(store.indices)
function newindex!(store::IndexStore{LT}) where LT <: Integer
    index = length(store) == 0 ? 1 : store.indices[end] + 1
    push!(store.indices, index)
    return index
end

function truncated_svd(tmat, dmax::Int, atol::Real)
    u, s, v = svd(tmat)
    dmax = min(searchsortedfirst(s, atol, rev=true), dmax, length(s))
    return u[:, 1:dmax], s[1:dmax], v'[1:dmax, :], sum(s[dmax+1:end])
end

function mps_dot_mpo(mps::MPS,mpo::MPO)
    @assert length(mps.tensors)==length(mpo.tensors)
    new_mps=deepcopy(mps)
    for i in 1:length(mps.tensors)
        mps_mpo=ein"mpjn,ijk->mipnk"(mpo.tensors[i],mps.tensors[i])     
        new_shape=(size(mps_mpo,1)*size(mps_mpo,2),
                   size(mps_mpo,3),
                   size(mps_mpo,4)*size(mps_mpo,5))
        new_mps.tensors[i]=reshape(mps_mpo,new_shape)
    end
    return new_mps  
end
  




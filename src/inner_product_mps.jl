struct MPS{T}
    tensors::Vector{T}
    function MPS(tensors::Vector{T}) where T<:AbstractArray
        new{T}(tensors)
    end
end

function dagger_mps(mps::MPS{T}) where T<:AbstractArray
    daggered_tensors = Vector{AbstractArray}(undef, nsites)
    daggered_tensors[1]=conj(transpose(mps.tensors[1]))
    for i in 2:length(mps.tensors)-1
        daggered_tensors[i]= conj(permutedims(mps.tensors[i], (3, 2, 1))) # Conjugate and permute (3D)
    end
    daggered_tensors[end]= conj(transpose(mps.tensors[end]))
    return MPS(daggered_tensors)
end


function inner_product(mps1::MPS{T}, mps2::MPS{T}) where T<:AbstractArray

    @assert length(mps1.tensors) == length(mps2.tensors)
    mps2_dagger = dagger_mps(mps2)
    n_sites = length(mps1.tensors)
    indices = Vector{Int}[]
    # generate indices for mps1
    for i in 1:n_sites
        mps11_indices=[]
        push!(mps11_indices,i)
        if i==1
            push!(mps11_indices,i+n_sites)
        
        elseif i==n_sites
            push!(mps11_indices,i+n_sites-1)
        
        else
            push!(mps11_indices,i+n_sites-1)
            push!(mps11_indices,i+n_sites)
            
        end
        push!(indices,mps11_indices)
    end

    # generate indices for mps2
  
    for i in 1:n_sites
        mps22_indices=[]
        push!(mps22_indices,i)
        if i==1
            push!(mps22_indices,i+2*n_sites-1)
        elseif i==n_sites
            push!(mps22_indices,i+2*n_sites-2)
        else
            push!(mps22_indices,i+2*n_sites-1)
            push!(mps22_indices,i+2*n_sites-2)
        end
        push!(indices,mps22_indices)
    end
    
    eincode = OMEinsum.DynamicEinCode(indices, Int[])
    nested_ein = optimize_code(eincode, uniformsize(eincode, 2), TreeSA())
    return nested_ein(mps1.tensors..., mps2_dagger.tensors...)[1]
end
#eincode(mps1.tensors..., mps2.tensors...)[1]

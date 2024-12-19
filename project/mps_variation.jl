using IsoPEPS

nsites=10
bond_dim=2
h=0.5

function optimize_groundstate()
    psi=generate_mps(bond_dim,nsites)
    @show size(psi.tensors)
    params=vcat(map(vec, psi.tensors)...)
    #params = code_mps2vec(psi)
    @show size(params)
    H=transverse_ising_mpo(nsites,h)
    @show size(H.tensors)
    energy = code_sandwich(psi, H, psi)
    @show energy
    #=
    function energy_fn(params)
       
        #psi = vectomps(params)
        
        energy = code_sandwich(psi, H, psi) / code_dot(psi, psi)
        return energy
    end
    =#
    function update_mps_from_params!(psi, params)
        idx = 1
        for i in 1:length(psi.tensors)
            size_tensor = size(psi.tensors[i])
            n_elements = prod(size_tensor)
            psi.tensors[i] .= reshape(params[idx:idx + n_elements - 1], size_tensor)
            idx += n_elements
        end
    end

    function f(params)
        update_mps_from_params!(psi, params)
        code1, energy = code_sandwich(psi, H, psi)
        @show code1
        code2,norm_factor = code_dot(psi, psi)
        cost = energy / norm_factor
        return cost
    end
    #G = zeros(length(params))
    function g!(G, params)
        update_mps_from_params!(psi, params)
        code1, energy = code_sandwich(psi, H, psi)
        cost1, mg1 = OMEinsum.cost_and_gradient(code1, (conj.(psi.tensors)..., H.tensors..., psi.tensors...))
        mg1 = vcat(map(vec, mg1[nsites+1:2*nsites])...)

        code2, norm_factor = code_dot(psi, psi)
        cost2, mg2 = OMEinsum.cost_and_gradient(code2, (conj.(psi.tensors)..., psi.tensors...))
        mg2 = vcat(map(vec, mg2)...)

        flattened_mg = (cost2[] * mg1 - cost1[] * mg2) / cost2[]^2
        G .= flattened_mg
    end
    result = optimize(
        f,g!,
        params,
        LBFGS()
    )
end
    

optimize_groundstate()  
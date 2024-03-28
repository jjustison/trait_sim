

function compute_leaf_gammas(net,hyb_sort)
    
    leaf_gammas = Vector{Vector}()
    for (hyb_nd_no,par_dict) in hyb_sort
        for (par_nd_no, leaves) in par_dict
            
            ##only do things if the sorting has more than 1 leaf
            if length(leaves)>1
                
                ##Find the edge that connects the parent node and child node 
                hyb_nd = net.node[(x-> x.number == hyb_nd_no).(net.node)]
                par_nd = net.node[(x-> x.number == par_nd_no).(net.node)]
                par_child_edges = PhyloNetworks.getchildrenedges(par_nd[1])
                par_child_nodes = getchild.(par_child_edges)

                the_gamma = par_child_edges[par_child_nodes .== hyb_nd][1].gamma
                push!(leaf_gammas,[the_gamma,collect(leaves)])
            end
        end

    end
    return leaf_gammas
end



function outsouRcing(pop,ngenes,rep_s,gene_nomial,pt_leaf_gammas=nothing)
    R_seed = rand(Int32)

    @rput(pop,ngenes,R_seed,gene_nomial)
    (!isnothing(pt_leaf_gammas)) && @rput(pt_leaf_gammas)
    R"source('../../baba.R')"
    Random.seed!(rep_s)
    @rget(hib_vcv)

    return hib_vcv


end

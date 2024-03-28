

preorder!(net)
tipnames=tipLabels(net)
ntips=length(net.leaf)

##Set gammas to ensure correct major and minor edges - I don't nessecarily follow the Newick convention that is used to denote major/minor edges, this accounts for it if an issue.
major_hyb_edges = net.edge[(x->(x.isMajor & x.hybrid)).(net.edge)]
(x->setGamma!(x,0.5)).(major_hyb_edges)

##Get parent Tree info sorted
pt_dict = Dict(map(reverse,PhyloNetworks.parentTrees(net;report_hybsorting=true))) ##dictionary where the hybrid_sorting points to the correct parent tree
ref_hyb_sorts = PhyloNetworks.parentTreeProbs(net; pop=1) ##Get the hybrid sortings in the ordering of that this function gives
ref_hyb_sorts = (x-> x[2]).(ref_hyb_sorts)

sorted_pts = (x->pt_dict[x]).(ref_hyb_sorts) ##Parent trees in the correct ordering
writeMultiTopology(sorted_pts,"parent_trees.net")

##get dummy array of the hyb_sortings
sorts = (x -> x[2]).(PhyloNetworks.parentTreeProbs(net; pop=1))

#####################################
#### Prep things for outsouRcing ####
#####################################

##for each lineage sorting we need to associate lineage sorting probabilities with each event
include("./leaf_gammas.jl")
pt_leaf_gammas=(x-> compute_leaf_gammas(net,x)).(ref_hyb_sorts)
@rput(pt_leaf_gammas,ntips,hyb_time,tipnames)
R"source('../../outsouRcing_prep.R')"


###Trait evolution params
sigma2=3
mu=zeros(length(net.leaf))


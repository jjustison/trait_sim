using Revise
using PhyloNetworks
using Distributions
using DataFrames
using CSV
using GLM
using Random


cd("C:/Users/justison/Documents/chapt3")

#cd("./16_tree1")
cd("./8_tree1")



## Read in Phylogeny
#net = readTopology("((((A:3,((B:1,C:1):1)#H1:1::0.5):1,(#H1:1::0.5,D:3):1):1,((E:3,((F:1,G:1):1)#H2:1::0.5):1,(#H2:1::0.5,H:3):1):1):1.0,(((I:3,((J:1,K:1):1)#H3:1::0.5):1,(#H3:1::0.5,L:3):1):1,((M:3,((N:1,O:1):1)#H4:1::0.5):1,(#H4:1::0.5,P:3):1):1):1.0);");
net = readTopology("(((A:3,((B:1,C:1):1)#H1:1::0.5):1,(#H1:1::0.5,D:3):1):1,((E:3,((F:1,G:1):1)#H2:1::0.5):1,(#H2:1::0.5,H:3):1):1);");
#net = readTopology("(((A:2,B:2):2,(((C:1.0,D:1.0):1.0,(E:1.0,F:1.0):1.0):1.0)#H1:1::0.5):1,(#H1:1::0.5,(G:2,H:2):2):1);")
preorder!(net)
tipnames=tipLabels(net)

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





############################################
##### Generate & fit Observed PT props #####
############################################

##These are the things that get iterated over 
no_pops=50
pop_sizes=tan.(range(atan(1/1000),atan(1000),length=no_pops)) ### generate population sizes to range "evenly" from 0 to Inf
pop_dict= Dict(zip(pop_sizes,1:no_pops))
ngene_sizes = 1:5:250 
nreps = 10000

###Trait evolution params
sigma2=3
mu=zeros(length(net.leaf))



Random.seed!(138937) ##The exponent of the largest known Wagstaff prime
pars=DataFrame(
    i      = repeat(1:length(pop_sizes)  , inner=length(ngene_sizes)),
    j      = repeat(1:length(ngene_sizes), outer=length(pop_sizes)),
    pop    = repeat(pop_sizes            , inner=length(ngene_sizes)),
    ngenes = repeat(ngene_sizes          , outer=length(pop_sizes)),
    seed   = rand!(zeros(UInt128,length(pop_sizes)*length(ngene_sizes))) 
    )
CSV.write("pars.csv",pars) 


##Bastide model
bas_vcv = vcv(net)

##Major Tree model
maj_vcv = vcv(majorTree(net))

##Min pop size model
s_p = (PhyloNetworks.parentTreeProbs(net; pop=0.000000000001)) ##should be zero but this causes some type funkiness. Results are effectively the same, however
min_vcv = obs_vcv = PhyloNetworks.vcvParent3(net,s_p)

##Max pop size model
s_p = (PhyloNetworks.parentTreeProbs(net; pop=Inf))
max_vcv = obs_vcv = PhyloNetworks.vcvParent3(net,s_p)

##unscaled pop model
s_p = (PhyloNetworks.parentTreeProbs(net; pop=1))
uns_vcv = obs_vcv = PhyloNetworks.vcvParent3(net,s_p)



for (i,j,pop,ngenes,seed) in eachrow(pars)
    println(i," ",j)
    
    ##Do population level things
    i=pop_dict[pop]
    sorts_n_props = PhyloNetworks.parentTreeProbs(net; pop=pop)
    props = (x -> x[1]).(sorts_n_props)
    exp_vcv= PhyloNetworks.vcvParent3(net,sorts_n_props)


    ##Make data frame for storing results
    dat_ests= DataFrame(i = Int64[], j = Int64[], 
    bas_lik=Float64[], bas_mu=Float64[], bas_sigma=Float64[],
    exp_lik=Float64[], exp_mu=Float64[], exp_sigma=Float64[],
    tru_lik=Float64[], tru_mu=Float64[], tru_sigma=Float64[],
    maj_lik=Float64[], maj_mu=Float64[], maj_sigma=Float64[],
    max_lik=Float64[], max_mu=Float64[], max_sigma=Float64[],
    min_lik=Float64[], min_mu=Float64[], min_sigma=Float64[],
    uns_lik=Float64[], uns_mu=Float64[], uns_sigma=Float64[])


    Random.seed!(seed)
    rep_seeds = rand!(zeros(UInt128,nreps))
    for rep in 1:nreps
        rep_s = rep_seeds[rep]
        Random.seed!(rep_s)


        ##Make empirical VCVs with changing number of controlling genes
        gene_nomial = rand(Multinomial(ngenes,props)) ##number of genes belonging to each pt
        CSV.write("gene_no.csv",DataFrame(gene_no=gene_nomial))

        ##Call R script to generate hibbins vcv 


        obs_props = (gene_nomial)/ngenes
        obs_sorts_n_props = collect(zip(obs_props,sorts))
        obs_vcv = PhyloNetworks.vcvParent3(net,obs_sorts_n_props)

        ##Simulate trait given  the obs gene props
        sigma= sigma2 .* Matrix(obs_vcv)
        mvdist = MvNormal(mu,sigma)
        traits = rand(mvdist)
        traits = DataFrame(b = Vector(traits), tipNames = names(obs_vcv))

        ##Fit the models 
        bas_fit = phylolm(@formula(b ~ 1), traits, net,model="BMVCV",VCV=Matrix(bas_vcv))
        exp_fit = phylolm(@formula(b ~ 1), traits, net,model="BMVCV",VCV=Matrix(exp_vcv))
        tru_fit = phylolm(@formula(b ~ 1), traits, net,model="BMVCV",VCV=Matrix(obs_vcv))
        maj_fit = phylolm(@formula(b ~ 1), traits, net,model="BMVCV",VCV=Matrix(maj_vcv))
        max_fit = phylolm(@formula(b ~ 1), traits, net,model="BMVCV",VCV=Matrix(max_vcv))
        min_fit = phylolm(@formula(b ~ 1), traits, net,model="BMVCV",VCV=Matrix(min_vcv))
        uns_fit = phylolm(@formula(b ~ 1), traits, net,model="BMVCV",VCV=Matrix(uns_vcv))

        push!(dat_ests,[
            i,j,
            loglikelihood(bas_fit), mu_phylo(bas_fit), sigma2_phylo(bas_fit),
            loglikelihood(exp_fit), mu_phylo(exp_fit), sigma2_phylo(exp_fit),
            loglikelihood(tru_fit), mu_phylo(tru_fit), sigma2_phylo(tru_fit),
            loglikelihood(maj_fit), mu_phylo(maj_fit), sigma2_phylo(maj_fit),
            loglikelihood(max_fit), mu_phylo(max_fit), sigma2_phylo(max_fit),
            loglikelihood(min_fit), mu_phylo(min_fit), sigma2_phylo(min_fit),
            loglikelihood(uns_fit), mu_phylo(uns_fit), sigma2_phylo(uns_fit)    
        ])
    end
    ##Save the trait data
    data_tag = string(i)*"_"*string(j)
    CSV.write("fit_ests_"*data_tag*".csv",dat_ests)


end
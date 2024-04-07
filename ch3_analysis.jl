using Revise
using PhyloNetworks
using Distributions
using DataFrames
using CSV
using GLM
using Random
using RCall


#cd("C:/Users/justison/Documents/chapt3/8_tree1")
#cd("./8_tree1")

############################################
############################################
##### Generate & fit Observed PT props #####
############################################
############################################

#cd("C:/Users/justison/Documents/chapt3/8_tree1/obs")
#cd("C:/Users/justison/Documents/chapt3/8_tree2/obs") ### Still need to finish this!
#cd("C:/Users/justison/Documents/chapt3/8_tree3/obs")

#Random.seed!(138937) ##The exponent of the largest known Wagstaff prime
#Random.seed!(15643) ## tree2
#Random.seed!(16) ###tree3

function obs_analysis(folder_path,starter_sd,rw_start=1,rw_end=nothing)


    cd(folder_path)

    include("../common_net.jl")
    include("../../common_pars.jl")

    Random.seed!(starter_sd)

    ##These are the things that get iterated over 
    no_pops=25
    pop_sizes=exp.(collect(range(log(0.001),log(1000),length=no_pops))) ### generate population sizes to range "evenly" from 0 to Inf
    #ngene_sizes = reduce(vcat,[1,(10:10:250)])
    ngene_sizes = 1:25
    nreps = 10000




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
    min_vcv  = PhyloNetworks.vcvParent3(net,s_p)

    ##Max pop size model
    s_p = (PhyloNetworks.parentTreeProbs(net; pop=Inf))
    max_vcv  = PhyloNetworks.vcvParent3(net,s_p)

    ##unscaled pop model
    s_p = (PhyloNetworks.parentTreeProbs(net; pop=1))
    uns_vcv = PhyloNetworks.vcvParent3(net,s_p)


    isnothing(rw_end) && (pars_iter=eachrow(pars)[rw_start:end])
    (!isnothing(rw_end)) && (pars_iter=eachrow(pars)[rw_start:rw_end])

    for (i,j,pop,ngenes,seed) in pars_iter
        println(i," ",j)
        
        ##Do population level things
        sorts_n_props = PhyloNetworks.parentTreeProbs(net; pop=pop)
        props = (x -> x[1]).(sorts_n_props)
        exp_vcv= PhyloNetworks.vcvParent3(net,sorts_n_props)


        ##Make data frame for storing results
        dat_ests= DataFrame(i = Int64[], j = Int64[], 
        bas_lik=Float64[], bas_mu=Float64[], bas_sigma=Float64[],
        #hib_lik=Float64[], hib_mu=Float64[], hib_sigma=Float64[],
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
        

            ##Call R script to generate gene trees and hibbins vcv
            #outsouRcing(pop,ngenes,rep_s,gene_nomial)

            ##Gernerate observed vcv
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
            obs_fit = phylolm(@formula(b ~ 1), traits, net,model="BMVCV",VCV=Matrix(obs_vcv))
            maj_fit = phylolm(@formula(b ~ 1), traits, net,model="BMVCV",VCV=Matrix(maj_vcv))
            max_fit = phylolm(@formula(b ~ 1), traits, net,model="BMVCV",VCV=Matrix(max_vcv))
            min_fit = phylolm(@formula(b ~ 1), traits, net,model="BMVCV",VCV=Matrix(min_vcv))
            uns_fit = phylolm(@formula(b ~ 1), traits, net,model="BMVCV",VCV=Matrix(uns_vcv))

            push!(dat_ests,[
                i,j,
                loglikelihood(bas_fit), mu_phylo(bas_fit), sigma2_phylo(bas_fit),
                loglikelihood(exp_fit), mu_phylo(exp_fit), sigma2_phylo(exp_fit),
                loglikelihood(obs_fit), mu_phylo(obs_fit), sigma2_phylo(obs_fit),
                loglikelihood(maj_fit), mu_phylo(maj_fit), sigma2_phylo(maj_fit),
                loglikelihood(max_fit), mu_phylo(max_fit), sigma2_phylo(max_fit),
                loglikelihood(min_fit), mu_phylo(min_fit), sigma2_phylo(min_fit),
                loglikelihood(uns_fit), mu_phylo(uns_fit), sigma2_phylo(uns_fit)    
            ])
        end
        ##Save the trait data
        data_tag = string(i)*"_"*string(j)
        CSV.write("fits/fit_ests_"*data_tag*".csv",dat_ests)
    end

end


####################################################
####################################################
##### Generate & fit Bastide Model over gammas #####
####################################################
####################################################


#cd("C:/Users/justison/Documents/chapt3/8_tree1/bas")
#cd("C:/Users/justison/Documents/chapt3/8_tree2/bas")
#cd("C:/Users/justison/Documents/chapt3/8_tree3/bas")
#cd("C:/Users/justison/Documents/chapt3/16_tree1/bas")

#Random.seed!(82589933) ##The exponent of the largest known Mersenne prime
#Random.seed!(6543) ## Tree 2
#Random.seed!(3784278)
#Random.seed!(8675309)




function bas_analysis(folder_path,starter_sd,rw_start=1,rw_end=nothing)

    cd(folder_path)

    include("../common_net.jl")
    include("../../common_pars.jl")

    Random.seed!(starter_sd)

    ##These are the things that get iterated over 
    no_pops=51
    gammas=range(0.5,1.0,length=no_pops) 
    nreps = 10000




    pars=DataFrame(
        i      = 1:no_pops,
        gammas = gammas,
        seed   = rand!(zeros(UInt128,no_pops))
    )
    CSV.write("pars.csv",pars) 

    ##Used if continuing a partially completed simulation set
    isnothing(rw_end) && (pars_iter=eachrow(pars)[rw_start:end])
    (!isnothing(rw_end)) && (pars_iter=eachrow(pars)[rw_start:rw_end])


    for (i,gamma,seed) in pars_iter
        println(i," gamma")
        
        ###Set gamma values
        major_hyb_edges = net.edge[(x->(x.isMajor & x.hybrid)).(net.edge)]
        (x->setGamma!(x,gamma)).(major_hyb_edges)

        ##Bastide model
        bas_vcv = vcv(net)

        ##Major Tree model
        maj_vcv = vcv(majorTree(net))

        ##Min pop size model
        s_p_min = (PhyloNetworks.parentTreeProbs(net; pop=0.000000000001)) ##should be zero but this causes some type funkiness. Results are effectively the same, however
        min_vcv = PhyloNetworks.vcvParent3(net,s_p_min)

        ##Max pop size model
        s_p_max = (PhyloNetworks.parentTreeProbs(net; pop=Inf))
        max_vcv = PhyloNetworks.vcvParent3(net,s_p_max)

        ##unscaled pop model
        s_p_uns = (PhyloNetworks.parentTreeProbs(net; pop=10))
        uns_vcv = PhyloNetworks.vcvParent3(net,s_p_uns)


        ##Make data frame for storing results
        dat_ests= DataFrame(i = Int64[], 
        bas_lik=Float64[], bas_mu=Float64[], bas_sigma=Float64[],
        maj_lik=Float64[], maj_mu=Float64[], maj_sigma=Float64[],
        max_lik=Float64[], max_mu=Float64[], max_sigma=Float64[],
        min_lik=Float64[], min_mu=Float64[], min_sigma=Float64[],
        uns_lik=Float64[], uns_mu=Float64[], uns_sigma=Float64[])


        Random.seed!(seed)
        rep_seeds = rand!(zeros(UInt128,nreps))
        for rep in 1:nreps
            rep_s = rep_seeds[rep]
            Random.seed!(rep_s)
        
            ##Simulate trait given  the obs gene props
            sigma= sigma2 .* Matrix(bas_vcv)
            mvdist = MvNormal(mu,sigma)
            traits = rand(mvdist)
            traits = DataFrame(b = Vector(traits), tipNames = names(bas_vcv))

            ##Fit the models 
            bas_fit = phylolm(@formula(b ~ 1), traits, net,model="BMVCV",VCV=Matrix(bas_vcv))
            maj_fit = phylolm(@formula(b ~ 1), traits, net,model="BMVCV",VCV=Matrix(maj_vcv))
            max_fit = phylolm(@formula(b ~ 1), traits, net,model="BMVCV",VCV=Matrix(max_vcv))
            min_fit = phylolm(@formula(b ~ 1), traits, net,model="BMVCV",VCV=Matrix(min_vcv))
            uns_fit = phylolm(@formula(b ~ 1), traits, net,model="BMVCV",VCV=Matrix(uns_vcv))

            push!(dat_ests,[
                i,
                loglikelihood(bas_fit), mu_phylo(bas_fit), sigma2_phylo(bas_fit),
                loglikelihood(maj_fit), mu_phylo(maj_fit), sigma2_phylo(maj_fit),
                loglikelihood(max_fit), mu_phylo(max_fit), sigma2_phylo(max_fit),
                loglikelihood(min_fit), mu_phylo(min_fit), sigma2_phylo(min_fit),
                loglikelihood(uns_fit), mu_phylo(uns_fit), sigma2_phylo(uns_fit)    
            ])
        end
        ##Save the trait data
        data_tag = string(i)
        CSV.write("fits/fit_ests_"*data_tag*".csv",dat_ests)
    end

end

##################################################
##################################################
##### Generate & fit exp over pop and gamma  #####
##################################################
##################################################


#cd("C:/Users/justison/Documents/chapt3/8_tree1/exp")
#cd("C:/Users/justison/Documents/chapt3/8_tree2/exp")
#cd("C:/Users/justison/Documents/chapt3/8_tree3/exp")

#Random.seed!(13) ##random prime
#Random.seed!(8675309) ##Tree 2

function exp_analysis(folder_path,starter_sd,rw_start=1,rw_end=nothing)
    cd(folder_path)
    include("../common_net.jl")
    include("../../common_pars.jl")
    Random.seed!(starter_sd)

    ##These are the things that get iterated over 
    no_pops=25
    pop_sizes=exp.(collect(range(log(0.001),log(1000),length=25))) ### generate population sizes to range "evenly" from 0 to Inf
    no_gammas=26
    gammas=range(0.5,1.0,length=no_gammas) ### generate population sizes to range "evenly" from 0 to Inf
    nreps = 10000


    pars=DataFrame(
        i      = repeat(1:no_pops  , inner=no_gammas),
        j      = repeat(1:no_gammas, outer=no_pops),
        pop    = repeat(pop_sizes  , inner=no_gammas),
        gammas = repeat(gammas     , outer=no_pops),
        seed   = rand!(zeros(UInt128,no_pops*no_gammas)) 
        )

    CSV.write("pars.csv",pars) 

    isnothing(rw_end) && (pars_iter=eachrow(pars)[rw_start:end])
    (!isnothing(rw_end)) && (pars_iter=eachrow(pars)[rw_start:rw_end])

    for (i,j,pop,gamma,seed) in pars_iter
        println(i," ",j)
        
        ###Set gamma values
        major_hyb_edges = net.edge[(x->(x.isMajor & x.hybrid)).(net.edge)]
        (x->setGamma!(x,gamma)).(major_hyb_edges)
        

        ##Do population level things
        sorts_n_props = PhyloNetworks.parentTreeProbs(net; pop=pop)
    
        ##Expected model    
        exp_vcv= PhyloNetworks.vcvParent3(net,sorts_n_props)
            
        ##Bastide model
        bas_vcv = vcv(net)

        ##Major Tree model
        maj_vcv = vcv(majorTree(net))

        ##Min pop size model
        s_p = (PhyloNetworks.parentTreeProbs(net; pop=0.000000000001)) ##should be zero but this causes some type funkiness. Results are effectively the same, however
        min_vcv = PhyloNetworks.vcvParent3(net,s_p)

        ##Max pop size model
        s_p = (PhyloNetworks.parentTreeProbs(net; pop=Inf))
        max_vcv = PhyloNetworks.vcvParent3(net,s_p)

        ##unscaled pop model
        s_p = (PhyloNetworks.parentTreeProbs(net; pop=1))
        uns_vcv = PhyloNetworks.vcvParent3(net,s_p)


        ##Make data frame for storing results
        dat_ests= DataFrame(i = Int64[], j = Int64[], 
        bas_lik=Float64[], bas_mu=Float64[], bas_sigma=Float64[],
        exp_lik=Float64[], exp_mu=Float64[], exp_sigma=Float64[],
        maj_lik=Float64[], maj_mu=Float64[], maj_sigma=Float64[],
        max_lik=Float64[], max_mu=Float64[], max_sigma=Float64[],
        min_lik=Float64[], min_mu=Float64[], min_sigma=Float64[],
        uns_lik=Float64[], uns_mu=Float64[], uns_sigma=Float64[])


        Random.seed!(seed)
        rep_seeds = rand!(zeros(UInt128,nreps))
        for rep in 1:nreps
            rep_s = rep_seeds[rep]
            Random.seed!(rep_s)

            ##Simulate trait given  the obs gene props
            sigma= sigma2 .* Matrix(exp_vcv)
            mvdist = MvNormal(mu,sigma)
            traits = rand(mvdist)
            traits = DataFrame(b = Vector(traits), tipNames = names(exp_vcv))

            ##Fit the models 
            bas_fit = phylolm(@formula(b ~ 1), traits, net,model="BMVCV",VCV=Matrix(bas_vcv))
            exp_fit = phylolm(@formula(b ~ 1), traits, net,model="BMVCV",VCV=Matrix(exp_vcv))
            maj_fit = phylolm(@formula(b ~ 1), traits, net,model="BMVCV",VCV=Matrix(maj_vcv))
            max_fit = phylolm(@formula(b ~ 1), traits, net,model="BMVCV",VCV=Matrix(max_vcv))
            min_fit = phylolm(@formula(b ~ 1), traits, net,model="BMVCV",VCV=Matrix(min_vcv))
            uns_fit = phylolm(@formula(b ~ 1), traits, net,model="BMVCV",VCV=Matrix(uns_vcv))

            push!(dat_ests,[
                i,j,
                loglikelihood(bas_fit), mu_phylo(bas_fit), sigma2_phylo(bas_fit),
                loglikelihood(exp_fit), mu_phylo(exp_fit), sigma2_phylo(exp_fit),
                loglikelihood(maj_fit), mu_phylo(maj_fit), sigma2_phylo(maj_fit),
                loglikelihood(max_fit), mu_phylo(max_fit), sigma2_phylo(max_fit),
                loglikelihood(min_fit), mu_phylo(min_fit), sigma2_phylo(min_fit),
                loglikelihood(uns_fit), mu_phylo(uns_fit), sigma2_phylo(uns_fit)    
            ])
        end
        ##Save the trait data
        data_tag = string(i)*"_"*string(j)
        CSV.write("fits/fit_ests_"*data_tag*".csv",dat_ests)
    end

end



############################################
############################################
#####   Generate & fit Hibbins model   #####
############################################
############################################

#cd("C:/Users/justison/Documents/chapt3/8_tree1/hib")
#cd("C:/Users/justison/Documents/chapt3/8_tree2/hib") 
#cd("C:/Users/justison/Documents/chapt3/8_tree3/hib")

#Random.seed!(76546) ##The exponent of the largest known Wagstaff prime
#Random.seed!(546) ## tree2
#Random.seed!(96030789) ###tree3

function hib_analysis(folder_path,starter_sd,rw_start=1,rw_end=nothing)

    cd(folder_path)

    include("../common_net.jl")
    include("../../common_pars.jl")

    Random.seed!(starter_sd)

    ##These are the things that get iterated over 
    no_pops=25
    pop_sizes=exp.(collect(range(log(0.001),log(1000),length=no_pops))) ### generate population sizes to range "evenly" from 0 to Inf
    #ngene_sizes = reduce(vcat,[1,(10:10:250)])
    ngene_sizes = 1:25
    nreps = 10000
    
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
    min_vcv = PhyloNetworks.vcvParent3(net,s_p)

    ##Max pop size model
    s_p = (PhyloNetworks.parentTreeProbs(net; pop=Inf))
    max_vcv = PhyloNetworks.vcvParent3(net,s_p)

    ##unscaled pop model
    s_p = (PhyloNetworks.parentTreeProbs(net; pop=1))
    uns_vcv = PhyloNetworks.vcvParent3(net,s_p)

    isnothing(rw_end) && (pars_iter=eachrow(pars)[rw_start:end])
    (!isnothing(rw_end)) && (pars_iter=eachrow(pars)[rw_start:rw_end])

    for (i,j,pop,ngenes,seed) in pars_iter
        println(i," ",j)
        
        ##Do population level things
        sorts_n_props = PhyloNetworks.parentTreeProbs(net; pop=pop)
        props = (x -> x[1]).(sorts_n_props)
        exp_vcv= PhyloNetworks.vcvParent3(net,sorts_n_props)


        ##Make data frame for storing results
        dat_ests= DataFrame(i = Int64[], j = Int64[], 
        bas_lik=Float64[], bas_mu=Float64[], bas_sigma=Float64[],
        hib_lik=Float64[], hib_mu=Float64[], hib_sigma=Float64[],
        exp_lik=Float64[], exp_mu=Float64[], exp_sigma=Float64[],
        tru_lik=Float64[], tru_mu=Float64[], tru_sigma=Float64[],
        maj_lik=Float64[], maj_mu=Float64[], maj_sigma=Float64[],
        max_lik=Float64[], max_mu=Float64[], max_sigma=Float64[],
        min_lik=Float64[], min_mu=Float64[], min_sigma=Float64[],
        uns_lik=Float64[], uns_mu=Float64[], uns_sigma=Float64[])


        Random.seed!(seed)
        rep_seeds = rand!(zeros(UInt128,nreps))
        for rep in 1:nreps
            ((rep % 100)==0) && println(i," ",j," ",rep)
            rep_s = rep_seeds[rep]
            Random.seed!(rep_s)


            ##Make empirical VCVs with changing number of controlling genes
            gene_nomial = rand(Multinomial(ngenes,props)) ##number of genes belonging to each pt
        

            ##Call R script to generate gene trees and hibbins vcv
            hib_vcv=outsouRcing(pop,ngenes,rep_s,gene_nomial)

            ##Gernerate observed vcv
            obs_props = (gene_nomial)/ngenes
            obs_sorts_n_props = collect(zip(obs_props,sorts))
            obs_vcv = PhyloNetworks.vcvParent3(net,obs_sorts_n_props)

            ##Simulate trait given  the obs gene props
            sigma= sigma2 .* hib_vcv
            mvdist = MvNormal(mu,sigma)
            traits = rand(mvdist)
            traits = DataFrame(b = Vector(traits), tipNames = names(obs_vcv))

            ##Fit the models 
            bas_fit = phylolm(@formula(b ~ 1), traits, net,model="BMVCV",VCV=Matrix(bas_vcv))
            hib_fit = phylolm(@formula(b ~ 1), traits, net,model="BMVCV",VCV=Matrix(hib_vcv))
            exp_fit = phylolm(@formula(b ~ 1), traits, net,model="BMVCV",VCV=Matrix(exp_vcv))
            tru_fit = phylolm(@formula(b ~ 1), traits, net,model="BMVCV",VCV=Matrix(obs_vcv))
            maj_fit = phylolm(@formula(b ~ 1), traits, net,model="BMVCV",VCV=Matrix(maj_vcv))
            max_fit = phylolm(@formula(b ~ 1), traits, net,model="BMVCV",VCV=Matrix(max_vcv))
            min_fit = phylolm(@formula(b ~ 1), traits, net,model="BMVCV",VCV=Matrix(min_vcv))
            uns_fit = phylolm(@formula(b ~ 1), traits, net,model="BMVCV",VCV=Matrix(uns_vcv))

            push!(dat_ests,[
                i,j,
                loglikelihood(bas_fit), mu_phylo(bas_fit), sigma2_phylo(bas_fit),
                loglikelihood(hib_fit), mu_phylo(hib_fit), sigma2_phylo(hib_fit),
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
        CSV.write("fits/fit_ests_"*data_tag*".csv",dat_ests)
    end

end





############################################################
############################################################
########### Covaraiance distance Analysis for VCV ##########
############################################################
############################################################

####Copied code from the obs_analysis function but eliminated superfulous code
function obs_dists(folder_path,starter_sd,rw_start=1,rw_end=nothing)



    cd(folder_path)

    include("../common_net.jl")
    include("../../common_pars.jl")

    Random.seed!(starter_sd)

    ##These are the things that get iterated over 
    no_pops=25
    pop_sizes=exp.(collect(range(log(0.001),log(1000),length=no_pops))) ### generate population sizes to range "evenly" from 0 to Inf
    #ngene_sizes = reduce(vcat,[1,(10:10:250)])
    ngene_sizes = 1:25
    nreps = 10000




    pars=DataFrame(
        i      = repeat(1:length(pop_sizes)  , inner=length(ngene_sizes)),
        j      = repeat(1:length(ngene_sizes), outer=length(pop_sizes)),
        pop    = repeat(pop_sizes            , inner=length(ngene_sizes)),
        ngenes = repeat(ngene_sizes          , outer=length(pop_sizes)),
        seed   = rand!(zeros(UInt128,length(pop_sizes)*length(ngene_sizes))) 
        )
    CSV.write("pars.csv",pars) 



    isnothing(rw_end) && (pars_iter=eachrow(pars)[rw_start:end])
    (!isnothing(rw_end)) && (pars_iter=eachrow(pars)[rw_start:rw_end])

    for (i,j,pop,ngenes,seed) in pars_iter
        println(i," ",j)
        
        ##Do population level things
        sorts_n_props = PhyloNetworks.parentTreeProbs(net; pop=pop)
        props = (x -> x[1]).(sorts_n_props)
        exp_vcv= PhyloNetworks.vcvParent3(net,sorts_n_props)


        ##Make data frame for storing results
        dat_dists= DataFrame(rho_dist = Float64[], vcv_dist = Float64[])


        Random.seed!(seed)
        rep_seeds = rand!(zeros(UInt128,nreps))
        for rep in 1:nreps
            rep_s = rep_seeds[rep]
            Random.seed!(rep_s)


            ##Make empirical VCVs with changing number of controlling genes
            gene_nomial = rand(Multinomial(ngenes,props)) ##number of genes belonging to each pt
        

            ##Gernerate observed vcv
            obs_props = (gene_nomial)/ngenes
            obs_sorts_n_props = collect(zip(obs_props,sorts))
            obs_vcv = PhyloNetworks.vcvParent3(net,obs_sorts_n_props)

            ###compute the distance between the proportions 
            obs_rho_dist = sqrt(sum((props .- obs_props).^2))

            ##Compute the distance between covariance matrices
            @rput(obs_vcv,exp_vcv)
            R"obs_vcv_dist <- relative.eigen(S1=obs_vcv, S2=exp_vcv)$distCov"
            @rget(obs_vcv_dist)




            push!(dat_dists,[obs_rho_dist,obs_vcv_dist])
        end
        ##Save the trait data
        data_tag = string(i)*"_"*string(j)
        CSV.write("fits/rho_vcv_dists_"*data_tag*".csv",dat_dists)
    end

end
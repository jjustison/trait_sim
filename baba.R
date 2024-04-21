parent_trees<-read.tree('./parent_trees.net')
set.seed(R_seed)## set seed for reproducibility. 'R_seed' is imported in the hib_analysis() function 
gene_trees<-list()
for(pt_ind in 1:length(parent_trees)){
  pt<-parent_trees[[pt_ind]]
  pt<-collapse.singles(pt)

  num_genes<-gene_nomial[pt_ind]
  leaf_gammas<-pt_leaf_gammas[[pt_ind]]
  
  completed_sims<-0
  while(completed_sims<num_genes){
    gene_tree<-generate_gene_tree_msc(
      pt,population_sizes = pop) ##All other parameters are assumed to be 1, what we want
    gt<-gene_tree$tree
    prob_sample<-1
    for(l_gamma in leaf_gammas){
      tip_nos<-match(l_gamma[[2]],gt$tip.label)

      tip_branches <- gt$edge.length[gt$edge[,2] %in% tip_nos]
      coal_over<-tip_branches>hyb_time ##Branches that are longer than hyb time coalesce after the hybridization
      gps<-0
      if(all(coal_over)){
        gps<-gps-1
      }
      gps<-gps+sum(coal_over)
      prob_sample<-prob_sample*(l_gamma[[1]]^gps)
    }
    
    if(prob_sample!=1 & (runif(1)>=prob_sample)){ ##sample fails
      next
    }
    gene_trees<-c(gene_trees,list(gt))
    completed_sims<-completed_sims+1
  }
}

gene_freqs<-rep(1/ngenes,ngenes)
gene_trees<-c(gene_trees,list(gene_freqs))
hib_vcv<-trees_to_vcv(tree_list = gene_trees)
matched_tips<- match(tipnames,colnames(hib_vcv)) ##match tip name ordering to the ones given
hib_vcv<-hib_vcv[matched_tips,matched_tips]



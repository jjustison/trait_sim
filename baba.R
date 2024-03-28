parent_trees<-read.tree('./parent_trees.net')

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
    
    prob_sample<-1
    
    for(l_gamma in leaf_gammas){
      coal_time<-5-gene_tree$gene_clade_times[getMRCA(gene_tree$tree,l_gamma[[2]])]
      if(coal_time>hyb_time){
        prob_sample<-prob_sample*l_gamma[[1]]
      }
    }
    if(prob_sample!=1 & (runif(1)>=prob_sample)){ ##sample fails
      next
    }
    gene_trees<-c(gene_trees,list(gene_tree$tree))
    completed_sims<-completed_sims+1
  }
}

gene_freqs<-rep(1/ngenes,ngenes)
gene_trees<-c(gene_trees,list(gene_freqs))
hib_vcv<-trees_to_vcv(tree_list = gene_trees)
matched_tips<- match(tipnames,colnames(hib_vcv)) ##match tip name ordering to the ones given
hib_vcv<-hib_vcv[matched_tips,matched_tips]



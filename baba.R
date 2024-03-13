library(treeducken)
library(SiPhyNetwork)

##From : http://blog.phytools.org/2012/01/function-to-get-descendant-node-numbers.html
getDescendants<-function(tree,node,curr=NULL){
  if(is.null(curr)) curr<-vector()
  daughters<-tree$edge[which(tree$edge[,1]==node),2]
  curr<-c(curr,daughters)
  w<-which(daughters>=length(tree$tip))
  if(length(w)>0) for(i in 1:length(w))
    curr<-getDescendants(tree,daughters[w[i]],curr)
  return(curr)
}


parent_trees<-read.tree('./8_tree1/parent_trees.net')
parent_trees[[1]]
tree<-parent_trees[[1]]

ntips<-8
get_hyb_nodes<-function(tree,ntips){
  hyb_nodes<-which(tree$node.label!="")+ntips
  kept<-list()
  for(nd in hyb_nodes){
    descs<-getDescendants(tree=tree,node=nd)
    tips<-descs[descs<=ntips]
    if(length(tips)>1){
      kept<-c(kept,list(list(nd,tips)))
    }
  }
  return(kept)
}




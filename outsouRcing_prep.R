library(castor)
library(SiPhyNetwork)
library(ape)
library(seastaR)
library(vcvComp)
library(partitions)

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

get_hyb_nodes<-function(tree,ntips){
  nd_times<-node.depth.edgelength(tree)
  hyb_nodes<-which(tree$node.label!="")+ntips
  kept<-list()
  for(nd in hyb_nodes){
    descs<-getDescendants(tree=tree,node=nd)
    tips<-descs[descs<=ntips]
    if(length(tips)>1){
      tip_names<-tree$tip.label[tips]
      kept<-c(kept,list(list(node=nd,
                             nd_time=nd_times[nd],
                             tip_nos=tips,
                             tip_names=tip_names)))
    }
  }
  return(kept)
}








summarize_data<-function(path,n_iters,models,true_mu=0,true_sig=3){
  nmodels<-length(models)
  pars_name<-paste(path,'/pars.csv',sep = '')
  pars <- read.csv(file =pars_name)
  pars$seed<-NULL ##Save memory, ride a cowboy
  

  
  ##extract n_iters columns from pars
  iter_names<-rep('',nrow(pars))
  for(iter_iter in 1:n_iters){
    iter_names<-paste(iter_names,'_',pars[,iter_iter],sep='')
  }
  par_names<-colnames(pars)
  # pars<-pars[,-(1:n_iters),drop=FALSE]

  

  ##Construct analyzed frame that stores the analyzed results
  col_names<-c('ave','var','mae','me','0.25','0.75','min','max')
  col_names<-c(par_names,'model',paste(col_names,rep(c('_mu','_sig'),each=8),sep=''))
  comp_dat<-data.frame(matrix(NA,nrow=nrow(pars)*nmodels,ncol=length(col_names)))
  colnames(comp_dat)<- col_names
  
  #####NOTE: this didn't really seem to save efficiency. I.e. I don't think copying the columns is the bottleneck but more likely the operations on the columns
  # ##set up pointer shenanigans, hopefully for efficiency
  # iter_name<-iter_names[1]  
  # raw_dat<-paste(path,'/fit_ests',iter_name,'.csv',sep='')
  # col_offset<-0
  # mu_col  %=% raw_dat[,(col_offset+2)]
  # sig_col %=% raw_dat[,(col_offset+3)]
  # 
  
  for(rw in 1:nrow(pars)){
    iter_name<-iter_names[rw]  
    print(iter_name)
    raw_dat<-paste(path,'/fits/fit_ests',iter_name,'.csv',sep='')
    raw_dat<-read.csv(file=raw_dat)
    raw_dat<-raw_dat[,-(1:n_iters)] ##get rid of columns specifying the iteration
    
    

    
    for(model_no in 1:nmodels){
        model_name<-models[model_no]
        
        col_offset<-(model_no-1)*3
        mu_col  <- raw_dat[,(col_offset+2)]
        sig_col <- raw_dat[,(col_offset+3)]

    
        mu_quant<-quantile(mu_col,probs=c(0.25,0.75))
        sig_quant<-quantile(sig_col,probs=c(0.25,0.75))
        comp_rw<- c(pars[rw,],
                    model_name,
                    mean(mu_col),var(mu_col),mean(abs(mu_col-true_mu)),mean(mu_col-true_mu),mu_quant[1],mu_quant[2],min(mu_col),max(mu_col),
                    mean(sig_col),var(sig_col),mean(abs(sig_col-true_sig)),mean(sig_col-true_sig),sig_quant[1],sig_quant[2],min(mu_col),max(mu_col))
        
        rw_offset<- ((rw-1)*nmodels)+model_no
        comp_dat[rw_offset,] <- comp_rw
    }
  }
  rownames(comp_dat)<-NULL
  write.csv(comp_dat,paste(path,'/compiled_data.csv',sep=''),row.names = F)
}


# path<-"./8_tree1/obs"
b_n_iters<-1
oe_n_iters<-2
bas_models<-c('bas','maj','max','min','uns')
exp_models<-c('bas','exp','maj','max','min','uns')
obs_models<-c('bas','exp','obs','maj','max','min','uns')

true_mu<-0
true_sig<-3


###Summarize all tree1 sims
summarize_data(path="./8_tree1/bas",n_iters=b_n_iters, models = bas_models)
summarize_data(path="./8_tree1/exp",n_iters=oe_n_iters,models = exp_models)
summarize_data(path="./8_tree1/obs",n_iters=oe_n_iters,models = obs_models)

###Summarize all tree2 sims
summarize_data(path="./8_tree2/bas",n_iters=b_n_iters, models = bas_models)
summarize_data(path="./8_tree2/exp",n_iters=oe_n_iters,models = exp_models)
summarize_data(path="./8_tree2/obs",n_iters=oe_n_iters,models = obs_models)

###Summarize all tree3 sims
summarize_data(path="./8_tree3/bas",n_iters=b_n_iters, models = bas_models)
summarize_data(path="./8_tree3/exp",n_iters=oe_n_iters,models = exp_models)
summarize_data(path="./8_tree3/obs",n_iters=oe_n_iters,models = obs_models)




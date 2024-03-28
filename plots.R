library(ggplot2)
library(cowplot)
library(tidyr)


#####################################################
#####################################################
########                                     ########
######## Section for plotting Tree Functions ########
########                                     ########
#####################################################
#####################################################

col_pal <- c('#CC6677','#332288','#DDCC77','#117733','#88CCEE','#44AA99','#999933','#AA4499')
all_models<-c('maj','bas','exp','obs','hib','min','max','uns')




plot_single_all_mod <- function(file_ext,dat,file_name,
                                xvar,yvar){
  col_pal <- c('#CC6677','#332288','#DDCC77','#117733','#88CCEE','#44AA99','#999933','#AA4499')
  all_models<-c('maj','bas','exp','obs','hib','min','max','uns')
  
  
  
  present_models<- all_models %in% dat$model
  present_cols<-col_pal[present_models]
  present_models<- all_models[present_models]
  

  name_dict<-list(
    gamma = expression('Inheritance Proportions '*gamma),
    ave_mu= expression('Average Initial Trait Value '*mu),
    sd_mu = expression('SD of Initial Trait Value '*mu),
    mae_mu= expression('MAE of Initial Trait Value '*mu),
    me_mu = expression('Mean Error of Initial Trait Value '*mu),
    ave_sig= expression('Average Evolutionary Rate '*sigma^2),
    sd_sig = expression('SD of Evolutionary Rate '*sigma^2),
    mae_sig= expression('MAE of '*sigma^2),
    me_sig = expression('Mean Error of '*sigma^2)
  )
  
  
  ggplot(data=dat,mapping = aes(x=!!sym(xvar),y=!!sym(yvar),color=factor(x=model,levels=present_models)))+
    geom_line(linewidth=3)+
    xlab(name_dict[[xvar]])+
    ylab(name_dict[[yvar]])+
    theme(
      axis.title=element_text(size=30,face="bold"),
      legend.position = "none",
      axis.text=element_text(size=15))+
    scale_color_manual(values=present_cols)
  
  ggsave(paste(file_ext,'/figs/',file_name,'.png',sep=''),width = 10,height = 10,dpi=600)
}


plot_iter_all_single<- function(file_ext,dat,
                                file_name,iter_loc,
                                xvar,yvar,
                                iter,iter_vals){
  ##Make single plots for all iterations
  for(val in iter_vals){
    dat_subset<- dat[dat[,iter]==val,] ##Only grab the data that has val 
    iter_name <- dat_subset[1,which(colnames(dat)==iter)-2] ##match the value to the respective i/j numbering in the pars.csv
    
    iter_folder<-paste(iter_loc,'/',iter,'_',iter_name,sep='')
    dir.create(paste(file_ext,'/figs/',iter_folder,sep=''),showWarnings = F)
    new_file_name<-paste('/',iter_folder,'/',y,sep='')
    plot_single_all_mod(file_ext,dat_subset,new_file_name,
                        xvar,yvar)
    
  }
}


make_model_legend <- function(dat,file_ext){
  col_pal <- c('#CC6677','#332288','#DDCC77','#117733','#88CCEE','#44AA99','#999933','#AA4499')
  all_models<-c('maj','bas','exp','obs','hib','min','max','uns')
  
  model_names<-expression('Major Tree',
                          'Network',
                          'Sorting, known '*N_e,
                          'Sorting, Observed Proportions',
                          'Quantitative Trait',
                          'Sorting, '*N[e]==0,
                          'Sorting, '*N[e]== infinity,
                          'Sorting, '*N[e]==1)
  
  present_models<- all_models %in% dat$model
  present_cols<-col_pal[present_models]
  model_names <-model_names[present_models]
  present_models<- all_models[present_models]
  
  
  leg_plot<-ggplot(data=dat, mapping = aes(x=ave_sig,y=ave_sig,color=factor(x=model,levels=present_models)))+
    geom_point()+
    theme(
      legend.title=element_text(size=15),
      legend.key.size = unit(0.1,'npc'),
      legend.text = element_text(size=15),
      legend.key.width = unit(0.1,'npc'))+
    scale_color_manual(
      values=present_cols,
      name="Fit Model",
      labels=model_names)
  leg <- ggdraw(get_legend(leg_plot))
  leg
  ggsave(paste(file_ext,'/figs/model_legend.png',sep=''),width = 10,height = 10,dpi=600,device=png)
  
}


#############################
####### Bastide Model #######
#############################
bas_plots<-function(file_loc){
  bas_dat<-read.csv(paste(file_loc,'/compiled_data.csv',sep=''))
  
  ##Get all columns that could be made into graphs
  y_vals<-colnames(bas_dat)[c(4:7,12:15)]
  for(y in y_vals){
    plot_single_all_mod(file_ext = file_loc,dat=bas_dat,file_name=y,
                        xvar='gamma',yvar=y)
  }

  ##Make a dummy thicc legend plot
  make_model_legend(bas_dat,file_loc)
}

bas_plots('./8_tree1/bas')
bas_plots('./8_tree2/bas')
bas_plots('./8_tree3/bas')


############################
####### Exp Model ##########
############################
exp_plots<-function(file_loc){
  exp_dat<-read.csv(paste(file_loc,'/compiled_data.csv',sep=''))
  
  ##Get all columns that could be made into graphs
  y_vals<-colnames(exp_dat[c(6:9,14:17)])
  
  gamma_iter_vals<-seq(0.5,1,0.1)
  
  iter<-'gammas'
  iter_ext<-paste(iter,"_iter",sep='')
  iter_loc<-paste(file_loc,'/figs/',iter_ext,sep = '')
  dir.create(iter_loc,showWarnings = F)
  
  for(y in y_vals){
    plot_iter_all_single(file_ext = file_ext,iter_loc=iter_ext,
                         dat=exp_dat,file_name=y,
                        xvar='pop',yvar=y,
                        iter=iter,iter_vals=gamma_iter_vals)
  }
  
  ##Make a dummy thicc legend plot
  make_model_legend(exp_dat,file_loc)
}

exp_plots('./8_tree1/exp')

exp_dat<-read.csv('./8_tree1/exp/compiled_data.csv')



file_loc<-'./8_tree1/exp'




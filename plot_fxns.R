library(ggplot2)
library(cowplot)
library(tidyr)
library(scales)



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
                                xvar,yvar,
                                ylims=NULL,
                                log_transform=F,
                                labs=F){
  col_pal <- c('#CC6677','#332288','#DDCC77','#117733','#88CCEE','#44AA99','#999933','#AA4499')
  all_models<-c('maj','bas','exp','obs','hib','min','max','uns')
  
  
  
  present_models<- all_models %in% dat$model
  present_cols<-col_pal[present_models]
  present_models<- all_models[present_models]
  

  name_dict<-list(
    pop   = expression('Population Size '*N[e]),
    ngenes= expression('Number of Genes'),
    gammas = expression('Inheritance Proportions '*gamma),
    ave_mu= expression('Average Initial Trait Value '*mu),
    sd_mu = expression('SD of Initial Trait Value '*mu),
    mae_mu= expression('MAE of Initial Trait Value '*mu),
    me_mu = expression('Mean Error of Initial Trait Value '*mu),
    ave_sig= expression('Average Evolutionary Rate '*sigma^2),
    sd_sig = expression('SD of Evolutionary Rate '*sigma^2),
    mae_sig= expression('MAE of '*sigma^2),
    me_sig = expression('Mean Error of '*sigma^2),
    rho_dist = expression('Distance Between Sorting weights '*Omega^{'*'}*' and '*Omega^{'obs'}),
    vcv_dist = expression('Distance Between Covariance Matrices '*C^{'*'}*' and '*C^{'obs'})
  )
  
  
  my_plot<-ggplot(data=dat,mapping = aes(x=!!sym(xvar),y=!!sym(yvar),color=factor(x=model,levels=present_models)))+
    geom_line(linewidth=3)+
    #geom_point(size=5)+
    theme(
      legend.position = "none",
      axis.text=element_text(size=25))+
    scale_color_manual(values=present_cols)
  if(labs){
    my_plot<-my_plot+theme(axis.title=element_text(size=30,face="bold"))+xlab(name_dict[[xvar]])+ylab(name_dict[[yvar]])
  }else{
    my_plot<-my_plot+theme(axis.title=element_blank())    
  }
   
  if(!is.null(ylims)){
    my_plot<-my_plot+ylim(ylims)
  }
  if(log_transform){
    my_plot<-my_plot+scale_x_continuous(transform = 'log',
                                        breaks = c(0.001,0.01,0.1,1,10,100,1000),
                                        labels = label_log())
  }
  
  my_plot
  ggsave(paste(file_ext,'/figs/',file_name,'.png',sep=''),width = 10,height = 10,dpi=600)
}


plot_iter_all_single<- function(file_ext,dat,
                                file_name,iter_loc,
                                xvar,yvar,
                                iter,iter_vals,
                                common_axis,
                                log_tranform=F,
                                ylims=NULL){
  
  if(is.null(ylims) && common_axis){
    y_dat<-dat[,yvar]
    ylims<-c(min(y_dat),max(y_dat))
  }
  
  
  ##Make single plots for all iterations
  for(val in iter_vals){
    dat_subset<- dat[dat[,iter]==val,] ##Only grab the data that has val 
    iter_name <- dat_subset[1,which(colnames(dat)==iter)-2] ##match the value to the respective i/j numbering in the pars.csv
    
    iter_folder<-paste(iter_loc,'/',iter,'_',iter_name,sep='')
    dir.create(paste(file_ext,'/figs/',iter_folder,sep=''),showWarnings = F)
    new_file_name<-paste('/',iter_folder,'/',file_name,sep='')
    plot_single_all_mod(file_ext,dat_subset,new_file_name,
                        xvar,yvar,
                        ylims=ylims,
                        log_tranform)
    
  }
}


make_model_legend <- function(dat,file_ext,simp=F){
  col_pal <- c('#CC6677','#332288','#DDCC77','#117733','#88CCEE','#44AA99','#999933','#AA4499')
  all_models<-c('maj','bas','exp','obs','hib','min','max','uns')
  
  model_names<-expression('Major Tree',
                          'Network',
                          'Sorting, known '*N[e],
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
      labels=model_names)+
    guides(color = guide_legend(override.aes = list(size = 15,shape=15))) 
  leg <- ggdraw(get_legend(leg_plot))
  leg
  
  if(simp){
    file_nm<-paste(file_ext,'/figs/model_legend.png',sep='')
  }else{
    file_nm<-paste(file_ext,'/figs/model_legend_simp.png',sep='')
  }
  ggsave(fil_nm,width = 10,height = 10,dpi=600,device=png)
  
}


###############################################
####### Bastide Model Plotting Function #######
###############################################
bas_plots<-function(tree_file,tree_exts,common_y=T){
  
  ##Get all columns that could be made into graphs
  bas_dat<-read.csv(paste(tree_file,'1','/bas/compiled_data.csv',sep=''))
  y_vals<-colnames(bas_dat)[c(4:7,12:15)]
  for(y in y_vals){
    print(y)
    ###read all the data if we want to find max and mins
    ylims<-NULL
    if(common_y){
      mins<-c()
      maxs<-c()
      for(ext in tree_exts){
        dat_file<-read.csv(paste(tree_file,ext,'/bas/compiled_data.csv',sep=''))
        mins<-c(mins,min(dat_file[,y]))
        maxs<-c(maxs,max(dat_file[,y]))
      }
      ylims<-c(min(mins),max(maxs))
    }
    
    for(ext in tree_exts){

      file_ext<-paste(tree_file,ext,'/bas',sep='')
      bas_dat<-read.csv(paste(file_ext,'/compiled_data.csv',sep=''))
      plot_single_all_mod(file_ext = file_ext,dat=bas_dat,file_name=y,
                          xvar='gammas',yvar=y,
                          ylims = ylims)
    }
    

  }

  ##Make a dummy thicc legend plot
  make_model_legend(bas_dat,file_ext)
}


##############################################
####### Exp Model Plotting Function ##########
##############################################
exp_plots<-function(tree_file,tree_exts,common_y=T,common_axis=F){
  
  exp_dat<-read.csv(paste(tree_file,'1','/exp/compiled_data.csv',sep=''))
  #exp_dat$pop <- log(exp_dat$pop)
  ##Get all columns that could be made into graphs
  y_vals<-colnames(exp_dat[c(6:9,14:17)])
  
  

  
  
  gamma_iter_vals<-seq(0.5,1,0.1)
  
  iter<-'gammas'
  iter_ext<-paste(iter,"_iter",sep='')

  for(y in y_vals){
    print(y)
    
    ###read all the data if we want to find max and mins
    ylims<-NULL
    if(common_y){
      mins<-c()
      maxs<-c()
      for(ext in tree_exts){
        dat_file<-read.csv(paste(tree_file,ext,'/exp/compiled_data.csv',sep=''))
        mins<-c(mins,min(dat_file[,y]))
        maxs<-c(maxs,max(dat_file[,y]))
      }
      ylims<-c(min(mins),max(maxs))
    }
    
    for(ext in tree_exts){
      file_ext<-paste(tree_file,ext,'/exp',sep='')
      dir.create(paste(file_ext,'/figs/',iter_ext,sep = ''),showWarnings = F)
      exp_dat<-read.csv(paste(file_ext,'/compiled_data.csv',sep=''))
      plot_iter_all_single(file_ext = file_ext,iter_loc=iter_ext,
                           dat=exp_dat,file_name=y,
                           xvar='pop',yvar=y,
                           iter=iter,iter_vals=gamma_iter_vals,
                           common_axis=common_axis,
                           ylims = ylims,
                           log_tranform = T)
    }
    
    
  }
  
  ##Make a dummy thicc legend plot
  make_model_legend(exp_dat,file_ext)
}

##############################################
####### Obs Model Plotting Function ##########
##############################################
obs_plots<-function(tree_file,tree_exts,
                    common_y=T,common_axis=F,
                    simple_models=F,
                    log_transform=T){
  
  obs_dat<-read.csv(paste(tree_file,'1','/obs/compiled_data.csv',sep=''))
  ##Get all columns that could be made into graphs
  y_vals<-colnames(obs_dat[c(6:9,14:17)])
  
  pop_iter_vals<-unique(obs_dat$pop)
  
  iter<-'pop'
  iter_ext<-paste(iter,"_iter",sep='')

  for(y in y_vals){
    print(y)
    f_name<-y
    if(simple_models){
      f_name<-paste(y,'_simp',sep='')
    }
    
    ###read all the data if we want to find max and mins
    ylims<-NULL
    if(common_y){
      mins<-c()
      maxs<-c()
      for(ext in tree_exts){
        dat_file<-read.csv(paste(tree_file,ext,'/obs/compiled_data.csv',sep=''))
        if(simple_models){
          dat_file<-dat_file[dat_file$model %in% c('obs','exp'),]
        }
        mins<-c(mins,min(dat_file[,y]))
        maxs<-c(maxs,max(dat_file[,y]))
      }
      ylims<-c(min(mins),max(maxs))
    }
    
    for(ext in tree_exts){
      file_ext<-paste(tree_file,ext,'/obs',sep='')
      dir.create(paste(file_ext,'/figs/',iter_ext,sep = ''),showWarnings = F)
      obs_dat<-read.csv(paste(file_ext,'/compiled_data.csv',sep=''))
      if(simple_models){
        obs_dat<-obs_dat[obs_dat$model %in% c('obs','exp'),]
      }
      plot_iter_all_single(file_ext = file_ext,iter_loc=iter_ext,
                           dat=obs_dat,file_name=f_name,
                           xvar='ngenes',yvar=y,
                           iter=iter,iter_vals=pop_iter_vals,
                           common_axis=common_axis,
                           ylims = ylims,
                           log_tranform = log_transform)
    }
  }
  
  if(simple_models){
    obs_dat<-obs_dat[obs_dat$model %in% c('obs','exp'),]
  }
  
  ##Make a dummy thicc legend plot
  make_model_legend(obs_dat,file_ext,simp=simple_models)
}



########################################
####### Matrix plotting function #######
########################################

##Note: these are called from Julia and not here 

matrix_plotting <- function(mat, mat_names,
                            file_loc,file_name){
  ##Reformat our data to x,y coords
  mat_ncol<-ncol(mat)
  mat2<-data.frame(matrix(NA,ncol = 3,nrow = mat_ncol^2))
  colnames(mat2)<-c('xs','ys','Cov')
  mat2_row<-1
  for(i in 1:mat_ncol){
    for(j in 1:mat_ncol){
      mat2[mat2_row,]<-c(mat_names[i],mat_names[j],mat[i,j])
      mat2_row<-mat2_row+1
    }
  }
  mat2$xs<-factor(mat2$xs,levels=mat_names)
  mat2$ys<-factor(mat2$ys,levels=mat_names)
  mat2$Cov<-as.numeric(mat2$Cov)
  mat2$Cov[mat2$Cov==0]<-NA
  
  
  ggplot(mat2, aes(x=xs, y=ys,fill=Cov)) +
    geom_tile() +
    geom_text(aes(label=Cov),size=0.36*25)+
    scale_fill_gradient(low="#ffd9d9",high="#FF0000",na.value='white')+
    theme(axis.title=element_blank(),
          axis.text=element_text(size=25),
          legend.title=element_text(size=15),
          legend.key.size = unit(0.1,'npc'),
          legend.text = element_text(size=15),
          legend.key.width = unit(0.1,'npc'))
  ggsave(paste(file_loc,'/',file_name,'.png',sep=''),width = 10,height = 10,dpi=600)
}




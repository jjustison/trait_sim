############################
############################
###### Make the Plots ######
############################
############################
source('./plot_fxns.R')


###############
## Bas Plots ##
###############
bas_plots('./8_tree',1:3)
#bas_dat<-read.csv('./8_tree1/bas/compiled_data.csv')
#make_model_legend(bas_dat,'./8_tree1/bas')


###############
## Exp Plots ##
###############
exp_plots('./8_tree',1:3)
exp_dat<-read.csv('./8_tree1/exp/compiled_data.csv')


###############
## Obs Plots ##
###############
obs_plots('./8_tree',1:3,log_transform = F)
obs_plots('./8_tree',1:3,simple_models=T,log_transform = F)


###############
## Hib Plots ##
###############
hib_plots('./8_tree',1:3,log_transform = F,common_axis = T)



##############################################
#### dist tables for bas and exp datasets ####
##############################################

make_table(file_loc = './8_tree1/bas/vcv_dists.csv',
           idvar = 'gamma',
           drop.cols<-1,
           idvar_subset=c(0.5,0.6,0.7,0.8,0.9,1.0),
           digits=4)

make_table(file_loc = './8_tree2/bas/vcv_dists.csv',
           idvar = 'gamma',
           drop.cols<-1,
           idvar_subset=c(0.5,0.6,0.7,0.8,0.9,1.0),
           digits=4)

make_table(file_loc = './8_tree3/bas/vcv_dists.csv',
           idvar = 'gamma',
           drop.cols<-1,
           idvar_subset=c(0.5,0.6,0.7,0.8,0.9,1.0),
           digits=4)






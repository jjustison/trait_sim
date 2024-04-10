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


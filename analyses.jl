####NOTE: You need to set the working directory to the root folder
cd("C:/Users/justison/Documents/chapt3")


##load in packages and analyses functions
include("./ch3_analysis.jl")



###########################
######### BAS FIT #########
###########################
bas_analysis("C:/Users/justison/Documents/chapt3/8_tree1/bas",82589933)
bas_analysis("C:/Users/justison/Documents/chapt3/8_tree2/bas",6543)
bas_analysis("C:/Users/justison/Documents/chapt3/8_tree3/bas",3784278)
bas_analysis("C:/Users/justison/Documents/chapt3/16_tree1/bas",8675309)


###########################
######### EXP FIT #########
###########################
exp_analysis("C:/Users/justison/Documents/chapt3/8_tree1/exp",13)
exp_analysis("C:/Users/justison/Documents/chapt3/8_tree2/exp",44444)
exp_analysis("C:/Users/justison/Documents/chapt3/8_tree3/exp",8675309)

###########################
######### OBS FIT #########
###########################
obs_analysis("C:/Users/justison/Documents/chapt3/8_tree1/obs",138937)
obs_analysis("C:/Users/justison/Documents/chapt3/8_tree2/obs",15643)
obs_analysis("C:/Users/justison/Documents/chapt3/8_tree3/obs",16)
obs_analysis("C:/Users/justison/Documents/chapt3/16_tree1/obs",572097346)


###############################
### Companion Dist Analyses ###
###############################
obs_dists("C:/Users/justison/Documents/chapt3/8_tree1/obs",138937)
obs_dists("C:/Users/justison/Documents/chapt3/8_tree2/obs",15643)
obs_dists("C:/Users/justison/Documents/chapt3/8_tree3/obs",16)
obs_dists("C:/Users/justison/Documents/chapt3/16_tree1/obs",572097346)

mat_dists("./8_tree1","bas")
mat_dists("./8_tree2","bas")
mat_dists("./8_tree3","bas")

mat_dists("./8_tree1","exp")
mat_dists("./8_tree2","exp")
mat_dists("./8_tree3","exp")



#############################
####### Plot Bas VCVs #######
#############################
plot_mat("./8_tree1","bas/figs/bas_vcv",model="bas")
plot_mat("./8_tree1","bas/figs/max_vcv",model="exp",pop=Inf)
plot_mat("./8_tree1","bas/figs/min_vcv",model="exp",pop=0.000000000001)

plot_mat("./8_tree2","bas/figs/bas_vcv",model="bas")
plot_mat("./8_tree2","bas/figs/max_vcv",model="exp",pop=Inf)
plot_mat("./8_tree2","bas/figs/min_vcv",model="exp",pop=0.000000000001)

plot_mat("./8_tree3","bas/figs/bas_vcv",model="bas")
plot_mat("./8_tree3","bas/figs/max_vcv",model="exp",pop=Inf)
plot_mat("./8_tree3","bas/figs/min_vcv",model="exp",pop=0.000000000001)


###########################
######### HIB FIT #########
###########################
## 1-51
##601-652
##1200-1250


hib_analysis("C:/Users/justison/Documents/chapt3/8_tree1/hib",58320,1200,1250)
hib_analysis("C:/Users/justison/Documents/chapt3/8_tree2/hib",509)
hib_analysis("C:/Users/justison/Documents/chapt3/8_tree3/hib",54332)
hib_analysis("C:/Users/justison/Documents/chapt3/16_tree1/hib",9532)

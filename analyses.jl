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
obs_dists("./8_tree1",138937)
obs_dists("./8_tree2",15643)
obs_dists("./8_tree3",16)
obs_dists("./16_tree1",572097346)

mat_dists("./8_tree1","bas")
mat_dists("./8_tree2","bas")
mat_dists("./8_tree3","bas")

mat_dists("./8_tree1","exp")
mat_dists("./8_tree2","exp")
mat_dists("./8_tree3","exp")


#############################
####### Plot Networks #######
#############################

h1 = readTopology("(((A:2,B:2):2,(((C:1.0,D:1.0):1.0,(E:1.0,F:1.0):1.0):1.0)#H1:1::0.5):1,(#H1:1::0.5,(G:2,H:2):2):1);")
plot(h1)
h2 = readTopology("(((A:3,((B:1,C:1):1)#H1:1::0.5):1,(#H1:1::0.5,D:3):1):1,((E:3,((F:1,G:1):1)#H2:1::0.5):1,(#H2:1::0.5,H:3):1):1);");
plot(h2)
h3 = readTopology("((((A:2,#H1:1::0.5):1,(B:1)#H1:2::0.5):1,((C:2,#H2:1::0.5):1,(D:1)#H2:2::0.5):1):1,(((E:2,#H3:1::0.5):1,(F:1)#H3:2::0.5):1,((G:2,#H4:1::0.5):1,(H:1)#H4:2::0.5):1):1);")
plot(h3)

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
## 1-50
##601-650
##1201-1250
#ind_rows = [(collect.([1:25,301:325,601:625])...)...]  ## half time. Only 25 diff genes
#ind_rows = [(collect.([1:50,601:650,1201:1250])...)...] ##full time. 50 diff genes
ind_rows1 = [(collect.([1:25,301:325,601:625])...)...]
ind_rows2 = [(collect.([1:25,301:325,601:625])...)...]
ind_rows3 = [(collect.([1:25,301:325,601:625])...)...]
#ind_rows2 = [(collect.([1:50,601:650,1201:1250])...)...]
#ind_rows3 = [(collect.([1:50,601:650,1201:1250])...)...]

hib_analysis("./8_tree1",58320,ind_rows1) 
hib_analysis("./8_tree2",509,ind_rows2) 
hib_analysis("./8_tree3",54332,ind_rows3)
hib_analysis("C:/Users/justison/Documents/chapt3/16_tree1/hib",9532)


##load in packages and analyses functions
cd("C:/Users/justison/Documents/chapt3")

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


##############################
### Companion OBS Analysis ###
##############################
obs_dists("C:/Users/justison/Documents/chapt3/8_tree1/obs",138937)
obs_dists("C:/Users/justison/Documents/chapt3/8_tree2/obs",15643)
obs_dists("C:/Users/justison/Documents/chapt3/8_tree3/obs",16)
obs_dists("C:/Users/justison/Documents/chapt3/16_tree1/obs",572097346)

###########################
######### HIB FIT #########
###########################
hib_analysis("C:/Users/justison/Documents/chapt3/8_tree1/hib",58320)
hib_analysis("C:/Users/justison/Documents/chapt3/8_tree2/hib",509)
hib_analysis("C:/Users/justison/Documents/chapt3/8_tree3/hib",54332)
hib_analysis("C:/Users/justison/Documents/chapt3/16_tree1/hib",9532)

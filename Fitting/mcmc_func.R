
#wrapper around model run

source("Models/transduction_mass-action.R")

run_mode

#data



#fix times, initial values and some paras

#mcmc function


#okay so fit based on mean but account for observation stochasticity byt resampling observations 
# from a distribution based on the confidence intervals estimated in the lab
# like when running det model but resampling from Poisson in MFIID course, but specify own distribution

#need flexible function to allow for fitting based on different curves (eg bacteria only, or bacteria+phage)
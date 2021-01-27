
mcmc_run = function(run_model, data, yinit,
                    fixed_theta = NULL, theta, theta_sd,
                    n_iter, likelihood) {
  
  #need error if model output length doesn't match data length, or better an autocorrect with notification
  times = seq(0,max(data$time),1)
  good_times = unique(data$time)
  
  curr_like = 0
  all_theta = c(fixed_theta, theta)
  obs = run_model(times, yinit, all_theta)
  obs = obs[obs$time %in% good_times,]

  #figure out which column to take from the model output based on data columns
  data_cols = colnames(data)[-1]
  model_cols = which(colnames(obs) %in% data_cols)
  
  if (length(model_cols) == 0 ) stop("Data and models do not match, check column names")
  
  # evaluate the function "likelihood" at initial parameter values
  for (j in model_cols) {
    
    added_like = likelihood(round(data[,colnames(obs)[j]]), round(obs[,j]))
    curr_like = curr_like + added_like
    
  }
  
  # initialise variables to store the current value of theta, the
  # vector of samples, and the number of accepted proposals
  samples = theta
  accepted = 0
  
  # repeat n.iterations times:
  
  for (i in 1:n_iter) {
    
    # - draw a new theta from the (Gaussian) proposal distribution
    #   with standard deviation sd.
    # new_theta = rnorm(length(theta), theta , theta_sd)
    
    #trick here is to force a lower band on growth rate variation
      #otherwise will get stuck going up because correlated
    
    #COULD also try fitting LOG instead of base
    
    new_theta = runif(length(theta), theta - theta_sd, theta + theta_sd)
    
    names(new_theta) = names(theta)
    
    # - evaluate the function "target" at the proposed theta
    all_theta = c(fixed_theta, new_theta)
    obs = run_model(times, yinit, all_theta)
    obs = obs[obs$time %in% good_times,]
    
    # evaluate the function "likelihood" at initial parameter values
    # exclude first column containing time
    new_like = 0
    for (j in model_cols) {
      
      added_like = likelihood(round(data[,colnames(obs)[j]]), round(obs[,j]))
      new_like = new_like + added_like
      
    }
    
    # - calculate the Metropolis-Hastings ratio
    log_accept = new_like - curr_like
    
    # - draw a random number between 0 and 1
    r = runif(1)
    
    # - accept or reject by comparing the random number to the
    #   Metropolis-Hastings ratio (acceptance probability); if accept,
    #   change the current value of theta to the proposed theta,
    #   update the current value of the target and keep track of the
    #   number of accepted proposals
    
    if (r < exp(log_accept)) {
      
      # add the current theta to the vector of samples
      theta = new_theta
      curr_like = new_like
      
      accepted = accepted + 1
      
    }
    
    samples = rbind(samples, theta, deparse.level = 0)
    
    if (i %% 100 == 0) {
      cat("Iteration:", i, "out of", n_iter,
          ", Chain:", theta,
          ", Acceptance rate:", accepted/i, "\n")
    }
  }
  
  # return the trace of the chain (i.e., the vector of samples)
  return(samples)
  
}


#parameter values:
fixed_theta = c(tre = 0,      #ery gene transduction probability
                trt = 0,
                Nmax = 1e9)     #maximum bacteria population size

parms = c(Ge = 5,         #growth parameter for ery resistant
          Gt = 2,
          Get = 2)       #growth parameter for tet resistant

#initial values:

yinit_mass = c(Be = 11000,    #resistant to ery
               Bt = 9000,    #resistant to tet
               Bet = 10000)    #resistant to ery and tet

trace = mcmc_run(run_mass, results, yinit_mass,
                 fixed_theta, parms, c(0.5,0.5,0.5),
                 5000, likelihood_func)

# results = run_mass(seq(0,24,1), c(Be = 11000, Bt = 9000, Bet = 10000),
#                    c(Ge = 1.4, Gt = 1.6, Get = 1.3, tre = 0, trt = 0, Nmax = 1e9))[1:10,]

par(mfrow = c(3,1))
plot(trace[,"Ge"], type = "l")
abline(mean(trace[,"Ge"]), 0, col = "red")
plot(trace[,"Gt"], type = "l")
abline(mean(trace[,"Gt"]), 0, col = "red")
plot(trace[,"Get"], type = "l")
abline(mean(trace[,"Get"]), 0, col = "red")
par(mfrow = c(1,1))

fitted_parms = c(fixed_theta, Ge = mean(trace[,"Ge"]), 
                 Gt = mean(trace[,"Gt"]), Get = mean(trace[,"Get"]))
fitted_parms

res_fit = run_mass(seq(0,24,1), yinit_mass, fitted_parms)

plot(log10(results$Be), type = "l", ylim = c(1, max(log10(results))))
lines(log10(results$Bt))
lines(log10(results$Bet))
lines(log10(res_fit$Be), col = "red")
lines(log10(res_fit$Bt), col = "red")
lines(log10(res_fit$Bet), col = "red")

#okay so fit based on mean with Poisson but account for observation stochasticity by resampling observations 
# from a distribution based on the confidence intervals estimated in the lab
# like when running det model but resampling from Poisson in MFIID course, but specify own distribution

#for resampling, sample from normal with mean and SD taken from lab estimates

#need flexible function to allow for fitting based on different curves (eg bacteria only, or bacteria+phage)
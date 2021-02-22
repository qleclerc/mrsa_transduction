model_name <- "transduction model"
model_state.names <- c("Be","Bt","Bet","Pl", "Pe", "Pt")
model_theta.names <- c("beta", "beta2", "L", "gamma", "alpha", "tau")

model_simulateDeterministic <- function(theta,init.state,times) {
  
}


## function to compute log-prior
model_prior <- function(theta, log = FALSE) {
  
  log.prior.L <- dunif(theta[["L"]], min = 1, max = 1000, log = TRUE)
  log.prior.beta <- dunif(theta[["beta"]], min = 1, max = 1e20, log = TRUE)
  log.prior.beta2 <- dunif(theta[["beta2"]], min = 0.001, max = 1, log = TRUE)
  log.prior.gamma <- dunif(theta[["gamma"]], min = 1, max = 1e20, log = TRUE)
  log.prior.alpha <- dunif(theta[["alpha"]], min = 1, max = 1e20, log = TRUE)
  
  log.prior.tau <- dunif(theta[["tau"]], min = 0.001, max = 1, log = TRUE)
  
  log.sum <- log.prior.L + log.prior.beta + log.prior.beta2 + log.prior.gamma + log.prior.alpha + log.prior.tau
  
  return(ifelse(log, log.sum, exp(log.sum)))
}

## function to compute the likelihood of one data point
model_pointLike <- function(data.point, model.point, theta, log = FALSE){
  
  # dpoisBe = dpois(x = data.point[["Be"]],
  #                 lambda = model.point[["Be"]],
  #                 log = log)
  # 
  # dpoisBt = dpois(x = data.point[["Bt"]],
  #                 lambda = model.point[["Bt"]],
  #                 log = log)
  dpoisBe = dpoisBt = dpoisPe = dpoisPt = 0
  
  dpoisBet = dpois(x = data.point[["Bet"]],
                   lambda = model.point[["Bet"]],
                   log = log)
  if(is.infinite(dpoisBet)) dpoisBet = -1e7
  
  
  dpoisPl = dpois(x = data.point[["P"]],
                 lambda = model.point[["Pl"]],
                 log = log)
  
  ## the prevalence is observed through a Poisson process
  return(sum(dpoisBe, dpoisBt, 10000*dpoisBet, dpoisPl, dpoisPe, dpoisPt))
}

## function to generate observation from a model simulation
# phagebac_genObsPoint <- function(model.point, theta){
#   
#   ## the prevalence is observed through a Poisson process
#   obs.point <- rpois(n = 1, lambda = model.point[["I"]])
#   
#   return(c(obs = obs.point))
# }

## create deterministic SIR fitmodel
model <- fitR::fitmodel(
  name = model_name,
  state.names = model_state.names,
  theta.names = model_theta.names,
  simulate = model_simulateDeterministic,
  dprior = model_prior,
  dPointObs = model_pointLike)

saveRDS(model, here::here("Model", "transduction_model.rds"))

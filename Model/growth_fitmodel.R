model_name <- "growth model"
model_state.names <- c("Be","Bt","Bet")
model_theta.names <- c("mu_e", "mu_t", "mu_et", "Nmax")

model_simulateDeterministic <- function(theta, init.state, times) {
  
  model_ode <- function(time, state, parameters) {
    
    ## parameters
    mu_e = parameters[["mu_e"]]
    mu_t = parameters[["mu_t"]] 
    mu_et = parameters[["mu_et"]]
    Nmax = parameters[["Nmax"]] 
    
    ## states
    Be = state[["Be"]]
    Bt = state[["Bt"]]
    Bet = state[["Bet"]]
    
    N = Be + Bt + Bet
    
    link = (1 - N/Nmax)
    
    #no link
    dBe = mu_e * link * Be 
    dBt = mu_t * link * Bt
    dBet = mu_et * link * Bet
    
    return(list(c(dBe, dBt, dBet)))
    
  }
  
  trajectory <- data.frame(ode(y = init.state,
                               times = times,
                               func = model_ode,
                               parms = theta,
                               method = "ode45"))
  
  return(trajectory)
  
}



## function to compute log-prior
model_prior <- function(theta, log = TRUE) {
  
  log.prior.mu_e <- dunif(theta[["mu_e"]], min = 0.1, max = 2, log = TRUE)
  log.prior.mu_t <- dunif(theta[["mu_t"]], min = 0.1, max = 2, log = TRUE)
  log.prior.mu_et <- dunif(theta[["mu_et"]], min = 0.1, max = 2, log = TRUE)
  log.prior.Nmax <- dunif(theta[["Nmax"]], min = 1, max = 1e12, log = TRUE)
  
  log.sum <- log.prior.mu_e + log.prior.mu_t + log.prior.mu_et + log.prior.Nmax
  
  return(ifelse(log, log.sum, exp(log.sum)))
}

## function to compute the likelihood of one data point
model_pointLike <- function(data.point, model.point, theta, log = FALSE){
  
  dpoisBe = dpois(x = round(data.point[["Be"]]/(10^(max(nchar(as.character(round(model.point[["Be"]]))),2)-2))),
                  lambda = model.point[["Be"]]/(10^(max(nchar(as.character(round(model.point[["Be"]]))),2)-2)),
                  log = log)
  
  dpoisBt = dpois(x = round(data.point[["Bt"]]/(10^(max(nchar(as.character(round(model.point[["Bt"]]))),2)-2))),
                  lambda = model.point[["Bt"]]/(10^(max(nchar(as.character(round(model.point[["Bt"]]))),2)-2)),
                  log = log)
  
  dpoisBet = dpois(x = round(data.point[["Bet"]]/(10^(max(nchar(as.character(round(model.point[["Bet"]]))),2)-2))),
                   lambda = model.point[["Bet"]]/(10^(max(nchar(as.character(round(model.point[["Bet"]]))),2)-2)),
                   log = log)
  
  ## the prevalence is observed through a Poisson process
  return(sum(dpoisBe, dpoisBt, dpoisBet))
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

saveRDS(model, here::here("Model", "growth_model.rds"))

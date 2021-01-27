mass_name <- "Mass-action model"
mass_state.names <- c("Be","Bt","Bet")
mass_theta.names <- c("Ge","Gt", "Get", "Nmax")

mass_simulateDeterministic <- function(theta,init.state,times) {
  
  mass_ode <- function(time, state, parameters) {
    
    ## parameters
    Ge = parameters[["Ge"]] 
    Gt = parameters[["Gt"]] 
    Get = parameters[["Get"]] 
    Nmax = parameters[["Nmax"]] 
    
    ## states
    Be = state[["Be"]]
    Bt = state[["Bt"]]
    Bet = state[["Bet"]]
    
    N = Be + Bt + Bet
    
    dBe = Ge*(1 - N/Nmax)*Be
    dBt = Gt*(1 - N/Nmax)*Bt
    dBet = Get*(1 - N/Nmax)*Bet
    
    return(list(c(dBe, dBt, dBet)))
  }
  
  trajectory <- data.frame(ode(y = init.state,
                               times = times,
                               func = mass_ode,
                               parms = theta,
                               method = "ode45"))
  
  return(trajectory)
}

## function to compute log-prior
mass_prior <- function(theta, log = FALSE) {
  
  log.prior.Nmax <- dunif(theta[["Nmax"]], min = 1e6, max = 1e11, log = TRUE)
  log.prior.Ge <- dunif(theta[["Ge"]], min = 0, max = 3, log = TRUE)
  log.prior.Gt <- dunif(theta[["Gt"]], min = 0, max = 3, log = TRUE)
  log.prior.Get <- dunif(theta[["Get"]], min = 0, max = 3, log = TRUE)
  
  log.sum <- log.prior.Nmax + log.prior.Ge + log.prior.Gt + log.prior.Get
  
  return(ifelse(log, log.sum, exp(log.sum)))
}

## function to compute the likelihood of one data point
mass_pointLike <- function(data.point, model.point, theta, log = FALSE){
  
  dpoisBe = dpois(x = data.point[["Be"]],
                  lambda = model.point[["Be"]],
                  log = log)
  
  dpoisBt = dpois(x = data.point[["Bt"]],
                  lambda = model.point[["Bt"]],
                  log = log)
  
  dpoisBet = dpois(x = data.point[["Bet"]],
                   lambda = model.point[["Bet"]],
                   log = log)
  
  ## the prevalence is observed through a Poisson process
  return(sum(dpoisBe, dpoisBt, dpoisBet))
}

## function to generate observation from a model simulation
# mass_genObsPoint <- function(model.point, theta){
#   
#   ## the prevalence is observed through a Poisson process
#   obs.point <- rpois(n = 1, lambda = model.point[["I"]])
#   
#   return(c(obs = obs.point))
# }

## create deterministic SIR fitmodel
mass_model <- fitR::fitmodel(
  name = mass_name,
  state.names = mass_state.names,
  theta.names = mass_theta.names,
  simulate = mass_simulateDeterministic,
  dprior = mass_prior,
  dPointObs = mass_pointLike)

saveRDS(mass_model, here::here("Fitting", "mass_model.rds"))
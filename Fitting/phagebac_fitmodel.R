phagebac_name <- "phagebac model"
phagebac_state.names <- c("Be","Bt","Bet","P")
phagebac_theta.names <- c("beta", "L", "phi", "tau", "tr")

phagebac_simulateDeterministic <- function(theta,init.state,times) {
  
  phagebac_ode <- function(time, state, parameters) {
    
    ## parameters
    Ge = 1.678177 
    Gt = 1.605017 
    Get = 0 #1.562914
    Nmax = 2.140317e+09 
    #ll for this set is -4.026931e+08 
    
    beta = 1/parameters[["beta"]]
    L = parameters[["L"]]
    phi = 1/parameters[["phi"]]
    
    ## states
    Be = state[["Be"]]
    Bt = state[["Bt"]]
    Bet = state[["Bet"]]
    P = state[["P"]]
    
    N = Be + Bt + Bet
    
    
    #no link
    dBe = Ge*(1 - N/Nmax)*Be - beta * Be * P
    dBt = Gt*(1 - N/Nmax)*Bt - beta * Bt * P
    dBet = Get*(1 - N/Nmax)*Bet - beta * Bet * P

    dP = beta * N * P * L - phi * P
    
    #link beta
    dBe = Ge*(1 - N/Nmax)*Be - beta*(1-N/Nmax) * Be * P
    dBt = Gt*(1 - N/Nmax)*Bt - beta*(1-N/Nmax) * Bt * P
    dBet = Get*(1 - N/Nmax)*Bet - beta*(1-N/Nmax) * Bet * P

    dP = beta*(1-N/Nmax) * N * P * L - phi * P
    
    #link L
    dBe = Ge*(1 - N/Nmax)*Be - beta * Be * P
    dBt = Gt*(1 - N/Nmax)*Bt - beta * Bt * P
    dBet = Get*(1 - N/Nmax)*Bet - beta * Bet * P

    dP = beta * N * P * L*(1-N/Nmax)
    
    #link both
    dBe = Ge*(1 - N/Nmax)*Be - beta*(1-N/Nmax) * Be * P
    dBt = Gt*(1 - N/Nmax)*Bt - beta*(1-N/Nmax) * Bt * P
    dBet = Get*(1 - N/Nmax)*Bet - beta*(1-N/Nmax) * Bet * P

    dP = beta*(1-N/Nmax) * N * P * L*(1-N/Nmax) - phi * P

    # #frequentist
    beta = beta*(1-N/Nmax)
    # P_inf = (1-exp(-beta*N))*P
    # B_inf = (1-exp(-P_inf/N))*N

    B_inf = (1-exp(-P*beta/N))*N


    dBe = Ge*(1 - N/Nmax)*Be - B_inf* (Be/N)
    dBt = Gt*(1 - N/Nmax)*Bt - B_inf* (Bt/N)
    dBet = Get*(1 - N/Nmax)*Bet - B_inf* (Bet/N)

    dP = B_inf * L * (1-N/Nmax) - phi*P
    
    #frequentist volkova
    # B_inf = (1-exp(-(P/N) * beta))
    # 
    # dBe = Ge*(1 - N/Nmax)*Be - Be * B_inf
    # dBt = Gt*(1 - N/Nmax)*Bt - Bt * B_inf
    # dBet = Get*(1 - N/Nmax)*Bet- Bet * B_inf
    # 
    # dP = B_inf * L
    
    
    
    
    return(list(c(dBe, dBt, dBet, dP)))
  }
  
  trajectory <- data.frame(ode(y = init.state,
                               times = times,
                               func = phagebac_ode,
                               parms = theta,
                               method = "ode45"))
  
  
  return(trajectory)
}

## function to compute log-prior
phagebac_prior <- function(theta, log = FALSE) {
  
  log.prior.L <- dunif(theta[["L"]], min = 1, max = 1000, log = TRUE)
  log.prior.beta <- dunif(theta[["beta"]], min = 0, max = 1e20, log = TRUE)
  log.prior.phi <- dunif(theta[["phi"]], min = 0, max = 1e20, log = TRUE)
  log.prior.tr <- dunif(theta[["tr"]], min = 0, max = 1e20, log = TRUE)
  
  #log.prior.tau <- dunif(theta[["tau"]], min = 0.001, max = 2, log = TRUE)
  
  log.sum <- log.prior.L + log.prior.beta + log.prior.phi + log.prior.tr#+ log.prior.tau
  
  return(ifelse(log, log.sum, exp(log.sum)))
}

## function to compute the likelihood of one data point
phagebac_pointLike <- function(data.point, model.point, theta, log = FALSE){
  
  # dpoisBe = dpois(x = data.point[["Be"]],
  #                 lambda = model.point[["Be"]],
  #                 log = log)
  # 
  # dpoisBt = dpois(x = data.point[["Bt"]],
  #                 lambda = model.point[["Bt"]],
  #                 log = log)
  dpoisBe = dpoisBt = 0
  
  dpoisBet = dpois(x = data.point[["Bet"]],
                   lambda = model.point[["Bet"]],
                   log = log)
  if(is.infinite(dpoisBet)) dpoisBet = 0

    
  dpoisP = dpois(x = data.point[["P"]],
                 lambda = model.point[["P"]],
                 log = log)
  
  ## the prevalence is observed through a Poisson process
  return(sum(dpoisBe, dpoisBt, dpoisBet, dpoisP))
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
phagebac_model <- fitR::fitmodel(
  name = phagebac_name,
  state.names = phagebac_state.names,
  theta.names = phagebac_theta.names,
  simulate = phagebac_simulateDeterministic,
  dprior = phagebac_prior,
  dPointObs = phagebac_pointLike)

saveRDS(phagebac_model, here::here("Fitting", "phagebac_model.rds"))


choose_model = function(model,
                        frequentist = FALSE,
                        second_beta = FALSE,
                        delay = FALSE, fixed_delay = NA,
                        decay = FALSE,
                        link_beta = FALSE, link_L = FALSE, link_delay = FALSE,
                        transduction = FALSE){
  
  if(delay){
    
    model_simulateDeterministic <- function(theta, init.state, times) {
      
      model_dde <- function(time, state, parameters) {
        
        ## parameters
        mu_e = 1.678177 
        mu_t = 1.605017 
        mu_et = 1.562914
        Nmax = 2.140317e+09 
        #ll for this set is -4.026931e+08 
        
        beta = 1/parameters[["beta"]]
        beta2 = ifelse(second_beta, parameters[["beta2"]], 1)
        L = parameters[["L"]]
        gamma = ifelse(decay, 1/parameters[["gamma"]], 0)
        alpha = ifelse(transduction, 1/parameters[["alpha"]], 0)
        tau = ifelse(is.na(fixed_delay), parameters[["tau"]], fixed_delay)
        
        ## states
        Be = state[["Be"]]
        Bt = state[["Bt"]]
        Bet = state[["Bet"]]
        Pl = state[["Pl"]]
        Pe = state[["Pe"]]
        Pt = state[["Pt"]]
        
        N = Be + Bt + Bet
        
        if(time <= tau){
          Be_past = 0
          Bt_past = 0
          Bet_past = 0
          Pl_past = 0
          Pe_past = 0
          Pt_past = 0
          N_past = 1
        } else {
          Be_past = lagvalue(time - tau, 1)
          Bt_past = lagvalue(time - tau, 2)
          Bet_past = lagvalue(time - tau, 3)
          Pl_past = lagvalue(time - tau, 4)
          Pe_past = lagvalue(time - tau, 5)
          Pt_past = lagvalue(time - tau, 6)
          N_past = Be_past + Bt_past + Bet_past
        }
        
        
        link = (1 - N/Nmax)
        link_past = (1 - N_past/Nmax)
        
        
        if(link_beta){
          beta = beta * link
          beta_past = beta * link_past
        } else beta_past = beta
        
        beta2_past = beta2
        
        if(link_L) L = L * link + 1
        
        if(link_delay) tau = tau * (N/Nmax)
        
        if(frequentist){
          
          lambda = (1 - exp(-beta * N))
          phi_Pl = (1 - exp(-lambda * Pl*beta2/N)) * N
          phi_Pe = (1 - exp(-lambda * Pe*beta2/N)) * N
          phi_Pt = (1 - exp(-lambda * Pt*beta2/N)) * N
          
          lambda_past = (1 - exp(-beta_past * N_past))
          phi_Pl_past = (1 - exp(-lambda_past * Pl_past*beta2_past/N_past)) * N_past
          phi_Pe_past = (1 - exp(-lambda_past * Pe_past*beta2_past/N_past)) * N_past
          phi_Pt_past = (1 - exp(-lambda_past * Pt_past*beta2_past/N_past)) * N_past
          
        } else {
          
          lambda = beta * N
          phi_Pl = lambda * Pl
          phi_Pe = lambda * Pe
          phi_Pt = lambda * Pt
          
          lambda_past = beta_past * N_past
          phi_Pl_past = lambda_past * Pl_past
          phi_Pe_past = lambda_past * Pe_past
          phi_Pt_past = lambda_past * Pt_past
          
        }
        
        #no link
        dBe = mu_e * link * (Be-(phi_Pl + phi_Pt) * Be/N) - (phi_Pl + phi_Pt) * Be/N
        dBt = mu_t * link * (Bt-(phi_Pl + phi_Pe) * Bt/N) - (phi_Pl + phi_Pe) * Bt/N
        dBet = mu_et * link * (Bet-phi_Pl*Bet/N) - phi_Pl * Bet/N +
          phi_Pe * Bt/N + phi_Pt * Be/N
        
        dPl = phi_Pl_past * L * (1 - alpha*(Be_past+Bt_past+2*Bet_past)/N_past) -
              lambda * Pl - gamma * Pl
        dPe = phi_Pl_past * L * alpha * (Be_past + Bet_past)/N_past - lambda * Pe - gamma * Pe
        dPt = phi_Pl_past * L * alpha * (Bt_past + Bet_past)/N_past - lambda * Pt - gamma * Pt
        
        return(list(c(dBe, dBt, dBet, dPl, dPe, dPt)))
        
      }
      
      trajectory <- data.frame(dede(y = init.state,
                                   times = times,
                                   func = model_dde,
                                   parms = theta))
      
      return(trajectory)
      
    }
    
  } else {
    
    model_simulateDeterministic <- function(theta, init.state, times) {
      
      model_ode <- function(time, state, parameters) {
        
        ## parameters
        mu_e = 1.678177 
        mu_t = 1.605017 
        mu_et = 1.562914
        Nmax = 2.140317e+09 
        #ll for this set is -4.026931e+08 
        
        beta = 1/parameters[["beta"]]
        L = parameters[["L"]]
        gamma = ifelse(decay, 1/parameters[["gamma"]], 0)
        alpha = ifelse(transduction, 1/parameters[["alpha"]], 0)

        ## states
        Be = state[["Be"]]
        Bt = state[["Bt"]]
        Bet = state[["Bet"]]
        Pl = state[["Pl"]]
        Pe = state[["Pe"]]
        Pt = state[["Pt"]]
        
        N = Be + Bt + Bet
        
        link = (1 - N/Nmax)
        
        if(link_beta) beta = beta * link
        if(link_L) L = L * link + 1
        
        if(frequentist){
          
          lambda = (1 - exp(-beta * N))
          phi_Pl = (1 - exp(-lambda * Pl/N)) * N
          phi_Pe = (1 - exp(-lambda * Pe/N)) * N
          phi_Pt = (1 - exp(-lambda * Pt/N)) * N
          
        } else {
          
          lambda = beta * N
          phi_Pl = lambda * Pl
          phi_Pe = lambda * Pe
          phi_Pt = lambda * Pt
          
        }
        
        #no link
        dBe = mu_e * link * Be - (phi_Pl + phi_Pt) * Be/N
        dBt = mu_t * link * Bt - (phi_Pl + phi_Pe) * Bt/N
        dBet = mu_et * link * Bet - phi_Pl * Bet/N +
          phi_Pe * Bt/N + phi_Pt * Be/N
        
        dPl = phi_Pl * L * (1 - alpha) - lambda * Pl - gamma * Pl
        dPe = phi_Pl * L * alpha * (Be + Bet)/N - lambda * Pe - gamma * Pe
        dPt = phi_Pl * L * alpha * (Bt + Bet)/N - lambda * Pt - gamma * Pt
        
        return(list(c(dBe, dBt, dBet, dPl, dPe, dPt)))
        
      }
      
      trajectory <- data.frame(ode(y = init.state,
                                   times = times,
                                   func = model_ode,
                                   parms = theta,
                                   method = "ode45"))
      
      return(trajectory)
      
    }
  }
  
  model$simulate = model_simulateDeterministic
  
  model
  
}


run_mcmc = function(model, lab_data,
                    init.theta = c(beta = 8e9, L = 20, gamma = 30000, alpha = 1e5),
                    proposal.sd = init.theta/50000,
                    n.iterations = 10000,
                    adapt.size.start = 200000,
                    adapt.size.cooling = 0.999,
                    adapt.shape.start = 200000){
  
  
  
  target_function = function(theta){
    
    my_init.state <- c(Be = lab_data$Be[1], Bt = lab_data$Bt[1], Bet = 0,
                       Pl = lab_data$P[1], Pt = 0, Pe = 0)
    
    return(dLogPosterior(fitmodel = model, theta = theta, init.state = my_init.state, 
                         data = lab_data, margLogLike = dTrajObs, log = TRUE))
    
    
  }
  
  
  mcmc_fit = mcmcMH(target = target_function,
                    init.theta = init.theta,
                    proposal.sd = proposal.sd, 
                    n.iterations = n.iterations,
                    adapt.size.start = adapt.size.start,
                    adapt.size.cooling = adapt.size.cooling,
                    adapt.shape.start = adapt.shape.start)
  
  mcmc_fit
  
}



multi_run = function(model, theta_trace, init.state, times = seq(0, 24, 1), nruns = 5000){
  
  theta_trace = theta_trace[,-ncol(theta_trace)]
  
  summary_runs = list()
  index = 1
  for (i in names(init.state)) {
    summary_runs[[index]] = matrix(0, length(times), nruns)
    index = index + 1
  }
  names(summary_runs) = names(init.state)
  
  
  for(i in 1:nruns){
    
    theta = apply(theta_trace, 2, FUN = function(x) sample(x, 1))
    
    traj = model$simulate(theta, init.state, times)
    
    for (name in names(init.state)) {
      summary_runs[[name]][,i] = traj[,name]
    }
    
  }
  
  
  #combine all results into a single dataframe with mean and sd
  summary_results = data.frame(time = times)
  
  for (name in names(init.state)) {
    summary_results = cbind(summary_results, 
                            rowMeans(summary_runs[[name]]), 
                            apply(summary_runs[[name]], 1, sd))
  }
  
  summary_colnames = c()
  for (name in names(init.state)) {
    summary_colnames = c(summary_colnames,
                         name,
                         paste0(name, "_sd"))
  }
  
  colnames(summary_results) = c("time", summary_colnames)
  
  summary_results
  
}

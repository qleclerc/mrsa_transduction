
## Prep ##########

library(fitR)
library(dplyr)
library(reshape2)
library(ggplot2)
library(cowplot)
library(scales)

#all_trace = readRDS(here::here("Fitting", "fitted_params.rds"))

phagebac_model = readRDS(here::here("Fitting", "phagebac_model.rds"))

lab_data_transM = read.csv(here::here("Lab", "Transduction", "summary_10_4.csv")) %>%
  select(Time, Bacteria, Mean, se) %>%
  #mutate(se = se*sqrt(3)) %>%
  filter(Bacteria != "Total")
lab_data_trans5M = read.csv(here::here("Lab", "Transduction", "summary_10_5.csv")) %>%
  select(Time, Bacteria, Mean, se)%>%
  #mutate(se = se*sqrt(3)) %>%
  filter(Bacteria != "Total")
lab_data_trans3M = read.csv(here::here("Lab", "Transduction", "summary_10_3.csv")) %>%
  select(Time, Bacteria, Mean, se) %>%
  #mutate(se = se*sqrt(3)) %>%
  filter(Bacteria != "Total")


lab_data_trans = read.csv(here::here("Lab", "Transduction", "summary_10_4.csv")) %>%
  select(Time, Bacteria, Mean) %>%
  dcast(Time~Bacteria) %>%
  select(-Total) %>%
  rename(time = Time, Bet = DRP, Be = EryR, Bt = TetR) %>%
  mutate(Be = round(Be),
         Bt = round(Bt),
         Bet = round(Bet),
         P = round(P)) 

lab_data_trans5 = read.csv(here::here("Lab", "Transduction", "summary_10_5.csv")) %>%
  select(Time, Bacteria, Mean) %>%
  dcast(Time~Bacteria) %>%
  select(-Total) %>%
  rename(time = Time, Bet = DRP, Be = EryR, Bt = TetR) %>%
  mutate(Be = round(Be),
         Bt = round(Bt),
         Bet = round(Bet),
         P = round(P)) 

lab_data_trans3 = read.csv(here::here("Lab", "Transduction", "summary_10_3.csv")) %>%
  select(Time, Bacteria, Mean) %>%
  dcast(Time~Bacteria) %>%
  select(-Total) %>%
  rename(time = Time, Bet = DRP, Be = EryR, Bt = TetR) %>%
  mutate(Be = round(Be),
         Bt = round(Bt),
         Bet = round(Bet),
         P = round(P))


choose_model = function(model = "mass"){
  #first set
  if(model == "mass"){
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
        
        dP = beta * N * P * L
        
        return(list(c(dBe, dBt, dBet, dP)))
      }
      trajectory <- data.frame(ode(y = init.state,
                                   times = times,
                                   func = phagebac_ode,
                                   parms = theta,
                                   method = "ode45"))
      return(trajectory)
    }
  }
  
  if(model == "mass_decay"){
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
        
        return(list(c(dBe, dBt, dBet, dP)))
      }
      trajectory <- data.frame(ode(y = init.state,
                                   times = times,
                                   func = phagebac_ode,
                                   parms = theta,
                                   method = "ode45"))
      return(trajectory)
    }
  }
  
  if(model == "dde_mass_decay"){
    phagebac_simulateDeterministic <- function(theta,init.state,times) {
      phagebac_dde <- function(time, state, parameters) {
        
        ## parameters
        Ge = 1.678177 
        Gt = 1.605017 
        Get = 0 #1.562914
        Nmax = 2.140317e+09 
        #ll for this set is -4.026931e+08 
        
        beta = 1/parameters[["beta"]]
        L = parameters[["L"]]
        phi = 1/parameters[["phi"]]
        tau = 0.5
        
        ## states
        Be = state[["Be"]]
        Bt = state[["Bt"]]
        Bet = state[["Bet"]]
        P = state[["P"]]
        
        N = Be + Bt + Bet
        
        if(time < tau){
          Be_past = 0
          Bt_past = 0
          Bet_past = 0
          P_past = 0
        } else {
          Be_past = lagvalue(time - tau, 1)
          Bt_past = lagvalue(time - tau, 2)
          Bet_past = lagvalue(time - tau, 3)
          P_past = lagvalue(time - tau, 4)
        }
        
        N_past = Be_past + Bt_past + Bet_past
        
        #link both
        dBe = Ge*(1 - N/Nmax)*Be - beta * Be * P
        dBt = Gt*(1 - N/Nmax)*Bt - beta * Bt * P
        dBet = Get*(1 - N/Nmax)*Bet - beta * Bet * P
        
        dP =  - beta * P * N + N_past * P_past * beta * L - phi * P
        
        return(list(c(dBe, dBt, dBet, dP)))
      }
      trajectory <- data.frame(dede(y = init.state,
                                    times = times,
                                    func = phagebac_dde,
                                    parms = theta))
      return(trajectory)
    }
  }
  
  #second set
  if(model == "mass_link_beta"){
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
        
        #link beta
        dBe = Ge*(1 - N/Nmax)*Be - beta*(1-N/Nmax) * Be * P
        dBt = Gt*(1 - N/Nmax)*Bt - beta*(1-N/Nmax) * Bt * P
        dBet = Get*(1 - N/Nmax)*Bet - beta*(1-N/Nmax) * Bet * P
        
        dP = beta*(1-N/Nmax) * N * P * L
        
        return(list(c(dBe, dBt, dBet, dP)))
      }
      trajectory <- data.frame(ode(y = init.state,
                                   times = times,
                                   func = phagebac_ode,
                                   parms = theta,
                                   method = "ode45"))
      return(trajectory)
    }
  }
  
  if(model == "mass_link_L"){
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
        
        #link L
        dBe = Ge*(1 - N/Nmax)*Be - beta * Be * P
        dBt = Gt*(1 - N/Nmax)*Bt - beta * Bt * P
        dBet = Get*(1 - N/Nmax)*Bet - beta * Bet * P
        
        dP = beta * N * P * L*(1-N/Nmax)
        
        return(list(c(dBe, dBt, dBet, dP)))
      }
      trajectory <- data.frame(ode(y = init.state,
                                   times = times,
                                   func = phagebac_ode,
                                   parms = theta,
                                   method = "ode45"))
      return(trajectory)
    }
  }
  
  if(model == "mass_link_both"){
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
        
        #link both
        dBe = Ge*(1 - N/Nmax)*Be - beta*(1-N/Nmax) * Be * P
        dBt = Gt*(1 - N/Nmax)*Bt - beta*(1-N/Nmax) * Bt * P
        dBet = Get*(1 - N/Nmax)*Bet - beta*(1-N/Nmax) * Bet * P
        
        dP = beta*(1-N/Nmax) * N * P * L*(1-N/Nmax)
        
        return(list(c(dBe, dBt, dBet, dP)))
      }
      trajectory <- data.frame(ode(y = init.state,
                                   times = times,
                                   func = phagebac_ode,
                                   parms = theta,
                                   method = "ode45"))
      return(trajectory)
    }
  }
  
  if(model == "mass_decay_link_beta"){
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
        
        #link beta
        dBe = Ge*(1 - N/Nmax)*Be - beta*(1-N/Nmax) * Be * P
        dBt = Gt*(1 - N/Nmax)*Bt - beta*(1-N/Nmax) * Bt * P
        dBet = Get*(1 - N/Nmax)*Bet - beta*(1-N/Nmax) * Bet * P
        
        dP = beta*(1-N/Nmax) * N * P * L - phi * P
        
        return(list(c(dBe, dBt, dBet, dP)))
      }
      trajectory <- data.frame(ode(y = init.state,
                                   times = times,
                                   func = phagebac_ode,
                                   parms = theta,
                                   method = "ode45"))
      return(trajectory)
    }
  }
  
  if(model == "mass_decay_link_L"){
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
        
        #link L
        dBe = Ge*(1 - N/Nmax)*Be - beta * Be * P
        dBt = Gt*(1 - N/Nmax)*Bt - beta * Bt * P
        dBet = Get*(1 - N/Nmax)*Bet - beta * Bet * P
        
        dP = beta * N * P * L*(1-N/Nmax) - phi * P
        
        return(list(c(dBe, dBt, dBet, dP)))
      }
      trajectory <- data.frame(ode(y = init.state,
                                   times = times,
                                   func = phagebac_ode,
                                   parms = theta,
                                   method = "ode45"))
      return(trajectory)
    }
  }
  
  if(model == "mass_decay_link_both"){
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
        
        #link both
        dBe = Ge*(1 - N/Nmax)*Be - beta*(1-N/Nmax) * Be * P
        dBt = Gt*(1 - N/Nmax)*Bt - beta*(1-N/Nmax) * Bt * P
        dBet = Get*(1 - N/Nmax)*Bet - beta*(1-N/Nmax) * Bet * P
        
        dP = beta*(1-N/Nmax) * N * P * L*(1-N/Nmax) - phi * P
        
        return(list(c(dBe, dBt, dBet, dP)))
      }
      trajectory <- data.frame(ode(y = init.state,
                                   times = times,
                                   func = phagebac_ode,
                                   parms = theta,
                                   method = "ode45"))
      return(trajectory)
    }
  }
  
  if(model == "dde_mass_link_both"){
    phagebac_simulateDeterministic <- function(theta,init.state,times) {
      phagebac_dde <- function(time, state, parameters) {
        
        ## parameters
        Ge = 1.678177 
        Gt = 1.605017 
        Get = 0 #1.562914
        Nmax = 2.140317e+09 
        #ll for this set is -4.026931e+08 
        
        beta = 1/parameters[["beta"]]
        L = parameters[["L"]]
        #phi = 1/parameters[["phi"]]
        tau = 0.5
        
        ## states
        Be = state[["Be"]]
        Bt = state[["Bt"]]
        Bet = state[["Bet"]]
        P = state[["P"]]
        
        N = Be + Bt + Bet
        
        if(time < 1){
          Be_past = 0
          Bt_past = 0
          Bet_past = 0
          P_past = 0
        } else {
          Be_past = lagvalue(time - tau, 1)
          Bt_past = lagvalue(time - tau, 2)
          Bet_past = lagvalue(time - tau, 3)
          P_past = lagvalue(time - tau, 4)
        }
        
        N_past = Be_past + Bt_past + Bet_past
        
        #link both
        dBe = Ge*(1 - N/Nmax)*Be - beta*(1-N/Nmax) * Be * P
        dBt = Gt*(1 - N/Nmax)*Bt - beta*(1-N/Nmax) * Bt * P
        dBet = Get*(1 - N/Nmax)*Bet - beta*(1-N/Nmax) * Bet * P
        
        dP =  - beta*(1-N/Nmax) * P * N + N_past * P_past * beta*(1-N_past/Nmax) * L*(1-N_past/Nmax)
        
        return(list(c(dBe, dBt, dBet, dP)))
      }
      trajectory <- data.frame(dede(y = init.state,
                                    times = times,
                                    func = phagebac_dde,
                                    parms = theta))
      return(trajectory)
    }
  }
  
  #third set
  if(model == "frequentist"){
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
        
        P_inf = (1-exp(-beta*N))*P
        B_inf = (1-exp(-P_inf/N))*N
        
        dBe = Ge*(1 - N/Nmax)*Be - B_inf* (Be/N)
        dBt = Gt*(1 - N/Nmax)*Bt - B_inf* (Bt/N)
        dBet = Get*(1 - N/Nmax)*Bet - B_inf* (Bet/N)
        
        dP = B_inf * L - P_inf
        
        return(list(c(dBe, dBt, dBet, dP)))
      }
      trajectory <- data.frame(ode(y = init.state,
                                   times = times,
                                   func = phagebac_ode,
                                   parms = theta,
                                   method = "ode45"))
      return(trajectory)
    }
  }
  
  if(model == "frequentist_link_both"){
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
        #phi = 1/parameters[["phi"]]
        
        ## states
        Be = state[["Be"]]
        Bt = state[["Bt"]]
        Bet = state[["Bet"]]
        P = state[["P"]]
        
        N = Be + Bt + Bet
        
        beta = beta * (1 - N/Nmax)
        L = L * (1 - N/Nmax)
        
        P_inf = (1-exp(-beta*N))*P
        B_inf = (1-exp(-P_inf/N))*N
        
        dBe = Ge*(1 - N/Nmax)*Be - B_inf* (Be/N)
        dBt = Gt*(1 - N/Nmax)*Bt - B_inf* (Bt/N)
        dBet = Get*(1 - N/Nmax)*Bet - B_inf* (Bet/N)
        
        dP = B_inf * L - P_inf
        
        return(list(c(dBe, dBt, dBet, dP)))
      }
      trajectory <- data.frame(ode(y = init.state,
                                   times = times,
                                   func = phagebac_ode,
                                   parms = theta,
                                   method = "ode45"))
      return(trajectory)
    }
  }
  
  if(model == "dde_frequentist_link_both"){
    phagebac_simulateDeterministic <- function(theta,init.state,times) {
      phagebac_dde <- function(time, state, parameters) {
        
        ## parameters
        Ge = 1.678177 
        Gt = 1.605017 
        Get = 0 #1.562914
        Nmax = 2.140317e+09 
        #ll for this set is -4.026931e+08 
        
        beta = 1/parameters[["beta"]]
        L = parameters[["L"]]
        #phi = 1/parameters[["phi"]]
        
        tau = 0.5#parameters[["tau"]]
        
        ## states
        Be = state[["Be"]]
        Bt = state[["Bt"]]
        Bet = state[["Bet"]]
        P = state[["P"]]
        
        N = Be + Bt + Bet
        
        #tau = tau * (N/Nmax)
        
        if(time < tau){
          Be_past = 0
          Bt_past = 0
          Bet_past = 0
          P_past = 0
        } else {
          Be_past = lagvalue(time - tau, 1)
          Bt_past = lagvalue(time - tau, 2)
          Bet_past = lagvalue(time - tau, 3)
          P_past = lagvalue(time - tau, 4)
        }
        
        N_past = ifelse(time < tau, 1, Be_past + Bt_past + Bet_past)
        
        L = L * (1 - N/Nmax)
        
        P_inf = (1-exp(-beta*(1 - N/Nmax)*N))*P
        B_inf = (1-exp(-P_inf/N))*N
        
        P_inf_past = ifelse(time < tau, 0, (1-exp(-beta*(1-N_past/Nmax) *N_past))*P_past)
        B_inf_past = ifelse(time < tau, 0, (1-exp(-P_inf_past/N_past))*N_past)
        
        dBe = Ge*(1 - N/Nmax)*Be - B_inf* (Be/N)
        dBt = Gt*(1 - N/Nmax)*Bt - B_inf* (Bt/N)
        dBet = Get*(1 - N/Nmax)*Bet - B_inf* (Bet/N)
        
        dP = B_inf_past * L - P_inf
        
        return(list(c(dBe, dBt, dBet, dP)))
      }
      trajectory <- data.frame(dede(y = init.state,
                                    times = times,
                                    func = phagebac_dde,
                                    parms = theta))
      return(trajectory)
    }
  }
  
  if(model == "dde_frequentist_decay_link_both"){
    phagebac_simulateDeterministic <- function(theta,init.state,times) {
      phagebac_dde <- function(time, state, parameters) {
        
        ## parameters
        Ge = 1.678177 
        Gt = 1.605017 
        Get = 0 #1.562914
        Nmax = 2.140317e+09 
        #ll for this set is -4.026931e+08 
        
        beta = 1/parameters[["beta"]]
        L = parameters[["L"]]
        phi = 1/parameters[["phi"]]
        
        tau = 0.5#parameters[["tau"]]
        
        ## states
        Be = state[["Be"]]
        Bt = state[["Bt"]]
        Bet = state[["Bet"]]
        P = state[["P"]]
        
        N = Be + Bt + Bet
        
        #tau = tau * (N/Nmax)
        
        if(time < tau){
          Be_past = 0
          Bt_past = 0
          Bet_past = 0
          P_past = 0
        } else {
          Be_past = lagvalue(time - tau, 1)
          Bt_past = lagvalue(time - tau, 2)
          Bet_past = lagvalue(time - tau, 3)
          P_past = lagvalue(time - tau, 4)
        }
        
        N_past = ifelse(time < tau, 1, Be_past + Bt_past + Bet_past)
        
        L = L * (1 - N/Nmax)
        
        P_inf = (1-exp(-beta*(1 - N/Nmax)*N))*P
        B_inf = (1-exp(-P_inf/N))*N
        
        P_inf_past = ifelse(time < tau, 0, (1-exp(-beta*(1-N_past/Nmax) *N_past))*P_past)
        B_inf_past = ifelse(time < tau, 0, (1-exp(-P_inf_past/N_past))*N_past)
        
        dBe = Ge*(1 - N/Nmax)*Be - B_inf* (Be/N)
        dBt = Gt*(1 - N/Nmax)*Bt - B_inf* (Bt/N)
        dBet = Get*(1 - N/Nmax)*Bet - B_inf* (Bet/N)
        
        dP = B_inf_past * L - P_inf - phi * P
        
        return(list(c(dBe, dBt, dBet, dP)))
      }
      trajectory <- data.frame(dede(y = init.state,
                                    times = times,
                                    func = phagebac_dde,
                                    parms = theta))
      return(trajectory)
    }
  }
  
  #run stuff
  phagebac_model$simulate = phagebac_simulateDeterministic
  
  phagebac_model
  
}


## Go ##########


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
    
    traj = phagebac_model$simulate(theta, init.state, times)
    
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



for(model_name in names(all_theta)){
  
  cat("\nWorking on", model_name)
  
  phagebac_model = choose_model(model = model_name)
  
  trace_model = all_trace[[model_name]]
  trace_model = trace_model[-c(1:5000),]
  
  #replicate 4
  init.state = c(Be = lab_data_trans$Be[1], Bt = lab_data_trans$Bt[1], Bet = 0,
                 P = lab_data_trans$P[1])
  
  traj = multi_run(phagebac_model, trace_model, init.state, times = seq(0, 30, 1))
  
  p4 = ggplot() +
    geom_line(data = traj, aes(time, Be, colour = "Model", linetype = "Bacteria")) +
    geom_ribbon(data = traj, aes(x = time, ymin = pmax(0.1, Be-1.96*Be_sd), ymax = Be+1.96*Be_sd),
                alpha = 0.3, fill = "darkturquoise") +
    geom_line(data = traj, aes(time, Bt, colour = "Model", linetype = "Bacteria")) +
    geom_ribbon(data = traj, aes(x = time, ymin = pmax(0.1,Bt-1.96*Bt_sd), ymax = Bt+1.96*Bt_sd),
                alpha = 0.3, fill = "darkturquoise") +
    geom_line(data = traj, aes(time, P, colour = "Model", linetype = "Phage")) +
    geom_ribbon(data = traj, aes(x = time, ymin =pmax(0.1, P-1.96*P_sd), ymax = P+1.96*P_sd),
                alpha = 0.3, fill = "darkturquoise") +
    geom_line(data = lab_data_transM %>% filter(Bacteria != "P" & Bacteria != "DRP"),
              aes(Time, Mean, group = Bacteria, colour = "Data", linetype = "Bacteria")) +
    geom_errorbar(data = lab_data_transM %>% filter(Bacteria != "P" & Bacteria != "DRP"), 
                  aes(x = Time, ymin = pmax(Mean - se, 0), ymax = Mean + se, group = Bacteria,
                      colour= "Data", linetype = "Bacteria")) +
    geom_line(data = lab_data_transM %>% filter(Bacteria == "P"), 
              aes(Time, Mean, group = Bacteria, colour = "Data", linetype = "Phage")) +
    geom_errorbar(data = lab_data_transM %>% filter(Bacteria == "P"), 
                  aes(x = Time, ymin = pmax(Mean - se, 0), ymax = Mean + se, group = Bacteria,
                      colour= "Data", linetype = "Phage")) +
    scale_y_continuous(trans=log10_trans(),
                       breaks=trans_breaks("log10", function(x) 10^x),
                       labels=trans_format("log10", math_format(10^.x))) +
    coord_cartesian(ylim = c(0.1, 10e10)) +
    labs(y = "cfu / pfu per mL", x = "Time (hours)", title = "10^4", 
         linetype = "Organism:", colour = "Source:") +
    guides(fill = F) +
    theme_bw()
  
  
  #replicate 5
  init.state = c(Be = lab_data_trans5$Be[1], Bt = lab_data_trans5$Bt[1], Bet = 0,
                 P = lab_data_trans5$P[1])
  
  traj = multi_run(phagebac_model, trace_model, init.state, times = seq(0, 30, 1))
  
  p5 = ggplot() +
    geom_line(data = traj, aes(time, Be, colour = "Model", linetype = "Bacteria")) +
    geom_ribbon(data = traj, aes(x = time, ymin = pmax(0.1,Be-1.96*Be_sd), ymax = Be+1.96*Be_sd),
                alpha = 0.3, fill = "darkturquoise") +
    geom_line(data = traj, aes(time, Bt, colour = "Model", linetype = "Bacteria")) +
    geom_ribbon(data = traj, aes(x = time, ymin = pmax(0.1,Bt-1.96*Bt_sd), ymax = Bt+1.96*Bt_sd),
                alpha = 0.3, fill = "darkturquoise") +
    geom_line(data = traj, aes(time, P, colour = "Model", linetype = "Phage")) +
    geom_ribbon(data = traj, aes(x = time, ymin = pmax(0.1,P-1.96*P_sd), ymax = P+1.96*P_sd),
                alpha = 0.3, fill = "darkturquoise") +
    geom_line(data = lab_data_trans5M %>% filter(Bacteria != "P" & Bacteria != "DRP"),
              aes(Time, Mean, group = Bacteria, colour = "Data", linetype = "Bacteria")) +
    geom_errorbar(data = lab_data_trans5M %>% filter(Bacteria != "P" & Bacteria != "DRP"), 
                  aes(x = Time, ymin = pmax(Mean - se, 0), ymax = Mean + se, group = Bacteria,
                      colour= "Data", linetype = "Bacteria")) +
    geom_line(data = lab_data_trans5M %>% filter(Bacteria == "P"), 
              aes(Time, Mean, group = Bacteria, colour = "Data", linetype = "Phage")) +
    geom_errorbar(data = lab_data_trans5M %>% filter(Bacteria == "P"), 
                  aes(x = Time, ymin = pmax(Mean - se, 0), ymax = Mean + se, group = Bacteria,
                      colour= "Data", linetype = "Phage")) +
    scale_y_continuous(trans=log10_trans(),
                       breaks=trans_breaks("log10", function(x) 10^x),
                       labels=trans_format("log10", math_format(10^.x))) +
    coord_cartesian(ylim = c(0.1, 10e10)) +
    labs(y = "cfu / pfu per mL", x = "Time (hours)", title = "10^5", 
         linetype = "Organism:", colour = "Source:") +
    theme_bw()
  
  
  #replicate 3
  init.state = c(Be = lab_data_trans3$Be[1], Bt = lab_data_trans3$Bt[1], Bet = 0,
                 P = lab_data_trans3$P[1])

  traj = multi_run(phagebac_model, trace_model, init.state, times = seq(0, 30, 1))
  
  p3 = ggplot() +
    geom_line(data = traj, aes(time, Be, colour = "Model", linetype = "Bacteria")) +
    geom_ribbon(data = traj, aes(x = time, ymin = pmax(0.1,Be-1.96*Be_sd), ymax = Be+1.96*Be_sd),
                alpha = 0.3, fill = "darkturquoise") +
    geom_line(data = traj, aes(time, Bt, colour = "Model", linetype = "Bacteria")) +
    geom_ribbon(data = traj, aes(x = time, ymin = pmax(0.1,Bt-1.96*Bt_sd), ymax = Bt+1.96*Bt_sd),
                alpha = 0.3, fill = "darkturquoise") +
    geom_line(data = traj, aes(time, P, colour = "Model", linetype = "Phage")) +
    geom_ribbon(data = traj, aes(x = time, ymin = pmax(0.1,P-1.96*P_sd), ymax = P+1.96*P_sd),
                alpha = 0.3, fill = "darkturquoise") +
    geom_line(data = lab_data_trans3M %>% filter(Bacteria != "P" & Bacteria != "DRP"),
              aes(Time, Mean, group = Bacteria, colour = "Data", linetype = "Bacteria")) +
    geom_errorbar(data = lab_data_trans3M %>% filter(Bacteria != "P" & Bacteria != "DRP"), 
                  aes(x = Time, ymin = pmax(Mean - se, 0), ymax = Mean + se, group = Bacteria,
                      colour= "Data", linetype = "Bacteria")) +
    geom_line(data = lab_data_trans3M %>% filter(Bacteria == "P"), 
              aes(Time, Mean, group = Bacteria, colour = "Data", linetype = "Phage")) +
    geom_errorbar(data = lab_data_trans3M %>% filter(Bacteria == "P"), 
                  aes(x = Time, ymin = pmax(Mean - se, 0), ymax = Mean + se, group = Bacteria,
                      colour= "Data", linetype = "Phage")) +
    scale_y_continuous(trans=log10_trans(),
                       breaks=trans_breaks("log10", function(x) 10^x),
                       labels=trans_format("log10", math_format(10^.x))) +
    coord_cartesian(ylim = c(0.1, 10e10)) +
    labs(y = "cfu / pfu per mL", x = "Time (hours)", title = "10^3", 
         linetype = "Organism:", colour = "Source:") +
    theme_bw()
  
  #final plot
  legend = get_legend(p4 + theme(legend.position = "bottom", legend.box = "vertical"))
  
  
  plot_grid(p3 + theme(legend.position = "none"), 
            p4 + theme(legend.position = "none"), 
            p5 + theme(legend.position = "none"),
            legend)
  filename = paste0("multi_runs_", model_name, ".png")
  ggsave(filename)
  
}


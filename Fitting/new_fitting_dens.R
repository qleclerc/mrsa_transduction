
library(dplyr)
library(openxlsx)
library(reshape2)
library(BayesianTools)
library(deSolve)
library(ggplot2)
library(scales)
library(cowplot)
library(openxlsx)

obs103 = read.csv(here::here("Data", "transduction_summary_10_3.csv")) %>%
  select(Time, Bacteria, Mean) %>%
  dcast(Time~Bacteria) %>%
  select(-Total) %>%
  rename(time = Time, Bet = DRP, Be = EryR, Bt = TetR, Pl = P) %>%
  mutate(Be = round(Be),
         Bt = round(Bt),
         Bet = round(Bet),
         Pl = round(Pl))

obs104 = read.csv(here::here("Data", "transduction_summary_10_4.csv")) %>%
  select(Time, Bacteria, Mean) %>%
  dcast(Time~Bacteria) %>%
  select(-Total) %>%
  rename(time = Time, Bet = DRP, Be = EryR, Bt = TetR, Pl = P) %>%
  mutate(Be = round(Be),
         Bt = round(Bt),
         Bet = round(Bet),
         Pl = round(Pl))

obs105 = read.csv(here::here("Data", "transduction_summary_10_5.csv")) %>%
  select(Time, Bacteria, Mean) %>%
  dcast(Time~Bacteria) %>%
  select(-Total) %>%
  rename(time = Time, Bet = DRP, Be = EryR, Bt = TetR, Pl = P) %>%
  mutate(Be = round(Be),
         Bt = round(Bt),
         Bet = round(Bet),
         Pl = round(Pl))



MOI_data = read.xlsx(here::here("Data", "varying_MOI_data.xlsx"))

MOI_data = apply(MOI_data, c(1,2), as.numeric)
MOI_data = MOI_data[-1, (-c(20:28))]


data_final = data.frame(init_bac = rep(4e6, nrow(MOI_data)))
data_final$init_pha = data_final$init_bac*MOI_data[,1]
data_final$Be = rowMeans(MOI_data[,c(2:10)])
data_final$Be_sd = apply(MOI_data[,c(2:10)], 1, function(x) sd(x)/sqrt(3))
data_final$Bt = rowMeans(MOI_data[,c(11:19)])
data_final$Bt_sd = apply(MOI_data[,c(11:19)], 1, function(x) sd(x)/sqrt(3))
data_final$Bet = rowMeans(MOI_data[,c(20:28)])
data_final$Bet_sd = apply(MOI_data[,c(20:28)], 1, function(x) sd(x)/sqrt(3))
data_final$Pl = rowMeans(MOI_data[,c(29:37)])
data_final$Pl_sd = apply(MOI_data[,c(29:37)], 1, function(x) sd(x)/sqrt(3))

MOI_data = round(data_final)


phage_tr_model = function(parameters,
                          init.state,
                          times,
                          mode = "dens",    #dens pow or hill
                          link_L = T){      #link burst size?
  
  
  model_dde = function(time, state, parameters) {
    
    mu_e = 1.679
    mu_t = 1.608
    mu_et = 1.566
    Nmax = 2.14e9
    
    beta = parameters[[1]]/1e10
    L = parameters[[2]]
    gamma = 0
    alpha = parameters[[3]]/1e8
    tau = parameters[[4]]
    pow = parameters[[5]]
    P_lim = parameters[[6]]*1e8
    
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
    
    if(link_L) L = L * link + 1
    
    if(mode == "pow"){
      
      F_PL = (beta * Pl)^pow
      F_PE = (beta * Pe)^pow
      F_PT = (beta * Pt)^pow
      
      F_PL_past = (beta * Pl_past)^pow
      
    } else if(mode == "dens"){
      
      F_PL = beta * Pl
      F_PE = beta * Pe
      F_PT = beta * Pt
      
      F_PL_past = beta * Pl_past
      
    } else if(mode == "hill"){
      
      F_PL = (beta * Pl)/(1+(Pl^pow)/(P_lim^pow))
      F_PE = (beta * Pe)/(1+(Pe^pow)/(P_lim^pow))
      F_PT = (beta * Pt)/(1+(Pt^pow)/(P_lim^pow))
      
      F_PL_past = (beta * Pl_past)/(1+(Pl_past^pow)/(P_lim^pow))
      
      
    } else {
      
      stop("Only pow, dens or hill are valid modes")
      
    }
    
    dBe = mu_e * link * Be - F_PL * Be - F_PT * Be
    dBt = mu_t * link * Bt - F_PL * Bt - F_PE * Bt
    dBet = mu_et * link * Bet  - F_PL * Bet + F_PT * Be + F_PE * Bt
    
    dPl = F_PL_past * N_past * L * (1 - alpha*(Be_past+Bt_past+2*Bet_past)/N_past) -
      F_PL * N - gamma * Pl
    dPe = F_PL_past * N_past * L * alpha * (Be_past + Bet_past)/N_past -
      F_PE * N - gamma * Pe
    dPt = F_PL_past * N_past * L * alpha * (Bt_past + Bet_past)/N_past -
      F_PT * N - gamma * Pt
    
    return(list(c(dBe, dBt, dBet, dPl, dPe, dPt)))
    
  }
  
  trajectory <- data.frame(dede(y = init.state,
                                times = times,
                                func = model_dde,
                                parms = parameters, method = "lsode"))
  
  return(trajectory)
  
}

#hill pars c(2.9, 60, 1.2, 0.67, 1, 90)
#pow pars c(0.31, 60, 1, 0.67, 0.72, 90)

# load reference parameter definition (upper, lower prior)
refPars <- data.frame(best = c(2.9, 60, 1.2, 0.67, 1, 90),
                      lower = c(0.1, 10, 0.1, 0.5, 0.6, 0.1),
                      upper = c(100, 120, 10, 0.8, 1, 1000))
rownames(refPars) = c("beta", "L", "alpha", "tau", "pow", "P_lim")

parSel = c(1:4)

# here is the likelihood 
likelihood <- function(par, mode = "dens"){
  # set parameters that are not calibrated on default values 
  x = refPars$best
  x[parSel] = par
  predicted103 <- phage_tr_model(x,
                                 c(Be = obs103$Be[1], Bt = obs103$Bt[1], Bet = 0,
                                   Pl = obs103$Pl[1], Pe = 0, Pt = 0),
                                 seq(0,24,1), mode = mode)[c(1:9, 17:25), -1] # replace here VSEM with your model 
  predicted103[predicted103<=0] = 0.00001
  
  predicted105 <- phage_tr_model(x,
                                 c(Be = obs105$Be[1], Bt = obs105$Bt[1], Bet = 0,
                                   Pl = obs105$Pl[1], Pe = 0, Pt = 0),
                                 seq(0,24,1), mode = mode)[c(1:9, 17:25), -1] # replace here VSEM with your model 
  predicted105[predicted105<=0] = 0.00001
  
  
  llValues1 = dpois(x = round(obs103$Pl/(10^(pmax(floor(log10(obs103$Pl)),1)-1))),
                    lambda = predicted103$Pl/(10^(pmax(floor(log10(obs103$Pl)),1)-1)),
                    log = T)
  llValues2 = 2*dpois(obs103$Bet, predicted103$Bet, log = TRUE)
  
  llValues3 = dpois(x = round(obs105$Pl/(10^(pmax(floor(log10(obs105$Pl)),1)-1))),
                    lambda = predicted105$Pl/(10^(pmax(floor(log10(obs105$Pl)),1)-1)),
                    log = T)
  llValues4 = 2*dpois(obs105$Bet, predicted105$Bet, log = TRUE)
  
  llValues5 = dpois(x = round(obs105$Be[12:18]/(10^(pmax(floor(log10(obs105$Be[12:18])),1)-1))),
                    lambda = predicted105$Be[12:18]/(10^(pmax(floor(log10(obs105$Be[12:18])),1)-1)),
                    log = T)
  llValues6 = dpois(x = round(obs105$Bt[12:18]/(10^(pmax(floor(log10(obs105$Bt[12:18])),1)-1))),
                    lambda = predicted105$Bt[12:18]/(10^(pmax(floor(log10(obs105$Bt[12:18])),1)-1)),
                    log = T)
  
  # #24h varying MOI
  llValues7 = rep(0, nrow(MOI_data))
  for(i in 1:nrow(MOI_data)){
    
    predictedi <- phage_tr_model(x,
                                 c(Be = MOI_data$init_bac[i], Bt = MOI_data$init_bac[i], Bet = 0,
                                   Pl = MOI_data$init_pha[i], Pe = 0, Pt = 0),
                                 seq(0,24,1), mode = mode)[25,-1]
    predictedi[predictedi<=0] = 0.00001
    
    llValues7[i] = (dpois(x = round(MOI_data$Pl[i]/(10^(max(floor(log10(MOI_data$Pl[i])),1)-1))),
                          lambda = predictedi$Pl/(10^(max(floor(log10(MOI_data$Pl[i])),1)-1)),
                          log = T)) + 2*dpois(MOI_data$Bet[i], predictedi$Bet, log = TRUE)
    
    
  }
  
  
  return(sum(llValues1,llValues2,llValues3,llValues4,llValues5,llValues6,llValues7))
}

# optional, you can also directly provide lower, upper in the createBayesianSetup, see help
prior <- createUniformPrior(lower = refPars$lower[parSel], 
                            upper = refPars$upper[parSel],
                            best = refPars$best[parSel])

bayesianSetup <- createBayesianSetup(likelihood, prior, names = rownames(refPars)[parSel])

# settings for the sampler, iterations should be increased for real applicatoin
settings <- list(iterations = 50000, nrChains = 2)

out <- runMCMC(bayesianSetup = bayesianSetup, settings = settings)


median_params = rbind(out[[1]][["chain"]][[1]],
                      out[[1]][["chain"]][[2]],
                      out[[1]][["chain"]][[3]],
                      out[[2]][["chain"]][[1]],
                      out[[2]][["chain"]][[2]],
                      out[[2]][["chain"]][[3]])

# median_params = rbind(out[[1]][["chain"]],
#                       out[[2]][["chain"]])
# colnames(median_params)[1:5] = c("beta", "L", "alpha", "tau", "pow")

median_params = apply(median_params, 2, median)[1:4]
names(median_params) = c("beta", "L", "alpha", "tau")
median_params["pow"] = 1
median_params["P_lim"] = 1

best103 = phage_tr_model(median_params,
                         c(Be = obs103$Be[1], Bt = obs103$Bt[1], Bet = 0,
                           Pl = obs103$Pl[1], Pe = 0, Pt = 0), seq(0,30,0.1),
                         mode = "dens")

best104 = phage_tr_model(median_params,
                         c(Be = obs104$Be[1], Bt = obs104$Bt[1], Bet = 0,
                           Pl = obs104$Pl[1], Pe = 0, Pt = 0), seq(0,30,0.1),
                         mode = "dens")


best105 = phage_tr_model(median_params,
                         c(Be = obs105$Be[1], Bt = obs105$Bt[1], Bet = 0,
                           Pl = obs105$Pl[1], Pe = 0, Pt = 0), seq(0,30,0.1),
                         mode = "dens")

p1 = ggplot()+
  geom_line(data = best105, aes(time, Be, colour = "Be")) +
  geom_line(data = best105, aes(time, Bt, colour = "Bt")) +
  geom_line(data = best105, aes(time, Bet, colour = "Bet")) +
  geom_line(data = best105, aes(time, Pl, colour = "Pl")) +
  geom_line(data = obs105, aes(time, Be, colour = "Be"), linetype = 2) +
  geom_line(data = obs105, aes(time, Bt, colour = "Bt"), linetype = 2) +
  geom_line(data = obs105, aes(time, Bet, colour = "Bet"), linetype = 2) +
  geom_line(data = obs105, aes(time, Pl, colour = "Pl"), linetype = 2) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  coord_cartesian(ylim = c(1,1e10)) +
  theme_bw() +
  labs(y = "Bacteria and phage", x = "Time", colour = "") +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.text = element_text(size=12))

p2 = ggplot()+
  geom_line(data = best104, aes(time, Be, colour = "Be")) +
  geom_line(data = best104, aes(time, Bt, colour = "Bt")) +
  geom_line(data = best104, aes(time, Bet, colour = "Bet")) +
  geom_line(data = best104, aes(time, Pl, colour = "Pl")) +
  geom_line(data = obs104, aes(time, Be, colour = "Be"), linetype = 2) +
  geom_line(data = obs104, aes(time, Bt, colour = "Bt"), linetype = 2) +
  geom_line(data = obs104, aes(time, Bet, colour = "Bet"), linetype = 2) +
  geom_line(data = obs104, aes(time, Pl, colour = "Pl"), linetype = 2) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  coord_cartesian(ylim = c(1,1e10)) +
  theme_bw() +
  labs(y = "Bacteria and phage", x = "Time", colour = "") +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.text = element_text(size=12))

p3 = ggplot()+
  geom_line(data = best103, aes(time, Be, colour = "Be")) +
  geom_line(data = best103, aes(time, Bt, colour = "Bt")) +
  geom_line(data = best103, aes(time, Bet, colour = "Bet")) +
  geom_line(data = best103, aes(time, Pl, colour = "Pl")) +
  geom_line(data = obs103, aes(time, Be, colour = "Be"), linetype = 2) +
  geom_line(data = obs103, aes(time, Bt, colour = "Bt"), linetype = 2) +
  geom_line(data = obs103, aes(time, Bet, colour = "Bet"), linetype = 2) +
  geom_line(data = obs103, aes(time, Pl, colour = "Pl"), linetype = 2) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  coord_cartesian(ylim = c(1,1e10)) +
  theme_bw() +
  labs(y = "Bacteria and phage", x = "Time", colour = "") +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.text = element_text(size=12))

pp=plot_grid(p1+theme(legend.position = "none"),
             p2+theme(legend.position = "none"),
             p3+theme(legend.position = "none"),
             nrow = 1)

ggsave(here::here("Fitting", "dens_fitted.png"), pp)
saveRDS(out, here::here("Fitting", "dens_fitted_out.rds"))


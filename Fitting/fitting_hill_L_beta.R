
library(dplyr)
library(openxlsx)
library(reshape2)
library(BayesianTools)
library(deSolve)
library(ggplot2)
library(scales)
library(cowplot)
library(openxlsx)

# obs103 = read.csv(here::here("Data", "transduction_summary_10_3.csv")) %>%
#   select(Time, Bacteria, Mean) %>%
#   dcast(Time~Bacteria) %>%
#   select(-Total) %>%
#   rename(time = Time, Bet = DRP, Be = EryR, Bt = TetR, Pl = P) %>%
#   mutate(Be = round(Be),
#          Bt = round(Bt),
#          Bet = round(Bet),
#          Pl = round(Pl))

obs104 = read.csv(here::here("Data", "transduction_summary_10_4.csv")) %>%
  select(Time, Bacteria, Mean) %>%
  dcast(Time~Bacteria) %>%
  select(-Total) %>%
  rename(time = Time, Bet = DRP, Be = EryR, Bt = TetR, Pl = P) %>%
  mutate(Be = round(Be),
         Bt = round(Bt),
         Bet = round(Bet),
         Pl = round(Pl))

# obs105 = read.csv(here::here("Data", "transduction_summary_10_5.csv")) %>%
#   select(Time, Bacteria, Mean) %>%
#   dcast(Time~Bacteria) %>%
#   select(-Total) %>%
#   rename(time = Time, Bet = DRP, Be = EryR, Bt = TetR, Pl = P) %>%
#   mutate(Be = round(Be),
#          Bt = round(Bt),
#          Bet = round(Bet),
#          Pl = round(Pl))

obs103 = read.csv(here::here("Fitting", "103res_L_beta.csv"))[,-8] %>% round()
obs105 = read.csv(here::here("Fitting", "105res_L_beta.csv"))[,-8] %>% round()

source(here::here("Model", "model.R"))

#hill pars c(2.9, 60, 1.2, 0.67, 1, 90)
#pow pars c(0.31, 60, 1, 0.67, 0.72, 90)

phage_tr_model(c(10, 60, 2, 0.67, 70),
               c(Be = obs105$Be[1], Bt = obs105$Bt[1], Bet = 0,
                 Pl = obs105$Pl[1], Pe = 0, Pt = 0), seq(0,24,0.1), mode = "hill",
               link_beta = T, link_L = T) %>%
  ggplot()+
  geom_line(aes(time, Be, colour = "Be")) +
  geom_line(aes(time, Bt, colour = "Bt")) +
  geom_line(aes(time, Bet, colour = "Bet")) +
  geom_line(aes(time, Pl, colour = "Pl")) +
  geom_line(data = obs105, aes(time, Be, colour = "Be"), linetype = 2) +
  geom_line(data = obs105, aes(time, Bt, colour = "Bt"), linetype = 2) +
  geom_line(data = obs105, aes(time, Bet, colour = "Bet"), linetype = 2) +
  geom_line(data = obs105, aes(time, Pl, colour = "Pl"), linetype = 2) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  coord_cartesian(ylim = c(1,1e11)) +
  theme_bw() +
  labs(y = "Bacteria and phage", x = "Time", colour = "") +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.text = element_text(size=12))
# 
# ggsave("better_hill.png")

# load reference parameter definition (upper, lower prior)
refPars <- data.frame(best = c(10, 60, 2, 0.67, 70),
                      lower = c(0.1, 10, 1, 0.6, 1),
                      upper = c(100, 90, 1000, 0.8, 1000))
rownames(refPars) = c("beta", "L", "alpha", "tau", "P_lim")

parSel = c(1:5)

# here is the likelihood 
likelihood <- function(par, mode = "hill"){
  # set parameters that are not calibrated on default values 
  x = refPars$best
  x[parSel] = par
  predicted103 <- phage_tr_model(x,
                                 c(Be = obs103$Be[1], Bt = obs103$Bt[1], Bet = 0,
                                   Pl = obs103$Pl[1], Pe = 0, Pt = 0),
                                 seq(0,30,1), mode = mode,
                                 link_L = T, link_beta = T) # replace here VSEM with your model 
  predicted103[predicted103<=0] = 0.00001
  
  predicted105 <- phage_tr_model(x,
                                 c(Be = obs105$Be[1], Bt = obs105$Bt[1], Bet = 0,
                                   Pl = obs105$Pl[1], Pe = 0, Pt = 0),
                                 seq(0,30,1), mode = mode,
                                 link_L = T, link_beta = T) # replace here VSEM with your model 
  predicted105[predicted105<=0] = 0.00001
  
  
  llValues1 = dpois(x = round(obs103$Pl/(10^(pmax(floor(log10(obs103$Pl)),1)-1))),
                    lambda = predicted103$Pl/(10^(pmax(floor(log10(obs103$Pl)),1)-1)),
                    log = T)
  llValues2 = dpois(obs103$Bet,
                    predicted103$Bet,
                    log = T)
  
  llValues3 = dpois(x = round(obs105$Pl/(10^(pmax(floor(log10(obs105$Pl)),1)-1))),
                    lambda = predicted105$Pl/(10^(pmax(floor(log10(obs105$Pl)),1)-1)),
                    log = T)
  llValues4 = dpois(obs105$Bet,
                    predicted105$Bet,
                    log = T)
  
  llValues5 = 2*dpois(x = round(obs105$Be/(10^(pmax(floor(log10(obs105$Be)),1)-1))),
                      lambda = predicted105$Be/(10^(pmax(floor(log10(obs105$Be)),1)-1)),
                      log = T)
  llValues6 = 2*dpois(x = round(obs105$Bt/(10^(pmax(floor(log10(obs105$Bt)),1)-1))),
                      lambda = predicted105$Bt/(10^(pmax(floor(log10(obs105$Bt)),1)-1)),
                      log = T)
  
  # #24h varying MOI
  # llValues7 = rep(0, nrow(MOI_data))
  # for(i in 1:nrow(MOI_data)){
  # 
  #   predictedi <- phage_tr_model(x,
  #                                c(Be = MOI_data$init_bac[i], Bt = MOI_data$init_bac[i], Bet = 0,
  #                                  Pl = MOI_data$init_pha[i], Pe = 0, Pt = 0),
  #                                seq(0,24,1), mode = mode)[25,-1]
  #   predictedi[predictedi<=0] = 0.00001
  #   
  #   llValues7[i] = (dpois(x = round(MOI_data$Pl[i]/(10^(max(floor(log10(MOI_data$Pl[i])),1)-1))),
  #                         lambda = predictedi$Pl/(10^(max(floor(log10(MOI_data$Pl[i])),1)-1)),
  #                         log = T)) + 2*dpois(MOI_data$Bet[i], predictedi$Bet, log = TRUE)
  # 
  # 
  # }
  
  
  return(sum(llValues1,llValues2,llValues3,llValues4,llValues5,llValues6))
}

# optional, you can also directly provide lower, upper in the createBayesianSetup, see help
prior <- createUniformPrior(lower = refPars$lower[parSel], 
                            upper = refPars$upper[parSel],
                            best = refPars$best[parSel])

bayesianSetup <- createBayesianSetup(likelihood, prior, names = rownames(refPars)[parSel])

# settings for the sampler, iterations should be increased for real applicatoin
settings <- list(iterations = 100000, nrChains = 2)
# settings <- list(iterations = 10000, adapt = T, DRlevels = 2,
#                  gibbsProbabilities = NULL, temperingFunction = NULL, optimize = T, message = FALSE)

out <- runMCMC(bayesianSetup = bayesianSetup, settings = settings)

#out = readRDS(here::here("Fitting", "hill_L_beta_fitted_out.rds"))

median_params = rbind(out[[1]][["chain"]][[1]],
                      out[[1]][["chain"]][[2]],
                      out[[1]][["chain"]][[3]],
                      out[[2]][["chain"]][[1]],
                      out[[2]][["chain"]][[2]],
                      out[[2]][["chain"]][[3]])

# median_params = rbind(out[[1]][["chain"]],
#                       out[[2]][["chain"]])
# colnames(median_params)[1:5] = c("beta", "L", "alpha", "tau", "pow")

median_params = apply(tail(median_params), 2, median)[1:5]

best103 = phage_tr_model(median_params,
                         c(Be = obs103$Be[1], Bt = obs103$Bt[1], Bet = 0,
                           Pl = obs103$Pl[1], Pe = 0, Pt = 0), seq(0,30,0.1),
                         mode = "hill",
                         link_L = T, link_beta = T)

best104 = phage_tr_model(median_params,
                         c(Be = obs104$Be[1], Bt = obs104$Bt[1], Bet = 0,
                           Pl = obs104$Pl[1], Pe = 0, Pt = 0), seq(0,30,0.1),
                         mode = "hill",
                         link_L = T, link_beta = T)


best105 = phage_tr_model(median_params,
                         c(Be = obs105$Be[1], Bt = obs105$Bt[1], Bet = 0,
                           Pl = obs105$Pl[1], Pe = 0, Pt = 0), seq(0,30,0.1),
                         mode = "hill",
                         link_L = T, link_beta = T)

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

ggsave(here::here("Fitting", "hill_L_beta_fitted.png"), pp)
saveRDS(out, here::here("Fitting", "hill_L_beta_fitted_out.rds"))


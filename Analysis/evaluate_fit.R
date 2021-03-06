# SETUP ####

library(fitR)
library(coda)
library(dplyr)
library(reshape2)

source(here::here("Model", "transduction_model_functions.R"))

# mass_model = readRDS(here::here("Fitting", "mass_model.rds"))
model = readRDS(here::here("Model", "transduction_model.rds"))

fitted_params4 = c(readRDS(here::here("Fitting", "Full_chains", "best_params_transduction.rds")),
                   readRDS(here::here("Fitting", "Full_chains", "best_params_transduction2.rds")),
                   readRDS(here::here("Fitting", "Full_chains", "best_params_transduction3.rds")))
fitted_params4b = c(readRDS(here::here("Fitting", "Full_chains", "best_params_transduction_b.rds")),
                   readRDS(here::here("Fitting", "Full_chains", "best_params_transduction2_b.rds")),
                   readRDS(here::here("Fitting", "Full_chains", "best_params_transduction3_b.rds")))


lab_data_trans = read.csv(here::here("Lab", "Transduction", "summary_10_4.csv")) %>%
  select(Time, Bacteria, Mean) %>%
  dcast(Time~Bacteria) %>%
  select(-Total) %>%
  rename(time = Time, Bet = DRP, Be = EryR, Bt = TetR) %>%
  mutate(Be = round(Be),
         Bt = round(Bt),
         Bet = round(Bet),
         P = round(P)) 


## GO #####


models_to_try = data.frame(model_name="mass_decay_link_L", frequentist=FALSE,
                           fixed_delay=NA, decay=TRUE,
                           link_beta=FALSE, link_L=TRUE, link_delay=FALSE, transduction=TRUE)
models_to_try = rbind(models_to_try,
                      data.frame(model_name="mass_decay_link_beta", frequentist=FALSE,
                                 fixed_delay=NA, decay=TRUE,
                                 link_beta=TRUE, link_L=FALSE, link_delay=FALSE, transduction=TRUE))
models_to_try = rbind(models_to_try,
                      data.frame(model_name="mass_decay_link_both", frequentist=FALSE,
                           fixed_delay=NA, decay=TRUE,
                           link_beta=TRUE, link_L=TRUE, link_delay=FALSE, transduction=TRUE))
models_to_try = rbind(models_to_try,
                      data.frame(model_name="frequentist_decay_link_beta", frequentist=TRUE,
                                 fixed_delay=NA, decay=TRUE,
                                 link_beta=TRUE, link_L=FALSE, link_delay=FALSE, transduction=TRUE))
models_to_try = rbind(models_to_try,
                      data.frame(model_name="frequentist_decay_link_both", frequentist=TRUE,
                           fixed_delay=NA, decay=TRUE,
                           link_beta=TRUE, link_L=TRUE, link_delay=FALSE, transduction=TRUE))
models_to_try = rbind(models_to_try,
                      data.frame(model_name="frequentist_decay_link_L", frequentist=TRUE,
                                 fixed_delay=NA, decay=TRUE,
                                 link_beta=FALSE, link_L=TRUE, link_delay=FALSE, transduction=TRUE))


DIC = function(trace, model, lab_data){
  
  theta.bar <- colMeans(trace[,model$theta.names])
  
  init.state <- c(Be = lab_data$Be[1], Bt = lab_data$Bt[1], Bet = 0,
                  Pl = lab_data$P[1], Pt = 0, Pe = 0)
  log.like.theta.bar <- dTrajObs(model, theta.bar, init.state, data = lab_data, 
                                 log = TRUE)
  
  D.theta.bar <- -2 * log.like.theta.bar
  
  p.D <- var(-2 * trace[,"log.density"])/2
  
  DIC <- D.theta.bar + 2 * p.D
  
  DIC
  
}

summary_dic = data.frame()

for(i in 1:nrow(models_to_try)){
  
  models_to_try$model_name[i]
  
  trace_model4 = fitted_params4[[models_to_try$model_name[i]]]
  trace_model4 = coda::mcmc(trace_model4)
  trace_model4b = fitted_params4b[[models_to_try$model_name[i]]][-c(300001:400000),]
  trace_model4b = coda::mcmc(trace_model4b)
  trace_model4 = mcmc.list(trace_model4, trace_model4b)
  trace_model4 = burnAndThin(trace_model4, burn = 100000, thin = 10)
  plot(trace_model4)
  gelman.diag(trace_model4)
  
  trace_model4 = fitted_params4[[models_to_try$model_name[i]]]
  trace_model4 = coda::mcmc(trace_model4)
  trace_model4 = burnAndThin(trace_model4, burn = 200000, thin = 1)
  
  trace_model4b = fitted_params4b[[models_to_try$model_name[i]]][-c(300001:400000),]
  trace_model4b = coda::mcmc(trace_model4b)
  trace_model4b = burnAndThin(trace_model4b, burn = 200000, thin = 1)
  
  trace_model4 = rbind(trace_model4, trace_model4b)
  write.csv(trace_model4, "params_dens_burst.csv", row.names = F)

  tt = readRDS(here::here("Fitting", "Full_chains", "best_params_transduction_b.rds"))
  tt$mass_decay_link_L = tt$mass_decay_link_L[-c(1:200000),]
  tt$frequentist_decay_link_L = tt$frequentist_decay_link_L[-c(1:200000),]
  saveRDS(tt, here::here("Fitting", "Full_chains", "best_params_transduction_b.rds"))
  
  model = choose_model(model,
                       frequentist = models_to_try$frequentist[i],
                       fixed_delay = models_to_try$fixed_delay[i],
                       decay = models_to_try$decay[i], 
                       link_beta = models_to_try$link_beta[i],
                       link_L = models_to_try$link_L[i], 
                       link_delay = models_to_try$link_delay[i],
                       transduction = models_to_try$transduction[i])
  
  dic4 = DIC(trace_model4, model, lab_data_trans)

  summary_dic = rbind(summary_dic, 
                      c(models_to_try$model_name[i], "10_4", dic4))
  
}

colnames(summary_dic) = c("model_name", "init_pha", "DIC")
summary_dic$DIC = as.numeric(summary_dic$DIC)
write.csv(summary_dic, here::here("Fitting", "DIC_values.csv"), row.names = F)

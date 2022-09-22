
# SETUP ####

library(BayesianTools)
library(dplyr)
library(reshape2)

source(here::here("Model", "model.R"))

files = list.files(here::here("Fitting", "Fitted_params"))

params = vector("list", length(files))
i=1

for(f in files){
  params[[i]] = read.csv(here::here("Fitting", "Fitted_params", f))
  i = i+1
}

names(params) = gsub(".csv", "", gsub("params_", "", files))


lab_data_transM = read.csv(here::here("Data", "transduction_summary_10_4.csv")) %>%
  select(Time, Bacteria, Mean, se) %>%
  filter(Bacteria != "Total")
lab_data_trans5M = read.csv(here::here("Data", "transduction_summary_10_5.csv")) %>%
  select(Time, Bacteria, Mean, se)%>%
  filter(Bacteria != "Total")
lab_data_trans3M = read.csv(here::here("Data", "transduction_summary_10_3.csv")) %>%
  select(Time, Bacteria, Mean, se) %>%
  filter(Bacteria != "Total")


lab_data_trans = read.csv(here::here("Data", "transduction_summary_10_4.csv")) %>%
  select(Time, Bacteria, Mean) %>%
  dcast(Time~Bacteria) %>%
  select(-Total) %>%
  rename(time = Time, Bet = DRP, Be = EryR, Bt = TetR) %>%
  mutate(Be = round(Be),
         Bt = round(Bt),
         Bet = round(Bet),
         Pl = round(P)) 

lab_data_trans5 = read.csv(here::here("Data", "transduction_summary_10_5.csv")) %>%
  select(Time, Bacteria, Mean) %>%
  dcast(Time~Bacteria) %>%
  select(-Total) %>%
  rename(time = Time, Bet = DRP, Be = EryR, Bt = TetR) %>%
  mutate(Be = round(Be),
         Bt = round(Bt),
         Bet = round(Bet),
         Pl = round(P)) 

lab_data_trans3 = read.csv(here::here("Data", "transduction_summary_10_3.csv")) %>%
  select(Time, Bacteria, Mean) %>%
  dcast(Time~Bacteria) %>%
  select(-Total) %>%
  rename(time = Time, Bet = DRP, Be = EryR, Bt = TetR) %>%
  mutate(Be = round(Be),
         Bt = round(Bt),
         Bet = round(Bet),
         Pl = round(P))

DIC = function(trace, mode, link_L, link_beta, lab_data, lab_data2){
  
  theta.bar <- colMeans(trace[,-which(colnames(trace) == "LL")])
  log.like.theta.bar <- likelihood(lab_data, theta.bar, mode, link_L, link_beta) +
    likelihood(lab_data2, theta.bar, mode, link_L, link_beta)
  
  Dhat = -2*log.like.theta.bar
  
  Dbar = mean(-2*trace[,"LL"])
  pD = Dbar - Dhat
  DIC = Dbar + pD
  
  DIC
  
}


## GO #####


trace_model = params[["hill_both"]]

for(j in 1:nrow(trace_model)){
  
  if(j %% 100 == 0) cat("Row", j, "\n")
  
  theta_j = trace_model[j,-6]
  
  trace_model$LL[j] = likelihood(lab_data_trans3, theta_j,
                                 mode = "hill",
                                 link_L = T,
                                 link_beta = T) +
    likelihood(lab_data_trans5, theta_j,
               mode = "hill",
               link_L = T,
               link_beta = T)
  
}

write.csv(trace_model, here::here("Fitting", "params_hill_both.csv"), row.names = F)

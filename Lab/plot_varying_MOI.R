

library(fitR)
library(dplyr)
library(reshape2)
library(ggplot2)
library(scales)
library(openxlsx)

source(here::here("Model", "transduction_model_functions.R"))

# mass_model = readRDS(here::here("Fitting", "mass_model.rds"))
model = readRDS(here::here("Model", "transduction_model.rds"))
fitted_params = c(readRDS(here::here("Fitting", "10_4", "best_params_transduction.rds")),
                  readRDS(here::here("Fitting", "10_4", "best_params_transduction2.rds")),
                  readRDS(here::here("Fitting", "10_4", "best_params_transduction3.rds")))
fitted_paramsb = c(readRDS(here::here("Fitting", "10_4", "best_params_transduction_b.rds")),
                  readRDS(here::here("Fitting", "10_4", "best_params_transduction2_b.rds")),
                  readRDS(here::here("Fitting", "10_4", "best_params_transduction3_b.rds")))


data = read.xlsx(here::here("Lab", "Varying_MOI", "jake_data3.xlsx"))
#data = read.xlsx(here::here("Lab", "Varying_MOI", "jake_data2.xlsx"))

data = data[c(2:25),c(1:29)]
data$se_bac = apply(data[,c(18,19,20)], 1, sd)

drp = data[,c(26:29)] %>%
  na.omit(.) %>%
  mutate(drp_se = apply(.[,c(1,2,3)], 1, sd))


data = data[,c(7,8,11,21,29,30)]

colnames(data) = c("init_bac", "init_pha", "plate", "bac_24h", "drp_24h", "bac_se")

for(i in 2:nrow(data)){
  if(is.na(data$init_pha[i])){
    data$init_pha[i] = data$init_pha[i-1]
  }
}

data$init_bac = data$init_bac[1]

data_se = data %>%
  dcast(., init_bac+init_pha~plate, value.var = "bac_se") %>%
  select("E", "T")
colnames(data_se) = c("Be_se", "Bt_se")

data = data %>%
  dcast(., init_bac+init_pha~plate, value.var = "bac_24h")

data$Bet = rev(drp$`Mean.DRP/mL`)
data$Bet_se = rev(drp$drp_se)

data = cbind(data, data_se)

data = as.data.frame(apply(data, c(1,2), as.numeric))[,-3]
colnames(data)[c(3,4)] = c("Be", "Bt")

data_pha = read.xlsx(here::here("Lab", "Varying_MOI", "jake_data3.xlsx"),sheet = 2)[c(2:9),]
#data_pha = read.xlsx(here::here("Lab", "Varying_MOI", "jake_data2.xlsx"),sheet = 2)[c(2:9),]

data_pha$se_pha = apply(data_pha[,c(8,9,10)], 1, sd)

data$Pl = rev(data_pha$`Mean.PFU/mL`)
data$Pl_se = rev(data_pha$se_pha)


models_to_try = data.frame(model_name="mass_decay_link_L", frequentist=FALSE,
                           delay=TRUE, 
                           fixed_delay=NA, decay=TRUE,
                           link_beta=FALSE, link_L=TRUE, link_delay=FALSE, transduction=TRUE)
models_to_try = rbind(models_to_try,
                      data.frame(model_name="mass_decay_link_beta", frequentist=FALSE,
                                 delay=TRUE,
                                 fixed_delay=NA, decay=TRUE,
                                 link_beta=TRUE, link_L=FALSE, link_delay=FALSE, transduction=TRUE))
models_to_try = rbind(models_to_try,
                      data.frame(model_name="mass_decay_link_both", frequentist=FALSE,
                                 delay=TRUE, 
                                 fixed_delay=NA, decay=TRUE,
                                 link_beta=TRUE, link_L=TRUE, link_delay=FALSE, transduction=TRUE))
models_to_try = rbind(models_to_try,
                      data.frame(model_name="frequentist_decay_link_beta", frequentist=TRUE,
                                 delay=TRUE,
                                 fixed_delay=NA, decay=TRUE,
                                 link_beta=TRUE, link_L=FALSE, link_delay=FALSE, transduction=TRUE))
models_to_try = rbind(models_to_try,
                      data.frame(model_name="frequentist_decay_link_both", frequentist=TRUE,
                                 delay=TRUE, 
                                 fixed_delay=NA, decay=TRUE,
                                 link_beta=TRUE, link_L=TRUE, link_delay=FALSE, transduction=TRUE))
models_to_try = rbind(models_to_try,
                      data.frame(model_name="frequentist_decay_link_L", frequentist=TRUE,
                                 delay=TRUE, 
                                 fixed_delay=NA, decay=TRUE,
                                 link_beta=FALSE, link_L=TRUE, link_delay=FALSE, transduction=TRUE))

## REPEAT WITH OTHER MODELS FOR SUPP MAT

data$model = "data"
all_results = data

for(i in 1:nrow(models_to_try)){
  
  cat("\nWorking on", models_to_try$model_name[i])
  
  model = choose_model(model,
                       frequentist = models_to_try$frequentist[i],
                       delay = models_to_try$delay[i],
                       fixed_delay = models_to_try$fixed_delay[i],
                       decay = models_to_try$decay[i], 
                       link_beta = models_to_try$link_beta[i],
                       link_L = models_to_try$link_L[i], 
                       link_delay = models_to_try$link_delay[i],
                       transduction = models_to_try$transduction[i])
  
  trace_model = fitted_params[[models_to_try$model_name[i]]]
  trace_model = coda::mcmc(trace_model)
  trace_model = burnAndThin(trace_model, burn = 20000, thin = 10)
  
  trace_modelb = fitted_paramsb[[models_to_try$model_name[i]]]
  trace_modelb = coda::mcmc(trace_modelb)
  trace_modelb = burnAndThin(trace_modelb, burn = 20000, thin = 10)
  
  trace_model = rbind(trace_model, trace_modelb)
  
  params = apply(trace_model, 2, median)

  data_model = data
  
  for(j in 1:nrow(data_model)){
    
    init.state = c(Be = data_model$init_bac[j]/2, Bt = data_model$init_bac[j]/2, Bet = 0,
                   Pl = data_model$init_pha[j], Pe = 0, Pt = 0)
    
    traj = multi_run2(model, params, init.state, times = seq(0, 24, 1), nruns = 10)
    
    data_model$Be[j] = traj$Be[25]
    data_model$Bt[j] = traj$Bt[25]
    data_model$Bet[j] = traj$Bet[25]
    data_model$Pl[j] = traj$Pl[25]
    
    data_model$Be_se[j] = traj$Be_sd[25]
    data_model$Bt_se[j] = traj$Bt_sd[25]
    data_model$Bet_se[j] = traj$Bet_sd[25]
    data_model$Pl_se[j] = traj$Pl_sd[25]
    

  }
  
  data_model$model = rep(models_to_try$model_name[i], nrow(data_model))
  
  all_results = rbind(all_results, data_model)

}

all_results[all_results<0.01] = 0.01
all_results[is.na(all_results)] = 0.01

all_results$init_pha = as.factor(format(unique(all_results$init_pha), scientific = T, digits = 2))
#all_results$init_pha = factor(all_results$init_pha, levels(all_results$init_pha)[c(1,3,5,7,2,4,6,8)])
all_results$init_pha = factor(all_results$init_pha, levels(all_results$init_pha)[c(6,1,3,5,7,2,4,8)])

all_results_L = all_results %>%
  filter(model %in% c("frequentist_decay_link_L", "data", "mass_decay_link_L"))
all_results_beta = all_results %>%
  filter(model %in% c("frequentist_decay_link_beta", "data", "mass_decay_link_beta"))
all_results_both = all_results %>%
  filter(model %in% c("frequentist_decay_link_both", "data", "mass_decay_link_both"))

all_results_L$model = as.factor(all_results_L$model)
all_results_L$model = factor(all_results_L$model, levels(all_results_L$model)[c(3,1,2)])
levels(all_results_L$model) = c("Density model", "Data", "Frequency model")

all_results_beta$model = as.factor(all_results_beta$model)
all_results_beta$model = factor(all_results_beta$model, levels(all_results_beta$model)[c(3,1,2)])
levels(all_results_beta$model) = c("Density model", "Data", "Frequency model")

all_results_both$model = as.factor(all_results_both$model)
all_results_both$model = factor(all_results_both$model, levels(all_results_both$model)[c(3,1,2)])
levels(all_results_both$model) = c("Density model", "Data", "Frequency model")


ggplot(all_results_L) + 
  geom_point(aes(init_pha, Be, colour = "Be"), size = 2.5) +
  geom_line(stat = "smooth", aes(as.numeric(init_pha), Be, colour = "Be"), alpha = 0.1, size = 5) +
  geom_errorbar(aes(init_pha, ymin = pmax(Be - Be_se, 0.01), ymax = Be + Be_se, colour = "Be"), size = 0.8) +
  geom_point(aes(init_pha, Bt, colour = "Bt"), size = 2.5) +
  geom_line(stat = "smooth", aes(as.numeric(init_pha), Bt, colour = "Bt"), alpha = 0.1, size = 5) +
  geom_errorbar(aes(init_pha, ymin = pmax(Bt - Bt_se, 0.01), ymax = Bt + Bt_se, colour = "Bt"), size = 0.8) +
  geom_point(aes(init_pha, Bet, colour = "Bet"), size = 2.5) +
  geom_line(stat = "smooth", aes(as.numeric(init_pha), Bet, colour = "Bet"), alpha = 0.1, size = 5) +
  geom_errorbar(aes(init_pha, ymin = pmax(Bet - Bet_se, 0.01), ymax = Bet + Bet_se, colour = "Bet"), size = 0.8) +
  geom_point(aes(init_pha, Pl, colour = "Pl"), size = 2.5, shape = 17) +
  geom_line(stat = "smooth", aes(as.numeric(init_pha), Pl, colour = "Pl"), alpha = 0.1, size = 5) +
  geom_errorbar(aes(init_pha, ymin = pmax(Pl - Pl_se, 0.01), ymax = Pl + Pl_se, colour = "Pl"), size = 0.8) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x, n = 6),
                     labels=c("x", "0", expression(10^0),
                              expression(10^2), expression(10^4), expression(10^6),
                              expression(10^8), expression(10^10), expression(10^12), "x")) +
  coord_cartesian(ylim = c(0.01, 1e12)) +
  facet_grid(~model) +
  labs(colour = "", x="Initial phage (pfu per mL)", y="cfu or pfu per mL") +
  theme_bw() +
  scale_colour_manual(breaks= c("Be", "Bt", "Bet", "Pl"),values=c("#685cc4","#6db356","#c2484d","#c88a33"),
                      labels = c(expression(B[E]),
                                 expression(B[T]),
                                 expression(B[ET]),
                                 expression(P[L])),
                      guide = guide_legend(override.aes = list(shape = c(19,19,19,17)))) +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text.x = element_text(size=12))

ggsave(here::here("Lab", "Plots", "varying_MOI_L.png"))

ggplot(all_results_beta) + 
  geom_point(aes(init_pha, Be, colour = "Be"), size = 2.5) +
  geom_line(stat = "smooth", aes(as.numeric(init_pha), Be, colour = "Be"), alpha = 0.1, size = 5) +
  geom_errorbar(aes(init_pha, ymin = pmax(Be - Be_se, 0.01), ymax = Be + Be_se, colour = "Be"), size = 0.8) +
  geom_point(aes(init_pha, Bt, colour = "Bt"), size = 2.5) +
  geom_line(stat = "smooth", aes(as.numeric(init_pha), Bt, colour = "Bt"), alpha = 0.1, size = 5) +
  geom_errorbar(aes(init_pha, ymin = pmax(Bt - Bt_se, 0.01), ymax = Bt + Bt_se, colour = "Bt"), size = 0.8) +
  geom_point(aes(init_pha, Bet, colour = "Bet"), size = 2.5) +
  geom_line(stat = "smooth", aes(as.numeric(init_pha), Bet, colour = "Bet"), alpha = 0.1, size = 5) +
  geom_errorbar(aes(init_pha, ymin = pmax(Bet - Bet_se, 0.01), ymax = Bet + Bet_se, colour = "Bet"), size = 0.8) +
  geom_point(aes(init_pha, Pl, colour = "Pl"), size = 2.5, shape = 17) +
  geom_line(stat = "smooth", aes(as.numeric(init_pha), Pl, colour = "Pl"), alpha = 0.1, size = 5) +
  geom_errorbar(aes(init_pha, ymin = pmax(Pl - Pl_se, 0.01), ymax = Pl + Pl_se, colour = "Pl"), size = 0.8) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x, n = 6),
                     labels=c("x", "0", expression(10^0),
                              expression(10^2), expression(10^4), expression(10^6),
                              expression(10^8), expression(10^10), expression(10^12), "x")) +
  coord_cartesian(ylim = c(0.01, 1e12)) +
  facet_grid(~model) +
  labs(colour = "", x="Initial phage (pfu per mL)", y="cfu or pfu per mL") +
  theme_bw() +
  scale_colour_manual(breaks= c("Be", "Bt", "Bet", "Pl"),values=c("#685cc4","#6db356","#c2484d","#c88a33"),
                      labels = c(expression(B[E]),
                                 expression(B[T]),
                                 expression(B[ET]),
                                 expression(P[L])),
                      guide = guide_legend(override.aes = list(shape = c(19,19,19,17)))) +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text.x = element_text(size=12))

ggsave(here::here("Lab", "Plots", "varying_MOI_beta.png"))

ggplot(all_results_both) + 
  geom_point(aes(init_pha, Be, colour = "Be"), size = 2.5) +
  geom_line(stat = "smooth", aes(as.numeric(init_pha), Be, colour = "Be"), alpha = 0.1, size = 5) +
  geom_errorbar(aes(init_pha, ymin = pmax(Be - Be_se, 0.01), ymax = Be + Be_se, colour = "Be"), size = 0.8) +
  geom_point(aes(init_pha, Bt, colour = "Bt"), size = 2.5) +
  geom_line(stat = "smooth", aes(as.numeric(init_pha), Bt, colour = "Bt"), alpha = 0.1, size = 5) +
  geom_errorbar(aes(init_pha, ymin = pmax(Bt - Bt_se, 0.01), ymax = Bt + Bt_se, colour = "Bt"), size = 0.8) +
  geom_point(aes(init_pha, Bet, colour = "Bet"), size = 2.5) +
  geom_line(stat = "smooth", aes(as.numeric(init_pha), Bet, colour = "Bet"), alpha = 0.1, size = 5) +
  geom_errorbar(aes(init_pha, ymin = pmax(Bet - Bet_se, 0.01), ymax = Bet + Bet_se, colour = "Bet"), size = 0.8) +
  geom_point(aes(init_pha, Pl, colour = "Pl"), size = 2.5, shape = 17) +
  geom_line(stat = "smooth", aes(as.numeric(init_pha), Pl, colour = "Pl"), alpha = 0.1, size = 5) +
  geom_errorbar(aes(init_pha, ymin = pmax(Pl - Pl_se, 0.01), ymax = Pl + Pl_se, colour = "Pl"), size = 0.8) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x, n = 6),
                     labels=c("x", "0", expression(10^0),
                              expression(10^2), expression(10^4), expression(10^6),
                              expression(10^8), expression(10^10), expression(10^12), "x")) +
  coord_cartesian(ylim = c(0.01, 1e12)) +
  facet_grid(~model) +
  labs(colour = "", x="Initial phage (pfu per mL)", y="cfu or pfu per mL") +
  theme_bw() +
  scale_colour_manual(breaks= c("Be", "Bt", "Bet", "Pl"),values=c("#685cc4","#6db356","#c2484d","#c88a33"),
                      labels = c(expression(B[E]),
                                 expression(B[T]),
                                 expression(B[ET]),
                                 expression(P[L])),
                      guide = guide_legend(override.aes = list(shape = c(19,19,19,17)))) +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text.x = element_text(size=12))

ggsave(here::here("Lab", "Plots", "varying_MOI_both.png"))

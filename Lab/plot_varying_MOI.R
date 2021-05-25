

library(fitR)
library(dplyr)
library(reshape2)
library(ggplot2)
library(scales)
library(openxlsx)

source(here::here("Model", "transduction_model_functions.R"))

# mass_model = readRDS(here::here("Fitting", "mass_model.rds"))
model = readRDS(here::here("Model", "transduction_model.rds"))

files = list.files(here::here("Fitting", "Fitted_params"))

all_params = vector("list", length(files))
i=1

for(f in files){
  all_params[[i]] = read.csv(here::here("Fitting", "Fitted_params", f))
  i = i+1
}

names(all_params) = gsub(".csv", "", gsub("params_", "", files))


data = read.xlsx(here::here("Lab", "Varying_MOI", "jake_data_final.xlsx"))

data = data[c(2:34),c(1:29)]
data$se_bac = apply(data[,c(18,19,20)], 1, function(x) sd(x)/sqrt(3))

drp = data[,c(26:29)] %>%
  na.omit(.) %>%
  mutate(drp_sd = apply(.[,c(1,2,3)], 1, function(x) sd(x)/sqrt(3)))


data = data[,c(7,8,11,21,29,30)]

colnames(data) = c("init_bac", "init_pha", "plate", "bac_24h", "drp_24h", "bac_sd")

for(i in 2:nrow(data)){
  if(is.na(data$init_pha[i])){
    data$init_pha[i] = data$init_pha[i-1]
  }
}

data$init_bac = data$init_bac[1]

data_sd = data %>%
  dcast(., init_bac+init_pha~plate, value.var = "bac_sd") %>%
  select("E", "T")
colnames(data_sd) = c("Be_sd", "Bt_sd")

data = data %>%
  dcast(., init_bac+init_pha~plate, value.var = "bac_24h")

data$Bet = rev(drp$`Mean.DRP/mL`)
data$Bet_sd = rev(drp$drp_sd)

data = cbind(data, data_sd)

data = as.data.frame(apply(data, c(1,2), as.numeric))[,-3]
colnames(data)[c(3,4)] = c("Be", "Bt")

data_pha = read.xlsx(here::here("Lab", "Varying_MOI", "jake_data_final.xlsx"),sheet = 2)[c(2:12),]

data_pha$se_pha = apply(data_pha[,c(8,9,10)], 1, function(x) sd(x)/sqrt(3))

data$Pl = rev(data_pha$`Mean.PFU/mL`)
data$Pl_sd = rev(data_pha$se_pha)


models_to_try = data.frame(model_name="dens_burst", frequentist=FALSE,
                           fixed_delay=NA, decay=TRUE,
                           link_beta=FALSE, link_L=TRUE, link_delay=FALSE, transduction=TRUE)
models_to_try = rbind(models_to_try,
                      data.frame(model_name="dens_beta", frequentist=FALSE,
                                 fixed_delay=NA, decay=TRUE,
                                 link_beta=TRUE, link_L=FALSE, link_delay=FALSE, transduction=TRUE))
models_to_try = rbind(models_to_try,
                      data.frame(model_name="dens_both", frequentist=FALSE,
                                 fixed_delay=NA, decay=TRUE,
                                 link_beta=TRUE, link_L=TRUE, link_delay=FALSE, transduction=TRUE))
models_to_try = rbind(models_to_try,
                      data.frame(model_name="freq_beta", frequentist=TRUE,
                                 fixed_delay=NA, decay=TRUE,
                                 link_beta=TRUE, link_L=FALSE, link_delay=FALSE, transduction=TRUE))
models_to_try = rbind(models_to_try,
                      data.frame(model_name="freq_both", frequentist=TRUE,
                                 fixed_delay=NA, decay=TRUE,
                                 link_beta=TRUE, link_L=TRUE, link_delay=FALSE, transduction=TRUE))
models_to_try = rbind(models_to_try,
                      data.frame(model_name="freq_burst", frequentist=TRUE,
                                 fixed_delay=NA, decay=TRUE,
                                 link_beta=FALSE, link_L=TRUE, link_delay=FALSE, transduction=TRUE))

## REPEAT WITH OTHER MODELS FOR SUPP MAT

data$model = "data"
all_results = data

for(i in 1:nrow(models_to_try)){
  
  cat("\nWorking on", models_to_try$model_name[i])
  
  model = choose_model(model,
                       frequentist = models_to_try$frequentist[i],
                       fixed_delay = models_to_try$fixed_delay[i],
                       decay = models_to_try$decay[i], 
                       link_beta = models_to_try$link_beta[i],
                       link_L = models_to_try$link_L[i], 
                       link_delay = models_to_try$link_delay[i],
                       transduction = models_to_try$transduction[i])
  
  trace_model = all_params[[models_to_try$model_name[i]]]
  
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
    
    data_model$Be_sd[j] = traj$Be_sd[25]
    data_model$Bt_sd[j] = traj$Bt_sd[25]
    data_model$Bet_sd[j] = traj$Bet_sd[25]
    data_model$Pl_sd[j] = traj$Pl_sd[25]
    
    
  }
  
  data_model$model = rep(models_to_try$model_name[i], nrow(data_model))
  
  all_results = rbind(all_results, data_model)
  
}

all_results$Be[all_results$Be<0.01] = 0.01
all_results$Bt[all_results$Bt<0.01] = 0.01
all_results$Bet[all_results$Bet<0.01] = 0.01
all_results$Pl[all_results$Pl<0.01] = 0.01
all_results$Be[is.na(all_results$Be)] = 0.01
all_results$Bt[is.na(all_results$Bt)] = 0.01
all_results$Bet[is.na(all_results$Bet)] = 0.01
all_results$Pl[is.na(all_results$Pl)] = 0.01


all_results$Be_sd[all_results$Be_sd<0.01] = 0
all_results$Bt_sd[all_results$Bt_sd<0.01] = 0
all_results$Bet_sd[all_results$Bet_sd<0.01] = 0
all_results$Pl_sd[all_results$Pl_sd<0.01] = 0
all_results$Be_sd[is.na(all_results$Be_sd)] = 0
all_results$Bt_sd[is.na(all_results$Bt_sd)] = 0
all_results$Bet_sd[is.na(all_results$Bet_sd)] = 0
all_results$Pl_sd[is.na(all_results$Pl_sd)] = 0

all_results$init_pha = as.factor(format(unique(all_results$init_pha), scientific = T, digits = 2))
all_results$init_pha = factor(all_results$init_pha, levels(all_results$init_pha)[c(9,2,5,8,1,4,7,10,3,6,11)])

all_results_L = all_results %>%
  filter(model %in% c("freq_burst", "data", "dens_burst"))
all_results_beta = all_results %>%
  filter(model %in% c("freq_beta", "data", "dens_beta"))
all_results_both = all_results %>%
  filter(model %in% c("freq_both", "data", "dens_both"))

all_results_L$model = as.factor(all_results_L$model)
all_results_L$model = factor(all_results_L$model, levels(all_results_L$model)[c(3,1,2)])
levels(all_results_L$model) = c("Frequency model", "Data", "Density model")

all_results_beta$model = as.factor(all_results_beta$model)
all_results_beta$model = factor(all_results_beta$model, levels(all_results_beta$model)[c(3,1,2)])
levels(all_results_beta$model) = c("Frequency model", "Data", "Density model")

all_results_both$model = as.factor(all_results_both$model)
all_results_both$model = factor(all_results_both$model, levels(all_results_both$model)[c(3,1,2)])
levels(all_results_both$model) = c("Frequency model", "Data", "Density model")


ggplot(all_results_L) + 
  geom_pointrange(aes(x = init_pha, y = Be, ymin = pmax(Be - Be_sd, 0.01), ymax = Be + Be_sd, colour = "Be"),
                  position = position_jitter(height=0, width=0.1), alpha = 1, size = 0.5, linetype = "solid") +
  geom_pointrange(aes(x = init_pha, y = Bt, ymin = pmax(Bt - Bt_sd, 0.01), ymax = Bt + Bt_sd, colour = "Bt"),
                  position = position_jitter(height=0, width=0.1), alpha = 1, size = 0.5, linetype = "solid") +
  geom_pointrange(aes(x = init_pha, y = Bet, ymin = pmax(Bet - Bet_sd, 0.01), ymax = Bet + Bet_sd, colour = "Bet"),
                  position = position_jitter(height=0, width=0.1), alpha = 1, size = 0.5, linetype = "solid") +
  geom_pointrange(aes(x = init_pha, y = Pl, ymin = pmax(Pl - Pl_sd, 0.01), ymax = Pl + Pl_sd, colour = "Pl"),
                  position = position_jitter(height=0, width=0.1), alpha = 1, size = 0.5, linetype = "solid") +
  geom_line(stat = "smooth", aes(as.numeric(init_pha), Be, colour = "Be"), alpha = 0.1, size = 5) +
  geom_line(stat = "smooth", aes(as.numeric(init_pha), Bt, colour = "Bt"), alpha = 0.1, size = 5) +
  geom_line(stat = "smooth", aes(as.numeric(init_pha), Bet, colour = "Bet"), alpha = 0.1, size = 5) +
  geom_line(stat = "smooth", aes(as.numeric(init_pha), Pl, colour = "Pl"), alpha = 0.1, size = 5) +
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
  geom_pointrange(aes(x = init_pha, y = Be, ymin = pmax(Be - Be_sd, 0.01), ymax = Be + Be_sd, colour = "Be"),
                  position = position_jitter(height=0, width=0.1), alpha = 1, size = 0.5, linetype = "solid") +
  geom_pointrange(aes(x = init_pha, y = Bt, ymin = pmax(Bt - Bt_sd, 0.01), ymax = Bt + Bt_sd, colour = "Bt"),
                  position = position_jitter(height=0, width=0.1), alpha = 1, size = 0.5, linetype = "solid") +
  geom_pointrange(aes(x = init_pha, y = Bet, ymin = pmax(Bet - Bet_sd, 0.01), ymax = Bet + Bet_sd, colour = "Bet"),
                  position = position_jitter(height=0, width=0.1), alpha = 1, size = 0.5, linetype = "solid") +
  geom_pointrange(aes(x = init_pha, y = Pl, ymin = pmax(Pl - Pl_sd, 0.01), ymax = Pl + Pl_sd, colour = "Pl"),
                  position = position_jitter(height=0, width=0.1), alpha = 1, size = 0.5, linetype = "solid") +
  geom_line(stat = "smooth", aes(as.numeric(init_pha), Be, colour = "Be"), alpha = 0.1, size = 5) +
  geom_line(stat = "smooth", aes(as.numeric(init_pha), Bt, colour = "Bt"), alpha = 0.1, size = 5) +
  geom_line(stat = "smooth", aes(as.numeric(init_pha), Bet, colour = "Bet"), alpha = 0.1, size = 5) +
  geom_line(stat = "smooth", aes(as.numeric(init_pha), Pl, colour = "Pl"), alpha = 0.1, size = 5) +
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
  geom_pointrange(aes(x = init_pha, y = Be, ymin = pmax(Be - Be_sd, 0.01), ymax = Be + Be_sd, colour = "Be"),
                  position = position_jitter(height=0, width=0.1), alpha = 1, size = 0.5, linetype = "solid") +
  geom_pointrange(aes(x = init_pha, y = Bt, ymin = pmax(Bt - Bt_sd, 0.01), ymax = Bt + Bt_sd, colour = "Bt"),
                  position = position_jitter(height=0, width=0.1), alpha = 1, size = 0.5, linetype = "solid") +
  geom_pointrange(aes(x = init_pha, y = Bet, ymin = pmax(Bet - Bet_sd, 0.01), ymax = Bet + Bet_sd, colour = "Bet"),
                  position = position_jitter(height=0, width=0.1), alpha = 1, size = 0.5, linetype = "solid") +
  geom_pointrange(aes(x = init_pha, y = Pl, ymin = pmax(Pl - Pl_sd, 0.01), ymax = Pl + Pl_sd, colour = "Pl"),
                  position = position_jitter(height=0, width=0.1), alpha = 1, size = 0.5, linetype = "solid") +
  geom_line(stat = "smooth", aes(as.numeric(init_pha), Be, colour = "Be"), alpha = 0.1, size = 5) +
  geom_line(stat = "smooth", aes(as.numeric(init_pha), Bt, colour = "Bt"), alpha = 0.1, size = 5) +
  geom_line(stat = "smooth", aes(as.numeric(init_pha), Bet, colour = "Bet"), alpha = 0.1, size = 5) +
  geom_line(stat = "smooth", aes(as.numeric(init_pha), Pl, colour = "Pl"), alpha = 0.1, size = 5) +
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

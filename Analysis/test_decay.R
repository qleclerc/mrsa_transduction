

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

all_results3 = c()

try_decay = c(10^-4, 10^-3.5, 10^-3, 10^-2.5, 10^-2, 10^-1.5, 10^-1)

for(i in 1:nrow(models_to_try)){
  
  cat("\nWorking on", models_to_try$model_name[i])
  model_result = data.frame(model = models_to_try$model_name[i], decay = try_decay, Pl = 0, Pl_sd = 0)
  
  for(j in 1:nrow(model_result)){
    
    model = choose_model(model,
                         frequentist = models_to_try$frequentist[i],
                         fixed_delay = models_to_try$fixed_delay[i],
                         decay = model_result$decay[j], 
                         link_beta = models_to_try$link_beta[i],
                         link_L = models_to_try$link_L[i], 
                         link_delay = models_to_try$link_delay[i],
                         transduction = models_to_try$transduction[i])
    
    trace_model = all_params[[models_to_try$model_name[i]]]
    
    params = apply(trace_model, 2, median)
    
    init.state = c(Be = 1e4, Bt = 1e4, Bet = 0,
                   Pl = 1e3, Pe = 0, Pt = 0)
    
    traj = multi_run2(model, params, init.state, times = seq(0, 24, 1), nruns = 10)
    
    model_result$Pl[j] = traj$Pl[25]
    model_result$Pl_sd[j] = traj$Pl_sd[25]
    
  }
  
  all_results3 = rbind(all_results3, model_result)
  
}  

all_results3$initial = "10^3"


all_results4 = c()

try_decay = c(10^-4, 10^-3.5, 10^-3, 10^-2.5, 10^-2, 10^-1.5, 10^-1)

for(i in 1:nrow(models_to_try)){
  
  cat("\nWorking on", models_to_try$model_name[i])
  model_result = data.frame(model = models_to_try$model_name[i], decay = try_decay, Pl = 0, Pl_sd = 0)
  
  for(j in 1:nrow(model_result)){
    
    model = choose_model(model,
                         frequentist = models_to_try$frequentist[i],
                         fixed_delay = models_to_try$fixed_delay[i],
                         decay = model_result$decay[j], 
                         link_beta = models_to_try$link_beta[i],
                         link_L = models_to_try$link_L[i], 
                         link_delay = models_to_try$link_delay[i],
                         transduction = models_to_try$transduction[i])
    
    trace_model = all_params[[models_to_try$model_name[i]]]
    
    params = apply(trace_model, 2, median)
    
    init.state = c(Be = 1e4, Bt = 1e4, Bet = 0,
                   Pl = 1e4, Pe = 0, Pt = 0)
    
    traj = multi_run2(model, params, init.state, times = seq(0, 24, 1), nruns = 10)
    
    model_result$Pl[j] = traj$Pl[25]
    model_result$Pl_sd[j] = traj$Pl_sd[25]
    
  }
  
  all_results4 = rbind(all_results4, model_result)
  
}  

all_results4$initial = "10^4"


all_results5 = c()

try_decay = c(10^-4, 10^-3.5, 10^-3, 10^-2.5, 10^-2, 10^-1.5, 10^-1)

for(i in 1:nrow(models_to_try)){
  
  cat("\nWorking on", models_to_try$model_name[i])
  model_result = data.frame(model = models_to_try$model_name[i], decay = try_decay, Pl = 0, Pl_sd = 0)
  
  for(j in 1:nrow(model_result)){
    
    model = choose_model(model,
                         frequentist = models_to_try$frequentist[i],
                         fixed_delay = models_to_try$fixed_delay[i],
                         decay = model_result$decay[j], 
                         link_beta = models_to_try$link_beta[i],
                         link_L = models_to_try$link_L[i], 
                         link_delay = models_to_try$link_delay[i],
                         transduction = models_to_try$transduction[i])
    
    trace_model = all_params[[models_to_try$model_name[i]]]
    
    params = apply(trace_model, 2, median)
    
    init.state = c(Be = 1e4, Bt = 1e4, Bet = 0,
                   Pl = 1e5, Pe = 0, Pt = 0)
    
    traj = multi_run2(model, params, init.state, times = seq(0, 24, 1), nruns = 10)
    
    model_result$Pl[j] = traj$Pl[25]
    model_result$Pl_sd[j] = traj$Pl_sd[25]
    
  }
  
  all_results5 = rbind(all_results5, model_result)
  
}  

all_results5$initial = "10^5"


all_results = rbind(all_results3, all_results4, all_results5)


ggplot(all_results) + 
  geom_pointrange(aes(x = log10(decay), y = Pl, ymin = pmax(Pl - Pl_sd, 0.01), ymax = Pl + Pl_sd, colour = model),
                  alpha = 0.6, size = 0.5, linetype = "solid", position = position_dodge(width=0.2)) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x, n = 6),
                     labels=trans_format("log10", math_format(10^.x))) +
  # scale_x_continuous(trans=log10_trans(),
  #                    breaks=trans_breaks("log10", function(x) 10^x, n = 6),
  #                    labels=trans_format("log10", math_format(10^.x))) +
  coord_cartesian(ylim = c(1e4, 1e8)) +
  facet_wrap(~initial) +
  labs(colour = "Model:", x="Log phage decay rate (phage per h)", y="pfu per mL after 24h") +
  theme_bw() +
  scale_colour_manual(breaks= c("dens_beta", "dens_burst", "dens_both", "freq_beta", "freq_burst", "freq_both"),
                      values = c("#579dd9", "#5383ea", "#2c56a7", "firebrick1", "firebrick3", "firebrick4"),
                      labels = c("Density - adsorption",
                                 "Density - burst",
                                 "Density - both",
                                 "Frequency - adsorption",
                                 "Frequency - burst",
                                 "Frequency - both")) +
  theme(axis.text.x = element_text(size=12, angle = 0, hjust = 1),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.x = element_text(size=12))

ggsave(here::here("Lab", "Plots", "decay_impact.png"))

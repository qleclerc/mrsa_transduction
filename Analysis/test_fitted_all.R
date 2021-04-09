# SETUP ####

library(fitR)
library(coda)
library(dplyr)
library(reshape2)
library(ggplot2)
library(cowplot)
library(scales)

source(here::here("Model", "transduction_model_functions.R"))

# mass_model = readRDS(here::here("Fitting", "mass_model.rds"))
model = readRDS(here::here("Model", "transduction_model.rds"))

fitted_params3 = c(readRDS(here::here("Fitting", "10_3", "best_params_transduction.rds")),
                   readRDS(here::here("Fitting", "10_3", "best_params_transduction2.rds")),
                   readRDS(here::here("Fitting", "10_3", "best_params_transduction3.rds")))
fitted_params4 = c(readRDS(here::here("Fitting", "10_4", "best_params_transduction.rds")),
                   readRDS(here::here("Fitting", "10_4", "best_params_transduction2.rds")),
                   readRDS(here::here("Fitting", "10_4", "best_params_transduction3.rds")))
fitted_params5 = c(readRDS(here::here("Fitting", "10_5", "best_params_transduction.rds")),
                   readRDS(here::here("Fitting", "10_5", "best_params_transduction2.rds")),
                   readRDS(here::here("Fitting", "10_5", "best_params_transduction3.rds")))

# lab_dataM = read.csv(here::here("Lab", "Triculture", "summary.csv")) %>%
#   select(Time, Bacteria, Mean, se) %>%
#   #mutate(se = se*sqrt(3)) %>%
#   filter(Bacteria != "Total")
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


# lab_data = read.csv(here::here("Lab", "Triculture", "summary.csv")) %>%
#   select(Time, Bacteria, Mean) %>%
#   dcast(Time~Bacteria) %>%
#   select(-Total) %>%
#   rename(time = Time, Bet = DRP, Be = EryR, Bt = TetR) %>%
#   mutate(Be = round(Be),
#          Bt = round(Bt),
#          Bet = round(Bet))

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


## GO #####


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


best_params = data.frame()

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
  
  # trace_model3 = fitted_params3[[models_to_try$model_name[i]]]
  # trace_model3 = trace_model3[which.max(trace_model3[,"log.density"]),]
  # trace_model3 = coda::mcmc(trace_model3)
  # trace_model3 = burnAndThin(trace_model3, burn = 20000, thin = 10)
  
  trace_model4 = fitted_params4[[models_to_try$model_name[i]]]
  #trace_model4 = trace_model4[which.max(trace_model4[,"log.density"]),]
  trace_model4 = coda::mcmc(trace_model4)
  trace_model4 = burnAndThin(trace_model4, burn = 20000, thin = 10)

  # trace_model5 = fitted_params5[[models_to_try$model_name[i]]]
  # trace_model5 = trace_model5[which.max(trace_model5[,"log.density"]),]
  # trace_model5 = coda::mcmc(trace_model5)
  # trace_model5 = burnAndThin(trace_model5, burn = 20000, thin = 10)
  
  # best_params = rbind(best_params, 
  #                     c(models_to_try$model_name[i], "10_3", trace_model3[which.max(trace_model3[,"log.density"]),]),
  #                     c(models_to_try$model_name[i], "10_4", trace_model4[which.max(trace_model4[,"log.density"]),]),
  #                     c(models_to_try$model_name[i], "10_5", trace_model5[which.max(trace_model5[,"log.density"]),]))
  
  best_params = rbind(best_params,
                      cbind(rep(models_to_try$model_name[i], 3),
                            c(0.025, 0.5, 0.975),
                            apply(trace_model4, 2, function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))))
  
  #plot(trace_model)
  #next
  
  #replicate 4
  init.state = c(Be = lab_data_trans$Be[1], Bt = lab_data_trans$Bt[1], Bet = 0,
                 Pl = lab_data_trans$P[1], Pe = 0, Pt = 0)
  
  traj = multi_run2(model, trace_model4, init.state,
                    times = seq(0, 24, 1), nruns = 100)
  
  p4 = ggplot() +
    geom_line(data = traj, aes(time, Be, colour = "Model", linetype = "Bacteria"), size = 0.8) +
    geom_ribbon(data = traj, 
                aes(x = time, ymin = pmax(0, Be-1.96*Be_sd), ymax = Be+1.96*Be_sd,
                    fill = "1"),
                alpha = 0.3) +
    geom_line(data = traj, aes(time, Bt, colour = "Model", linetype = "Bacteria"), size = 0.8) +
    geom_ribbon(data = traj,
                aes(x = time, ymin = pmax(0,Bt-1.96*Bt_sd), ymax = Bt+1.96*Bt_sd,
                    fill = "1"),
                alpha = 0.3) +
    geom_line(data = traj, aes(time, Bet, colour = "Model", linetype = "Bacteria"), size = 0.8) +
    geom_ribbon(data = traj,
                aes(x = time, ymin = pmax(0, Bet-1.96*Bet_sd), ymax = Bet+1.96*Bet_sd,
                    fill = "1"),
                alpha = 0.3) +
    geom_line(data = traj, aes(time, Pl, colour = "Model", linetype = "Phage"), size = 0.8) +
    geom_ribbon(data = traj,
                aes(x = time, ymin = pmax(0, Pl-1.96*Pl_sd), ymax = Pl+1.96*Pl_sd,
                    fill = "1"),
                alpha = 0.3) +
    geom_line(data = lab_data_transM %>% filter(Bacteria != "P"),
              aes(Time, Mean, group = Bacteria, colour = "Data", linetype = "Bacteria"), size = 0.8) +
    geom_errorbar(data = lab_data_transM %>% filter(Bacteria != "P"), 
                  aes(x = Time, ymin = pmax(Mean - se, 0), ymax = Mean + se, group = Bacteria,
                      colour= "Data", linetype = "Bacteria"), size = 0.8) +
    geom_line(data = lab_data_transM %>% filter(Bacteria == "P"), 
              aes(Time, Mean, group = Bacteria, colour = "Data", linetype = "Phage"), size = 0.8) +
    geom_errorbar(data = lab_data_transM %>% filter(Bacteria == "P"), 
                  aes(x = Time, ymin = pmax(Mean - se, 0), ymax = Mean + se, group = Bacteria,
                      colour= "Data", linetype = "Phage"), size = 0.8) +
    scale_y_continuous(trans=log10_trans(),
                       breaks=trans_breaks("log10", function(x) 10^x),
                       labels=trans_format("log10", math_format(10^.x))) +
    coord_cartesian(ylim = c(0.1, 3e10)) +
    scale_x_continuous(breaks = seq(0,24,4)) +
    labs(y = "cfu or pfu per mL", x = "Time (hours)", title = bquote(paste(10^4, " phage (fitting data)")),
         linetype = "Organism:", colour = "Source:") +
    theme_bw() +
    scale_colour_manual(values=c("#685cc4","#6db356")) +
    scale_fill_manual(values = c("#6db356")) +
    guides(fill = F) +
    theme(axis.text.x = element_text(size=12),
          axis.title.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=12),
          legend.text = element_text(size=12),
          legend.title = element_text(size=12),
          strip.text.x = element_text(size=12))
  
  
  
  
  #replicate 5
  init.state = c(Be = lab_data_trans5$Be[1], Bt = lab_data_trans5$Bt[1], Bet = 0,
                 Pl = lab_data_trans5$P[1], Pe = 0, Pt = 0)
  
  traj = multi_run2(model, trace_model4, init.state,
                    times = seq(0, 24, 1), nruns = 100)
  
  p5 = ggplot() +
    geom_line(data = traj, aes(time, Be, colour = "Model", linetype = "Bacteria"), size = 0.8) +
    geom_ribbon(data = traj, 
                aes(x = time, ymin = pmax(0, Be-1.96*Be_sd), ymax = Be+1.96*Be_sd,
                    fill = "1"),
                alpha = 0.3) +
    geom_line(data = traj, aes(time, Bt, colour = "Model", linetype = "Bacteria"), size = 0.8) +
    geom_ribbon(data = traj,
                aes(x = time, ymin = pmax(0,Bt-1.96*Bt_sd), ymax = Bt+1.96*Bt_sd,
                    fill = "1"),
                alpha = 0.3) +
    geom_line(data = traj, aes(time, Bet, colour = "Model", linetype = "Bacteria"), size = 0.8) +
    geom_ribbon(data = traj,
                aes(x = time, ymin =pmax(0, Bet-1.96*Bet_sd), ymax = Bet+1.96*Bet_sd,
                    fill = "1"),
                alpha = 0.3) +
    geom_line(data = traj, aes(time, Pl, colour = "Model", linetype = "Phage"), size = 0.8) +
    geom_ribbon(data = traj,
                aes(x = time, ymin =pmax(0, Pl-1.96*Pl_sd), ymax = Pl+1.96*Pl_sd,
                    fill = "1"),
                alpha = 0.3) +
    geom_line(data = lab_data_trans5M %>% filter(Bacteria != "P"),
              aes(Time, Mean, group = Bacteria, colour = "Data", linetype = "Bacteria"), size = 0.8) +
    geom_errorbar(data = lab_data_trans5M %>% filter(Bacteria != "P"), 
                  aes(x = Time, ymin = pmax(Mean - se, 0), ymax = Mean + se, group = Bacteria,
                      colour= "Data", linetype = "Bacteria"), size = 0.8) +
    geom_line(data = lab_data_trans5M %>% filter(Bacteria == "P"), 
              aes(Time, Mean, group = Bacteria, colour = "Data", linetype = "Phage"), size = 0.8) +
    geom_errorbar(data = lab_data_trans5M %>% filter(Bacteria == "P"), 
                  aes(x = Time, ymin = pmax(Mean - se, 0), ymax = Mean + se, group = Bacteria,
                      colour= "Data", linetype = "Phage"), size = 0.8) +
    scale_y_continuous(trans=log10_trans(),
                       breaks=trans_breaks("log10", function(x) 10^x),
                       labels=trans_format("log10", math_format(10^.x))) +
    coord_cartesian(ylim = c(0.1, 3e10)) +
    scale_x_continuous(breaks = seq(0,24,4)) +
    labs(y = "cfu or pfu per mL", x = "Time (hours)", title = bquote(paste(10^5, " phage")),
         linetype = "Organism:", colour = "Source:") +
    theme_bw() +
    scale_colour_manual(values=c("#685cc4","#6db356")) +
    scale_fill_manual(values = c("#6db356")) +
    guides(fill = F) +
    theme(axis.text.x = element_text(size=12),
          axis.title.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=12),
          legend.text = element_text(size=12),
          legend.title = element_text(size=12),
          strip.text.x = element_text(size=12))
  
  
  #replicate 3
  init.state = c(Be = lab_data_trans3$Be[1], Bt = lab_data_trans3$Bt[1], Bet = 0,
                 Pl = lab_data_trans3$P[1], Pe = 0, Pt = 0)
  
  traj = multi_run2(model, trace_model4, init.state,
                    times = seq(0, 24, 1), nruns = 100)
  
  p3 = ggplot() +
    geom_line(data = traj, aes(time, Be, colour = "Model", linetype = "Bacteria"), size = 0.8) +
    geom_ribbon(data = traj, 
                aes(x = time, ymin = pmax(0, Be-1.96*Be_sd), ymax = Be+1.96*Be_sd,
                    fill = "1"),
                alpha = 0.3) +
    geom_line(data = traj, aes(time, Bt, colour = "Model", linetype = "Bacteria"), size = 0.8) +
    geom_ribbon(data = traj,
                aes(x = time, ymin = pmax(0,Bt-1.96*Bt_sd), ymax = Bt+1.96*Bt_sd,
                    fill = "1"),
                alpha = 0.3) +
    geom_line(data = traj, aes(time, Bet, colour = "Model", linetype = "Bacteria"), size = 0.8) +
    geom_ribbon(data = traj,
                aes(x = time, ymin =pmax(0, Bet-1.96*Bet_sd), ymax = Bet+1.96*Bet_sd,
                    fill = "1"),
                alpha = 0.3) +
    geom_line(data = traj, aes(time, Pl, colour = "Model", linetype = "Phage"), size = 0.8) +
    geom_ribbon(data = traj,
                aes(x = time, ymin =pmax(0, Pl-1.96*Pl_sd), ymax = Pl+1.96*Pl_sd,
                    fill = "1"),
                alpha = 0.3) +
    geom_line(data = lab_data_trans3M %>% filter(Bacteria != "P"),
              aes(Time, Mean, group = Bacteria, colour = "Data", linetype = "Bacteria"), size = 0.8) +
    geom_errorbar(data = lab_data_trans3M %>% filter(Bacteria != "P"), 
                  aes(x = Time, ymin = pmax(Mean - se, 0), ymax = Mean + se, group = Bacteria,
                      colour= "Data", linetype = "Bacteria"), size = 0.8) +
    geom_line(data = lab_data_trans3M %>% filter(Bacteria == "P"), 
              aes(Time, Mean, group = Bacteria, colour = "Data", linetype = "Phage"), size = 0.8) +
    geom_errorbar(data = lab_data_trans3M %>% filter(Bacteria == "P"), 
                  aes(x = Time, ymin = pmax(Mean - se, 0), ymax = Mean + se, group = Bacteria,
                      colour= "Data", linetype = "Phage"), size = 0.8) +
    scale_y_continuous(trans=log10_trans(),
                       breaks=trans_breaks("log10", function(x) 10^x),
                       labels=trans_format("log10", math_format(10^.x))) +
    coord_cartesian(ylim = c(0.1, 3e10)) +
    scale_x_continuous(breaks = seq(0,24,4)) +
    labs(y = "cfu or pfu per mL", x = "Time (hours)", title = bquote(paste(10^3, " phage")), 
         linetype = "Organism:", colour = "Source:") +
    theme_bw() +
    scale_colour_manual(values=c("#685cc4","#6db356")) +
    scale_fill_manual(values = c("#6db356")) +
    guides(fill = F) +
    theme(axis.text.x = element_text(size=12),
          axis.title.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=12),
          legend.text = element_text(size=12),
          legend.title = element_text(size=12),
          strip.text.x = element_text(size=12))
  
  
  #final plot
  legend = get_legend(p4 + theme(legend.position = "bottom", legend.box = "vertical"))
  
  plot_grid(p3 + theme(legend.position = "none"), 
            p4 + theme(legend.position = "none"), 
            p5 + theme(legend.position = "none"),
            legend)
  
  filename = paste0("multi_", models_to_try$model_name[i], ".png")
  ggsave(here::here("Fitting", "Multi_runs", filename))
  
}

colnames(best_params) = c("model_name", "value", colnames(trace_model4))
best_params[,-1] = apply(best_params[,-1], c(1,2), as.numeric)
best_params$beta = 1/best_params$beta
best_params$gamma = 1/best_params$gamma
best_params$alpha = 1/best_params$alpha

write.csv(best_params, here::here("Fitting", "best_params.csv"), row.names = F)

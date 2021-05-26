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

files = list.files(here::here("Fitting", "Fitted_params"))

params = vector("list", length(files))
i=1

for(f in files){
  params[[i]] = read.csv(here::here("Fitting", "Fitted_params", f))
  i = i+1
}

names(params) = gsub(".csv", "", gsub("params_", "", files))


lab_data_transM = read.csv(here::here("Lab", "Transduction", "summary_10_4.csv")) %>%
  select(Time, Bacteria, Mean, se) %>%
  filter(Bacteria != "Total")
lab_data_trans5M = read.csv(here::here("Lab", "Transduction", "summary_10_5.csv")) %>%
  select(Time, Bacteria, Mean, se)%>%
  filter(Bacteria != "Total")
lab_data_trans3M = read.csv(here::here("Lab", "Transduction", "summary_10_3.csv")) %>%
  select(Time, Bacteria, Mean, se) %>%
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


## GO #####


models_to_try = data.frame(model_name="dens_beta", frequentist=FALSE,
                           fixed_delay=NA, decay=TRUE,
                           link_beta=TRUE, link_L=FALSE, link_delay=FALSE, transduction=TRUE)
models_to_try = rbind(models_to_try,
                      data.frame(model_name="dens_burst", frequentist=FALSE,
                                 fixed_delay=NA, decay=TRUE,
                                 link_beta=FALSE, link_L=TRUE, link_delay=FALSE, transduction=TRUE))
models_to_try = rbind(models_to_try,
                      data.frame(model_name="dens_both", frequentist=FALSE,
                                 fixed_delay=NA, decay=TRUE,
                                 link_beta=TRUE, link_L=TRUE, link_delay=FALSE, transduction=TRUE))
models_to_try = rbind(models_to_try,
                      data.frame(model_name="freq_beta", frequentist=TRUE,
                                 fixed_delay=NA, decay=TRUE,
                                 link_beta=TRUE, link_L=FALSE, link_delay=FALSE, transduction=TRUE))
models_to_try = rbind(models_to_try,
                      data.frame(model_name="freq_burst", frequentist=TRUE,
                                 fixed_delay=NA, decay=TRUE,
                                 link_beta=FALSE, link_L=TRUE, link_delay=FALSE, transduction=TRUE))
models_to_try = rbind(models_to_try,
                      data.frame(model_name="freq_both", frequentist=TRUE,
                                 fixed_delay=NA, decay=TRUE,
                                 link_beta=TRUE, link_L=TRUE, link_delay=FALSE, transduction=TRUE))


best_params = data.frame()

all_traj_4 = data.frame()
all_traj_3 = data.frame()
all_traj_5 = data.frame()

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
  
  trace_model = params[[models_to_try$model_name[i]]]
  
  quants = apply(trace_model, 2, function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
  quants_s = c(models_to_try$model_name[i])
  for(j in 1:(ncol(quants)-1)){
    quants_s = c(quants_s, quants[2,j], quants[1,j], quants[3,j])
  }
  
  best_params = rbind(best_params,
                      quants_s)
  
  
  #replicate 4
  init.state = c(Be = lab_data_trans$Be[1], Bt = lab_data_trans$Bt[1], Bet = 0,
                 Pl = lab_data_trans$P[1], Pe = 0, Pt = 0)
  
  traj = multi_run2(model, trace_model, init.state,
                    times = seq(0, 24, 1), nruns = 500)
  
  traj$model = models_to_try$model_name[i]
  all_traj_4 = rbind(all_traj_4, traj)
  
  
  #replicate 5
  init.state = c(Be = lab_data_trans5$Be[1], Bt = lab_data_trans5$Bt[1], Bet = 0,
                 Pl = lab_data_trans5$P[1], Pe = 0, Pt = 0)
  
  traj = multi_run2(model, trace_model, init.state,
                    times = seq(0, 24, 1), nruns = 500)
  
  traj$model = models_to_try$model_name[i]
  all_traj_5 = rbind(all_traj_5, traj)
  
  #replicate 3
  init.state = c(Be = lab_data_trans3$Be[1], Bt = lab_data_trans3$Bt[1], Bet = 0,
                 Pl = lab_data_trans3$P[1], Pe = 0, Pt = 0)
  
  traj = multi_run2(model, trace_model, init.state,
                    times = seq(0, 24, 1), nruns = 500)
  
  traj$model = models_to_try$model_name[i]
  all_traj_3 = rbind(all_traj_3, traj)
  
  
}

quants_names = c("model_name")
for(j in 1:(ncol(quants)-1)){
  quants_names = c(quants_names, colnames(quants)[j],
                   paste0(colnames(quants)[j], "_0.025"),
                   paste0(colnames(quants)[j], "_0.975"))
}

colnames(best_params) = quants_names

best_params[,-1] = apply(best_params[,-1], c(1,2), as.numeric)
best_params[,c(2:4,8:10)] = apply(best_params[,c(2:4,8:10)], c(1,2), function(x) 1/x)

write.csv(best_params, here::here("Fitting", "best_params.csv"), row.names = F)




all_traj_4[all_traj_4<0.01] = 0.01
all_traj_5[all_traj_5<0.01] = 0.01
all_traj_3[all_traj_3<0.01] = 0.01

all_traj_4$group = all_traj_4$model
all_traj_4 = all_traj_4 %>% 
  mutate(model = replace(model, group %in% c("dens_beta",
                                             "dens_both",
                                             "freq_beta",
                                             "freq_both"), "other"))

all_traj_5$group = all_traj_5$model
all_traj_5 = all_traj_5 %>% 
  mutate(model = replace(model, group %in% c("dens_beta",
                                             "dens_both",
                                             "freq_beta",
                                             "freq_both"), "other"))
all_traj_3$group = all_traj_3$model
all_traj_3 = all_traj_3 %>% 
  mutate(model = replace(model, group %in% c("dens_beta",
                                             "dens_both",
                                             "freq_beta",
                                             "freq_both"), "other"))

p4 = ggplot() +
  geom_line(data = all_traj_4, aes(time, Pl, group = group, colour = model,
                                   linetype = "Model", size = model)) +
  geom_ribbon(data = all_traj_4, aes(x = time, ymin = pmax(Pl - 1.96*Pl_sd, 0), ymax = Pl + 1.96*Pl_sd,
                                     group= group, fill=model), alpha = 0.1) +
  geom_line(data = lab_data_transM %>% filter(Bacteria == "P"), 
            aes(Time, Mean, linetype = "Data", colour = "PL"), size = 0.8) +
  geom_errorbar(data = lab_data_transM %>% filter(Bacteria == "P"), 
                aes(x = Time, ymin = pmax(Mean - se, 0), ymax = Mean + se,
                    colour= "PL"), size = 0.8) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(breaks = seq(0,24,4)) +
  coord_cartesian(ylim = c(1e3, 1e11)) +
  scale_linetype_manual(breaks = c("Data", "Model"),
                        labels = c("Data", "Model"),
                        values = c(1, 2)) +
  scale_colour_manual(breaks = c("PL"),
                      labels = c("Lytic phage"),
                      values = c("darkgreen", "#5ad3ec", "grey", "#c88a33")) +
  scale_fill_manual(breaks = c("dens_burst", "freq_burst", "other"),
                    labels = c("Density, link burst size", "Frequency, link burst size", "Other"),
                    values = c("darkgreen", "#5ad3ec", "grey"),
                    guide = guide_legend(override.aes = list(size = c(0.8,0.8,0.5),
                                                             linetype = c(2,2,2)))) +
  scale_size_manual(breaks = c("dens_burst", "freq_burst", "other"),
                    labels = c("Density, link burst size", "Frequency, link burst size", "Other"),
                    values = c(0.8,0.8,0.5)) +
  labs(x = "Time (hours)", y = "cfu or pfu per mL", linetype = "Source:",
       colour = "Organism:", fill = "Model type:", size = "Model type:",
       title = bquote(paste("B) ", 10^4, " dataset, phage"))) +
  guides(linetype = guide_legend(order = 1)) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.x = element_text(size=12))



p5 = ggplot() +
  geom_line(data = all_traj_5, aes(time, Pl, group = group, colour = model,
                                   linetype = "Model", size = model)) +
  geom_ribbon(data = all_traj_5, aes(x = time, ymin = pmax(Pl - 1.96*Pl_sd, 0), ymax = Pl + 1.96*Pl_sd,
                                     group= group, fill=model), alpha = 0.1) +
  geom_line(data = lab_data_trans5M %>% filter(Bacteria == "P"), 
            aes(Time, Mean, linetype = "Data", colour = "PL"), size = 0.8) +
  geom_errorbar(data = lab_data_trans5M %>% filter(Bacteria == "P"), 
                aes(x = Time, ymin = pmax(Mean - se, 0), ymax = Mean + se,
                    colour= "PL"), size = 0.8) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(breaks = seq(0,24,4)) +
  coord_cartesian(ylim = c(1e3, 1e11)) +
  scale_linetype_manual(breaks = c("Data", "Model"),
                        labels = c("Data", "Model"),
                        values = c(1, 2)) +
  scale_colour_manual(breaks = c("PL"),
                      labels = c("Lytic phage"),
                      values = c("darkgreen", "#5ad3ec", "grey", "#c88a33")) +
  scale_fill_manual(breaks = c("dens_burst", "freq_burst", "other"),
                    labels = c("Density, link burst size", "Frequency, link burst size", "Other"),
                    values = c("darkgreen", "#5ad3ec", "grey"),
                    guide = guide_legend(override.aes = list(size = c(0.8,0.8,0.5),
                                                             linetype = c(2,2,2)))) +
  scale_size_manual(breaks = c("dens_burst", "freq_burst", "other"),
                    labels = c("Density, link burst size", "Frequency, link burst size", "Other"),
                    values = c(0.8,0.8,0.5)) +
  labs(x = "Time (hours)", y = "cfu or pfu per mL", linetype = "Source:",
       colour = "Organism:", fill = "Model type:", size = "Model type:",
       title = bquote(paste("C) ", 10^5, " dataset, phage"))) +
  guides(linetype = guide_legend(order = 1)) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.x = element_text(size=12))


p3 = ggplot() +
  geom_line(data = all_traj_3, aes(time, Pl, group = group, colour = model,
                                   linetype = "Model", size = model)) +
  geom_ribbon(data = all_traj_3, aes(x = time, ymin = pmax(Pl - 1.96*Pl_sd, 0), ymax = Pl + 1.96*Pl_sd,
                                     group= group, fill=model), alpha = 0.1) +
  geom_line(data = lab_data_trans3M %>% filter(Bacteria == "P"), 
            aes(Time, Mean, linetype = "Data", colour = "PL"), size = 0.8) +
  geom_errorbar(data = lab_data_trans3M %>% filter(Bacteria == "P"), 
                aes(x = Time, ymin = pmax(Mean - se, 0), ymax = Mean + se,
                    colour= "PL"), size = 0.8) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(breaks = seq(0,24,4)) +
  coord_cartesian(ylim = c(1e3, 1e11)) +
  scale_linetype_manual(breaks = c("Data", "Model"),
                        labels = c("Data", "Model"),
                        values = c(1, 2)) +
  scale_colour_manual(breaks = c("PL"),
                      labels = c("Lytic phage"),
                      values = c("darkgreen", "#5ad3ec", "grey", "#c88a33")) +
  scale_fill_manual(breaks = c("dens_burst", "freq_burst", "other"),
                    labels = c("Density, link burst size", "Frequency, link burst size", "Other"),
                    values = c("darkgreen", "#5ad3ec", "grey"),
                    guide = guide_legend(override.aes = list(size = c(0.8,0.8,0.5),
                                                             linetype = c(2,2,2)))) +
  scale_size_manual(breaks = c("dens_burst", "freq_burst", "other"),
                    labels = c("Density, link burst size", "Frequency, link burst size", "Other"),
                    values = c(0.8,0.8,0.5)) +
  labs(x = "Time (hours)", y = "cfu or pfu per mL", linetype = "Source:",
       colour = "Organism:", fill = "Model type:", size = "Model type:",
       title = bquote(paste("A) ", 10^3, " dataset, phage"))) +
  guides(linetype = guide_legend(order = 1)) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.x = element_text(size=12))


p5_extra = ggplot() +
  geom_line(data = all_traj_5, aes(time, Pl, group = group, colour = model,
                                   linetype = "Model", size = model)) +
  geom_ribbon(data = all_traj_5, aes(x = time, ymin = pmax(Pl - 1.96*Pl_sd, 0), ymax = Pl + 1.96*Pl_sd,
                                     group= group, fill=model), alpha = 0.1) +
  geom_line(data = all_traj_5, aes(time, Bet, group = group, colour = model,
                                   linetype = "Model", size = model)) +
  geom_ribbon(data = all_traj_5, aes(x = time, ymin = pmax(Bet - 1.96*Bet_sd, 0), ymax = Bet + 1.96*Bet_sd,
                                     group = group, fill = model), alpha = 0.1) +
  geom_line(data = lab_data_trans5M %>% filter(Bacteria == "P"), 
            aes(Time, Mean, linetype = "Data", colour = "PL"), size = 0.8) +
  geom_errorbar(data = lab_data_trans5M %>% filter(Bacteria == "P"), 
                aes(x = Time, ymin = pmax(Mean - se, 0), ymax = Mean + se,
                    colour= "PL"), size = 0.8) +
  geom_line(data = lab_data_trans5M %>% filter(Bacteria == "DRP"), 
            aes(Time, Mean, linetype = "Data", colour = "BET"), size = 0.8) +
  geom_errorbar(data = lab_data_trans5M %>% filter(Bacteria == "DRP"), 
                aes(x = Time, ymin = pmax(Mean - se, 0), ymax = Mean + se,
                    colour= "BET"), size = 0.8) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(breaks = seq(0,24,4)) +
  coord_cartesian(ylim = c(1e5, 1e11), xlim = c(16, 24)) +
  scale_linetype_manual(breaks = c("Data", "Model"),
                        labels = c("Data", "Model"),
                        values = c(1, 2)) +
  scale_colour_manual(breaks = c("PL", "BET"),
                      labels = c("Lytic phage", "Double resistant bacteria (DRP)"),
                      values = c("#c2484d", "darkgreen", "#5ad3ec", "grey", "#c88a33")) +
  scale_fill_manual(breaks = c("dens_burst", "freq_burst", "other"),
                    labels = c("Density, link burst size", "Frequency, link burst size", "Other"),
                    values = c("darkgreen", "#5ad3ec", "grey"),
                    guide = guide_legend(override.aes = list(size = c(0.8,0.8,0.5),
                                                             linetype = c(2,2,2)))) +
  scale_size_manual(breaks = c("dens_burst", "freq_burst", "other"),
                    labels = c("Density, link burst size", "Frequency, link burst size", "Other"),
                    values = c(0.8,0.8,0.5)) +
  labs(x = "Time (hours)", y = "cfu or pfu per mL", linetype = "Source:",
       colour = "Organism:", fill = "Model type:", size = "Model type:",
       title = bquote(paste("D) ", 10^5, " phage, zoomed in"))) +
  guides(linetype = guide_legend(order = 1)) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.x = element_text(size=12))


p4_bet = ggplot() +
  geom_line(data = all_traj_4, aes(time, Bet, group = group, colour = model,
                                   linetype = "Model", size = model)) +
  geom_ribbon(data = all_traj_4, aes(x = time, ymin = pmax(Bet - 1.96*Bet_sd, 0), ymax = Bet + 1.96*Bet_sd,
                                     group = group, fill = model), alpha = 0.1) +
  geom_line(data = lab_data_transM %>% filter(Bacteria == "DRP"), 
            aes(Time, Mean, linetype = "Data", colour = "BET"), size = 0.8) +
  geom_errorbar(data = lab_data_transM %>% filter(Bacteria == "DRP"), 
                aes(x = Time, ymin = pmax(Mean - se, 0), ymax = Mean + se,
                    colour= "BET"), size = 0.8) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x, n = 4),
                     labels=trans_format("log10", math_format(10^.x)),) +
  scale_x_continuous(breaks = seq(0,24,4)) +
  coord_cartesian(ylim = c(0.1, 1e2)) +
  scale_linetype_manual(breaks = c("Data", "Model"),
                        labels = c("Data", "Model"),
                        values = c(1, 2)) +
  scale_colour_manual(breaks = c("BET"),
                      labels = c("Double resistant bacteria"),
                      values = c("#c2484d", "darkgreen", "#5ad3ec", "grey")) +
  scale_fill_manual(breaks = c("dens_burst", "freq_burst", "other"),
                    labels = c("Density, link burst size", "Frequency, link burst size", "Other"),
                    values = c("darkgreen", "#5ad3ec", "grey"),
                    guide = guide_legend(override.aes = list(size = c(0.8,0.8,0.5),
                                                             linetype = c(2,2,2)))) +
  scale_size_manual(breaks = c("dens_burst", "freq_burst", "other"),
                    labels = c("Density, link burst size", "Frequency, link burst size", "Other"),
                    values = c(0.8,0.8,0.5)) +
  labs(x = "Time (hours)", y = "cfu or pfu per mL", linetype = "Source:",
       colour = "Organism:", fill = "Model type:", size = "Model type:",
       title = bquote(paste("E) ", 10^4, " dataset, DRP"))) +
  guides(linetype = guide_legend(order = 1)) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.x = element_text(size=12))



p5_bet = ggplot() +
  geom_line(data = all_traj_5, aes(time, Bet, group = group, colour = model,
                                   linetype = "Model", size = model)) +
  geom_ribbon(data = all_traj_5, aes(x = time, ymin = pmax(Bet - 1.96*Bet_sd, 0), ymax = Bet + 1.96*Bet_sd,
                                     group = group, fill = model), alpha = 0.1) +
  geom_line(data = lab_data_trans5M %>% filter(Bacteria == "DRP"), 
            aes(Time, Mean, linetype = "Data", colour = "BET"), size = 0.8) +
  geom_errorbar(data = lab_data_trans5M %>% filter(Bacteria == "DRP"), 
                aes(x = Time, ymin = pmax(Mean - se, 0), ymax = Mean + se,
                    colour= "BET"), size = 0.8) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x, n = 4),
                     labels=trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(breaks = seq(0,24,4)) +
  coord_cartesian(ylim = c(0.1, 1e2)) +
  scale_linetype_manual(breaks = c("Data", "Model"),
                        labels = c("Data", "Model"),
                        values = c(1, 2)) +
  scale_colour_manual(breaks = c("BET"),
                      labels = c("Double resistant bacteria"),
                      values = c("#c2484d", "darkgreen", "#5ad3ec", "grey")) +
  scale_fill_manual(breaks = c("dens_burst", "freq_burst", "other"),
                    labels = c("Density, link burst size", "Frequency, link burst size", "Other"),
                    values = c("darkgreen", "#5ad3ec", "grey"),
                    guide = guide_legend(override.aes = list(size = c(0.8,0.8,0.5),
                                                             linetype = c(2,2,2)))) +
  scale_size_manual(breaks = c("dens_burst", "freq_burst", "other"),
                    labels = c("Density, link burst size", "Frequency, link burst size", "Other"),
                    values = c(0.8,0.8,0.5)) +
  labs(x = "Time (hours)", y = "cfu or pfu per mL", linetype = "Source:",
       colour = "Organism:", fill = "Model type:", size = "Model type:",
       title = bquote(paste("F) ", 10^5, " dataset, DRP"))) +
  guides(linetype = guide_legend(order = 1)) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.x = element_text(size=12))


p3_bet = ggplot() +
  geom_line(data = all_traj_3, aes(time, Bet, group = group, colour = model,
                                   linetype = "Model", size = model)) +
  geom_ribbon(data = all_traj_3, aes(x = time, ymin = pmax(Bet - 1.96*Bet_sd, 0), ymax = Bet + 1.96*Bet_sd,
                                     group = group, fill = model), alpha = 0.1) +
  geom_line(data = lab_data_trans3M %>% filter(Bacteria == "DRP"), 
            aes(Time, Mean, linetype = "Data", colour = "BET"), size = 0.8) +
  geom_errorbar(data = lab_data_trans3M %>% filter(Bacteria == "DRP"), 
                aes(x = Time, ymin = pmax(Mean - se, 0), ymax = Mean + se,
                    colour= "BET"), size = 0.8) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x, n = 4),
                     labels=trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(breaks = seq(0,24,4)) +
  coord_cartesian(ylim = c(0.1, 1e2)) +
  scale_linetype_manual(breaks = c("Data", "Model"),
                        labels = c("Data", "Model"),
                        values = c(1, 2)) +
  scale_colour_manual(breaks = c("BET"),
                      labels = c("Double resistant bacteria"),
                      values = c("#c2484d", "darkgreen", "#5ad3ec", "grey")) +
  scale_fill_manual(breaks = c("dens_burst", "freq_burst", "other"),
                    labels = c("Density, link burst size", "Frequency, link burst size", "Other"),
                    values = c("darkgreen", "#5ad3ec", "grey"),
                    guide = guide_legend(override.aes = list(size = c(0.8,0.8,0.5),
                                                             linetype = c(2,2,2)))) +
  scale_size_manual(breaks = c("dens_burst", "freq_burst", "other"),
                    labels = c("Density, link burst size", "Frequency, link burst size", "Other"),
                    values = c(0.8,0.8,0.5)) +
  labs(x = "Time (hours)", y = "cfu or pfu per mL", linetype = "Source:",
       colour = "Organism:", fill = "Model type:", size = "Model type:",
       title = bquote(paste("D) ", 10^3, " dataset, DRP"))) +
  guides(linetype = guide_legend(order = 1)) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.x = element_text(size=12))

#final plot
legend = get_legend(p5_extra + theme(legend.position = "right", legend.box = "horizontal"))
# plot_grid(plot_grid(p3 + theme(legend.position = "none"), 
#                     p4 + theme(legend.position = "none"), 
#                     p5 + theme(legend.position = "none"),
#                     p5_extra + theme(legend.position = "none"),
#                     p3_bet + theme(legend.position = "none"),
#                     p4_bet + theme(legend.position = "none"),
#                     p5_bet + theme(legend.position = "none"), ncol = 2),
#           NULL,
#           legend,
#           nrow = 3, rel_heights = c(1,0.05,0.2))

plot_grid(plot_grid(p3 + theme(legend.position = "none"), 
                    p4 + theme(legend.position = "none"), 
                    p5 + theme(legend.position = "none"),
                    p3_bet + theme(legend.position = "none"),
                    p4_bet + theme(legend.position = "none"),
                    p5_bet + theme(legend.position = "none"), nrow = 2),
          NULL,
          legend,
          nrow = 3, rel_heights = c(1,0.05,0.2))

ggsave("plot_compare_all.png", height = 8, width = 18)


all_traj_4 %>%
  filter(group == "freq_burst") %>%
  filter(time %in% lab_data_trans$time) %>%
  mutate(data_Pl = lab_data_trans$P[-1]) %>%
  mutate(diff = data_Pl/Pl)


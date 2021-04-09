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

fitted_params4 = c(readRDS(here::here("Fitting", "10_4", "best_params_transduction.rds")),
                   readRDS(here::here("Fitting", "10_4", "best_params_transduction2.rds")),
                   readRDS(here::here("Fitting", "10_4", "best_params_transduction3.rds")))

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

all_traj_4 = data.frame()
all_traj_3 = data.frame()
all_traj_5 = data.frame()

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
  
  trace_model4 = fitted_params4[[models_to_try$model_name[i]]]
  #trace_model4 = trace_model4[which.max(trace_model4[,"log.density"]),]
  trace_model4 = coda::mcmc(trace_model4)
  trace_model4 = burnAndThin(trace_model4, burn = 20000, thin = 10)
  
  #plot(trace_model)
  #next
  
  #replicate 4
  init.state = c(Be = lab_data_trans$Be[1], Bt = lab_data_trans$Bt[1], Bet = 0,
                 Pl = lab_data_trans$P[1], Pe = 0, Pt = 0)
  
  traj = multi_run2(model, trace_model4, init.state,
                    times = seq(0, 24, 1), nruns = 500)
  
  traj$model = models_to_try$model_name[i]
  all_traj_4 = rbind(all_traj_4, traj)
  
  
  #replicate 5
  init.state = c(Be = lab_data_trans5$Be[1], Bt = lab_data_trans5$Bt[1], Bet = 0,
                 Pl = lab_data_trans5$P[1], Pe = 0, Pt = 0)
  
  traj = multi_run2(model, trace_model4, init.state,
                    times = seq(0, 24, 1), nruns = 500)
  
  traj$model = models_to_try$model_name[i]
  all_traj_5 = rbind(all_traj_5, traj)
  
  #replicate 3
  init.state = c(Be = lab_data_trans3$Be[1], Bt = lab_data_trans3$Bt[1], Bet = 0,
                 Pl = lab_data_trans3$P[1], Pe = 0, Pt = 0)
  
  traj = multi_run2(model, trace_model4, init.state,
                    times = seq(0, 24, 1), nruns = 500)
  
  traj$model = models_to_try$model_name[i]
  all_traj_3 = rbind(all_traj_3, traj)
  

}

all_traj_4$group = all_traj_4$model
all_traj_4 = all_traj_4 %>% 
  mutate(model = replace(model, group %in% c("mass_decay_link_beta",
                                             "mass_decay_link_both",
                                             "frequentist_decay_link_beta",
                                             "frequentist_decay_link_both"), "other"))
  
all_traj_5$group = all_traj_5$model
all_traj_5 = all_traj_5 %>% 
  mutate(model = replace(model, group %in% c("mass_decay_link_beta",
                                             "mass_decay_link_both",
                                             "frequentist_decay_link_beta",
                                             "frequentist_decay_link_both"), "other"))
all_traj_3$group = all_traj_3$model
all_traj_3 = all_traj_3 %>% 
  mutate(model = replace(model, group %in% c("mass_decay_link_beta",
                                             "mass_decay_link_both",
                                             "frequentist_decay_link_beta",
                                             "frequentist_decay_link_both"), "other"))

p4 = ggplot() +
  geom_line(data = all_traj_4, aes(time, Pl, group = group, colour = model, linetype = "Model")) +
  geom_ribbon(data = all_traj_4, aes(x = time, ymin = pmax(Pl - 1.96*Pl_sd, 0), ymax = Pl + 1.96*Pl_sd,
                                     group= group, fill=model), alpha = 0.1) +
  geom_line(data = all_traj_4, aes(time, Bet, group = group, colour = model,linetype = "Model")) +
  geom_ribbon(data = all_traj_4, aes(x = time, ymin = pmax(Bet - 1.96*Bet_sd, 0), ymax = Bet + 1.96*Bet_sd,
                                     group = group, fill = model), alpha = 0.1) +
  geom_line(data = lab_data_transM %>% filter(Bacteria == "P"), 
            aes(Time, Mean, linetype = "Data", colour = "PL"), size = 0.8) +
  geom_errorbar(data = lab_data_transM %>% filter(Bacteria == "P"), 
                aes(x = Time, ymin = pmax(Mean - se, 0), ymax = Mean + se,
                    colour= "PL"), size = 0.8) +
  geom_line(data = lab_data_transM %>% filter(Bacteria == "DRP"), 
            aes(Time, Mean, linetype = "Data", colour = "BET"), size = 0.8) +
  geom_errorbar(data = lab_data_transM %>% filter(Bacteria == "DRP"), 
                aes(x = Time, ymin = pmax(Mean - se, 0), ymax = Mean + se,
                    colour= "BET"), size = 0.8) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(breaks = seq(0,24,4)) +
  coord_cartesian(ylim = c(0.1, 1e11), xlim = c(-0.1, 24.1)) +
  scale_linetype_manual(breaks = c("Model", "Data"),
                        labels = c("Model", "Data"),
                        values = c("dashed", "solid")) +
  scale_colour_manual(breaks = c("PL", "BET"),
                      labels = c(bquote(P[L]), bquote(B[ET])),
                      values = c("#c2484d", "grey", "grey", "grey", "#c88a33")) +
  scale_fill_manual(breaks = c("mass_decay_link_L", "frequentist_decay_link_L", "other"),
                    labels = c("Density, link burst size", "Frequency, link burst size", "Other"),
                    values = c("red", "blue", "grey"),
                    guide = guide_legend(override.aes = list(linetype = c(1,1,1)))) +
  labs(x = "Time (hours)", y = "cfu or pfu per mL", linetype = "Source:",
       colour = "Organism:", fill = "Model type:",
       title = bquote(paste("B) ", 10^4, " phage"))) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.x = element_text(size=12))



p5 = ggplot() +
  geom_line(data = all_traj_5, aes(time, Pl, group = group, colour = model, linetype = "Model")) +
  geom_ribbon(data = all_traj_5, aes(x = time, ymin = pmax(Pl - 1.96*Pl_sd, 0), ymax = Pl + 1.96*Pl_sd,
                                     group= group, fill=model), alpha = 0.1) +
  geom_line(data = all_traj_5, aes(time, Bet, group = group, colour = model,linetype = "Model")) +
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
  coord_cartesian(ylim = c(0.1, 1e11), xlim = c(-0.1, 24.1)) +
  scale_linetype_manual(breaks = c("Model", "Data"),
                        labels = c("Model", "Data"),
                        values = c("dashed", "solid")) +
  scale_colour_manual(breaks = c("PL", "BET"),
                      labels = c(bquote(P[L]), bquote(B[ET])),
                      values = c("#c2484d", "grey", "grey", "grey", "#c88a33")) +
  scale_fill_manual(breaks = c("mass_decay_link_L", "frequentist_decay_link_L", "other"),
                    labels = c("Density, link burst size", "Frequency, link burst size", "Other"),
                    values = c("red", "blue", "grey"),
                    guide = guide_legend(override.aes = list(linetype = c(1,1,1)))) +
  labs(x = "Time (hours)", y = "cfu or pfu per mL", linetype = "Source:",
       colour = "Organism:", fill = "Model type:",
       title = bquote(paste("C) ", 10^5, " phage"))) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.x = element_text(size=12))


p3 = ggplot() +
  geom_line(data = all_traj_3, aes(time, Pl, group = group, colour = model, linetype = "Model")) +
  geom_ribbon(data = all_traj_3, aes(x = time, ymin = pmax(Pl - 1.96*Pl_sd, 0), ymax = Pl + 1.96*Pl_sd,
                                     group= group, fill=model), alpha = 0.1) +
  geom_line(data = all_traj_3, aes(time, Bet, group = group, colour = model,linetype = "Model")) +
  geom_ribbon(data = all_traj_3, aes(x = time, ymin = pmax(Bet - 1.96*Bet_sd, 0), ymax = Bet + 1.96*Bet_sd,
                                     group = group, fill = model), alpha = 0.1) +
  geom_line(data = lab_data_trans3M %>% filter(Bacteria == "P"), 
            aes(Time, Mean, linetype = "Data", colour = "PL"), size = 0.8) +
  geom_errorbar(data = lab_data_trans3M %>% filter(Bacteria == "P"), 
                aes(x = Time, ymin = pmax(Mean - se, 0), ymax = Mean + se,
                    colour= "PL"), size = 0.8) +
  geom_line(data = lab_data_trans3M %>% filter(Bacteria == "DRP"), 
            aes(Time, Mean, linetype = "Data", colour = "BET"), size = 0.8) +
  geom_errorbar(data = lab_data_trans3M %>% filter(Bacteria == "DRP"), 
                aes(x = Time, ymin = pmax(Mean - se, 0), ymax = Mean + se,
                    colour= "BET"), size = 0.8) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(breaks = seq(0,24,4)) +
  coord_cartesian(ylim = c(0.1, 1e11), xlim = c(-0.1, 24.1)) +
  scale_linetype_manual(breaks = c("Model", "Data"),
                        labels = c("Model", "Data"),
                        values = c("dashed", "solid")) +
  scale_colour_manual(breaks = c("PL", "BET"),
                      labels = c(bquote(P[L]), bquote(B[ET])),
                      values = c("#c2484d", "grey", "grey", "grey", "#c88a33")) +
  scale_fill_manual(breaks = c("mass_decay_link_L", "frequentist_decay_link_L", "other"),
                    labels = c("Density, link burst size", "Frequency, link burst size", "Other"),
                    values = c("red", "blue", "grey"),
                    guide = guide_legend(override.aes = list(linetype = c(1,1,1)))) +
  labs(x = "Time (hours)", y = "cfu or pfu per mL", linetype = "Source:",
       colour = "Organism:", fill = "Model type:",
       title = bquote(paste("A) ", 10^3, " phage"))) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.x = element_text(size=12))


#final plot
legend = get_legend(p4 + theme(legend.position = "right", legend.box = "horizontal"))

plot_grid(p3 + theme(legend.position = "none"), 
          p4 + theme(legend.position = "none"), 
          p5 + theme(legend.position = "none"),
          legend)

ggsave("plot_compare_all.png")

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

DIC = function(trace, model, lab_data, lab_data2){
  
  theta.bar <- colMeans(trace[,model$theta.names])
  init.state <- c(Be = lab_data$Be[1], Bt = lab_data$Bt[1], Bet = 0,
                  Pl = lab_data$P[1], Pt = 0, Pe = 0)
  init.state2 <- c(Be = lab_data2$Be[1], Bt = lab_data2$Bt[1], Bet = 0,
                   Pl = lab_data2$P[1], Pt = 0, Pe = 0)
  log.like.theta.bar <- dTrajObs(model, theta.bar, init.state, data = lab_data, 
                                 log = TRUE) +
    dTrajObs(model, theta.bar, init.state2, data = lab_data2, 
             log = TRUE)
  
  Dhat = -2*log.like.theta.bar
  
  Dbar = mean(-2*trace[,"log.density"])
  pD = Dbar - Dhat
  DIC = Dbar + pD
  
  DIC
  
}


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
  
  dic = DIC(trace_model, model, lab_data_trans3, lab_data_trans5)
  
  quants = apply(trace_model, 2, function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
  quants_s = c(models_to_try$model_name[i])
  for(j in 1:(ncol(quants)-1)){
    quants_s = c(quants_s, quants[2,j], quants[1,j], quants[3,j])
  }
  
  best_params = rbind(best_params,
                      c(quants_s, dic))
  
  
  #replicate 4
  init.state = c(Be = lab_data_trans$Be[1], Bt = lab_data_trans$Bt[1], Bet = 0,
                 Pl = lab_data_trans$P[1], Pe = 0, Pt = 0)
  
  traj = multi_run2(model, trace_model, init.state,
                    times = seq(0, 24, 1), nruns = 300, median = T, sampling_error = T)
  
  traj$model = models_to_try$model_name[i]
  all_traj_4 = rbind(all_traj_4, traj)
  
  
  #replicate 5
  init.state = c(Be = lab_data_trans5$Be[1], Bt = lab_data_trans5$Bt[1], Bet = 0,
                 Pl = lab_data_trans5$P[1], Pe = 0, Pt = 0)
  
  traj = multi_run2(model, trace_model, init.state,
                    times = seq(0, 24, 1), nruns = 300, median = T, sampling_error = T)
  
  traj$model = models_to_try$model_name[i]
  all_traj_5 = rbind(all_traj_5, traj)
  
  #replicate 3
  init.state = c(Be = lab_data_trans3$Be[1], Bt = lab_data_trans3$Bt[1], Bet = 0,
                 Pl = lab_data_trans3$P[1], Pe = 0, Pt = 0)
  
  traj = multi_run2(model, trace_model, init.state,
                    times = seq(0, 24, 1), nruns = 300, median = T, sampling_error = T)
  
  traj$model = models_to_try$model_name[i]
  all_traj_3 = rbind(all_traj_3, traj)
  
}

quants_names = c("model_name")
for(j in 1:(ncol(quants)-1)){
  quants_names = c(quants_names, colnames(quants)[j],
                   paste0(colnames(quants)[j], "_0.025"),
                   paste0(colnames(quants)[j], "_0.975"))
}

colnames(best_params) = c(quants_names, "DIC")

best_params[,-1] = apply(best_params[,-1], c(1,2), as.numeric)
best_params[,c(2:4,8:10)] = apply(best_params[,c(2:4,8:10)], c(1,2), function(x) 1/x)

best_params$DIC = round(best_params$DIC-(min(best_params$DIC)))

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

all_traj_3$dataset = "103"
all_traj_4$dataset = "104"
all_traj_5$dataset = "105"
all_traj_all = rbind(all_traj_3, all_traj_4, all_traj_5)
all_traj_all$dataset = as.factor(all_traj_all$dataset)
levels(all_traj_all$dataset) = c(expression(paste(10^3, " phage")),
                                expression(paste(10^4, " phage")),
                                expression(paste(10^5, " phage")))


lab_data_trans3M$dataset = "103"
lab_data_transM$dataset = "104"
lab_data_trans5M$dataset = "105"
lab_data_transMall = rbind(lab_data_trans3M, lab_data_transM, lab_data_trans5M)
lab_data_transMall$dataset = as.factor(lab_data_transMall$dataset)
levels(lab_data_transMall$dataset) = c(expression(paste(10^3, " phage")),
                                 expression(paste(10^4, " phage")),
                                 expression(paste(10^5, " phage")))

pp = ggplot() +
  geom_line(data = all_traj_all, aes(time, Pl, group = group, colour = model,
                                   linetype = "Model", size = model)) +
  geom_ribbon(data = all_traj_all %>% filter(model == "other"),
              aes(x = time, ymin = pmax(Pl - 1.96*Pl_sd, 0), ymax = Pl + 1.96*Pl_sd,
                  group= group, fill=model), alpha = 0.05) +
  geom_ribbon(data = all_traj_all %>% filter(model != "other"),
              aes(x = time, ymin = pmax(Pl - 1.96*Pl_sd, 0), ymax = Pl + 1.96*Pl_sd,
                  group= group, fill=model), alpha = 0.2) +
  geom_line(data = lab_data_transMall %>% filter(Bacteria == "P"), 
            aes(Time, Mean, linetype = "Data", colour = "PL"), size = 0.8) +
  geom_errorbar(data = lab_data_transMall %>% filter(Bacteria == "P"), 
                aes(x = Time, ymin = pmax(Mean - se, 0), ymax = Mean + se,
                    colour= "PL"), size = 0.8) +
  facet_wrap(~dataset, labeller = label_parsed)+
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
                      values = c("purple", "green", "black", "#c88a33")) +
  scale_fill_manual(breaks = c("dens_burst", "freq_burst", "other"),
                    labels = c("Density, link burst size to bacterial growth",
                               "Frequency, link burst size to bacterial growth", "Other"),
                    values = c("purple", "green", "black"),
                    guide = guide_legend(override.aes = list(size = c(0.8,0.8,0.5),
                                                             colour = c("purple", "green", "black")))) +
  scale_size_manual(breaks = c("dens_burst", "freq_burst", "other"),
                    labels = c("Density, link burst size to bacterial growth",
                               "Frequency, link burst size to bacterial growth", "Other"),
                    values = c(0.8,0.8,0.5)) +
  labs(x = "Time (hours)", y = "pfu per mL", linetype = "Source:",
       colour = "Organism:", fill = "Model type:", size = "Model type:") +
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
  geom_ribbon(data = all_traj_5 %>% filter(model == "other"),
              aes(x = time, ymin = pmax(Pl - 1.96*Pl_sd, 0), ymax = Pl + 1.96*Pl_sd,
                  group= group, fill=model), alpha = 0.05) +
  geom_ribbon(data = all_traj_5 %>% filter(model != "other"),
              aes(x = time, ymin = pmax(Pl - 1.96*Pl_sd, 0), ymax = Pl + 1.96*Pl_sd,
                  group= group, fill=model), alpha = 0.2) +
  geom_line(data = all_traj_5, aes(time, Bet, group = group, colour = model,
                                   linetype = "Model", size = model)) +
  geom_ribbon(data = all_traj_5 %>% filter(model == "other"),
              aes(x = time, ymin = pmax(Bet - 1.96*Bet_sd, 0), ymax = Bet + 1.96*Bet_sd,
                  group= group, fill=model), alpha = 0.05) +
  geom_ribbon(data = all_traj_5 %>% filter(model != "other"),
              aes(x = time, ymin = pmax(Pl - 1.96*Pl_sd, 0), ymax = Pl + 1.96*Pl_sd,
                  group= group, fill=model), alpha = 0.2) +
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
                      values = c("#c2484d", "purple", "green", "black", "#c88a33")) +
  scale_fill_manual(breaks = c("dens_burst", "freq_burst", "other"),
                    labels = c("Density, link burst size to bacterial growth",
                               "Frequency, link burst size to bacterial growth", "Other"),
                    values = c("purple", "green", "black"),
                    guide = guide_legend(override.aes = list(size = c(0.8,0.8,0.5),
                                                             linetype = c(2,2,2),
                                                             colour = c("purple", "green", "black")))) +
  scale_size_manual(breaks = c("dens_burst", "freq_burst", "other"),
                    labels = c("Density, link burst size to bacterial growth",
                               "Frequency, link burst size to bacterial growth", "Other"),
                    values = c(0.8,0.8,0.5)) +
  labs(x = "Time (hours)", y = "pfu per mL", linetype = "Source:",
       colour = "Organism:", fill = "Model type:", size = "Model type:") +
  guides(linetype = guide_legend(order = 1), colour = guide_legend(order = 2)) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.x = element_text(size=12))


pbet = ggplot() +
  geom_line(data = all_traj_all, aes(time, Bet, group = group, colour = model,
                                   linetype = "Model", size = model)) +
  geom_ribbon(data = all_traj_all, aes(x = time, ymin = pmax(Bet - 1.96*Bet_sd, 0), ymax = Bet + 1.96*Bet_sd,
                                     group = group, fill = model), alpha = 0.1) +
  geom_line(data = lab_data_transMall %>% filter(Bacteria == "DRP"), 
            aes(Time, Mean, linetype = "Data", colour = "BET"), size = 0.8) +
  geom_errorbar(data = lab_data_transMall %>% filter(Bacteria == "DRP"), 
                aes(x = Time, ymin = pmax(Mean - se, 0), ymax = Mean + se,
                    colour= "BET"), size = 0.8) +
  scale_x_continuous(breaks = seq(0,24,4)) +
  scale_linetype_manual(breaks = c("Data", "Model"),
                        labels = c("Data", "Model"),
                        values = c(1, 2)) +
  scale_colour_manual(breaks = c("BET"),
                      labels = c("Double resistant bacteria"),
                      values = c("#c2484d", "purple", "green", "black")) +
  scale_fill_manual(breaks = c("dens_burst", "freq_burst", "other"),
                    labels = c("Density, link burst size to bacterial growth",
                               "Frequency, link burst size to bacterial growth", "Other"),
                    values = c("purple", "green", "black"),
                    guide = guide_legend(override.aes = list(size = c(0.8,0.8,0.5),
                                                             colour = c("purple", "green", "black")))) +
  scale_size_manual(breaks = c("dens_burst", "freq_burst", "other"),
                    labels = c("Density, link burst size to bacterial growth",
                               "Frequency, link burst size to bacterial growth", "Other"),
                    values = c(0.8,0.8,0.5)) +
  facet_wrap(~dataset, labeller = label_parsed)+
  labs(x = "Time (hours)", y = "cfu per mL", linetype = "Source:",
       colour = "Organism:", fill = "Model type:", size = "Model type:") +
  guides(linetype = guide_legend(order = 1)) +
  theme_bw() +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x, n = 4),
                     labels=trans_format("log10", math_format(10^.x))) +
  coord_cartesian(ylim = c(1e-1, 1e2)) +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.x = element_text(size=12))

#final plot
legend = get_legend(p5_extra + theme(legend.position = "right", legend.box = "horizontal"))

fig4ab = plot_grid(pp + theme(legend.position = "none"),
          NULL,
          pbet + theme(legend.position = "none"),
          NULL, legend, nrow = 5, rel_heights = c(1,0.05,0.8,0.05,0.2),
          labels = c("a","", "b", "", ""))

fig4ab

ggsave("plot_compare_all.png", height = 10, width = 18)


all_traj_4 %>%
  filter(group == "freq_burst") %>%
  filter(time %in% lab_data_trans$time) %>%
  mutate(data_Pl = lab_data_trans$P[-1]) %>%
  mutate(diff = data_Pl/Pl)


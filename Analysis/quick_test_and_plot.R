
library(fitR)
library(dplyr)
library(reshape2)
library(ggplot2)
library(scales)
library(cowplot)
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
lab_data_transM = read.csv(here::here("Lab", "Transduction", "summary_10_4.csv")) %>%
  select(Time, Bacteria, Mean, se) %>%
  filter(Bacteria != "Total")
lab_data_trans5M = read.csv(here::here("Lab", "Transduction", "summary_10_5.csv")) %>%
  select(Time, Bacteria, Mean, se)%>%
  filter(Bacteria != "Total")
lab_data_trans3M = read.csv(here::here("Lab", "Transduction", "summary_10_3.csv")) %>%
  select(Time, Bacteria, Mean, se) %>%
  filter(Bacteria != "Total")

model = choose_model(model,
                     frequentist = F,
                     fixed_delay = NA,
                     decay = F,
                     link_beta = T,
                     link_L = T,
                     link_delay = F,
                     transduction = T, prop_Bet = T)

trace_model = all_params[["dens_both"]]
params = apply(trace_model, 2, median)
#params = c(beta = 7646387086.68, L = 91.19, gamma = 11.09, alpha = 46604726.39, tau = 0.41)


#params[["L"]] = 80
# params[["beta"]] = 5.9e9
# params[["alpha"]] = 7e6

#replicate 4
init.state = c(Be = lab_data_trans$Be[1], Bt = lab_data_trans$Bt[1], Bet = 0,
               Pl = lab_data_trans$P[1], Pe = 0, Pt = 0)

traj = model$simulate(params, init.state, seq(0,24,1))

traj4 = multi_run2(model, params, init.state,
                   times = seq(0, 24, 1), nruns = 10)

#replicate 5
init.state = c(Be = lab_data_trans5$Be[1], Bt = lab_data_trans5$Bt[1], Bet = 0,
               Pl = lab_data_trans5$P[1], Pe = 0, Pt = 0)

traj = model$simulate(params, init.state, seq(0,24,1))

traj5 = multi_run2(model, params, init.state,
                   times = seq(0, 24, 1), nruns = 10)

#replicate 3
init.state = c(Be = lab_data_trans3$Be[1], Bt = lab_data_trans3$Bt[1], Bet = 0,
               Pl = lab_data_trans3$P[1], Pe = 0, Pt = 0)

traj = model$simulate(params, init.state, seq(0,24,1))

traj3 = multi_run2(model, params, init.state,
                   times = seq(0, 24, 1), nruns = 10)

traj4[traj4<0.01] = 0.01
traj5[traj5<0.01] = 0.01
traj3[traj3<0.01] = 0.01


p4 = ggplot() +
  geom_line(data = traj4, aes(time, Pl, colour = "PL", linetype = "Model"), size = 0.8) +
  geom_ribbon(data = traj4, aes(x = time, ymin = pmax(Pl - 1.96*Pl_sd, 0), ymax = Pl + 1.96*Pl_sd,
                                fill = "PL"), alpha = 0.1) +
  geom_line(data = traj4, aes(time, Bet, colour = "BET", linetype = "Model"), size = 0.8) +
  geom_ribbon(data = traj4, aes(x = time, ymin = pmax(Bet - 1.96*Bet_sd, 0), ymax = Bet + 1.96*Bet_sd,
                                fill = "BET"), alpha = 0.1) +
  geom_line(data = lab_data_transM %>% filter(Bacteria == "P"), 
            aes(Time, Mean, linetype = "Data", colour = "PL")) +
  geom_errorbar(data = lab_data_transM %>% filter(Bacteria == "P"), 
                aes(x = Time, ymin = pmax(Mean - se, 0), ymax = Mean + se,
                    colour= "PL"), size = 0.8) +
  geom_line(data = lab_data_transM %>% filter(Bacteria == "DRP"), 
            aes(Time, Mean, linetype = "Data", colour = "BET")) +
  geom_errorbar(data = lab_data_transM %>% filter(Bacteria == "DRP"), 
                aes(x = Time, ymin = pmax(Mean - se, 0), ymax = Mean + se,
                    colour= "BET"), size = 0.8) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(breaks = seq(0,24,4)) +
  coord_cartesian(ylim = c(0.1, 1e11)) +
  scale_linetype_manual(breaks = c("Data", "Model"),
                        labels = c("Data", "Model"),
                        values = c(1, 2)) +
  scale_colour_manual(breaks = c("PL", "BET"),
                      labels = c("Lytic phage", "Double resistant progeny"),
                      values = c("#c88a33", "#c2484d")) +
  scale_fill_manual(breaks = c("PL", "BET"),
                    labels = c("Lytic phage", "Double resistant progeny"),
                    values = c("#c88a33", "#c2484d")) +
  labs(x = "Time (hours)", y = "cfu or pfu per mL", linetype = "Source:",
       colour = "Organism:", fill = "Organism:",
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
  geom_line(data = traj5, aes(time, Pl, colour = "PL", linetype = "Model"), size = 0.8) +
  geom_ribbon(data = traj5, aes(x = time, ymin = pmax(Pl - 1.96*Pl_sd, 0), ymax = Pl + 1.96*Pl_sd,
                                fill = "PL"), alpha = 0.1) +
  geom_line(data = traj5, aes(time, Bet, colour = "BET", linetype = "Model"), size = 0.8) +
  geom_ribbon(data = traj5, aes(x = time, ymin = pmax(Bet - 1.96*Bet_sd, 0), ymax = Bet + 1.96*Bet_sd,
                                fill = "BET"), alpha = 0.1) +
  geom_line(data = lab_data_trans5M %>% filter(Bacteria == "P"), 
            aes(Time, Mean, linetype = "Data", colour = "PL")) +
  geom_errorbar(data = lab_data_trans5M %>% filter(Bacteria == "P"), 
                aes(x = Time, ymin = pmax(Mean - se, 0), ymax = Mean + se,
                    colour= "PL"), size = 0.8) +
  geom_line(data = lab_data_trans5M %>% filter(Bacteria == "DRP"), 
            aes(Time, Mean, linetype = "Data", colour = "BET")) +
  geom_errorbar(data = lab_data_trans5M %>% filter(Bacteria == "DRP"), 
                aes(x = Time, ymin = pmax(Mean - se, 0), ymax = Mean + se,
                    colour= "BET"), size = 0.8) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(breaks = seq(0,24,4)) +
  coord_cartesian(ylim = c(0.1, 1e11)) +
  scale_linetype_manual(breaks = c("Data", "Model"),
                        labels = c("Data", "Model"),
                        values = c(1, 2)) +
  scale_colour_manual(breaks = c("PL", "BET"),
                      labels = c("Lytic phage", "Double resistant progeny"),
                      values = c("#c88a33", "#c2484d")) +
  scale_fill_manual(breaks = c("PL", "BET"),
                    labels = c("Lytic phage", "Double resistant progeny"),
                    values = c("#c88a33", "#c2484d")) +
  labs(x = "Time (hours)", y = "cfu or pfu per mL", linetype = "Source:",
       colour = "Organism:", fill = "Organism:",
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
  geom_line(data = lab_data_trans3M %>% filter(Bacteria == "P"), 
            aes(Time, Mean, linetype = "Data", colour = "PL")) +
  geom_errorbar(data = lab_data_trans3M %>% filter(Bacteria == "P"), 
                aes(x = Time, ymin = pmax(Mean - se, 0), ymax = Mean + se,
                    colour= "PL"), size = 0.8) +
  geom_line(data = lab_data_trans3M %>% filter(Bacteria == "DRP"), 
            aes(Time, Mean, linetype = "Data", colour = "BET")) +
  geom_errorbar(data = lab_data_trans3M %>% filter(Bacteria == "DRP"), 
                aes(x = Time, ymin = pmax(Mean - se, 0), ymax = Mean + se,
                    colour= "BET"), size = 0.8) +
  geom_line(data = traj3, aes(time, Pl, colour = "PL", linetype = "Model"), size = 0.8) +
  geom_ribbon(data = traj3, aes(x = time, ymin = pmax(Pl - 1.96*Pl_sd, 0), ymax = Pl + 1.96*Pl_sd,
                                fill = "PL"), alpha = 0.1) +
  geom_line(data = traj3, aes(time, Bet, colour = "BET", linetype = "Model"), size = 0.8) +
  geom_ribbon(data = traj3, aes(x = time, ymin = pmax(Bet - 1.96*Bet_sd, 0), ymax = Bet + 1.96*Bet_sd,
                                fill = "BET"), alpha = 0.1) +scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(breaks = seq(0,24,4)) +
  coord_cartesian(ylim = c(0.1, 1e11)) +
  scale_linetype_manual(breaks = c("Data", "Model"),
                        labels = c("Data", "Model"),
                        values = c(1, 2)) +
  scale_colour_manual(breaks = c("PL", "BET"),
                      labels = c("Lytic phage", "Double resistant progeny"),
                      values = c("#c88a33", "#c2484d")) +
  scale_fill_manual(breaks = c("PL", "BET"),
                    labels = c("Lytic phage", "Double resistant progeny"),
                    values = c("#c88a33", "#c2484d")) +
  labs(x = "Time (hours)", y = "cfu or pfu per mL", linetype = "Source:",
       colour = "Organism:", fill = "Organism:",
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

#ggsave("nice_model_output.png")


ggplot() +
  geom_line(data = traj5, aes(time, Pl, colour = "PL", linetype = "Model"), size = 0.8) +
  geom_ribbon(data = traj5, aes(x = time, ymin = pmax(Pl - 1.96*Pl_sd, 0), ymax = Pl + 1.96*Pl_sd,
                                fill = "PL"), alpha = 0.1) +
  geom_line(data = traj5, aes(time, Bet, colour = "BET", linetype = "Model"), size = 0.8) +
  geom_ribbon(data = traj5, aes(x = time, ymin = pmax(Bet - 1.96*Bet_sd, 0), ymax = Bet + 1.96*Bet_sd,
                                fill = "BET"), alpha = 0.1) +
  geom_line(data = lab_data_trans5M %>% filter(Bacteria == "P"), 
            aes(Time, Mean, linetype = "Data", colour = "PL")) +
  geom_errorbar(data = lab_data_trans5M %>% filter(Bacteria == "P"), 
                aes(x = Time, ymin = pmax(Mean - se, 0), ymax = Mean + se,
                    colour= "PL"), size = 0.8) +
  geom_line(data = lab_data_trans5M %>% filter(Bacteria == "DRP"), 
            aes(Time, Mean, linetype = "Data", colour = "BET")) +
  geom_errorbar(data = lab_data_trans5M %>% filter(Bacteria == "DRP"), 
                aes(x = Time, ymin = pmax(Mean - se, 0), ymax = Mean + se,
                    colour= "BET"), size = 0.8) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(breaks = seq(0,24,4)) +
  coord_cartesian(ylim = c(1e7, 1e11), xlim = c(6,24)) +
  scale_linetype_manual(breaks = c("Data", "Model"),
                        labels = c("Data", "Model"),
                        values = c(1, 2)) +
  scale_colour_manual(breaks = c("PL", "BET"),
                      labels = c("Lytic phage", ""),
                      values = c("#c88a33", NA)) +
  scale_fill_manual(breaks = c("PL", "BET"),
                    labels = c("Lytic phage",""),
                    values = c("#c88a33", NA)) +
  labs(x = "Time (hours)", y = "cfu or pfu per mL", linetype = "Source:",
       colour = "Organism:", fill = "Organism:",
       title = bquote(paste("C) ", 10^5, " phage"))) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.x = element_text(size=12))

#ggsave("zoomed_in.png")

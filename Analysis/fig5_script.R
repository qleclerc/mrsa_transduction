
library(fitR)
library(dplyr)
library(reshape2)
library(ggplot2)
library(scales)
library(cowplot)
library(openxlsx)

source(here::here("Model", "transduction_model_functions.R"))

model = readRDS(here::here("Model", "transduction_model.rds"))

files = list.files(here::here("Fitting", "Fitted_params"))

all_params = vector("list", length(files))
i=1

for(f in files){
  all_params[[i]] = read.csv(here::here("Fitting", "Fitted_params", f))
  i = i+1
}

names(all_params) = gsub(".csv", "", gsub("params_", "", files))

lab_data_trans = read.csv(here::here("Data", "transduction_summary_10_4.csv")) %>%
  select(Time, Bacteria, Mean) %>%
  dcast(Time~Bacteria) %>%
  select(-Total) %>%
  rename(time = Time, Bet = DRP, Be = EryR, Bt = TetR) %>%
  mutate(Be = round(Be),
         Bt = round(Bt),
         Bet = round(Bet),
         P = round(P))
lab_data_trans5 = read.csv(here::here("Data", "transduction_summary_10_5.csv")) %>%
  select(Time, Bacteria, Mean) %>%
  dcast(Time~Bacteria) %>%
  select(-Total) %>%
  rename(time = Time, Bet = DRP, Be = EryR, Bt = TetR) %>%
  mutate(Be = round(Be),
         Bt = round(Bt),
         Bet = round(Bet),
         P = round(P))
lab_data_trans3 = read.csv(here::here("Data", "transduction_summary_10_3.csv")) %>%
  select(Time, Bacteria, Mean) %>%
  dcast(Time~Bacteria) %>%
  select(-Total) %>%
  rename(time = Time, Bet = DRP, Be = EryR, Bt = TetR) %>%
  mutate(Be = round(Be),
         Bt = round(Bt),
         Bet = round(Bet),
         P = round(P))
lab_data_transM = read.csv(here::here("Data", "transduction_summary_10_4.csv")) %>%
  select(Time, Bacteria, Mean, se) %>%
  filter(Bacteria != "Total")
lab_data_trans5M = read.csv(here::here("Data", "transduction_summary_10_5.csv")) %>%
  select(Time, Bacteria, Mean, se)%>%
  filter(Bacteria != "Total")
lab_data_trans3M = read.csv(here::here("Data", "transduction_summary_10_3.csv")) %>%
  select(Time, Bacteria, Mean, se) %>%
  filter(Bacteria != "Total")

model = choose_model(model,
                     frequentist = T,
                     fixed_delay = NA,
                     decay = F,
                     link_beta = F,
                     link_L = T,
                     link_delay = F,
                     transduction = T, fig5 = T)

trace_model = all_params[["freq_burst"]]
params = apply(trace_model, 2, median)

#replicate 4
init.state = c(Be = lab_data_trans$Be[1], Bt = lab_data_trans$Bt[1], Bet = 0,
               Pl = lab_data_trans$P[1], Pe = 0, Pt = 0)

traj4 = model$simulate(params, init.state, seq(0,24,1))
traj4$init_phage = "104"
traj4$rel_Pl = c(0,diff(traj4$Pl)/traj4$Pl[-25])
traj4$rel_Pt = c(0,diff(traj4$Pt+traj4$Pe)/(traj4$Pt+traj4$Pe)[-25])
traj4$B_tot = traj4$Be+traj4$Bt+traj4$Bet
traj4$rel_B = c(0,diff(traj4$B_tot)/(traj4$B_tot)[-25])


#replicate 5
init.state = c(Be = lab_data_trans5$Be[1], Bt = lab_data_trans5$Bt[1], Bet = 0,
               Pl = lab_data_trans5$P[1], Pe = 0, Pt = 0)

traj5 = model$simulate(params, init.state, seq(0,24,1))
traj5$init_phage = "105"
traj5$rel_Pl = c(0,diff(traj5$Pl)/traj5$Pl[-25])
traj5$rel_Pt = c(0,diff(traj5$Pt+traj5$Pe)/(traj5$Pt+traj5$Pe)[-25])
traj5$B_tot = traj5$Be+traj5$Bt+traj4$Bet
traj5$rel_B = c(0,diff(traj5$B_tot)/(traj5$B_tot)[-25])


#replicate 3
init.state = c(Be = lab_data_trans3$Be[1], Bt = lab_data_trans3$Bt[1], Bet = 0,
               Pl = lab_data_trans3$P[1], Pe = 0, Pt = 0)

traj3 = model$simulate(params, init.state, seq(0,24,1))
traj3$init_phage = "103"
traj3$rel_Pl = c(0,diff(traj3$Pl)/traj3$Pl[-25])
traj3$rel_Pt = c(0,diff(traj3$Pt+traj3$Pe)/(traj3$Pt+traj3$Pe)[-25])
traj3$B_tot = traj3$Be+traj3$Bt+traj3$Bet
traj3$rel_B = c(0,diff(traj3$B_tot)/(traj3$B_tot)[-25])


#bind together
traj_all = rbind(traj3, traj4, traj5)
traj_all$rel_Pt[c(1,2,26,27,51,52)] = NA

traj_all$init_phage = as.factor(traj_all$init_phage)
levels(traj_all$init_phage) = c(expression(paste(10^3, " phage")),
                                expression(paste(10^4, " phage")),
                                expression(paste(10^5, " phage")))




pa = ggplot(traj_all) +
  geom_line(aes(time, 75*(1-(Be+Bt+Bet)/2.763e+09)+1), size = 1) +
  facet_wrap(~init_phage, labeller = label_parsed) +
  scale_y_continuous(breaks = seq(0,75,15)) +
  scale_x_continuous(breaks = seq(0,24,4), limits = c(0,24)) +
  labs(x = "Time (hours)", y = "Phage burst size",
       title = "a") +
  theme_bw() +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.x = element_text(size=12),
        plot.title = element_text(face = "bold"))

pb = ggplot(traj_all) +
  geom_line(aes(time, rel_Pl, colour = "PL"), size = 1) +
  geom_line(aes(time, rel_Pt, colour = "PT"), size = 1) +
  geom_line(aes(time, rel_B, colour = "B"), size = 1) +
  geom_hline(aes(yintercept = 0), size = 1, linetype = "dashed", alpha = 0.5) +
  facet_wrap(~init_phage, labeller = label_parsed) +
  scale_x_continuous(breaks = seq(0,24,4), limits = c(1,24)) +
  scale_y_continuous(breaks = c(-2,0,2,4,6,8,10,12), limits = c(-2,12)) +
  scale_colour_manual(breaks = c("PL", "PT", "B"),
                      labels = c("Lytic phage", "Transducing phage", "Bacteria"),
                      values = c("#c88a33", "green4", "blue")) +
  labs(x = "Time (hours)", y = "Relative change in pfu or cfu per mL",
       colour = "Organism:", title = "b") +
  theme_bw() +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.x = element_text(size=12),
        plot.title = element_text(face = "bold"))

pc = ggplot(traj_all)+
  geom_line(aes(time, Pl_inc, linetype = init_phage, colour = "PL"), size = 1) +
  geom_line(aes(time, Pe_inc+Pt_inc, linetype = init_phage, colour = "PT"), size = 1) +
  scale_x_continuous(breaks = seq(0,24,4), limits = c(1,24)) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x, n=8),
                     labels=trans_format("log10", math_format(10^.x))) +
  labs(x = "Time (hours)", y = "Phage incidence (pfu per mL)", colour = "Organism:",
       linetype = "Starting phage\nconcentration:", title = "c") +
  scale_colour_manual(breaks = c("PL", "PT"),
                      labels = c("Lytic phage", "Transducing phage"),
                      values = c("#c88a33", "green4")) +
  scale_linetype_manual(values = c(1,2,3),
                        labels = c(bquote(10^3),
                                   bquote(10^4),
                                   bquote(10^5))) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.x = element_text(size=12),
        plot.title = element_text(face = "bold")) +
  guides(colour = guide_legend(order = 1))

pd = ggplot(traj_all)+
  geom_line(aes(time, Bet_tr/(Bet_tr+Bet_gr), linetype = init_phage), size = 1) +
  geom_line(aes(time, Bet_tr/(Bet_tr+Bet_gr), linetype = init_phage), size = 1) +
  geom_line(aes(time, Bet_tr/(Bet_tr+Bet_gr), linetype = init_phage), size = 1) +
  scale_x_continuous(breaks = seq(0,24,4), limits = c(1,24)) +
  scale_y_continuous(limits = c(0.5,1)) +
  labs(x = "Time (hours)", y = "Fraction of DRP\ngenerated by transduction",
       linetype = "Starting phage\nconcentration:", title = "d") +
  scale_linetype_manual(values = c(1,2,3),
                        labels = c(bquote(10^3),
                                   bquote(10^4),
                                   bquote(10^5))) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.x = element_text(size=12),
        plot.title = element_text(face = "bold"))


legendc = get_legend(pc + theme(legend.position = "bottom", legend.box = "horizontal"))

plot_grid(pa,
          pb + theme(legend.position = "bottom", legend.box = "horizontal"),
          plot_grid(pc + theme(legend.position = "none"),
                    pd + theme(legend.position = "none"),
                    nrow = 1,
                    rel_widths = c(1,1)),
          NULL,
          legendc,
          nrow = 5, rel_heights = c(0.8,1.2,1,0.05,0.1))


ggsave(here::here("Figures","fig5.png"), height = 10, width = 9)


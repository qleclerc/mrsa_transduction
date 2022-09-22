
#fig4ab required below comes from "analyse_fitted_models.R" script!

library(fitR)
library(dplyr)
library(reshape2)
library(ggplot2)
library(cowplot)
library(scales)
library(openxlsx)

source(here::here("Model", "model.R"))


files = list.files(here::here("Fitting", "Fitted_params"))

all_params = vector("list", length(files))
i=1

for(f in files){
  all_params[[i]] = read.csv(here::here("Fitting", "Fitted_params", f))
  i = i+1
}

names(all_params) = gsub(".csv", "", gsub("params_", "", files))


data = read.xlsx(here::here("Data", "varying_MOI_data.xlsx"))

data = apply(data, c(1,2), as.numeric)
data = data[-1, (-c(20:28))]


data_final = data.frame(init_bac = rep(4e6, nrow(data)))
data_final$init_pha = data_final$init_bac*data[,1]
data_final$Be = rowMeans(data[,c(2:10)])
data_final$Be_sd = apply(data[,c(2:10)], 1, function(x) sd(x)/sqrt(3))
data_final$Bt = rowMeans(data[,c(11:19)])
data_final$Bt_sd = apply(data[,c(11:19)], 1, function(x) sd(x)/sqrt(3))
data_final$Bet = rowMeans(data[,c(20:28)])
data_final$Bet_sd = apply(data[,c(20:28)], 1, function(x) sd(x)/sqrt(3))
data_final$Pl = rowMeans(data[,c(29:37)])
data_final$Pl_sd = apply(data[,c(29:37)], 1, function(x) sd(x)/sqrt(3))

data = data_final

models_to_try = data.frame(model_name="dens_beta",
                           mode = "dens", link_L=FALSE, link_beta=TRUE)
models_to_try = rbind(models_to_try,
                      data.frame(model_name="dens_burst",
                                 mode = "dens", link_L=T, link_beta=F))
models_to_try = rbind(models_to_try,
                      data.frame(model_name="dens_both",
                                 mode = "dens", link_L=T, link_beta=TRUE))
models_to_try = rbind(models_to_try,
                      data.frame(model_name="hill_beta",
                                 mode = "hill", link_L=FALSE, link_beta=TRUE))
models_to_try = rbind(models_to_try,
                      data.frame(model_name="hill_burst",
                                 mode = "hill", link_L=T, link_beta=F))
models_to_try = rbind(models_to_try,
                      data.frame(model_name="hill_both",
                                 mode = "hill", link_L=T, link_beta=TRUE))

## REPEAT WITH OTHER MODELS FOR SUPP MAT

data$model = "data"
all_results = data

for(i in 1:nrow(models_to_try)){
  
  cat("\nWorking on", models_to_try$model_name[i])
  
  trace_model = all_params[[models_to_try$model_name[i]]]
  
  params = apply(trace_model, 2, median)
  
  data_model = data
  
  for(j in 1:nrow(data_model)){
    
    init.state = c(Be = data_model$init_bac[j]/2, Bt = data_model$init_bac[j]/2, Bet = 0,
                   Pl = data_model$init_pha[j], Pe = 0, Pt = 0)
    
    traj = multi_run2(trace_model, init.state, mode = models_to_try$mode[i],
                      link_L = models_to_try$link_L[i], link_beta = models_to_try$link_beta[i],
                      times = seq(0, 24, 1), nruns = 300, median = T, sampling_error = T)
    
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


all_results_L = all_results %>%
  filter(model %in% c("hill_burst", "data", "dens_burst"))
all_results_beta = all_results %>%
  filter(model %in% c("hill_beta", "data", "dens_beta"))
all_results_both = all_results %>%
  filter(model %in% c("hill_both", "data", "dens_both"))

all_results_L$model = as.factor(all_results_L$model)
all_results_L$model = factor(all_results_L$model, levels(all_results_L$model)[c(2,1,3)])
levels(all_results_L$model) = c("Density model", "Data", "Saturated model")

all_results_beta$model = as.factor(all_results_beta$model)
all_results_beta$model = factor(all_results_beta$model, levels(all_results_beta$model)[c(2,1,3)])
levels(all_results_beta$model) = c("Density model", "Data", "Saturated model")

all_results_both$model = as.factor(all_results_both$model)
all_results_both$model = factor(all_results_both$model, levels(all_results_both$model)[c(2,1,3)])
levels(all_results_both$model) = c("Density model", "Data", "Saturated model")


d = ggplot(all_results_L %>% filter(model == "Data")) + 
  geom_point(aes(x=10^4, y=10^13, colour = "Bacteria:"))+
  geom_point(aes(x=10^4, y=10^13, colour = "Phage:"))+
  geom_pointrange(aes(x = init_pha, y = Be, ymin = pmax(Be - Be_sd, 0.01), ymax = Be + Be_sd, colour = "Be"),
                  position = position_jitter(height=0, width=0.05),
                  shape = 19, alpha = 1, size = 0.5, linetype = "solid") +
  geom_pointrange(aes(x = init_pha, y = Bt, ymin = pmax(Bt - Bt_sd, 0.01), ymax = Bt + Bt_sd, colour = "Bt"),
                  position = position_jitter(height=0, width=0.05),
                  shape = 19, alpha = 1, size = 0.5, linetype = "solid") +
  geom_pointrange(aes(x = init_pha, y = Bet, ymin = pmax(Bet - Bet_sd, 0.01), ymax = Bet + Bet_sd, colour = "Bet"),
                  position = position_jitter(height=0, width=0.05),
                  shape = 19, alpha = 1, size = 0.5, linetype = "solid") +
  geom_pointrange(aes(x = init_pha, y = Pl, ymin = pmax(Pl - Pl_sd, 0.01), ymax = Pl + Pl_sd, colour = "Pl"),
                  position = position_jitter(height=0, width=0.05),
                  shape = 17, alpha = 1, size = 0.5, linetype = "solid") +
  scale_x_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x, n=4),
                     labels=trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x, n = 6),
                     labels=c("x", "0", expression(10^0),
                              expression(10^2), expression(10^4), expression(10^6),
                              expression(10^8), expression(10^10), expression(10^12), "x")) +
  coord_cartesian(ylim = c(0.01, 1e12), xlim = c(4.5e3, 1.5e7)) +
  facet_grid(~model) +
  labs(colour = "", x="Starting phage (pfu per mL)", y="cfu or pfu per mL after 24h") +
  theme_bw() +
  scale_colour_manual(breaks= c("Bacteria:","Be", "Bt", "Bet", "Phage:", "Pl"),
                      values=c("white","#685cc4","#6db356","#c2484d","white","#c88a33"),
                      labels = c("Bacteria:",
                                 bquote(B[E]),
                                 bquote(B[T]),
                                 bquote(B[ET]),
                                 "Phage:",
                                 bquote(P[L])),
                      guide = guide_legend(override.aes = list(shape = c(19,19,19,19,19,17)))) +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text.x = element_text(size=12, color = "white", face = "bold"),
        strip.background = element_rect(fill="gray48"))


burst = ggplot(all_results_L %>% filter(model != "Data")) + 
  geom_point(aes(x=10^4, y=10^13, colour = "Bacteria:"))+
  geom_point(aes(x=10^4, y=10^13, colour = "Phage:"))+
  geom_pointrange(aes(x = init_pha, y = Be, ymin = pmax(Be - Be_sd, 0.01), ymax = Be + Be_sd, colour = "Be"),
                  position = position_jitter(height=0, width=0.05),
                  shape = 19, alpha = 1, size = 0.5, linetype = "solid") +
  geom_pointrange(aes(x = init_pha, y = Bt, ymin = pmax(Bt - Bt_sd, 0.01), ymax = Bt + Bt_sd, colour = "Bt"),
                  position = position_jitter(height=0, width=0.05),
                  shape = 19, alpha = 1, size = 0.5, linetype = "solid") +
  geom_pointrange(aes(x = init_pha, y = Bet, ymin = pmax(Bet - Bet_sd, 0.01), ymax = Bet + Bet_sd, colour = "Bet"),
                  position = position_jitter(height=0, width=0.05),
                  shape = 19, alpha = 1, size = 0.5, linetype = "solid") +
  geom_pointrange(aes(x = init_pha, y = Pl, ymin = pmax(Pl - Pl_sd, 0.01), ymax = Pl + Pl_sd, colour = "Pl"),
                  position = position_jitter(height=0, width=0.05),
                  shape = 17, alpha = 1, size = 0.5, linetype = "solid") +
  scale_x_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x, n = 4),
                     labels=trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x, n = 6),
                     labels=c("x", "0", expression(10^0),
                              expression(10^2), expression(10^4), expression(10^6),
                              expression(10^8), expression(10^10), expression(10^12), "x")) +
  coord_cartesian(ylim = c(0.01, 1e12), xlim = c(4.5e3, 1.5e7)) +
  facet_grid(~model) +
  labs(colour = "", x="Starting phage (pfu per mL)", y="cfu or pfu per mL after 24h") +
  theme_bw() +
  scale_colour_manual(breaks= c("Bacteria:","Be", "Bt", "Bet", "Phage:", "Pl"),
                      values=c("white","#685cc4","#6db356","#c2484d","white","#c88a33"),
                      labels = c("Bacteria:",
                                 bquote(B[E]),
                                 bquote(B[T]),
                                 bquote(B[ET]),
                                 "Phage:",
                                 bquote(P[L])),
                      guide = guide_legend(override.aes = list(shape = c(19,19,19,19,19,17)))) +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text.x = element_text(size=12, face = "bold"))

legend = get_legend(d + theme(legend.position = "right", legend.box = "vertical"))
fig4c = plot_grid(burst + theme(legend.position = "none"),
                  d + theme(legend.position = "none"),
                  legend,
                  rel_widths = c(1,0.56,0.2),
                  labels = c("c", "", ""),
                  nrow = 1)

#reminder, fig4ab comes from "analyse_fitted_models.R" script!
plot_grid(fig4ab, NULL, fig4c, nrow = 3, rel_heights = c(2,0.05,1))
ggsave(here::here("Figures", "fig4.png"), height = 12, width = 15, dpi=900)



beta = ggplot(all_results_beta %>% filter(model != "Data")) + 
  geom_point(aes(x=10^4, y=10^13, colour = "Bacteria:"))+
  geom_point(aes(x=10^4, y=10^13, colour = "Phage:"))+
  geom_pointrange(aes(x = init_pha, y = Be, ymin = pmax(Be - Be_sd, 0.01), ymax = Be + Be_sd, colour = "Be"),
                  position = position_jitter(height=0, width=0.05),
                  shape = 19, alpha = 1, size = 0.5, linetype = "solid") +
  geom_pointrange(aes(x = init_pha, y = Bt, ymin = pmax(Bt - Bt_sd, 0.01), ymax = Bt + Bt_sd, colour = "Bt"),
                  position = position_jitter(height=0, width=0.05),
                  shape = 19, alpha = 1, size = 0.5, linetype = "solid") +
  geom_pointrange(aes(x = init_pha, y = Bet, ymin = pmax(Bet - Bet_sd, 0.01), ymax = Bet + Bet_sd, colour = "Bet"),
                  position = position_jitter(height=0, width=0.05),
                  shape = 19, alpha = 1, size = 0.5, linetype = "solid") +
  geom_pointrange(aes(x = init_pha, y = Pl, ymin = pmax(Pl - Pl_sd, 0.01), ymax = Pl + Pl_sd, colour = "Pl"),
                  position = position_jitter(height=0, width=0.05),
                  shape = 17, alpha = 1, size = 0.5, linetype = "solid") +
  scale_x_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x, n = 4),
                     labels=trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x, n = 6),
                     labels=c("x", "0", expression(10^0),
                              expression(10^2), expression(10^4), expression(10^6),
                              expression(10^8), expression(10^10), expression(10^12), "x")) +
  coord_cartesian(ylim = c(0.01, 1e12), xlim = c(4.5e3, 1.5e7)) +
  facet_grid(~model) +
  labs(colour = "", x="Starting phage (pfu per mL)", y="cfu or pfu per mL after 24h") +
  theme_bw() +
  scale_colour_manual(breaks= c("Bacteria:","Be", "Bt", "Bet", "Phage:", "Pl"),
                      values=c("white","#685cc4","#6db356","#c2484d","white","#c88a33"),
                      labels = c("Bacteria:",
                                 bquote(B[E]),
                                 bquote(B[T]),
                                 bquote(B[ET]),
                                 "Phage:",
                                 bquote(P[L])),
                      guide = guide_legend(override.aes = list(shape = c(19,19,19,19,19,17)))) +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text.x = element_text(size=12, face = "bold"))

legend = get_legend(d + theme(legend.position = "right", legend.box = "vertical"))
pbeta = plot_grid(beta + theme(legend.position = "none"),
                  d + theme(legend.position = "none"),
                  legend,
                  rel_widths = c(1,0.56,0.2),
                  nrow = 1)

both = ggplot(all_results_both %>% filter(model != "Data")) + 
  geom_point(aes(x=10^4, y=10^13, colour = "Bacteria:"))+
  geom_point(aes(x=10^4, y=10^13, colour = "Phage:"))+
  geom_pointrange(aes(x = init_pha, y = Be, ymin = pmax(Be - Be_sd, 0.01), ymax = Be + Be_sd, colour = "Be"),
                  position = position_jitter(height=0, width=0.05),
                  shape = 19, alpha = 1, size = 0.5, linetype = "solid") +
  geom_pointrange(aes(x = init_pha, y = Bt, ymin = pmax(Bt - Bt_sd, 0.01), ymax = Bt + Bt_sd, colour = "Bt"),
                  position = position_jitter(height=0, width=0.05),
                  shape = 19, alpha = 1, size = 0.5, linetype = "solid") +
  geom_pointrange(aes(x = init_pha, y = Bet, ymin = pmax(Bet - Bet_sd, 0.01), ymax = Bet + Bet_sd, colour = "Bet"),
                  position = position_jitter(height=0, width=0.05),
                  shape = 19, alpha = 1, size = 0.5, linetype = "solid") +
  geom_pointrange(aes(x = init_pha, y = Pl, ymin = pmax(Pl - Pl_sd, 0.01), ymax = Pl + Pl_sd, colour = "Pl"),
                  position = position_jitter(height=0, width=0.05),
                  shape = 17, alpha = 1, size = 0.5, linetype = "solid") +
  scale_x_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x, n= 4),
                     labels=trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x, n = 6),
                     labels=c("x", "0", expression(10^0),
                              expression(10^2), expression(10^4), expression(10^6),
                              expression(10^8), expression(10^10), expression(10^12), "x")) +
  coord_cartesian(ylim = c(0.01, 1e12), xlim = c(4.5e3, 1.5e7)) +
  facet_grid(~model) +
  labs(colour = "", x="Starting phage (pfu per mL)", y="cfu or pfu per mL after 24h") +
  theme_bw() +
  scale_colour_manual(breaks= c("Bacteria:","Be", "Bt", "Bet", "Phage:", "Pl"),
                      values=c("white","#685cc4","#6db356","#c2484d","white","#c88a33"),
                      labels = c("Bacteria:",
                                 bquote(B[E]),
                                 bquote(B[T]),
                                 bquote(B[ET]),
                                 "Phage:",
                                 bquote(P[L])),
                      guide = guide_legend(override.aes = list(shape = c(19,19,19,19,19,17)))) +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text.x = element_text(size=12, face = "bold"))

legend = get_legend(d + theme(legend.position = "right", legend.box = "vertical"))
pboth = plot_grid(both + theme(legend.position = "none"),
                  d + theme(legend.position = "none"),
                  legend,
                  rel_widths = c(1,0.56,0.2),
                  nrow = 1)

plot_grid(pbeta, pboth, nrow = 2, labels = c("a", "b"))
ggsave(here::here("Figures", "supp_fig4.png"), height = 10, width = 9, dpi=600)


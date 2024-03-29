
# SETUP ####

library(fitR)
library(coda)
library(dplyr)
library(reshape2)
library(ggplot2)
library(cowplot)
library(scales)

source(here::here("Model", "transduction_model_functions.R"))

model = readRDS(here::here("Model", "transduction_model.rds"))

lab_data_transM = read.csv(here::here("Data", "transduction_summary_10_4.csv")) %>%
  select(Time, Bacteria, Mean, se) %>%
  filter(Bacteria != "Total")
lab_data_trans5M = read.csv(here::here("Data", "transduction_summary_10_5.csv")) %>%
  select(Time, Bacteria, Mean, se)%>%
  filter(Bacteria != "Total")
lab_data_trans3M = read.csv(here::here("Data", "transduction_summary_10_3.csv")) %>%
  select(Time, Bacteria, Mean, se) %>%
  filter(Bacteria != "Total")


lab_data_trans4 = read.csv(here::here("Data", "transduction_summary_10_4.csv")) %>%
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


# FIT PHAGE #####

models_to_try = data.frame(model_name="mass_decay_link_beta", frequentist=FALSE,
                           fixed_delay=NA, decay=TRUE,
                           link_beta=TRUE, link_L=FALSE, link_delay=FALSE, transduction=TRUE)
models_to_try = rbind(models_to_try,
                      data.frame(model_name="frequentist_decay_link_beta", frequentist=TRUE,
                                 fixed_delay=NA, decay=TRUE,
                                 link_beta=TRUE, link_L=FALSE, link_delay=FALSE, transduction=TRUE))

all_theta = vector("list", nrow(models_to_try))

for(i in 1:nrow(models_to_try)){

  cat("Working on", models_to_try$model_name[i], "\n")

  model = choose_model(model,
                       frequentist = models_to_try$frequentist[i],
                       fixed_delay = models_to_try$fixed_delay[i],
                       decay = models_to_try$decay[i],
                       link_beta = models_to_try$link_beta[i],
                       link_L = models_to_try$link_L[i],
                       link_delay = models_to_try$link_delay[i],
                       transduction = models_to_try$transduction[i])

  init.theta = c(beta = 4e9, L = 20, alpha = 9e5, tau = 0.4)
  mcmc_fit = run_mcmc(model, lab_data_trans3, lab_data_trans5,
                      init.theta = init.theta,
                      proposal.sd = c(init.theta[1]/300,
                                      init.theta[2]/300,
                                      init.theta[3]/300,
                                      init.theta[4]/300),
                      n.iterations = 100000,
                      adapt.size.start = 1000,
                      adapt.shape.start = 1100,
                      adapt.size.cooling = 0.999)


  #replicate 4
  init.state = c(Be = lab_data_trans4$Be[1], Bt = lab_data_trans4$Bt[1], Bet = 0,
                 Pl = lab_data_trans4$P[1], Pe = 0, Pt = 0)
  theta = mcmc_fit$trace[which.max(mcmc_fit$trace[,"log.density"]),]

  traj = model$simulate(theta, init.state, times = seq(0, 30, 1))

  p4 = ggplot() +
    geom_line(data = traj, aes(time, Be, colour = "Model", linetype = "Bacteria"), size = 0.8) +
    geom_line(data = traj, aes(time, Bt, colour = "Model", linetype = "Bacteria"), size = 0.8) +
    geom_line(data = traj, aes(time, Bet, colour = "Model", linetype = "Bacteria"), size = 0.8) +
    geom_line(data = traj, aes(time, Pl, colour = "Model", linetype = "Phage"), size = 0.8) +
    geom_line(data = lab_data_transM %>% filter(Bacteria != "P"),
              aes(Time, Mean, group = Bacteria, colour = "Data", linetype = "Bacteria"), size = 0.8) +
    geom_errorbar(data = lab_data_transM %>% filter(Bacteria != "P"),
                  aes(x = Time, ymin = pmax(Mean - se, 0), ymax = Mean + se, group = Bacteria,
                      colour= "Data", linetype = "Bacteria"), size = 0.7) +
    geom_line(data = lab_data_transM %>% filter(Bacteria == "P"),
              aes(Time, Mean, group = Bacteria, colour = "Data", linetype = "Phage"), size = 0.8) +
    geom_errorbar(data = lab_data_transM %>% filter(Bacteria == "P"),
                  aes(x = Time, ymin = pmax(Mean - se, 0), ymax = Mean + se, group = Bacteria,
                      colour= "Data", linetype = "Phage"), size = 0.7) +
    scale_y_continuous(trans=log10_trans(),
                       breaks=trans_breaks("log10", function(x) 10^x),
                       labels=trans_format("log10", math_format(10^.x)),
                       limits = c(0.1, 3e10)) +
    scale_x_continuous(breaks = seq(0,30,5)) +
    labs(y = "cfu or pfu per mL", x = "Time (hours)", title = "10^4",
         linetype = "Organism:", colour = "Source:") +
    theme_bw() +
    theme(axis.text.x = element_text(size=12),
          axis.title.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=12),
          legend.text = element_text(size=12),
          legend.title = element_text(size=12),
          strip.text.x = element_text(size=12)) +
    scale_colour_manual(values=c("#685cc4","#6db356"))



  #replicate 5
  init.state = c(Be = lab_data_trans5$Be[1], Bt = lab_data_trans5$Bt[1], Bet = 0,
                 Pl = lab_data_trans5$P[1], Pe = 0, Pt = 0)
  theta = mcmc_fit$trace[which.max(mcmc_fit$trace[,"log.density"]),]

  traj = model$simulate(theta, init.state, times = seq(0, 30, 1))

  p5 = ggplot() +
    geom_line(data = traj, aes(time, Be, colour = "Model", linetype = "Bacteria"), size = 0.8) +
    geom_line(data = traj, aes(time, Bt, colour = "Model", linetype = "Bacteria"), size = 0.8) +
    geom_line(data = traj, aes(time, Bet, colour = "Model", linetype = "Bacteria"), size = 0.8) +
    geom_line(data = traj, aes(time, Pl, colour = "Model", linetype = "Phage"), size = 0.8) +
    geom_line(data = lab_data_trans5M %>% filter(Bacteria != "P"),
              aes(Time, Mean, group = Bacteria, colour = "Data", linetype = "Bacteria"), size = 0.8) +
    geom_errorbar(data = lab_data_trans5M %>% filter(Bacteria != "P"),
                  aes(x = Time, ymin = pmax(Mean - se, 0), ymax = Mean + se, group = Bacteria,
                      colour= "Data", linetype = "Bacteria"), size = 0.7) +
    geom_line(data = lab_data_trans5M %>% filter(Bacteria == "P"),
              aes(Time, Mean, group = Bacteria, colour = "Data", linetype = "Phage"), size = 0.8) +
    geom_errorbar(data = lab_data_trans5M %>% filter(Bacteria == "P"),
                  aes(x = Time, ymin = pmax(Mean - se, 0), ymax = Mean + se, group = Bacteria,
                      colour= "Data", linetype = "Phage"), size = 0.7) +
    scale_y_continuous(trans=log10_trans(),
                       breaks=trans_breaks("log10", function(x) 10^x),
                       labels=trans_format("log10", math_format(10^.x)),
                       limits = c(0.1, 3e10)) +
    scale_x_continuous(breaks = seq(0,30,5)) +
    labs(y = "cfu or pfu per mL", x = "Time (hours)", title = "10^5",
         linetype = "Organism:", colour = "Source:") +
    theme_bw() +
    theme(axis.text.x = element_text(size=12),
          axis.title.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=12),
          legend.text = element_text(size=12),
          legend.title = element_text(size=12),
          strip.text.x = element_text(size=12)) +
    scale_colour_manual(values=c("#685cc4","#6db356"))


  #replicate 3
  init.state = c(Be = lab_data_trans3$Be[1], Bt = lab_data_trans3$Bt[1], Bet = 0,
                 Pl = lab_data_trans3$P[1], Pe = 0, Pt = 0)
  theta = mcmc_fit$trace[which.max(mcmc_fit$trace[,"log.density"]),]

  traj = model$simulate(theta, init.state, times = seq(0, 30, 1))

  p3 = ggplot() +
    geom_line(data = traj, aes(time, Be, colour = "Model", linetype = "Bacteria"), size = 0.8) +
    geom_line(data = traj, aes(time, Bt, colour = "Model", linetype = "Bacteria"), size = 0.8) +
    geom_line(data = traj, aes(time, Bet, colour = "Model", linetype = "Bacteria"), size = 0.8) +
    geom_line(data = traj, aes(time, Pl, colour = "Model", linetype = "Phage"), size = 0.8) +
    geom_line(data = lab_data_trans3M %>% filter(Bacteria != "P"),
              aes(Time, Mean, group = Bacteria, colour = "Data", linetype = "Bacteria"), size = 0.8) +
    geom_errorbar(data = lab_data_trans3M %>% filter(Bacteria != "P"),
                  aes(x = Time, ymin = pmax(Mean - se, 0), ymax = Mean + se, group = Bacteria,
                      colour= "Data", linetype = "Bacteria"), size = 0.7) +
    geom_line(data = lab_data_trans3M %>% filter(Bacteria == "P"),
              aes(Time, Mean, group = Bacteria, colour = "Data", linetype = "Phage"), size = 0.8) +
    geom_errorbar(data = lab_data_trans3M %>% filter(Bacteria == "P"),
                  aes(x = Time, ymin = pmax(Mean - se, 0), ymax = Mean + se, group = Bacteria,
                      colour= "Data", linetype = "Phage"), size = 0.7) +
    scale_y_continuous(trans=log10_trans(),
                       breaks=trans_breaks("log10", function(x) 10^x),
                       labels=trans_format("log10", math_format(10^.x)),
                       limits = c(0.1, 3e10)) +
    scale_x_continuous(breaks = seq(0,30,5)) +
    labs(y = "cfu or pfu per mL", x = "Time (hours)", title = "10^3",
         linetype = "Organism:", colour = "Source:") +
    theme_bw() +
    theme(axis.text.x = element_text(size=12),
          axis.title.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=12),
          legend.text = element_text(size=12),
          legend.title = element_text(size=12),
          strip.text.x = element_text(size=12)) +
    scale_colour_manual(values=c("#685cc4","#6db356"))


  #final plot
  legend = get_legend(p4 + theme(legend.position = "bottom", legend.box = "vertical"))

  plot_grid(p3 + theme(legend.position = "none"),
            p4 + theme(legend.position = "none"),
            p5 + theme(legend.position = "none"),
            legend)

  filename = paste0(models_to_try$model_name[i], "_b.png")
  ggsave(here::here("Fitting", "Full_chains", "Best_fits", filename))

  all_theta[[i]] = mcmc_fit$trace
  names(all_theta)[i] = models_to_try$model_name[i]

}

saveRDS(all_theta, here::here("Fitting", "Full_chains", "best_params_transduction2_b.rds"))



# SETUP ####

library(fitR)
library(coda)
library(dplyr)
library(reshape2)
library(ggplot2)
library(scales)

source(here::here("Model", "transduction_model_functions.R"))

model = readRDS(here::here("Model", "growth_model.rds"))

lab_dataM = read.csv(here::here("Data", "growth_summary.csv")) %>%
  select(Time, Bacteria, Mean, se) %>%
  filter(Bacteria != "Total")


lab_data = read.csv(here::here("Data", "growth_summary.csv")) %>%
  select(Time, Bacteria, Mean) %>%
  dcast(Time~Bacteria) %>%
  select(-Total) %>%
  rename(time = Time, Bet = DRP, Be = EryR, Bt = TetR) %>%
  mutate(Be = round(Be),
         Bt = round(Bt),
         Bet = round(Bet))


# FIT PHAGE #####


init.theta = c(mu_e = 1, mu_t = 1, mu_et = 1, Nmax = 1e9)
mcmc_fit = run_mcmc(model, lab_data,
                    init.theta = init.theta,
                    proposal.sd = c(init.theta[1]/100,
                                    init.theta[2]/100,
                                    init.theta[3]/100,
                                    init.theta[4]/100),
                    n.iterations = 100000,
                    adapt.size.start = NULL,
                    adapt.shape.start = NULL,
                    adapt.size.cooling = 0.999,
                    growth = TRUE)


trace = coda::mcmc(mcmc_fit$trace)
trace = burnAndThin(trace, 2000, 2)
plot(trace)


params = apply(trace, 2, function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
write.csv(params, here::here("Fitting", "growth_params.csv"), row.names = F)

#replicate 4
init.state = c(Be = lab_data$Be[1], Bt = lab_data$Bt[1], Bet = lab_data$Bet[1])

traj = multi_run2(model, trace, init.state, nruns = 100)


ggplot() +
  geom_line(data = traj, aes(time, Be, colour = "Be", linetype = "Model")) +
  geom_ribbon(data = traj, aes(x = time, ymin = pmax(Be - 1.96*Be_sd, 0), ymax = Be + 1.96*Be_sd,
                               fill = "Be"), alpha = 0.1) +
  geom_line(data = traj, aes(time, Bt, colour = "Bt", linetype = "Model")) +
  geom_ribbon(data = traj, aes(x = time, ymin = pmax(Bt - 1.96*Bt_sd, 0), ymax = Bt + 1.96*Bt_sd,
                               fill = "Bt"), alpha = 0.1) +
  geom_line(data = traj, aes(time, Bet, colour = "Bet", linetype = "Model")) +
  geom_ribbon(data = traj, aes(x = time, ymin = pmax(Bet - 1.96*Bet_sd, 0), ymax = Bet + 1.96*Bet_sd,
                               fill = "Bet"), alpha = 0.1) +
  geom_line(data = lab_dataM %>% filter(Bacteria == "EryR"), 
            aes(Time, Mean, linetype = "Data", colour = "Be"), size = 0.8) +
  geom_errorbar(data = lab_dataM %>% filter(Bacteria == "EryR"), 
                aes(x = Time, ymin = pmax(Mean - se, 0), ymax = Mean + se,
                    colour= "Be"), size = 0.8) +
  geom_line(data = lab_dataM %>% filter(Bacteria == "TetR"), 
            aes(Time, Mean, linetype = "Data", colour = "Bt"), size = 0.8) +
  geom_errorbar(data = lab_dataM %>% filter(Bacteria == "TetR"), 
                aes(x = Time, ymin = pmax(Mean - se, 0), ymax = Mean + se,
                    colour= "Bt"), size = 0.8) +
  geom_line(data = lab_dataM %>% filter(Bacteria == "DRP"), 
            aes(Time, Mean, linetype = "Data", colour = "Bet"), size = 0.8) +
  geom_errorbar(data = lab_dataM %>% filter(Bacteria == "DRP"), 
                aes(x = Time, ymin = pmax(Mean - se, 0), ymax = Mean + se,
                    colour= "Bet"), size = 0.8) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(breaks = seq(0,24,4)) +
  coord_cartesian(ylim = c(1e3, 1e10)) +
  scale_linetype_manual(breaks = c("Data", "Model"),
                        labels = c("Data", "Model"),
                        values = c(1, 2)) +
  scale_colour_manual(breaks = c("Be", "Bt", "Bet"),
                      labels = c(bquote(B[E]), bquote(B[T]), bquote(B[ET])),
                      values = c("#685cc4","#6db356","#c2484d")) +
  scale_fill_manual(breaks = c("Be", "Bt", "Bet"),
                    labels = c(bquote(B[E]), bquote(B[T]), bquote(B[ET])),
                    values = c("#685cc4","#6db356","#c2484d")) +
  labs(x = "Time (hours)", y = "cfu per mL", linetype = "Source:",
       colour = "Bacteria:", fill = "Bacteria:") +
  theme_bw() +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.x = element_text(size=12))


ggsave(here::here("Figures", "supp_fig1.png"))




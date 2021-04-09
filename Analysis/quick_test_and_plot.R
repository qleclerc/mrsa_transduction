
library(fitR)
library(dplyr)
library(reshape2)
library(ggplot2)
library(scales)
library(openxlsx)

source(here::here("Model", "transduction_model_functions.R"))

# mass_model = readRDS(here::here("Fitting", "mass_model.rds"))
model = readRDS(here::here("Model", "transduction_model.rds"))
fitted_params = readRDS(here::here("Fitting", "10_4", "best_params_transduction2.rds"))
lab_data_trans = read.csv(here::here("Lab", "Transduction", "summary_10_4.csv")) %>%
  select(Time, Bacteria, Mean) %>%
  dcast(Time~Bacteria) %>%
  select(-Total) %>%
  rename(time = Time, Bet = DRP, Be = EryR, Bt = TetR) %>%
  mutate(Be = round(Be),
         Bt = round(Bt),
         Bet = round(Bet),
         P = round(P))

model = choose_model(model,
                     frequentist = T,
                     delay = T,
                     fixed_delay = NA,
                     decay = T,
                     link_beta = T,
                     link_L = F,
                     link_delay = F,
                     transduction = T)

trace_model = fitted_params[["frequentist_decay_link_beta"]]
params = trace_model[which.max(trace_model[,"log.density"]),]
params[["gamma"]] = 1/1e-2
init.state = c(Be = 1e6, Bt = 1e6, Bet = 0,
               Pl = 4e5, Pe = 0, Pt = 0)

traj = model$simulate(params, init.state, seq(0, 24, 1))


ggplot() +
  geom_line(data = traj, aes(time, Be, colour = "Be")) +
  geom_line(data = traj, aes(time, Bt, colour = "Bt")) +
  geom_line(data = traj, aes(time, Bet, colour = "Bet")) +
  geom_line(data = traj, aes(time, Pl, colour = "Pl")) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x, n = 6),
                     labels=trans_format("log10", math_format(10^.x))) +
  coord_cartesian(ylim = c(0.1, 1e12))
  
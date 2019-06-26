#A quick script to compare the 3 model outputs with equal parameter values

library(ggplot2)
library(scales)
library(deSolve)

#get the full model function:
source("Models/transduction_full.R")

#get the intermediate model function:
source("Models/transduction_intermediate.R")

#get the mass-action model function:
source("Models/transduction_mass-action.R")

#define parameters:
#time (hours):
times = seq(0,24,0.2)

#parameter values:
parms = c(Ge = 1.9,         #growth parameter for ery resistant
          Gt = 1.9,       #growth parameter for tet resistant
          Get = 1.9,      #growth parameter for ery+tet double resistant
          lysis = 0.08,    #lysis rate
          tre = 1e-7,      #ery gene transduction probability
          trt = 1e-7,      #tet gene transduction probability
          decay = 0.7,    #phage decay rate
          L = 20,         #burst size
          Nmax = 1e9)     #maximum bacteria population size


#initial values:
yinit_full = c(Be = 10000,    #resistant to ery
               Bt = 10000,    #resistant to tet
               Bet = 0,    #resistant to ery and tet
               p = 3300000,         #normal phage
               pe = 0,        #transducing phage (ery)
               pt = 0)        #transducing phage (tet)

yinit_inter = c(Be = 10000,   #resistant to ery
                Bt = 10000,   #resistant to tet
                Bet = 0,   #resistant to ery and tet
                p = 3300000)        #normal phage

yinit_mass = c(Be = 50000,    #resistant to ery
               Bt = 50000,    #resistant to tet
               Bet = 0)    #resistant to ery and tet

#run models:
res_full = as.data.frame(deSolve::ode(func = full_model, times = times, y = yinit_full, parms = parms))
res_inter = as.data.frame(deSolve::ode(func = inter_model, times = times, y = yinit_inter, parms = parms))
res_mass = as.data.frame(deSolve::ode(func = mass_model, times = times, y = yinit_mass, parms = parms))


#plot outputs together:
ggplot() +
  geom_line(data = res_full, aes(time, Bet, colour = "1",linetype = "2"), size = 1) +
  geom_line(data = res_inter, aes(time, Bet, colour = "2",linetype = "2"), size = 1) +
  geom_line(data = res_mass, aes(time, Bet, colour = "3",linetype = "2"), size = 1) +
  geom_line(data = res_full, aes(time, (Bet + Be + Bt), colour = "1",linetype = "1"), size = 1) +
  geom_line(data = res_inter, aes(time, (Bet + Be + Bt), colour = "2",linetype = "1"), size = 1) +
  geom_line(data = res_mass, aes(time, (Bet + Be + Bt), colour = "3",linetype = "1"), size = 1) +
  scale_x_continuous(limits = c(0,24), breaks = seq(0,24,2)) +
  # scale_y_continuous(breaks = seq(0,1e9,2e8)) +
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  xlab("Time (hours)") +
  ylab("Number of bacteria") +
  scale_color_manual(name = "Model:", breaks = c("1","2","3"), values = c("royalblue", "red3", "green3"), labels = c("Full", "Intermediate", "Mass-action")) +
  scale_linetype_manual(name = "Bacteria:", breaks = c("1","2"), values = c("solid","dotted"), labels = c("Total", "DRP"))


#save image:
#ggsave("compare_models.png")

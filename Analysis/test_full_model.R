#A quick script to produce figure 7 from my upgrade report comparing the 3 model outputs

library(ggplot2)
library(scales)
library(deSolve)

#get the full model function:
source("Models/transduction_full.R")

#define parameters:
#time (hours):
times = seq(0,24,0.2)

#parameter values:
parms = c(Ge = 1.9,         #growth parameter for ery resistant
          Gt = 1.9,       #growth parameter for tet resistant
          Get = 1.9,      #growth parameter for ery+tet double resistant
          lysis = 0.08,    #lysis rate
          tre = 0,      #ery gene transduction probability
          trt = 0,      #tet gene transduction probability
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

#run models:
res_full = as.data.frame(deSolve::ode(func = full_model, times = times, y = yinit_full, parms = parms))

ggplot(res_full) +
  geom_line(aes(time, (Bet + Be + Bt), colour = "1"), size = 1) +
  geom_line(aes(time, p, colour = "2"), size = 1) +
  scale_x_continuous(limits = c(0, 24), breaks = seq(0, 24, 2)) +
  scale_y_continuous(trans = log10_trans(),
                     breaks =  trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  xlab("Time (hours)") +
  ylab("Counts") +
  scale_color_manual(name = "", breaks = c("1", "2"), values = c("royalblue", "red3"), labels = c("Bacteria", "Phage"))


#save plot:
#ggsave("test_model.png")

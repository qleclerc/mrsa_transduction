#A quick script to test the output of the full model

library(ggplot2)
library(scales)
library(deSolve)

#get the full model function:
source("Models/transduction_full_alt.R")

res_full = run_full()

#define parameters:
#time (hours):
times = seq(0,24,0.2)

#parameter values:
parms = c(Ge = 1.4,         #growth parameter for ery resistant
          Gt = 1.4,       #growth parameter for tet resistant
          Get = 1.2,      #growth parameter for ery+tet double resistant
          lysis = 10,    #lysis rate
          tre = 0.00001,      #ery gene transduction probability
          trt = 0.00001,      #tet gene transduction probability
          decay = 0.001,    #phage decay rate
          L = 1000,         #burst size
          Nmax = 1*10^9)     #maximum bacteria population size


#initial values:
yinit_full = c(Be = 20*10*10^1,    #resistant to ery
               Bt = 21*10*10^1,    #resistant to tet
               Bet = 0,    #resistant to ery and tet
               p = 50*10^2,         #normal phage
               pe = 0,        #transducing phage (ery)
               pt = 0)        #transducing phage (tet)

#run models:
res_full = run_full(times, yinit_full, parms)

res_full$Bet[which(res_full$Bet < 1)] = NA

#plot just total bacteria and phage:
# ggplot(res_full) +
#   geom_line(aes(time, (Bet + Be + Bt), colour = "1"), size = 1) +
#   geom_line(aes(time, p, colour = "2"), size = 1) +
#   scale_x_continuous(limits = c(0, max(times)), breaks = seq(0, max(times), 2)) +
#   scale_y_continuous(trans = log10_trans(),
#                      breaks =  trans_breaks("log10", function(x) 10^x),
#                      labels = trans_format("log10", math_format(10^.x))) +
#   xlab("Time (hours)") +
#   ylab("Counts") +
#   scale_color_manual(name = "", breaks = c("1", "2"), values = c("royalblue", "red3"),
#                      labels = c("Total bacteria", "Phage"))


#plot everything:
ggplot(res_full) +
  geom_line(aes(time, Bet, colour = "1"), size = 1) +
  geom_line(aes(time, Be, colour = "2"), size = 1) +
  geom_line(aes(time, Bt, colour = "3"), size = 1) +
  geom_line(aes(time, (Bet + Be + Bt), colour = "4"), size = 1) +
  geom_line(aes(time, p, colour = "5"), size = 1) +
  scale_x_continuous(limits = c(0, max(times)), breaks = seq(0, max(times), 2)) +
  scale_y_continuous(trans = log10_trans(),
                     breaks =  trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  xlab("Time (hours)") +
  ylab("Counts") +
  scale_color_manual(name = "", breaks = c("1", "2", "3", "4", "5"), 
                     values = c("royalblue", "green3", "purple3", "black", "red3"),
                     labels = c("DRP", "EryR", "TetR", "Total", "Phage"))

#save plot:
#ggsave("Plots/test_model.png")

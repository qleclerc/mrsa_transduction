

library(deSolve)
library(dplyr)
library(ggplot2)
library(scales)
library(cowplot)


## Load the lab data: ####

#adjust to what your WD is!!
lab_data4 = read.csv(here::here("Lab", "Transduction", "summary_10_4.csv"), stringsAsFactors = F)
lab_data3 = read.csv(here::here("Lab", "Transduction", "summary_10_3.csv"), stringsAsFactors = F)
lab_data5 = read.csv(here::here("Lab", "Transduction", "summary_10_5.csv"), stringsAsFactors = F)

#also you might have some issues on mac due to compatibility - fix column names in the datasets in that case


## MODEL ############################################

#I have a wrapper around the desolve model for simplicity
#parameter "def_model" allows you to choose which model to run
test_model = function(times = seq(0,24,1),
                      params = c(Nmax = 1e9,
                                 mu1 = 1.9,
                                 mu2 = 1.9,
                                 beta = 0.000000001,
                                 L = 10),
                      y = c(B1 = 1e2,
                            B2 = 1e2,
                            P = 1e4),
                      def_model = "mass"){
  
  if(!(def_model %in% c("mass", "mass_link", "mass_link_burst", "mass_link_both", "frequentist"))){
    stop("Function parameter \"def_model\" must be \"mass\", \"mass_link\", \"mass_link_burst\",
         \"mass_link_both\" or \"frequentist\"")
  }
  
  #MASS ACTION
  mass_model = function(times, y, params){
    
    with(as.list(c(y,params)),{
      
      
      N = B1 + B2
      
      
      dB1 = B1 * mu1 * (1 - N/Nmax) - B1 * beta * P   
      dB2 = B2 * mu2 * (1 - N/Nmax) - B2 * beta * P
      dP = N * beta * P * L
      
      
      return(list(c(dB1, dB2,dP)))
      
    }
    )
  }
  
  
  ##MASS ACTION WITH LINK
  mass_link_model = function(times, y, params){
    
    with(as.list(c(y,params)),{
      
      
      N = B1 + B2
      
      beta_t = beta * (1 - N/Nmax)
      
      dB1 = B1 * mu1 * (1 - N/Nmax) - B1 * beta_t * P 
      dB2 = B2 * mu2 * (1 - N/Nmax) - B2 * beta_t * P
      dP = N * beta_t * L * P
      
      
      return(list(c(dB1, dB2,dP)))
      
    }
    )
  }
  
  
  ##MASS ACTION WITH LINK, BURST
  mass_link_burst_model = function(times, y, params){
    
    with(as.list(c(y,params)),{
      
      
      N = B1 + B2
      
      L_t = L * (1 - N/Nmax)
      
      dB1 = B1 * mu1 * (1 - N/Nmax) - B1 * beta * P 
      dB2 = B2 * mu2 * (1 - N/Nmax) - B2 * beta * P
      dP = N * beta * L_t * P
      
      
      return(list(c(dB1, dB2,dP)))
      
    }
    )
  }
  
  
  
  ##MASS ACTION WITH LINK, BOTH
  mass_link_both_model = function(times, y, params){
    
    with(as.list(c(y,params)),{
      
      
      N = B1 + B2
      
      beta_t = beta * (1 - N/Nmax)
      L_t = L * (1 - N/Nmax)
      
      dB1 = B1 * mu1 * (1 - N/Nmax) - B1 * beta_t * P 
      dB2 = B2 * mu2 * (1 - N/Nmax) - B2 * beta_t * P
      dP = N * beta_t * L_t * P
      
      
      return(list(c(dB1, dB2,dP)))
      
    }
    )
  }
  
  
  ##FREQUENTIST MODEL
  frequentist_model = function(times, y, params){
    
    with(as.list(c(y,params)),{
      
      
      N = B1 + B2
      
      beta_t = beta * (1 - N/Nmax)
      L_t = L * (1 - N/Nmax)
      
      P_inf = (1-exp(-beta_t*N)) * P
      B_inf = (1-exp(-P_inf/N)) * N
      
      dB1 = B1 * mu1 * (1 - N/Nmax) - B1/N * B_inf 
      dB2 = B2 * mu2 * (1 - N/Nmax) - B2/N * B_inf 
      dP = B_inf * L_t - P_inf
      
      
      return(list(c(dB1, dB2,dP)))
      
    }
    )
  }
  
  
  
  
  if(def_model == "mass") {
    res_full = as.data.frame(deSolve::ode(func = mass_model, times = times, y = y, parms = params))
  }
  if(def_model == "mass_link") {
    res_full = as.data.frame(deSolve::ode(func = mass_link_model, times = times, y = y, parms = params))
  }
  if(def_model == "mass_link_burst") {
    res_full = as.data.frame(deSolve::ode(func = mass_link_burst_model, times = times, y = y, parms = params))
  }
  if(def_model == "mass_link_both") {
    res_full = as.data.frame(deSolve::ode(func = mass_link_both_model, times = times, y = y, parms = params))
  }
  if(def_model == "frequentist") {
    res_full = as.data.frame(deSolve::ode(func = frequentist_model, times = times, y = y, parms = params))
  }
  
  
  #replace values < 1 by 0, except "times" column ofc
  res_full[,-1][res_full[,-1] < 1] = 0
  
  
  return(res_full)
  
  
}



## PARAMETERS #####################################

#play around with this!

params = c(Nmax = 1e9, #max number of bacteria possible in the system
           mu1 = 1.9, #max growth rate of bacteria 1
           mu2 = 1.8, #max growth rate of bacteria 2
           beta = 0.000000005, #this corresponds to individual phage adsorption rate
           L = 7) #burst size (number of new phage released after infection)


#which model do you want to run?
#mass, mass_link, mass_link_burst, mass_link_both, or frequentist

model = "frequentist"




## RUN ALL MODELS & PLOT ################################

#fixed parameters, no touching!!
#times
times = seq(0,24,0.2)

#initial number of bacteria and phage (matches average from lab data)
yinit = c(B1 = 1e4, 
          B2 = 1e4,
          P = 5e4)

#run model starting phage at 1e4
res4 = test_model(times, params, yinit, def_model = model)

#run model starting phage at 1e3
yinit["P"] = 1.5e3
res3 = test_model(times, params, yinit, def_model = model)

#run model starting phage at 1e5
yinit["P"] = 2e5
res5 = test_model(times, params, yinit, def_model = model)

#plot for Pstart 1e4
pp4 = ggplot(res4) +
  geom_line(aes(time, B1, colour = "EryR")) +
  geom_line(aes(time, B2, colour = "TetR")) +
  geom_line(aes(time, P, colour = "P")) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  geom_line(data = lab_data4 %>% filter(Bacteria != "Total"),
            aes(Time, Mean, colour = Bacteria), linetype = "dashed") +
  theme_bw()

#plot for Pstart 1e3
pp3 = ggplot(res3) +
  geom_line(aes(time, B1, colour = "EryR")) +
  geom_line(aes(time, B2, colour = "TetR")) +
  geom_line(aes(time, P, colour = "P")) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  geom_line(data = lab_data3 %>% filter(Bacteria != "Total"),
            aes(Time, Mean, colour = Bacteria), linetype = "dashed") +
  theme_bw()

#plot for Pstart 1e5
pp5 = ggplot(res5) +
  geom_line(aes(time, B1, colour = "EryR")) +
  geom_line(aes(time, B2, colour = "TetR")) +
  geom_line(aes(time, P, colour = "P")) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  geom_line(data = lab_data5 %>% filter(Bacteria != "Total"),
            aes(Time, Mean, colour = Bacteria), linetype = "dashed") +
  theme_bw()

#plot all together
plot_grid(pp3, pp4, pp5,
          labels = c("10^3", "10^4", "10^5"))
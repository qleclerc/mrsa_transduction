
library(deSolve)
library(reshape2)
library(dplyr)
library(ggplot2)
library(scales)
library(cowplot)



## Load the lab data: ####

#adjust to what your WD is!!
lab_data4 = read.csv(here::here("Lab", "Transduction", "summary_10_4.csv"), stringsAsFactors = F)
lab_data3 = read.csv(here::here("Lab", "Transduction", "summary_10_3.csv"), stringsAsFactors = F)
lab_data5 = read.csv(here::here("Lab", "Transduction", "summary_10_5.csv"), stringsAsFactors = F)


## Model: ####

#https://jvi.asm.org/content/jvi/76/11/5557.full.pdf
#basically consider that sort of link between phage and bacteria

#I have a wrapper around the desolve model for simplicity
#some default values included in function definition, which correspond to best "eyball" fit to lab data
test_model = function(times = seq(0,24,1),
                      params = c(Nmax = 6e8,
                                 mu1 = 1.9,
                                 mu2 = 1.9,
                                 mu3 = 1.9,
                                 beta = 0.0000002,
                                 pr1 = 0.01,
                                 pr2 = 0.01,
                                 L = 5,
                                 decay = 0),
                      y = c(B1 = 1e2,
                            B2 = 1e2,
                            B3 = 0,
                            P = 1e4,
                            P1 = 0,
                            P2 = 0),
                      quick_plot = F,
                      lab_data = NULL, 
                      def_model = T){
  
  #the original model I designed
  full_model = function(times, y, parms){
    
    with(as.list(c(y,parms)),{
      
      #total bacteria and phage populations:
      N = B1 + B2 + B3
      p_tot = P + P1 + P2
      if (p_tot == 0) p_tot = 1 #safety to prevent model breaking
      
      #probability for one bacterium to be infected by one phage:
      pr_inf_tot = beta
      
      #bacteria resistant to ery:
      dB1 = mu1*(1 - N/Nmax)*B1 - pr_inf_tot*p_tot/B1
      #growth - infected by normal phage or transducing phage (tet)
      
      #bacteria resistant to tet:
      dB2 = mu2*(1 - N/Nmax)*B2 - pr_inf_tot*p_tot/B2
      #growth - infected by normal phage or transducing phage (ery)
      
      #bacteria resistant to ery and tet:
      dB3 = mu3*(1 - N/Nmax)*B3 - pr_inf_tot*p_tot/B3 + pr_inf_tot*p_tot*(P2/p_tot)/B1 + pr_inf_tot*p_tot*(P1/p_tot)/B2
      #growth - infected by normal phage + tet transduction event + ery transduction event
      
      
      #normal phage:
      dP = pr_inf_tot*p_tot*(L - 1)*(1 - pr1)/B1 + pr_inf_tot*p_tot*(L - 1)*(1 - pr2)/B2 + pr_inf_tot*p_tot*(L - 1)*(1 - pr1 - pr2)/B3 - decay*P
      #new phage produced by lysis in each bacteria group - decay
      
      #transducing phage (ery):
      dP1 = pr_inf_tot*p_tot*pr1/B1 + pr_inf_tot*p_tot*pr1/B3 - pr_inf_tot*P1/N - decay*P1
      #new phage by transduction in ery resistant + new phage by transduction in ery+tet resistant - infections - decay
      
      #transducing phage (tet):
      dP2 = pr_inf_tot*p_tot*pr2/B2 + pr_inf_tot*p_tot*pr2/B3 - pr_inf_tot*P2/N - decay*P2
      #new phage by transduction in tet resistant + new phage by transduction in ery+tet resistant - infections - decay
      
      
      return(list(c(dB1, dB2, dB3, dP, dP1, dP2)))
      
    }
    )
  }
  
  #more complicated, converts to probas to estimate n phage binding and n bacteria bound
  run_model = function(times, y, params){
    
    with(as.list(c(y,params)),{
      
      #a little "cheating"; lab data shows bacteria lags for 1h, but growth functions will systematically fail to reproduce this...
      #so instead, I'm forcing the bacteria to not grow for 1h
      #if(times < 1) mu1 = mu2 = mu3 = 0
      
      N = B1 + B2 + B3
      
      #Number of phage who will infect a bacteria at time t
      P_inf = (1 - exp(-beta*N)) * P
      P_inf1 = (1 - exp(-beta*N)) * P1
      P_inf2 = (1 - exp(-beta*N)) * P2
      
      #Proportion of bacteria who will be infected at time t
      lambda = 1 - exp(-P_inf/N)
      lambda1 = 1 - exp(-P_inf1/N)
      lambda2 = 1 - exp(-P_inf2/N)
      
      
      dB1 = B1 * mu1 * (1 - N/Nmax) - B1 * lambda / N * B1 - B1 * lambda2 / N * B1 
      dB2 = B2 * mu2 * (1 - N/Nmax) - B2 * lambda / N * B2 - B2 * lambda1 / N * B2
      dB3 = B3 * mu3 * (1 - N/Nmax) - B3 * lambda / N * B3 +
        B2 * lambda1 / N * B2 + B1 * lambda2 / N * B1 
      dP = N * lambda * L - P_inf - decay * P
      
      dP1 = (B1+B3) * lambda * pr1 / N * (B1+B3) - P_inf1
      dP2 = (B2+B3) * lambda * pr2 / N * (B2+B3) - P_inf2
      
      
      return(list(c(dB1, dB2, dB3, dP, dP1, dP2)))
      
    }
    )
  }
  
  
  ##MASS ACTION WITH LINK
  mass_link_model = function(times, y, params){
    
    with(as.list(c(y,params)),{
      
      
      N = B1 + B2
      
      #Proportion of bacteria who will be infected at time t
      lambda = beta * (1 - N/Nmax)

      dB1 = B1 * mu1 * (1 - N/Nmax) - B1 * lambda * P 
      dB2 = B2 * mu2 * (1 - N/Nmax) - B2 * lambda * P
      dB3 = 0 
      dP = N * lambda * (L-1) * P
      
      dP1 = 0
      dP2 = 0
      
      
      return(list(c(dB1, dB2, dB3, dP, dP1, dP2)))
      
    }
    )
  }
  
  
  
  #MASS ACTION
  mass_model = function(times, y, params){
    
    with(as.list(c(y,params)),{
      
      
      N = B1 + B2
      
    
      dB1 = B1 * mu1 * (1 - N/Nmax) - B1 * beta * P   
      dB2 = B2 * mu2 * (1 - N/Nmax) - B2 * beta * P
      dB3 = 0
      dP = N * beta * P * (L-1)
      
      dP1 = 0
      dP2 = 0
      
      
      return(list(c(dB1, dB2, dB3, dP, dP1, dP2)))
      
    }
    )
  }
  
  link_model = function(times, y, params){
    
    with(as.list(c(y,params)),{
      

      N = B1 + B2 + B3

      #Proportion of bacteria who will be infected at time t
      P_inf = beta * (1 - N/Nmax)

      
      #Proportion of bacteria who will be infected at time t
      lambda = P_inf*P #1 - exp(-P_inf)# * P/N)
      lambda1 = 1 - exp(-P_inf * P1/N)
      lambda2 = 1 - exp(-P_inf * P2/N)
      
      print(paste0("P_inf: ", P_inf, ", lambda: ", lambda))
      
      dB1 = B1 * mu1 * (1 - N/Nmax) - B1 * lambda  - B1 * lambda2  
      dB2 = B2 * mu2 * (1 - N/Nmax) - B2 * lambda  - B2 * lambda1 
      dB3 = B3 * mu3 * (1 - N/Nmax) - B3 * lambda  +
        B2 * lambda1 + B1 * lambda2 
      dP = N * lambda * L - P_inf * P - decay * P
      
      dP1 = (B1+B3) * lambda * pr1 - N * lambda1 * P1
      dP2 = (B2+B3) * lambda * pr2 - N * lambda2 * P2
      
      
      return(list(c(dB1, dB2, dB3, dP, dP1, dP2)))
      
    }
    )
  }
  
  if(def_model == "alt" ) {
    res_full = as.data.frame(deSolve::ode(func = run_model, times = times, y = y, parms = params))
  }
  if(def_model == "full") {
    res_full = as.data.frame(deSolve::ode(func = full_model, times = times, y = y, parms = params))
  }
  if(def_model == "link") {
    res_full = as.data.frame(deSolve::ode(func = link_model, times = times, y = y, parms = params))
  }
  if(def_model == "mass") {
    res_full = as.data.frame(deSolve::ode(func = mass_model, times = times, y = y, parms = params))
  }
  if(def_model == "mass_link") {
    res_full = as.data.frame(deSolve::ode(func = mass_link_model, times = times, y = y, parms = params))
  }
  
  #replace values < 1 by 0, except "times" column ofc
  res_full[,-1][res_full[,-1] < 1] = 0
  
  if (quick_plot == T) {
    
    #Quick plotting if you're playing around with the model:
    pp = ggplot(res_full) +
      geom_line(aes(time, B1, colour = "EryR")) +
      geom_line(aes(time, B2, colour = "TetR")) +
      geom_line(aes(time, B3, colour = "DRP")) +
      geom_line(aes(time, P, colour = "P")) +
      scale_y_continuous(trans=log10_trans(),
                         breaks=trans_breaks("log10", function(x) 10^x),
                         labels=trans_format("log10", math_format(10^.x)))
    
    
    if(!(is.null(lab_data))){
      
      lab_data = lab_data %>% 
        filter(Bacteria != "Total")
      
      pp = pp +
        geom_line(data = lab_data, aes(Time, Mean, colour = Bacteria), linetype = "dashed")
      
      
    }
    
    plot(pp)
    
  }
  
  return(res_full)
  
  
}

times = seq(0,24,0.2)

params = c(Nmax = 1e9, #max number of bacteria possible in the system
           mu1 = 1.9, #max growth rate of bacteria
           mu2 = 1.8,
           mu3 = 1.5,
           beta = 0.000000001, #this corresponds to individual phage adsorption rate
           pr1 = 0,
           pr2 = 0,
           L = 10, #burst size (number of new phage released after infection)
           decay = 0) #phage decay rate in the system

#initial number of bacteria and phage (matches average from lab data)
yinit = c(B1 = 1e4, 
          B2 = 1e4,
          B3 = 0,
          P = 5e4,
          P1 = 0,
          P2 = 0)

#quick_plot option generates a plot of the results, useful when playing around with parameter values (set FALSE otherwise)
#lab_data option adds mean lab data points to plot, useful when attempting to eyeball parameter values (leave NULL otherwise)
res4 = test_model(times, params, yinit, quick_plot = F, def_model = "mass_link")
yinit["P"] = 1.5e3
res3 = test_model(times, params, yinit, quick_plot = F, def_model = "mass_link")
yinit["P"] = 2e5
res5 = test_model(times, params, yinit, quick_plot = F, def_model = "mass_link")

pp4 = ggplot(res4) +
  geom_line(aes(time, B1, colour = "EryR")) +
  geom_line(aes(time, B2, colour = "TetR")) +
  geom_line(aes(time, B3, colour = "DRP")) +
  geom_line(aes(time, P, colour = "P")) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  geom_line(data = lab_data4 %>% filter(Bacteria != "Total"),
            aes(Time, Mean, colour = Bacteria), linetype = "dashed") +
  theme_bw()

pp3 = ggplot(res3) +
  geom_line(aes(time, B1, colour = "EryR")) +
  geom_line(aes(time, B2, colour = "TetR")) +
  geom_line(aes(time, B3, colour = "DRP")) +
  geom_line(aes(time, P, colour = "P")) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  geom_line(data = lab_data3 %>% filter(Bacteria != "Total"),
            aes(Time, Mean, colour = Bacteria), linetype = "dashed") +
  theme_bw()

pp5 = ggplot(res5) +
  geom_line(aes(time, B1, colour = "EryR")) +
  geom_line(aes(time, B2, colour = "TetR")) +
  geom_line(aes(time, B3, colour = "DRP")) +
  geom_line(aes(time, P, colour = "P")) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  geom_line(data = lab_data5 %>% filter(Bacteria != "Total"),
            aes(Time, Mean, colour = Bacteria), linetype = "dashed") +
  theme_bw()


plot_grid(pp3, pp4, pp5,
          labels = c("10^3", "10^4", "10^5"))
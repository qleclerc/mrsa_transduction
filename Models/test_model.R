
library(deSolve)
library(reshape2)
library(dplyr)
library(ggplot2)
library(scales)



## Load the lab data: ####

#adjust to what your WD is!!
lab_data = read.csv(here::here("Lab", "Transduction", "04_09_20", "data1.csv"), stringsAsFactors = F)


## Model: ####

#https://jvi.asm.org/content/jvi/76/11/5557.full.pdf
#basically consider that sort of link between phage and bacteria

#I have a wrapper around the desolve model for simplicity
#some default values included in function definition, which correspond to best "eyball" fit to lab data
test_model = function(times = seq(0,24,1),
                      params = c(Nmax = 6e8,
                                 mu = 1.9,
                                 beta = 0.0000002, 
                                 L = 5,
                                 decay = 0),
                      y = c(B = 1e2,
                            P = 1e4),
                      quick_plot = F,
                      lab_data = NULL,
                      mass = F){
  
  #complete MOI model
  run_model = function(times, y, params){
    
    with(as.list(c(y,params)),{
      
      #a little "cheating"; lab data shows bacteria lags for 1h, but growth functions will systematically fail to reproduce this...
      #so instead, I'm forcing the bacteria to not grow for 1h
      if(times < 1) mu = 0
      
      #Number of phage who will infect a bacteria at time t
      P_inf = (1 - exp(-beta*B)) * P
      
      #Proportion of bacteria who will be infected at time t
      lambda = 1 - exp(-P_inf/B)

      dB = B * mu * (1 - B/Nmax) - B * lambda
      dP = B * lambda * L - P_inf - decay * P
      

      return(list(c(dB, dP)))
      
    }
    )
  }
  
  #simple mass action model
  run_model_mass = function(times, y, params){
    
    with(as.list(c(y,params)),{
      
      #a little "cheating"; lab data shows bacteria lags for 1h, but growth functions will systematically fail to reproduce this...
      #so instead, I'm forcing the bacteria to not grow for 1h
      if(times < 1) mu = 0
      
      #Proportion of bacteria who will be infected at time t
      lambda = beta * P
      
      dB = B * mu * (1 - B/Nmax) - B * lambda
      dP = B * lambda * L - decay * P
      
      
      return(list(c(dB, dP)))
      
    }
    )
  }
  
  if(mass == T){
    res_full = as.data.frame(deSolve::ode(func = run_model_mass, times = times, y = y, parms = params))
  }
  else{
    res_full = as.data.frame(deSolve::ode(func = run_model, times = times, y = y, parms = params))
  }

  
  if (quick_plot == T) {
    
    #Quick plotting if you're playing around with the model:
    plot(res_full$time, log10(res_full$B), col="black", ylim = c(2, 12), type = "l", 
         xlab = "Time (hours)", ylab = "Bacteria & phage (log10(cfu/mL))", xaxt="n")
    axis(1, at = seq(0,24,2))
    lines(res_full$time, log10(res_full$P), col="red")
    legend("bottomright",lty=1, legend=c("Bacteria", "Phage"), col=c("black", "red"))
    
    #Also plot lab data alongside, if passed to function
    if (!is.null(lab_data)) {
      
      points(unique(lab_data$Time),
             lab_data %>% filter(Bacteria == "Total") %>% select("Mean") %>% log10() %>% unlist(.),
             col = "black",
             pch = 19)
      points(unique(lab_data$Time),
             lab_data %>% filter(Bacteria == "P") %>% select("Mean") %>% log10(.) %>% unlist(.),
             col = "red",
             pch = 19)
      
    }
  }
  
  return(res_full)
  
  
}

times = seq(0,24,1)

params = c(Nmax = 6e8, #max number of bacteria possible in the system
           mu = 1.9, #max growth rate of bacteria
           beta = 0.00000001, #this corresponds to individual phage adsorption rate
           L = 2, #burst size (number of new phage released after infection)
           decay = 0) #phage decay rate in the system

#initial number of bacteria and phage (matches average from lab data)
yinit = c(B = 4.4e3, 
          P = 9e3)

#quick_plot option generates a plot of the results, useful when playing around with parameter values (set FALSE otherwise)
#lab_data option adds mean lab data points to plot, useful when attempting to eyeball parameter values (leave NULL otherwise)
res = test_model(times, params, yinit, quick_plot = T, lab_data = lab_data,T)


lab_data_s = data.frame(lab_data %>% filter(Bacteria == "Total") %>% select("Mean"),
                        lab_data %>% filter(Bacteria == "P") %>% select("Mean"))
colnames(lab_data_s) = c("Bacteria_d", "Phage_d") 


#some optim play
opt_fun = function(x){
  params = c(Nmax = 6e8, #max number of bacteria possible in the system
             mu = 1.9, #max growth rate of bacteria
             beta = x[1], #this corresponds to individual phage adsorption rate
             L = x[2], #burst size (number of new phage released after infection)
             decay = 0) #phage decay rate in the system

  res = test_model(seq(0,8,1), params, yinit, F, NULL, F)
  colnames(res) = c("Time", "Bacteria", "Phage")
  
  res = cbind(res, lab_data_s)
  
  rss = res %>%
    transmute(RSS_bac = (Bacteria-Bacteria_d)^2, RSS_phage = (Phage-Phage_d)^2) %>%
    sum(.)
    
  return(rss)
}

optim_res = optim(c(0.0000002, 5), opt_fun, lower=c(0.0000000001, 1), upper=c(2,Inf), method = "L-BFGS-B")

params = c(Nmax = 6e8, #max number of bacteria possible in the system
           mu = 1.9, #max growth rate of bacteria
           beta = optim_res$par[1], #this corresponds to individual phage adsorption rate
           L = optim_res$par[2], #burst size (number of new phage released after infection)
           decay = 0) #phage decay rate in the system

res = test_model(times, params, yinit, quick_plot = T, lab_data = lab_data)

#reshape model output for simplicity in ggplot plotting below
res = reshape2::melt(res, id.vars = "time")

#This line will resample the model results assuming a normal distibution with SD = 0.25*mean
#This value for SD is roughly what I'm seeing in the lab
res$value = sapply(res$value, FUN = function(x) rnorm(1, x, x*0.5))

#Rename columns
colnames(res) = c("Time", "Bacteria", "Mean")

#generate summary res
res = test_model(times, params, yinit, quick_plot = T, lab_data = lab_data)
res = reshape2::melt(res, id.vars = "time")
res$value = sapply(res$value, FUN = function(x) rnorm(1, x, x*0.5))
colnames(res) = c("Time", "Bacteria", "Mean")
res2 = test_model(times, params, yinit, quick_plot = T, lab_data = lab_data)
res2 = reshape2::melt(res2, id.vars = "time")
res2$value = sapply(res2$value, FUN = function(x) rnorm(1, x, x*0.5))
colnames(res2) = c("Time", "Bacteria", "Mean")
res3 = test_model(times, params, yinit, quick_plot = T, lab_data = lab_data)
res3 = reshape2::melt(res3, id.vars = "time")
res3$value = sapply(res3$value, FUN = function(x) rnorm(1, x, x*0.5))
colnames(res3) = c("Time", "Bacteria", "Mean")

summary_res = res
colnames(summary_res)[3] = "data1"
summary_res$data2 = res2$Mean
summary_res$data3 = res3$Mean
summary_res$Mean = rowMeans(summary_res[,c(3:5)])
write.csv(summary_res, "summary_res.csv", row.names = F)

#for res mass
res_mass = test_model(times, params, yinit, quick_plot = T, lab_data = lab_data, T)
res_mass = reshape2::melt(res_mass, id.vars = "time")
res_mass$value = sapply(res_mass$value, FUN = function(x) rnorm(1, x, x*0.5))
colnames(res_mass) = c("Time", "Bacteria", "Mean")
res_mass2 = test_model(times, params, yinit, quick_plot = T, lab_data = lab_data, T)
res_mass2 = reshape2::melt(res_mass2, id.vars = "time")
res_mass2$value = sapply(res_mass2$value, FUN = function(x) rnorm(1, x, x*0.5))
colnames(res_mass2) = c("Time", "Bacteria", "Mean")
res_mass3 = test_model(times, params, yinit, quick_plot = T, lab_data = lab_data,T)
res_mass3 = reshape2::melt(res_mass3, id.vars = "time")
res_mass3$value = sapply(res_mass3$value, FUN = function(x) rnorm(1, x, x*0.5))
colnames(res_mass3) = c("Time", "Bacteria", "Mean")

summary_res_mass = res_mass
colnames(summary_res_mass)[3] = "data1"
summary_res_mass$data2 = res_mass2$Mean
summary_res_mass$data3 = res_mass3$Mean
summary_res_mass$Mean = rowMeans(summary_res_mass[,c(3:5)])
write.csv(summary_res_mass, "summary_res_mass.csv", row.names = F)


## Plotting: ####


#Plot lab data and model output together:
ggplot() +
  geom_line(data = res %>% filter(Bacteria == "B"), aes(Time, Mean, colour = "Bacteria", linetype="full")) +
  geom_line(data = res %>% filter(Bacteria == "P"), aes(Time, Mean, colour = "Phage", linetype="full")) +
  geom_line(data = res_mass %>% filter(Bacteria == "B"), aes(Time, Mean, colour = "Bacteria", linetype="mass")) +
  geom_line(data = res_mass %>% filter(Bacteria == "P"), aes(Time, Mean, colour = "Phage", linetype="mass")) +
  geom_point(data = lab_data %>% filter(Bacteria == "Total"), aes(Time, Mean, colour = "Bacteria"), size=3, shape=18) +
  geom_point(data = lab_data %>% filter(Bacteria == "P"), aes(Time, Mean, colour = "Phage"), size=3, shape=18) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  scale_color_manual(name="",breaks=c("Bacteria", "Phage"),values=c("red4", "royalblue"),labels=c("Bacteria", "Phage"))+
  scale_linetype_manual(name="", breaks=c("full", "mass"), values=c("solid", "dashed"),labels=c("MOI interaction", "Mass action interaction"))+
  ylab("log10(units/mL)")+
  xlab("Time (hours)")

ggsave("longterm_models.png")


#Load a dataset containing 3 replicates from the model (with resampling for variability)
summary_res = read.csv("summary_res.csv")
summary_res = summary_res %>% filter(Time < 9)
summary_res$SE = apply(summary_res[,c(3:5)], 1, FUN = function(x) sd(x)/sqrt(3))

summary_res_mass = read.csv("summary_res_mass.csv")
summary_res_mass = summary_res_mass %>% filter(Time < 9)
summary_res_mass$SE = apply(summary_res_mass[,c(3:5)], 1, FUN = function(x) sd(x)/sqrt(3))


#Plot 3 replicates of the model and the average of the model data:
ggplot() +
  # geom_point(data = summary_res %>% filter(Bacteria == "B"), aes(Time, data1, colour = "Bacteria")) +
  # geom_point(data = summary_res %>% filter(Bacteria == "P"), aes(Time, data1, colour = "Phage")) +
  # geom_point(data = summary_res %>% filter(Bacteria == "B"), aes(Time, data2, colour = "Bacteria")) +
  # geom_point(data = summary_res %>% filter(Bacteria == "P"), aes(Time, data2, colour = "Phage")) +
  # geom_point(data = summary_res %>% filter(Bacteria == "B"), aes(Time, data3, colour = "Bacteria")) +
  # geom_point(data = summary_res %>% filter(Bacteria == "P"), aes(Time, data3, colour = "Phage")) +
  geom_line(data = summary_res %>% filter(Bacteria == "B"), aes(Time, Mean, colour = "Bacteria")) +
  geom_line(data = summary_res %>% filter(Bacteria == "P"), aes(Time, Mean, colour = "Phage")) +
  # geom_point(data = summary_res_mass %>% filter(Bacteria == "B"), aes(Time, data1, colour = "Bacteria")) +
  # geom_point(data = summary_res_mass %>% filter(Bacteria == "P"), aes(Time, data1, colour = "Phage")) +
  # geom_point(data = summary_res_mass %>% filter(Bacteria == "B"), aes(Time, data2, colour = "Bacteria")) +
  # geom_point(data = summary_res_mass %>% filter(Bacteria == "P"), aes(Time, data2, colour = "Phage")) +
  # geom_point(data = summary_res_mass %>% filter(Bacteria == "B"), aes(Time, data3, colour = "Bacteria")) +
  # geom_point(data = summary_res_mass %>% filter(Bacteria == "P"), aes(Time, data3, colour = "Phage")) +
  geom_line(data = summary_res_mass %>% filter(Bacteria == "B"), aes(Time, Mean, colour = "Bacteria")) +
  geom_line(data = summary_res_mass %>% filter(Bacteria == "P"), aes(Time, Mean, colour = "Phage")) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x)))+
  scale_color_manual(name="",breaks=c("Bacteria", "Phage"),values=c("red4", "royalblue"),labels=c("Bacteria", "Phage"))+
  ylab("cfu/mL")



#Plot the average of 3 model replicates, and the average of 3 lab replicates
ggplot() +
  geom_line(data = summary_res %>% filter(Bacteria == "B"), aes(Time, Mean, colour = "Bacteria", linetype="full")) +
  geom_line(data = summary_res %>% filter(Bacteria == "P"), aes(Time, Mean, colour = "Phage", linetype="full")) +
  geom_errorbar(data = summary_res %>% filter(Bacteria == "B"), aes(Time, ymin=Mean-SE, ymax=Mean+SE, colour = "Bacteria")) +
  geom_errorbar(data = summary_res %>% filter(Bacteria == "P"), aes(Time, ymin=Mean-SE, ymax=Mean+SE, colour = "Phage")) +
  geom_line(data = summary_res_mass %>% filter(Bacteria == "B"), aes(Time, Mean, colour = "Bacteria", linetype="mass")) +
  geom_line(data = summary_res_mass %>% filter(Bacteria == "P"), aes(Time, Mean, colour = "Phage", linetype="mass")) +
  geom_errorbar(data = summary_res_mass %>% filter(Bacteria == "B"), aes(Time, ymin=Mean-SE, ymax=Mean+SE, colour = "Bacteria")) +
  geom_errorbar(data = summary_res_mass %>% filter(Bacteria == "P"), aes(Time, ymin=Mean-SE, ymax=Mean+SE, colour = "Phage")) +
  geom_point(data = lab_data %>% filter(Bacteria == "Total"), aes(Time, Mean, colour = "Bacteria"), size=3, shape=18) +
  geom_point(data = lab_data %>% filter(Bacteria == "P"), aes(Time, Mean, colour = "Phage"), size=3, shape=18) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  scale_color_manual(name="",breaks=c("Bacteria", "Phage"),values=c("red4", "royalblue"),labels=c("Bacteria", "Phage"))+
  scale_linetype_manual(name="", breaks=c("full", "mass"), values=c("solid", "dashed"),labels=c("MOI interaction", "Mass action interaction"))+
  ylab("log10(units/mL)") +
  xlab("Time (hours)")


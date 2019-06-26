
#This is the mass-action model, with no inclusion of phage at all, and the assumption that transduction occurs based on the density of bacteria


library(deSolve)


mass_model = function(times, y, parms){
  
  with(as.list(c(y,parms)),{
    
    #total bacteria population:
    N = Be+Bt+Bet

    #bacteria resistant to ery:
    dBe = Ge*(1-N/Nmax)*Be - trt*Be*(Bt+Bet)/N
    #growth - transduced by tet
    
    #bacteria resistant to tet:
    dBt = Gt*(1-N/Nmax)*Bt - tre*Bt*(Be+Bet)/N
    #growth - transduced by ery
    
    #bacteria resistant to ery and tet:
    dBet = Get*(1-N/Nmax)*Bet + trt*Be*(Bt+Bet)/N + tre*Bt*(Be+Bet)/N
    #growth + ery transduction event + tet transduction event
    
    
    return(list(c(dBe, dBt, dBet)))
    
  }
  )
}



#time (hours):
times = seq(1,25,0.2)

#parameter values:
#(dummy values that give a somewhat realistic result)
parms=c(Ge=2,         #growth parameter for ery resistant
        Gt=2,         #growth parameter for tet resistant
        Get=2,        #growth parameter for ery+tet double resistant
        tre=0.1,      #ery gene transduction rate
        trt=0.1,      #tet gene transduction rate
        Nmax=1e9)     #maximum bacteria population size


#initial values
yinit = c(Be=8000,    #resistant to ery
          Bt=7000,    #resistant to tet
          Bet=0)      #resistant to ery and tet


#event trigger to add erythromycin at later time point:
#(simple method is to remove all susceptible bacteria in model i.e. Bt)
#time = -1: turned off
eventdat = data.frame(var="Bt",
                       time=-1,
                       value=0,
                       method="rep")

#event trigger to add tetracycline at later time point:
#(simple method is to remove all susceptible bacteria in model i.e. Be)
#time = -1: turned off
eventdat2 = data.frame(var="Be",
                       time=-1,
                       value=0,
                       method="rep")



#run and plot:
results = as.data.frame(deSolve::ode(func=mass_model, times=times, y=yinit, parms=parms, events=list(data=rbind(eventdat,eventdat2))))
results$time = results$time-1   #adjust model time (in model, t=1 is the start, so effectively 0h)

plot(results$time, (results$Be+results$Bt+results$Bet), col="black", type="l", 
     ylab="Number of bacteria", xlab="Time (h)", lwd=3, 
     ylim=c(0, max((results$Be+results$Bt+results$Bet))), xaxt="n")
axis(1, at=seq(0,max(results$time),4))
lines(results$time, results$Be, col="green", lwd=3)
lines(results$time, results$Bt, col="red", lwd=3)
lines(results$time, results$Bet, col="blue", lwd=3)
legend("topleft", cex=0.8, legend=c("Total bacteria","Ery mono-resistant","Tet mono-resistant","Double resistant"), col=c("black","green","red","blue"), lty=1, lwd=3)



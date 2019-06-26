
#This is the intermediate model, with just lytic phage, and the assumption that transduction occurs based on the density of phage and bacteria


library(deSolve)


intermediate_model = function(times, y, parms){
  
  with(as.list(c(y,parms)),{
    
    #total bacteria population:
    N = Be+Bt+Bet

    #probability for one bacterium to be infected by one phage:
    pr_inf = 1 - exp(-(p/N)*lysis)

    #bacteria resistant to ery:
    dBe = Ge*(1-N/Nmax)*Be - pr_inf*Be
    #growth - infected by phage
    
    #bacteria resistant to tet:
    dBt = Gt*(1-N/Nmax)*Bt - pr_inf*Bt
    #growth - infected by phage
    
    #bacteria resistant to ery and tet:
    dBet = Get*(1-N/Nmax)*Bet - pr_inf*Bet + pr_inf*Be*trt*(Bt+Bet)/N + pr_inf*Bt*tre*(Be+Bet)/N
    #growth - infected by phage + tet transduction event + ery transduction event
    
    
    #normal phage:
    dp = pr_inf*Be*(1-tre)*(L-1) + pr_inf*Bt*(1-trt)*(L-1) + pr_inf*Bet*(1-tre-trt)*(L-1) - decay*p
    #new phage produced by lysis in each bacteria group - decay
    
    
    return(list(c(dBe, dBt, dBet, dp)))
    
  }
  )
}



#time (hours):
times = seq(1,25,0.2)

#parameter values:
#(dummy values that give a (somewhat?) realistic result)
parms=c(Ge=2,         #growth parameter for ery resistant
        Gt=2,         #growth parameter for tet resistant
        Get=2,        #growth parameter for ery+tet double resistant
        lysis=0.2,    #lysis rate
        tre=0.1,      #ery gene transduction probability
        trt=0.1,      #tet gene transduction probability
        decay=0.8,    #phage decay rate
        L=10,         #burst size
        Nmax=1e9)     #maximum bacteria population size


#initial values
yinit = c(Be=8000,    #resistant to ery
          Bt=7000,    #resistant to tet
          Bet=0,      #resistant to ery and tet
          p=1000)     #normal phage


#event trigger to add phage at later time point:
#time = -1: turned off
eventdat = data.frame(var="p",
                      time=-1,
                      value=1000,
                      method="add")

#event trigger to remove transduction at later time point:
#(analogous to removing all phage in model, sodium citrate effectively blocks all of them)
#time = -1: turned off
eventdat2 = data.frame(var="p",
                       time=-1,
                       value=0,
                       method="rep")

#event trigger to add erythromycin at later time point:
#(simple method is to remove all susceptible bacteria in model i.e. Bt)
#time = -1: turned off
eventdat3 = data.frame(var="Bt",
                       time=-1,
                       value=0,
                       method="rep")

#event trigger to add tetracycline at later time point:
#(simple method is to remove all susceptible bacteria in model i.e. Be)
#time = -1: turned off
eventdat4 = data.frame(var="Be",
                       time=-1,
                       value=0,
                       method="rep")


#run and plot:
results = as.data.frame(deSolve::ode(func=intermediate_model, times=times, y=yinit, parms=parms, events=list(data=rbind(eventdat,eventdat2,eventdat3,eventdat4))))
results$time = results$time-1   #adjust model time (in model, t=1 is the start, so effectively 0h)

par(mfrow=c(1,2), mar=c(5.1,4.1,5.5,2.1), xpd=T)

plot(results$time, (results$Be+results$Bt+results$Bet), col="black", type="l", 
     ylab="Number of bacteria", xlab="Time (h)", lwd=3, 
     ylim=c(0, max((results$Be+results$Bt+results$Bet))), xaxt="n")
axis(1, at=seq(0,max(results$time),4))
lines(results$time, results$Be, col="green", lwd=3)
lines(results$time, results$Bt, col="red", lwd=3)
lines(results$time, results$Bet, col="blue", lwd=3)
legend("topright", cex=0.8, inset=c(0,-0.19), legend=c("Total bacteria","Ery mono-resistant","Tet mono-resistant","Double resistant"), col=c("black","green","red","blue"), lty=1, lwd=3)

plot(results$time, results$p, col="black", type="l", ylim=c(0, max(results$p)), lwd=3, 
     ylab="Number of phages", xlab="Time (h)", xaxt="n")
axis(1, at=seq(0,24,4))

par(mfrow=c(1,1))


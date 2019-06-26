#this is inspired by previous work from payne&jansen (for simple bacteria-phage dynamics) and volkova (for transduction probabilities) models:
#https://www.sciencedirect.com/science/article/pii/S0022519300921982?via%3Dihub
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4068684/ 

#This is the full model, with normal and transducing phage


library(deSolve)


full_model = function(times, y, parms){
  
  with(as.list(c(y,parms)),{
    
    #total bacteria and phage populations:
    N = Be+Bt+Bet
    p_tot = p+pe+pt
    if(p_tot == 0) p_tot = 1 #safety to prevent model breaking
    
    #probability for one bacterium to be infected by one phage:
    pr_inf_tot = 1 - exp(-(p_tot/N)*lysis)
    if(p_tot == 0) p_inf_tot=0 #safety to prevent model breaking
    
    #bacteria resistant to ery:
    dBe = Ge*(1-N/Nmax)*Be - pr_inf_tot*Be*((p+pt)/p_tot)
    #growth - infected by normal phage or transducing phage (tet)
    
    #bacteria resistant to tet:
    dBt = Gt*(1-N/Nmax)*Bt - pr_inf_tot*Bt*((p+pe)/p_tot)
    #growth - infected by normal phage or transducing phage (ery)
    
    #bacteria resistant to ery and tet:
    dBet = Get*(1-N/Nmax)*Bet - pr_inf_tot*Bet*(p/p_tot) + pr_inf_tot*Be*(pt/p_tot) + pr_inf_tot*Bt*(pe/p_tot)
    #growth - infected by normal phage + tet transduction event + ery transduction event
    
    
    #normal phage:
    dp = pr_inf_tot*Be*(p/p_tot)*(L-1)*(1-tre) + pr_inf_tot*Bt*(p/p_tot)*(L-1)*(1-trt) + pr_inf_tot*Bet*(p/p_tot)*(L-1)*(1-tre-trt) - decay*p
    #new phage produced by lysis in each bacteria group - decay
    
    #transducing phage (ery):
    dpe = pr_inf_tot*Be*(p/p_tot)*tre + pr_inf_tot*Bet*(p/p_tot)*tre - pr_inf_tot*N*(pe/p_tot) - decay*pe
    #new phage by transduction in ery resistant + new phage by transduction in ery+tet resistant - infections - decay
    
    #transducing phage (tet):
    dpt = pr_inf_tot*Bt*(p/p_tot)*trt + pr_inf_tot*Bet*(p/p_tot)*trt - pr_inf_tot*N*(pt/p_tot) - decay*pt
    #new phage by transduction in tet resistant + new phage by transduction in ery+tet resistant - infections - decay
    
    
    return(list(c(dBe, dBt, dBet, dp, dpe, dpt)))
    
  }
  )
}



#time (hours):
times = seq(1,25,0.2)

#parameter values:
#(dummy values that give a (somewhat?) realistic result)
parms=c(Ge=2,         #growth parameter for ery resistant
        Gt=1.9,         #growth parameter for tet resistant
        Get=1.9,        #growth parameter for ery+tet double resistant
        lysis=0.2,    #lysis rate
        tre=0.1,      #ery gene transduction probability
        trt=0.1,      #tet gene transduction probability
        decay=0.8,    #phage decay rate
        L=10,         #burst size
        Nmax=1e9)     #maximum bacteria population size


#initial values:
yinit = c(Be=11000,    #resistant to ery
          Bt=16000,    #resistant to tet
          Bet=9100,      #resistant to ery and tet
          p=0,     #normal phage
          pe=0,       #transducing phage (ery)
          pt=0)       #transducing phage (tet)


#event trigger to add phage at later time point:
#time = -1: turned off
eventdat = data.frame(var="p",
                      time=-1,
                      value=1000,
                      method="add")

#event trigger to remove transduction at later time point:
#(analogous to removing all phage in model, sodium citrate effectively blocks all of them)
#time = -1: turned off
eventdat2 = data.frame(var=c("p","pe","pt"),
                       time=rep(-1,3),
                       value=rep(0,3),
                       method=rep("rep",3))

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
results = as.data.frame(deSolve::ode(func=full_model, times=times, y=yinit, parms=parms, events=list(data=rbind(eventdat,eventdat2,eventdat3,eventdat4))))
results$time = results$time-1   #adjust model time (in model, t=1 is the start, so effectively 0h)


#par(mfrow=c(1,2), mar=c(5.1,4.1,5.5,2.1), xpd=T)

plot(results$time, (results$Be+results$Bt+results$Bet),col="black", type="l", 
     ylab="Number of bacteria", xlab="Time (h)", lwd=3, 
     ylim=c(0, max((results$Be+results$Bt+results$Bet))), xaxt="n")
axis(1, at=seq(0,max(results$time),4))
lines(results$time, results$Be, col="green", lwd=3)
lines(results$time, results$Bt, col="red", lwd=3)
lines(results$time, results$Bet, col="blue", lwd=3)
legend("topright", cex=0.8, inset=c(0,-0.19), legend=c("Total bacteria","Ery mono-resistant","Tet mono-resistant","Double resistant"), col=c("black","green","red","blue"), lty=1, lwd=3)

# plot(results$time, results$p, col="black", type="l", ylim=c(0, max(results$p)), lwd=3, 
#      ylab="Number of phages", xlab="Time (h)", xaxt="n")
# axis(1, at=seq(0,24,4))
# lines(results$time, results$pe, col="green", lwd=3)
# lines(results$time, results$pt, col="red", lwd=3)
# legend("topright", cex=0.8, inset=c(0,-0.15), legend=c("Normal phage","Transducing phage (ery)","Transducing phage (tet)"), col=c("black","green","red"), lty=1, lwd=3)
# 
# par(mfrow=c(1,1))


# #code for a prettier figure:
# 
# library(ggplot2)
# library(scales)
# 
# # results$time = results$time-1
# # results = rbind(results, tail(results)[6,])
# # results$time[117] = 24
# ggplot(results)+
#   geom_line(aes(x=time, y=Be, colour="EryR"), size=1)+
#   geom_line(aes(x=time, y=Bet, colour="DRP"), size=1)+
#   geom_line(aes(x=time, y=Bt, colour="TetR"), size=1)+
#   geom_line(aes(x=time, y=(Be+Bt+Bet), colour="Total"), size=1)+
#   scale_x_continuous(limits=c(0,24), breaks=seq(0,24,2))+
#   scale_y_continuous(limits=c(0, 1e9),
#                      breaks=seq(0, 1e9, 2e8))+
#   labs(color="Bacteria:")+
#   ylab("cfu/ml")+
#   xlab("Time (hours)")
# 
# ggsave("expected_lab.png")
# 
# 
# 
# ggplot(results)+
#   geom_line(aes(x=time, y=p, colour="Lytic phage"), size=1)+
#   scale_x_continuous(limits=c(0,24), breaks=seq(0,24,2))+
#   scale_y_continuous(limits=c(0, 5.5e9),
#                      breaks=seq(0, 5.5e9, 5e8))+
#   labs(color="")+
#   ylab("pfu/ml")+
#   xlab("Time (hours)")
# 
# ggsave("expected_lab_phage.png")

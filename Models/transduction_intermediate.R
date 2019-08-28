#This is the intermediate model, with just lytic phage, and the assumption that transduction occurs based on the density of phage and bacteria

run_inter = function(times = seq(0,24,0.2),
                     yinit = c(Be = 11000, 
                               Bt = 9000, 
                               Bet = 0,  
                               p = 3300000,    
                               pe = 0,   
                               pt = 0),
                     parms = c(Ge = 1.5,       
                               Gt = 1.6,     
                               Get = 1.3,      
                               lysis = 0.08,   
                               tre = 0.00001,      
                               trt = 0.00001,     
                               decay = 0.7,    
                               L = 20,         
                               Nmax = 1e9)){
  
  inter_model = function(times, y, parms){
    
    with(as.list(c(y,parms)),{
      
      #total bacteria population:
      N = Be + Bt + Bet
      
      #probability for one bacterium to be infected by one phage:
      pr_inf = 1 - exp(-(p/N)*lysis)
      
      #bacteria resistant to ery:
      dBe = Ge*(1 - N/Nmax)*Be - pr_inf*Be
      #growth - infected by phage
      
      #bacteria resistant to tet:
      dBt = Gt*(1 - N/Nmax)*Bt - pr_inf*Bt
      #growth - infected by phage
      
      #bacteria resistant to ery and tet:
      dBet = Get*(1 - N/Nmax)*Bet - pr_inf*Bet + pr_inf*Be*trt*(Bt + Bet)/N + pr_inf*Bt*tre*(Be + Bet)/N
      #growth - infected by phage + tet transduction event + ery transduction event
      
      
      #normal phage:
      dp = pr_inf*Be*(1 - tre)*(L - 1) + pr_inf*Bt*(1 - trt)*(L - 1) + pr_inf*Bet*(1 - tre - trt)*(L - 1) - decay*p
      #new phage produced by lysis in each bacteria group - decay
      
      
      return(list(c(dBe, dBt, dBet, dp)))
      
    }
    )
  }
  
  res_inter = as.data.frame(deSolve::ode(func = inter_model, times = times, y = yinit, parms = parms))
  
  return(res_inter)
  
}

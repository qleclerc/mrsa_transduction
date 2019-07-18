#This is the mass-action model, with no inclusion of phage at all, and the assumption that transduction occurs based on the density of bacteria

run_mass = function(times, yinit, parms){
  
  mass_model = function(times, y, parms){
    
    with(as.list(c(y,parms)),{
      
      #total bacteria population:
      N = Be + Bt + Bet
      
      #bacteria resistant to ery:
      dBe = Ge*(1 - N/Nmax)*Be - trt*Be*(Bt + Bet)/N
      #growth - transduced by tet
      
      #bacteria resistant to tet:
      dBt = Gt*(1 - N/Nmax)*Bt - tre*Bt*(Be + Bet)/N
      #growth - transduced by ery
      
      #bacteria resistant to ery and tet:
      dBet = Get*(1 - N/Nmax)*Bet + trt*Be*(Bt + Bet)/N + tre*Bt*(Be + Bet)/N
      #growth + ery transduction event + tet transduction event
      
      
      return(list(c(dBe, dBt, dBet)))
      
    }
    )
  }
  
  res_mass = as.data.frame(deSolve::ode(func = mass_model, times = times, y = yinit, parms = parms))
  
  return(res_mass)
  
}
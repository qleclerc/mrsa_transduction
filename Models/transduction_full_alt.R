#The logic here is that the rate of infection must be tied to the ratio of bacteria and phage
#not sure yet how best to implement that, but seems on track
#must see an initial lag of 4h until phage start multiplying, could enforce hard threshold of bacteria, but more likely that this is some sort of threshold ratio

#currently problem is that phage numbers stay too low as the threshold catches up, could fix by tweaking formula for bacteria growth, or some how calculating threshold differently such that could end up with more phage

run_full = function(times = seq(0,24,0.2),
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
  
  library(deSolve)
  
  full_model = function(times, y, parms){
    
    with(as.list(c(y,parms)),{
      
      #total bacteria and phage populations:
      N = Be + Bt + Bet
      p_tot = p + pe + pt
      if (p_tot == 0) p_tot = 1 #safety to prevent model breaking
      
      #probability for one bacterium to be infected by one phage:
      pr_inf_tot = lysis

      #bacteria resistant to ery:
      dBe = Ge*(1 - N/Nmax)*Be - pr_inf_tot*Be/p_tot
      #growth - infected by normal phage or transducing phage (tet)
      
      #bacteria resistant to tet:
      dBt = Gt*(1 - N/Nmax)*Bt - pr_inf_tot*Bt/p_tot
      #growth - infected by normal phage or transducing phage (ery)
      
      #bacteria resistant to ery and tet:
      dBet = Get*(1 - N/Nmax)*Bet - pr_inf_tot*Bet/p_tot + pr_inf_tot*Be/p_tot*(pt/p_tot) + pr_inf_tot*Bt/p_tot*(pe/p_tot)
      #growth - infected by normal phage + tet transduction event + ery transduction event
      
      
      #normal phage:
      dp = pr_inf_tot*Be/p_tot*(L - 1)*(1 - tre) + pr_inf_tot*Bt/p_tot*(L - 1)*(1 - trt) + pr_inf_tot*Bet/p_tot*(L - 1)*(1 - tre - trt) - decay*p
      #new phage produced by lysis in each bacteria group - decay
      
      #transducing phage (ery):
      dpe = pr_inf_tot*Be/p_tot*tre + pr_inf_tot*Bet/p_tot*tre - pr_inf_tot*N*pe - decay*pe
      #new phage by transduction in ery resistant + new phage by transduction in ery+tet resistant - infections - decay
      
      #transducing phage (tet):
      dpt = pr_inf_tot*Bt/p_tot*trt + pr_inf_tot*Bet/p_tot*trt - pr_inf_tot*N*pt - decay*pt
      #new phage by transduction in tet resistant + new phage by transduction in ery+tet resistant - infections - decay
      
      
      return(list(c(dBe, dBt, dBet, dp, dpe, dpt)))
      
    }
    )
  }
  
  res_full = as.data.frame(deSolve::ode(func = full_model, times = times, y = yinit, parms = parms))
  
  return(res_full)
  
}
#this is inspired by previous work from payne&jansen (for simple bacteria-phage dynamics) and volkova (for transduction probabilities) models:
#https://www.sciencedirect.com/science/article/pii/S0022519300921982?via%3Dihub
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4068684/ 

#This is the full model, with normal and transducing phage


library(deSolve)

#redefine as wrapper around the model, so when source creates a function that runs the model

full_model = function(times, y, parms){
  
  with(as.list(c(y,parms)),{
    
    #total bacteria and phage populations:
    N = Be + Bt + Bet
    p_tot = p + pe + pt
    if (p_tot == 0) p_tot = 1 #safety to prevent model breaking
    
    #probability for one bacterium to be infected by one phage:
    pr_inf_tot = 1 - exp(-(p_tot/N)*lysis)
    if (p_tot == 0) pr_inf_tot = 0 #safety to prevent model breaking
    
    #bacteria resistant to ery:
    dBe = Ge*(1 - N/Nmax)*Be - pr_inf_tot*Be*((p + pt)/p_tot)
    #growth - infected by normal phage or transducing phage (tet)
    
    #bacteria resistant to tet:
    dBt = Gt*(1 - N/Nmax)*Bt - pr_inf_tot*Bt*((p + pe)/p_tot)
    #growth - infected by normal phage or transducing phage (ery)
    
    #bacteria resistant to ery and tet:
    dBet = Get*(1 - N/Nmax)*Bet - pr_inf_tot*Bet*(p/p_tot) + pr_inf_tot*Be*(pt/p_tot) + pr_inf_tot*Bt*(pe/p_tot)
    #growth - infected by normal phage + tet transduction event + ery transduction event
    
    
    #normal phage:
    dp = pr_inf_tot*Be*(p/p_tot)*(L - 1)*(1 - tre) + pr_inf_tot*Bt*(p/p_tot)*(L - 1)*(1 - trt) + pr_inf_tot*Bet*(p/p_tot)*(L - 1)*(1 - tre - trt) - decay*p
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
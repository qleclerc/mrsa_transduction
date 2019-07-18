
likelihood_func = function(data, observation){
  
  #value to store sum of log likelihood
  dens = 0
  
  #loop through data length
  for (i in length(data)) {
    
    pdens = dpois(data[i], observation[i], log = T)
    dens = dens + pdens
    
  }
  
  return(pdens)
  
}

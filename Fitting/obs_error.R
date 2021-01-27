
#resamples data with a given observation SD

obs_error = function(data, error_sd=1){
  
  #resample each point with corresponding SD
  #data_error = rnorm(data, data, error_sd)
  
  data_error = rpois(1, data)
  
  #this trick overcomes NA error if integer is past .Machine$integer.max
  if (is.na(data_error)) data_error = round(rnorm(1,mean=data,sd=sqrt(data)))
  
  return(data_error)
  
}

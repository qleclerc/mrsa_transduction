
#resamples data with a given observation SD

obs_error = function(data, error_sd){
  
  #resample each point with corresponding SD
  data_error = rnorm(data, data, error_sd)
  
  return(data_error)
  
}

#generates a replicate of the data, assuming Poisson distribution

generate_replicate = function(data){
  
  data[,2:7] = apply(res_full[,2:7], c(1,2), function(x) obs_error(x))
  
  data
}

res_full1 = generate_replicate(res_full)
res_full2 = generate_replicate(res_full)
res_full3 = generate_replicate(res_full)

plot_output(res_full1)
plot_output(res_full2)
plot_output(res_full3)

res_full1
rowMeans(cbind(res_full1$Be,res_full2$Be,res_full3$Be))

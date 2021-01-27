
plot_multi_output = function(results1, results2, results3, log_t = TRUE, title = NULL){
  
  mean_results = results1
  
  #combine data in means
  for (column in colnames(mean_results)[-1]) {
    
    colvals = rowMeans(cbind(res_full1[,column],res_full2[,column],res_full3[,column]))
    mean_results[,column] = colvals
  }
  
  
  #plot data
  
  p = ggplot(mean_results) +
    geom_line(aes(time, rowSums(cbind(mean_results$Bet, mean_results$Be, mean_results$Bt)), colour = "1"), size = 1) +
    scale_x_continuous(limits = c(0, 24), breaks = seq(0, 24, 2)) +
    xlab("Time (hours)") +
    ylab("Counts") +
    ggtitle(title) +
    scale_color_manual(name = "", breaks = c("1", "2", "3", "4", "5"), 
                       values = c("black", "royalblue", "green3", "purple3", "red3"),
                       labels = c("Total", "DRP", "EryR", "TetR", "Phage"))
  
  
  if (!is.null(mean_results$Bet)) p = p + 
    geom_line(aes(time, mean_results$Bet, colour = "2"), size = 1) +
    geom_point(data = results1, aes(time, Bet, colour = "2")) +
    geom_point(data = results2, aes(time, Bet, colour = "2")) +
    geom_point(data = results3, aes(time, Bet, colour = "2"))

  if (!is.null(mean_results$Be)) p = p + 
    geom_line(aes(time, mean_results$Be, colour = "3"), size = 1) +
    geom_point(data = results1, aes(time, Be, colour = "3")) +
    geom_point(data = results2, aes(time, Be, colour = "3")) +
    geom_point(data = results3, aes(time, Be, colour = "3"))
  
  if (!is.null(mean_results$Bt)) p = p + 
    geom_line(aes(time, mean_results$Bt, colour = "4"), size = 1) +
    geom_point(data = results1, aes(time, Bt, colour = "4")) +
    geom_point(data = results2, aes(time, Bt, colour = "4")) +
    geom_point(data = results3, aes(time, Bt, colour = "4"))

  if (!is.null(mean_results$p)) p = p + 
    geom_line(aes(x = time, y = mean_results$p, colour = "5"), size = 1) +
    geom_point(data = results1, aes(time, p, colour = "5")) +
    geom_point(data = results2, aes(time, p, colour = "5")) +
    geom_point(data = results3, aes(time, p, colour = "5"))
  
  
  if (log_t == T) p = p + scale_y_continuous(trans = log10_trans(),
                                             breaks =  trans_breaks("log10", function(x) 10^x, n = 6),
                                             labels = trans_format("log10", math_format(10^.x)))
  
  if (error_bar == T) {
   all_data = cbind(results1)
     p = p + geom_errorbar(aes(x=Time, ymin=1, ymax=2)) 
  
  }
  print(p)
  
  return(mean_results)
  
  
}

plot_multi_output(res_full1, res_full2, res_full3)

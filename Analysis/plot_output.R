
plot_output = function(results, log_t = TRUE, title = NULL){
  
  p = ggplot(results) +
    geom_line(aes(time, rowSums(cbind(results$Bet, results$Be, results$Bt)), colour = "1"), size = 1) +
    scale_x_continuous(limits = c(0, 24), breaks = seq(0, 24, 2)) +
    xlab("Time (hours)") +
    ylab("Counts") +
    ggtitle(title) +
    scale_color_manual(name = "", breaks = c("1", "2", "3", "4", "5"), 
                       values = c("black", "royalblue", "green3", "purple3", "red3"),
                       labels = c("Total", "DRP", "EryR", "TetR", "Phage"))
  
  
  if (!is.null(results$Bet)) p = p + geom_line(aes(time, results$Bet, colour = "2"), size = 1)
  if (!is.null(results$Be)) p = p + geom_line(aes(time, results$Be, colour = "3"), size = 1)
  if (!is.null(results$Bt)) p = p + geom_line(aes(time, results$Bt, colour = "4"), size = 1)
  if (!is.null(results$p)) p = p + geom_line(aes(x = time, y = results$p, colour = "5"), size = 1)
  
  if (log_t == T) p = p + scale_y_continuous(trans = log10_trans(),
                                             breaks =  trans_breaks("log10", function(x) 10^x, n = 6),
                                             labels = trans_format("log10", math_format(10^.x)))
  
  
  plot(p)
  
  
}

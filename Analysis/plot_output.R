
plot_output = function(results, log_t = TRUE, title = NULL){
  
  if (log_t == TRUE) {
    ggplot(results) +
      geom_line(aes(time, results$Bet, colour = "1"), size = 1) +
      geom_line(aes(time, results$Be, colour = "2"), size = 1) +
      geom_line(aes(time, results$Bt, colour = "3"), size = 1) +
      geom_line(aes(time, rowSums(cbind(results$Bet, results$Be, results$Bt)), colour = "4"), size = 1) +
      geom_line(aes(time, results$p, colour = "5"), size = 1) +
      scale_x_continuous(limits = c(0, 24), breaks = seq(0, 24, 2)) +
      scale_y_continuous(trans = log10_trans(),
                         breaks =  trans_breaks("log10", function(x) 10^x, n = 6),
                         labels = trans_format("log10", math_format(10^.x))) +
      xlab("Time (hours)") +
      ylab("Counts") +
      ggtitle(title) +
      scale_color_manual(name = "", breaks = c("1", "2", "3", "4", "5"), 
                         values = c("royalblue", "green3", "purple3", "black", "red3"),
                         labels = c("DRP", "EryR", "TetR", "Total", "Phage"))
  } else {
    
    ggplot(results) +
      geom_line(aes(time, results$Bet, colour = "1"), size = 1) +
      geom_line(aes(time, results$Be, colour = "2"), size = 1) +
      geom_line(aes(time, results$Bt, colour = "3"), size = 1) +
      geom_line(aes(time, rowSums(cbind(results$Bet, results$Be, results$Bt)), colour = "4"), size = 1) +
      geom_line(aes(time, results$p, colour = "5"), size = 1) +
      scale_x_continuous(limits = c(0, 24), breaks = seq(0, 24, 2)) +
      xlab("Time (hours)") +
      ylab("Counts") +
      ggtitle(title) +
      scale_color_manual(name = "", breaks = c("1", "2", "3", "4", "5"), 
                         values = c("royalblue", "green3", "purple3", "black", "red3"),
                         labels = c("DRP", "EryR", "TetR", "Total", "Phage"))
    
  }
  
}

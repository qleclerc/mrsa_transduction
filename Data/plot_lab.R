
plot_lab = function(data, errorbar = T, all_points = F, save_plot = F){
  
  filename = "lab_plot"
  
  gg = ggplot(data, aes(x=Time, y=Mean, colour=Bacteria))+
    geom_point()+
    geom_line(linetype="dashed")+
    scale_x_continuous(limits=c(0,max(data$Time)), breaks=seq(0,max(data$Time),2))+
    scale_y_continuous(trans=log10_trans(),
                       breaks=trans_breaks("log10", function(x) 10^x),
                       labels=trans_format("log10", math_format(10^.x)))+
    ylab("cfu/ml")+
    xlab("Time (hours)")
  
  if(all_points == T){
    gg = gg +
      geom_point(aes(Time, data1, colour=Bacteria))+
      geom_point(aes(Time, data2, colour=Bacteria))+
      geom_point(aes(Time, data3, colour=Bacteria))
    filename = paste0(filename, "_all")
  }
  
  if(errorbar == T){
    gg = gg + geom_errorbar(aes(x=Time, ymin=se_min, ymax=se_max))
    filename = paste0(filename, "_error")
  }
  
  print(gg)
  
  if(save_plot == T) ggsave(paste0(filename,".png"))
  
}

plot_lab(data, F, T, T)

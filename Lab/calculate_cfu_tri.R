library(ggplot2)
library(dplyr)
library(scales)


calculate_cfu_tri = function(data1, data2=NULL, data3=NULL, logged=TRUE){
  
  extra_label = ""

  count = 0
  for(i in 1:3){
    
    if(!is.null(get(paste0("data",i)))) count = count+1
    
  }
  
  timepoints = unique(data1$Time)
  
  all_data = data.frame(Time=rep(timepoints, 4))
  all_data$Time = all_data$Time[order(all_data$Time)]
  all_data$Bacteria = c("Total", "EryR", "TetR", "DRP")
  
  
  for(i in 1:count){
    
    d = get(paste0("data",i))
    
    d$cfu = d$Count*(10^(d$Dilution))*20
    
    d = d %>%
      group_by(Time, Plate) %>%
      summarise(cfu = mean(cfu)) %>%
      ungroup
    
    for(times in unique(d$Time)){
      
      good_rows = which(d$Time == times)
      d$cfu[good_rows[2]] = d$cfu[good_rows[2]] - d$cfu[good_rows[4]]
      d$cfu[good_rows[3]] = d$cfu[good_rows[3]] - d$cfu[good_rows[4]]
      
    }
    
    d$Bacteria = rep(c("Total", "EryR", "TetR", "DRP"), nrow(d)/4)
    
    #clean up:
    d = d[,-2]
    
    
    #plot:
    if(logged==T){
      
      gg = ggplot(d, aes(Time, cfu, colour=Bacteria))+
        geom_point()+
        geom_line()+
        scale_x_continuous(limits=c(0,24), breaks=seq(0,24,2))+
        scale_y_continuous(trans=log10_trans(),
                           breaks=trans_breaks("log10", function(x) 10^x),
                           labels=trans_format("log10", math_format(10^.x)))+
        ylab("cfu/ml")+
        xlab("Time (hours)")+
        theme_bw()
      
      print(gg)
      ggsave(paste0(extra_label,"data",i,"_log.png"), plot=gg)
      
    } else {
      
      gg = ggplot(d, aes(Time, cfu, colour=Bacteria))+
        geom_point()+
        geom_line()+
        scale_x_continuous(limits=c(0,24), breaks=seq(0,24,2))+
        ylab("cfu/ml")+
        xlab("Time (hours)")+
        theme_bw()
      
      print(gg)
      ggsave(paste0(extra_label,"data",i,".png"), plot=gg)
      
    }
    
    write.csv(d, paste0(extra_label,"data",i,".csv"), row.names = F)
    
    all_data = cbind(all_data, d$cfu)
    colnames(all_data)[2+i] = paste0("data",i)
    
  }
  
  #summary if more than 1 dataset
  if(count>1){
    
    all_data$Mean = rowMeans(all_data[,3:dim(all_data)[2]])
    all_data$se = apply(all_data[,3:(dim(all_data)[2]-1)], 1,
                        FUN = function(x){sd(x)/sqrt(count)})
    
    all_data$se_min = all_data$Mean-all_data$se
    #correction to avoid error when logging:
    all_data$se_min[which(all_data$se_min < 0)] = 1
    
    all_data$se_max = all_data$Mean+all_data$se
    
    if(logged==T){
      
      gg = ggplot(all_data, aes(x=Time, y=Mean, colour=Bacteria))+
        geom_point()+
        geom_line()+
        geom_errorbar(aes(x=Time, ymin=se_min, ymax=se_max))+
        scale_x_continuous(limits=c(0,24), breaks=seq(0,24,2))+
        scale_y_continuous(trans=log10_trans(),
                           breaks=trans_breaks("log10", function(x) 10^x),
                           labels=trans_format("log10", math_format(10^.x)))+
        ylab("cfu/ml")+
        xlab("Time (hours)")+
        theme_bw()
      
      print(gg)
      ggsave(paste0(extra_label,"summary_log.png"), plot=gg)
      
    } else {
      
      gg = ggplot(all_data, aes(x=Time, y=Mean, colour=Bacteria))+
        geom_point()+
        geom_line()+
        geom_errorbar(aes(x=Time, ymin=se_min, ymax=se_max))+
        scale_x_continuous(limits=c(0,24), breaks=seq(0,24,2))+
        ylab("cfu/ml")+
        xlab("Time (hours)")+
        theme_bw()
      
      print(gg)
      ggsave(paste0(extra_label,"summary.png"), plot=gg)
      
    }
    
    
    write.csv(all_data, paste0(extra_label,"summary.csv"), row.names = F)
    
  }
}

data1 = read.csv('Lab/Triculture/050619/050619_1.csv')[-81,]
data2 = read.csv('Lab/Triculture/050619/050619_2.csv')[-81,]
data3 = read.csv('Lab/Triculture/050619/050619_3.csv')[-81,]

calculate_cfu_tri(data1, data2, data3, logged=T)

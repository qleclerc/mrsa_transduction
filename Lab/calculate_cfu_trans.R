setwd("../R/Lab/Transduction/200619")
library(ggplot2)
library(dplyr)
library(scales)


calculate_cfu_trans = function(data1, data2=NULL, data3=NULL, data4=NULL, data5=NULL, data6=NUL, logged=TRUE){
  
  
  count = 0
  for(i in 1:5){
    
    if(!is.null(get(paste0("data",i)))) count = count+1
    
  }
  
  timepoints = unique(data1$Time)
  
  all_data = data.frame(Time=rep(timepoints, 5))
  all_data$Time = all_data$Time[order(all_data$Time)]
  all_data$Bacteria = c("Total", "EryR", "P", "TetR", "DRP")
  
  
  for(i in 1:count){
    
    d = get(paste0("data",i))
    
    d$cfu = d$Count*(10^(d$Dilution))*20
    
    d$cfu[which(d$Plate == "P")] = d$cfu[which(d$Plate == "P")]/2
    
    d = d %>%
      group_by(Time, Plate) %>%
      summarise(cfu = mean(cfu))
    
    #loop through each time and deduct number of DRP from Ery and Tet to work out number of 327 and 201kt7
    for(times in unique(d$Time)){
      
      #order is: BHIA, E, P, T, T+E
      good_rows = which(d$Time == times)
      d$cfu[good_rows[2]] = d$cfu[good_rows[2]] - d$cfu[good_rows[5]]
      d$cfu[good_rows[4]] = d$cfu[good_rows[4]] - d$cfu[good_rows[5]]
      
    }
    
    d$Bacteria = c("Total", "EryR", "P", "TetR", "DRP")
    
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
        xlab("Time (hours)")
      
      print(gg)
      ggsave(paste0("data",i,"_log.png"), plot=gg)
      
    } else {
      
      gg = ggplot(d, aes(Time, cfu, colour=Bacteria))+
        geom_point()+
        geom_line()+
        scale_x_continuous(limits=c(0,24), breaks=seq(0,24,2))+
        ylab("cfu/ml")+
        xlab("Time (hours)")
      
      print(gg)
      ggsave(paste0("data",i,".png"), plot=gg)
      
    }
    
    write.csv(d, paste0("data",i,".csv"), row.names = F)
    
    all_data = cbind(all_data, d$cfu)
    colnames(all_data)[2+i] = paste0("data",i)
    
  }
  
  #summary if more than 1 dataset
  if(count>1){
    
    all_data$Mean = rowMeans(all_data[,3:dim(all_data)[2]])
    all_data$se = apply(all_data[,3:(dim(all_data)[2]-1)], 1,
                        FUN = function(x){1.96*(sd(x)/sqrt(dim(all_data)[2]-2))})
    
    all_data$se_min = all_data$Mean-all_data$se
    #correction to avoid error when logging:
    all_data$se_min[which(all_data$se_min < 1)] = 1
    
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
        xlab("Time (hours)")
      
      print(gg)
      ggsave("summary_log.png", plot=gg)
      
    } else {
      
      gg = ggplot(all_data, aes(x=Time, y=Mean, colour=Bacteria))+
        geom_point()+
        geom_line()+
        geom_errorbar(aes(x=Time, ymin=se_min, ymax=se_max))+
        scale_x_continuous(limits=c(0,24), breaks=seq(0,24,2))+
        ylab("cfu/ml")+
        xlab("Time (hours)")
      
      print(gg)
      ggsave("summary.png", plot=gg)
      
    }
    
    
    write.csv(all_data, "summary.csv", row.names = F)
    
  }
}

data1 = read.csv('200619.csv')[-101,-(5:8)]
data2 = read.csv('050619_2.csv')[-101,]
data3 = read.csv('050619_3.csv')[-101,]

calculate_cfu_trans(data1, logged=T)

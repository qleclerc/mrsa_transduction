library(ggplot2)
library(dplyr)
library(scales)
library(reshape2)


calculate_cfu_trans = function(data1, data2=NULL, data3=NULL, logged=TRUE, volume = 30){
  
  #check how many datasets were inputed
  count = 0
  for(i in 1:3){
    
    if(!is.null(get(paste0("data",i)))) count = count+1
    
  }
  
  #extract timepoints
  timepoints = unique(data1$Time)
  
  #create summary dataframe to populate later
  all_data = data.frame(Time=rep(timepoints, 5))
  all_data$Time = all_data$Time[order(all_data$Time)]
  all_data$Bacteria = c("Total", "EryR", "TetR", "DRP", "P")
  
  #loop through all inputed datasets
  for(i in 1:count){
    
    #recover name of data and access as "d"
    d = get(paste0("data",i))
    
    #calculate cfu based on plate count and dilution factor
    d[,2] = d[,2]*(10^(d[,3]))*(1000/volume)
    d[,4] = d[,4]*(10^(d[,5]))*(1000/volume)
    d[,6] = d[,6]*(10^(d[,7]))*(1000/volume)
    
    #correct for T+E plate since 500uL are plated
    d[,8] = d[,8]*(10^(d[,9]))
    #correct for P plate since 200uL are plated
    d[,10] = d[,10]*(10^(d[,11]))*(1000/200)
    
    #remove dil
    d = d %>%
      select(Time, contains("count")) 
    
    colnames(d) = gsub("_count", "", colnames(d))
    
    d = melt(d, id.vars = "Time")
    
    #merge together the double timepoints (0 and 24)
    d = d %>%
      rename(Plate = variable,
             cfu = value) %>%
      group_by(Time, Plate) %>%
      summarise(cfu = mean(cfu)) %>%
      ungroup
    
    #loop through each time and deduct number of DRP from T+E to work out number of 327 and 201kt7
    for(times in unique(d$Time)){
      
      #order is: BHIA, E, T, T+E, P
      good_rows = which(d$Time == times)
      d$cfu[good_rows[2]] = d$cfu[good_rows[2]] - d$cfu[good_rows[4]]
      d$cfu[good_rows[3]] = d$cfu[good_rows[3]] - d$cfu[good_rows[4]]
      
    }
    
    d$Bacteria = rep(c("Total", "EryR", "TetR", "DRP", "P"), nrow(d)/5)
    
    #clean up:
    d = d[,-2]
    
    
    
    #plot:
    if(logged==T){
      
      gg = ggplot(d, aes(Time, cfu, colour=Bacteria))+
        geom_point()+
        geom_line()+
        scale_x_continuous(limits=c(0,max(timepoints)), breaks=seq(0,max(timepoints),2))+
        scale_y_continuous(trans=log10_trans(),
                           breaks=trans_breaks("log10", function(x) 10^x),
                           labels=trans_format("log10", math_format(10^.x)))+
        ylab("cfu/ml")+
        xlab("Time (hours)")+
        theme_bw()
      
      print(gg)
      ggsave(paste0("data",i,"_log.png"), plot=gg)
      
    } else {###
      
      gg = ggplot(d, aes(Time, cfu, colour=Bacteria))+
        geom_point()+
        geom_line()+
        scale_x_continuous(limits=c(0,max(timepoints)), breaks=seq(0,max(timepoints),2))+
        ylab("cfu/ml")+
        xlab("Time (hours)")+
        theme_bw()
      
      print(gg)
      ggsave(paste0("data",i,".png"), plot=gg)
      
    }
    
    #save dataset
    write.csv(d, paste0("data",i,".csv"), row.names = F)
    
    #bind together (then used below if inputed more than one dataset)
    all_data = cbind(all_data, d$cfu)
    colnames(all_data)[2+i] = paste0("data",i)
    
  }
  
  #summary if more than 1 dataset
  if(count>1){
    
    #calculate mean and standard error at 95% confidence
    all_data$Mean = rowMeans(all_data[,3:dim(all_data)[2]])
    all_data$se = apply(all_data[,3:(dim(all_data)[2]-1)], 1,
                        FUN = function(x){sd(x)/sqrt(count)})
    
    all_data$se_min = all_data$Mean-all_data$se
    #correction to avoid error when logging:
    all_data$se_min[which(all_data$se_min < 0)] = 0
    
    all_data$se_max = all_data$Mean+all_data$se
    
    if(logged==T){
      
      gg = ggplot(all_data, aes(x=Time, y=Mean, colour=Bacteria))+
        geom_point()+
        geom_line()+
        geom_errorbar(aes(x=Time, ymin=se_min, ymax=se_max))+
        scale_x_continuous(limits=c(0,max(timepoints)), breaks=seq(0,max(timepoints),2))+
        scale_y_continuous(trans=log10_trans(),
                           breaks=trans_breaks("log10", function(x) 10^x),
                           labels=trans_format("log10", math_format(10^.x)))+
        ylab("cfu/ml")+
        xlab("Time (hours)")+
        theme_bw()
      
      print(gg)
      ggsave("summary_log.png", plot=gg)
      
    } else {
      
      gg = ggplot(all_data, aes(x=Time, y=Mean, colour=Bacteria))+
        geom_point()+
        geom_line()+
        geom_errorbar(aes(x=Time, ymin=se_min, ymax=se_max))+
        scale_x_continuous(limits=c(0,max(timepoints)), breaks=seq(0,max(timepoints),2))+
        ylab("cfu/ml")+
        xlab("Time (hours)")+
        theme_bw()
      
      print(gg)
      ggsave("summary.png", plot=gg)
      
    }
    
    #save summary dataset
    write.csv(all_data, "summary.csv", row.names = F)
    
  }
}


#load datasets to analyse
data1 = rio::import(here::here("Lab", "Transduction", "07_10_20", "071020.xlsx"))
data2 = rio::import(here::here("Lab", "Transduction", "14_10_20", "141020_1.xlsx"))
data3 = rio::import(here::here("Lab", "Transduction", "14_10_20", "141020_2.xlsx"))

#run function
calculate_cfu_trans(data1, data2, data3, logged=T)

#A quick script to produce figure 7 from my upgrade report comparing the 3 model outputs
#This assumes that the 3 code files for the model are located in the current working directory

library(ggplot2)
library(deSolve)

#run the full model:
source("transduction_full.R")
res_full = results

#run the intermediate model:
source("transduction_intermediate.R")
res_inter = results

#run the mass-action model:
source("transduction_mass-action.R")
res_mass = results


#plot outputs together:
ggplot()+
  geom_line(data=res_full, aes(time, Bet, colour="1",linetype="2"), size=1)+
  geom_line(data=res_inter, aes(time, Bet, colour="2",linetype="2"), size=1)+
  geom_line(data=res_mass, aes(time, Bet, colour="3",linetype="2"), size=1)+
  geom_line(data=res_full, aes(time, (Bet+Be+Bt), colour="1",linetype="1"), size=1)+
  geom_line(data=res_inter, aes(time, (Bet+Be+Bt), colour="2",linetype="1"), size=1)+
  geom_line(data=res_mass, aes(time, (Bet+Be+Bt), colour="3",linetype="1"), size=1)+
  scale_x_continuous(limits=c(0,24), breaks=seq(0,24,2))+
  scale_y_continuous(breaks=seq(0,1e9,2e8))+
  xlab("Time (hours)")+
  ylab("Number of DRP")+
  scale_color_manual(name="Model:", breaks=c("1","2","3"), values=c("royalblue", "red3", "green3"), labels=c("Full", "Intermediate", "Mass-action"))+
  scale_linetype_manual(name="Bacteria:", breaks=c("1","2"), values=c("solid","dotted"), labels=c("Total", "DRP"))


#save image:
ggsave("compare_models.png")

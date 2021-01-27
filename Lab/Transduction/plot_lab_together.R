library(ggplot2)
library(scales)
library(dplyr)

data103 = read.csv(here::here("Lab", "Transduction", "summary_10_3.csv"))
data103$Data = "103"

data104 = read.csv(here::here("Lab", "Transduction", "summary_10_4.csv"))
data104$Data = "104"

data105 = read.csv(here::here("Lab", "Transduction", "summary_10_5.csv"))
data105$Data = "105"

data = rbind(data103, data104, data105)

data.labels = c("10^3", "10^4", "10^5")
names(data.labels) = c("103", "104", "105")

ggplot(data, aes(x=Time, y=Mean, colour=Bacteria))+
  geom_point()+
  geom_line()+
  geom_errorbar(aes(x=Time, ymin=se_min, ymax=se_max)) +
  scale_x_continuous(limits=c(0,max(data$Time)), breaks=seq(0,max(data$Time),2))+
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x)))+
  ylab("cfu/ml")+
  xlab("Time (hours)")+
  theme_bw() +
  facet_wrap(~Data,
             labeller = labeller(Data = data.labels))

ggsave("all_concentrations_wrap.png")


small_data = data %>%
  filter(Bacteria %in% c("Total", "P", "DRP"))


ggplot(small_data, aes(x=Time, y=Mean, colour=Data, linetype = Bacteria))+
  geom_point()+
  geom_line(size = 1)+
  geom_errorbar(aes(x=Time, ymin=se_min, ymax=se_max), linetype = "solid") +
  scale_x_continuous(limits=c(0,max(small_data$Time)), breaks=seq(0,max(small_data$Time),2))+
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x)))+
  ylab("cfu/ml")+
  xlab("Time (hours)")+
  theme_bw() +
  scale_colour_discrete(name = "Phage start", labels = c("10^3", "10^4", "10^5"))

ggsave("all_concentrations_together.png")

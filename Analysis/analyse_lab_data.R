
library(ggplot2)
library(scales)
library(dplyr)

data = read.csv(here::here("Data", "growth_summary.csv"))
data3 = read.csv(here::here("Data", "transduction_summary_10_3.csv"))
data4 = read.csv(here::here("Data", "transduction_summary_10_4.csv"))
data5 = read.csv(here::here("Data", "transduction_summary_10_5.csv"))

data$Type = "No phage"
data3$Type = "10^3"
data4$Type = "10^4"
data5$Type = "10^5"

data_all = rbind(data, data3, data4, data5) %>%
  filter(Bacteria != "Total") %>%
  mutate(Bacteria=replace(Bacteria, Bacteria=="EryR", "Be")) %>%
  mutate(Bacteria=replace(Bacteria, Bacteria=="TetR", "Bt")) %>%
  mutate(Bacteria=replace(Bacteria, Bacteria=="DRP", "Bet")) %>%
  mutate(Bacteria=replace(Bacteria, Bacteria=="P", "Pl"))

data_all$group = data_all$Bacteria
data_all = data_all %>% 
  mutate(group = replace(group, Bacteria %in% c("Be","Bt","Bet"), "other"))


data_all$Type = as.factor(data_all$Type)
data_all$Type = factor(data_all$Type, levels(data_all$Type)[c(4,1,2,3)])
levels(data_all$Type) = c("No phage",
                          expression(paste(10^3, " phage")),
                          expression(paste(10^4, " phage")),
                          expression(paste(10^5, " phage")))

data_all$Bacteria = as.factor(data_all$Bacteria)
data_all$Bacteria = factor(data_all$Bacteria, c("Bacteria:", levels(data_all$Bacteria)[c(1,3,2)],
                                                " ", "Phage:", levels(data_all$Bacteria)[4]))

ggplot(data_all %>% filter(Type != "No phage"), aes(x=Time, y=Mean, colour=Bacteria))+
  geom_point(aes(shape = group), size = 2)+
  geom_line(alpha=0.3)+
  geom_line(data = data_all %>% filter(Type != "No phage") %>% filter(Time < 9), size = 0.8) +
  geom_line(data = data_all %>% filter(Type != "No phage") %>% filter(Time > 15), size = 0.8) +
  geom_errorbar(aes(x=Time, ymin=se_min, ymax=se_max), size = 0.7) +
  scale_x_continuous(breaks=seq(0,max(data$Time),4))+
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x)))+
  labs(y = "cfu or pfu per mL", x = "Time (hours)", colour = "")+
  theme_bw() +
  facet_wrap(~Type, labeller = label_parsed) +
  scale_colour_manual(values=c("white","#685cc4","#6db356","#c2484d","white","white","#c88a33"),
                      drop = FALSE,
                      labels = c("Bacteria:",
                                 bquote(B[E]),
                                 bquote(B[T]),
                                 bquote(B[ET]),
                                 " ",
                                 "Phage:",
                                 bquote(P[L])),
                      guide = guide_legend(override.aes = list(shape = c(19,19,19,19,19,17,17)))) +
  guides(linetype = F, shape = F) +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text.x = element_text(size=12))


ggsave(here::here("Figures", "fig2.png"))

ggplot(data_all %>% filter(Type != "No phage"), aes(x=Time, y=Mean,
                                                    colour=Bacteria, linetype= Type))+
  geom_point(aes(shape = group), size = 2)+
  geom_line(alpha=0.3)+
  geom_line(data = data_all %>% filter(Type != "No phage") %>% filter(Time < 9), size = 0.8) +
  geom_line(data = data_all %>% filter(Type != "No phage") %>% filter(Time > 15), size = 0.8) +
  geom_errorbar(aes(x=Time, ymin=se_min, ymax=se_max), size = 0.7) +
  scale_x_continuous(breaks=seq(0,max(data$Time),4))+
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x)))+
  labs(y = "cfu or pfu per mL", x = "Time (hours)", colour = "Organism:", linetype = "Initial phage:")+
  theme_bw() +
  scale_colour_manual(values=c("white","#685cc4","#6db356","#c2484d","white","white","#c88a33"),
                      drop = FALSE,
                      labels = c("Bacteria:",
                                 bquote(B[E]),
                                 bquote(B[T]),
                                 bquote(B[ET]),
                                 " ",
                                 "Phage:",
                                 bquote(P[L])),
                      guide = guide_legend(override.aes = list(shape = c(19,19,19,19,19,19,17)))) +
  scale_linetype_manual(values = c(1,2,3),
                        labels = c(expression(10^3),
                                   expression(10^4),
                                   expression(10^5))) +
  guides(shape = F) +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.x = element_text(size=12))

ggsave(here::here("Figures", "supp_fig2.png"))


fitness = function(data, str1, str2){
  
  str1_0 = data %>%
    filter(Time == 0 & Bacteria == str1) %>%
    select(c("data1", "data2", "data3"))
  
  str1_24 = data %>%
    filter(Time == 24 & Bacteria == str1) %>%
    select(c("data1", "data2", "data3"))
  
  str2_0 = data %>%
    filter(Time == 0 & Bacteria == str2) %>%
    select(c("data1", "data2", "data3"))
  
  str2_24 = data %>%
    filter(Time == 24 & Bacteria == str2) %>%
    select(c("data1", "data2", "data3"))
  
  str1_w = log(str1_24/str1_0)
  str2_w = log(str2_24/str2_0)
  
  tt = str1_w/str2_w
  tt$Mean = rowMeans(tt)
  tt$sd = sd(tt[1,1:3])
  tt$sd_min = tt$Mean - tt$sd
  tt$sd_max = tt$Mean + tt$sd
  
  tt
  
}

fitness(data, "EryR", "TetR")
fitness(data, "DRP", "EryR")
fitness(data, "DRP", "TetR")


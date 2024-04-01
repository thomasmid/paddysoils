library(ggplot2)
library(readxl)
library(multcompView)
library(gridExtra)
library(grid)
library(tidyverse)
library(writexl)
library(ggthemes)    
library(ggpubr)

#function for labelling
generate_label <- function(tukey){
  tukey.labels <- data.frame(multcompLetters(tukey[[1]][,4])['Letters'])
  tukey.labels$group=rownames(tukey.labels)
  tukey.labels=tukey.labels[order(tukey.labels$group) , ]
  names(tukey.labels)<-c("letter",'group')
  return(tukey.labels)}

#ph plot
pH_data <- read_excel("soilplots_pH.xlsx")

pH_tukey <- TukeyHSD(aov(pH_data$pH ~ pH_data$group), ordered = FALSE, conf.level = 0.95)
pH_label <-generate_label(pH_tukey)

pH_plot <- ggplot(pH_data, aes(x = group, y = as.numeric(pH)))+
    theme_economist()+
    theme(text = element_text(size = 15))+
    ylim(cbind(6.4,7.5))+
    labs(caption="",x="",y="")+
    scale_x_discrete(limits=c("C","SORn","GCn","WOn","WHEn","SORc","SORf","GCc","GCf","WOc","WOf","WHEc","WHEf"))+
    geom_boxplot(color="#A2061C",fill="#DCD8E0")+
    geom_text(data = pH_label, aes(x = group, y = 6.4, label = letter), vjust=-.5,hjust=0, size=6.5) 
    
ggsave(filename="pH_plot.png",pH_plot, width = 15,height =8,dpi = 400)


#C_change plot
CN_data <- read_excel("soilplots_CN.xlsx")

C_tukey <- TukeyHSD(aov(CN_data$C_change ~ CN_data$group), ordered = FALSE, conf.level = 0.95)
C_label <-generate_label(C_tukey)

C_plot <- ggplot(CN_data, aes(x = group, y = as.numeric(C_change)))+
  theme_economist()+
  theme(text = element_text(size = 15))+
  coord_flip()+
  ylim(cbind(-16,3))+
  labs(caption="",x="",y="")+
  scale_x_discrete(limits=c("WHEf","WHEc","WOf","WOc","GCf","GCc","SORf","SORc","C"))+
  stat_summary(fun=mean, geom="point",shape=3,size=7,stroke=2,color="black")+
  geom_dotplot(binaxis = "y", stackdir = "center",dotsize=0.4,stroke = 3,color="#A2061C", stackratio = 0.5,fill="transparent",binwidth = 0.8)+
  geom_text(data = C_label, aes(x = group, y = -16, label = letter), vjust=-0.4,hjust=0, size=6.5) 

ggsave(filename="C_plot.png",C_plot, width = 15,height =8,dpi = 400)

#N_change plot
N_tukey <- TukeyHSD(aov(CN_data$N_change~ CN_data$group), ordered = FALSE, conf.level = 0.95)
N_label <-generate_label(N_tukey)

N_plot <- ggplot(CN_data, aes(x = group, y = as.numeric(N_change)))+
  theme_economist()+
  theme(text = element_text(size = 15))+
  coord_flip()+
  ylim(cbind(-400,60))+
  labs(caption="",x="",y="")+
  scale_x_discrete(limits=c("WHEf","WHEc","WOf","WOc","GCf","GCc","SORf","SORc","C"))+
  stat_summary(fun=mean, geom="point",shape=3,size=7,stroke=2,color="black")+
  geom_dotplot(binaxis = "y", stackdir = "center",dotsize=0.4,stroke = 3,color="#A2061C", stackratio = 0.5,fill="transparent",binwidth = 20)+
  geom_text(data = N_label, aes(x = group, y = -400, label = letter), vjust=-0.4,hjust=0, size=6.5) 

ggsave(filename="N_plot.png",N_plot, width = 15,height =8,dpi = 400)


C_loss_plot <- ggplot(CN_data, aes(x = as.numeric(total_miner_mean_2_3), y = as.numeric(C_change)))+
  theme_economist()+
  theme(text = element_text(size = 15))+
  #coord_flip()+
  #ylim(cbind(-400,60))+
  labs(caption="",x="",y="")+
  #stat_summary(fun=mean, geom="point",shape=3,size=7,stroke=2,color="black")+
  geom_point(shape= 1, size=3,stroke = 2.5,color="#136699",fill="transparent")+
  geom_smooth(method=lm, linewidth= 1.8,color="#A2061C", fill="#FCA9B5")+
  stat_regline_equation(size=7,label.y = 2,label.x=55, aes(label = after_stat(eq.label))) +
  stat_regline_equation(size=7,label.y = 1,label.x=55, aes(label = after_stat(rr.label)))
ggsave(filename="C_loss_plot.png",C_loss_plot, width = 15,height =8,dpi = 400)

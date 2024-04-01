#packages
library(readxl)
library(writexl)
library(multcompView)
library(ggplot2)
library(ggthemes)    
library(ggpubr)
library(tidyverse)

library(car) #D# for levene test
library(FSA) #D# for Dunn's test (posthoc for Kurskal Wallis)
# define target parameter
parameters <- c("exp","a","b","k","a_k","m","m_SE","m_SE_rel","d","d_SE","d_SE_rel","pmCO2","omCO2r","fm_CO2","fm_CH4","fm_Cmin","fm_respratio","am_t","nmCO2")

####fitting exp 1####
data_CO2 <- read_excel("1_CO2_full.xlsx")
data_CH4 <- read_excel("1_CH4_full.xlsx")

# model parameter start values
start_CO2 <- list(a = 0.3, b = 55, k = 0.005)
start_CH4 <- list(m = 0.05, d = 200)

# loop initiation
x <- colnames(subset(data_CO2, select = - c(day, GCn1, GCn2, GCn3, GCn4)))#remove columns if necessary 
exp_1 <- data.frame(parameters)

for (i in seq(x)) {
fit_CO2 <- nls(data_CO2[[paste0(x[[i]])]] ~ (-a/k) * exp(-k*day) + b,
  data = data_CO2,
  start = start_CO2)
fit_CH4 <- nls(data_CH4[[c(paste0(x[i]))]] ~ (m*(data_CH4$day-d))/(1+exp(-1000*(day - d))),
  data = data_CH4,
  start = start_CH4)

# model summary
sum_CO2 <- summary(fit_CO2)
sum_CH4 <- summary(fit_CH4)

# parameters CO2-fit
a <- sum_CO2$coefficients[[1]]
b <- sum_CO2$coefficients[[2]]
k <- sum_CO2$coefficients[[3]]
# methanogenic productivity
m <- sum_CH4$coefficients[[1]]
# delay of onset of methanogenesis
d <- sum_CH4$coefficients[[2]]

#errors of parameters
m_SE <- sum_CH4$coefficients[[3]]
d_SE <- sum_CH4$coefficients[[4]]
m_SE_rel <- sum_CH4$coefficients[[3]]/sum_CH4$coefficients[[1]]
d_SE_rel <- sum_CH4$coefficients[[4]]/sum_CH4$coefficients[[2]]

# define functions
f_CO2 <- 
    function(day){(-a/k)*exp(-k*day)+b}
f_CO2_deriv <-
    function(day){a*exp(-k*day)}
f_CH4 <-
    function(day){m*(day-d)/(1+exp(-10*(day- d)))}
  
# define results
m <- m
a_k <- a/k

# onset of methanogenesis
pmCO2 <- f_CO2(d)
omCO2r <- f_CO2_deriv(d)

# dominant methanogenesis (2 = CO2_deriv/m)
#earliest_t <-(-log(2*m/a)/k)
#dm_t <- max(c(om_t,earliest_t))
#dm_CO2 <- f_CO2(dm_t)
#dm_CO2_deriv <- f_CO2_deriv(dm_t)
#dm_CH4 <- f_CH4(dm_t)
#dm_Cmin <- f_CO2(dm_t)+f_CH4(dm_t)
#dm_NMCO2 <- f_CO2(dm_t)-f_CH4(dm_t)
#dm_respratio <- dm_CO2_deriv/CH4_deriv

# considerable methanogenesis (5 = CH4)
#cm_t <- (5/m)+d
#cm_CO2 <- f_CO2(cm_t)
#cm_CO2_deriv <- f_CO2_deriv(cm_t)
#cm_CH4 <- f_CH4(cm_t)
#cm_Cmin <- f_CO2(cm_t)+f_CH4(cm_t)
#cm_NMCO2 <- f_CO2(cm_t)-f_CH4(cm_t)
#cm_respratio <- cm_CO2_deriv/CH4_deriv

# final measurement
fm_t <- tail(data_CO2$day[], n=1)
fm_CO2 <- f_CO2(fm_t)
fm_CH4 <- f_CH4(fm_t)
fm_Cmin <- f_CO2(fm_t)+f_CH4(fm_t)
fm_CO2r <- f_CO2_deriv(fm_t)
fm_respratio <- fm_CO2r/m

# absolute methanogenesis
am_t <- (-log(m/a)/k)
am_NMCO2 <- f_CO2(am_t)-f_CH4(am_t)
nmCO2 <- am_NMCO2 - pmCO2

# get results
exp_1[,i+1] <-
    c("Common set-up",a,b,k,a_k,m,m_SE,m_SE_rel,d,d_SE,d_SE_rel,pmCO2,omCO2r,fm_CO2,fm_CH4,fm_Cmin,fm_respratio,am_t,nmCO2)
# close loop
}
# finish data frame
names(exp_1) <- c("parameter", x)

####fitting exp 2####
data_CO2_2 <- read_excel("2_CO2_full.xlsx")
data_CH4_2 <- read_excel("2_CH4_full.xlsx")

# model parameter start values
start_CO2 <- list(a = 0.5, b = 100, k = 0.005)
start_CH4 <- list(m = 0.05, d = 100)

# loop initiation
x <- colnames(subset(data_CO2_2, select = - c(day, GCc1, GCc2, GCc3, GCc4, GCf1, GCf2, GCf3, GCf4)))#remove columns if necessary 
exp_2 <- data.frame(parameters)

for (i in seq(x)) {
fit_CO2 <- nls(data_CO2_2[[paste0(x[[i]])]] ~ (-a/k) * exp(-k*day) + b,
  data = data_CO2_2,
  start = start_CO2)
fit_CH4 <- nls(data_CH4_2[[c(paste0(x[i]))]] ~ (m*(day-d))/(1+exp(-1000*(day - d))),
  data = data_CH4_2,
  start = start_CH4)

# model summary
sum_CO2 <- summary(fit_CO2)
sum_CH4 <- summary(fit_CH4)

# parameters CO2-fit
a <- sum_CO2$coefficients[[1]]
b <- sum_CO2$coefficients[[2]]
k <- sum_CO2$coefficients[[3]]
# methanogenic productivity
m <- sum_CH4$coefficients[[1]]
# delay of onset of methanogenesis
d <- sum_CH4$coefficients[[2]]

#errors of parameters
m_SE <- sum_CH4$coefficients[[3]]
d_SE <- sum_CH4$coefficients[[4]]
m_SE_rel <- sum_CH4$coefficients[[3]]/sum_CH4$coefficients[[1]]
d_SE_rel <- sum_CH4$coefficients[[4]]/sum_CH4$coefficients[[2]]

# define functions
f_CO2 <- 
  function(day){(-a/k)*exp(-k*day)+b}
f_CO2_deriv <-
  function(day){a*exp(-k*day)}
f_CH4 <-
  function(day){m*(day-d)/(1+exp(-10*(day- d)))}

# define results
m <- m
a_k <- a/k

# onset of methanogenesis
pmCO2 <- f_CO2(d)
omCO2r <- f_CO2_deriv(d)

# final measurement
fm_t <- tail(data_CO2$day[], n=1)
fm_CO2 <- f_CO2(fm_t)
fm_CH4 <- f_CH4(fm_t)
fm_Cmin <- f_CO2(fm_t)+f_CH4(fm_t)
fm_CO2r <- f_CO2_deriv(fm_t)
fm_respratio <- fm_CO2r/m

# absolute methanogenesis
am_t <- (-log(m/a)/k)
am_NMCO2 <- f_CO2(am_t)-f_CH4(am_t)
nmCO2 <- am_NMCO2 - pmCO2

# get results
exp_2[,i+1] <-
  c("Modified set-up",a,b,k,a_k,m,m_SE,m_SE_rel,d,d_SE,d_SE_rel,pmCO2,omCO2r,fm_CO2,fm_CH4,fm_Cmin,fm_respratio,am_t,nmCO2)
# close loop
}
#finish data frame
names(exp_2) <- c("parameter", x)

####transform data frames####
exp_1 <- t(exp_1)
exp_1 <- as.data.frame(exp_1)
names(exp_1) <- as.matrix(exp_1[1,])
exp_1 <- exp_1[-1,]
group_1 <- c("C","C","C","C","SORn","SORn","SORn","SORn","WOn","WOn","WOn","WOn","WHEn","WHEn","WHEn","WHEn")
exp_1$group <- group_1

exp_2 <- t(exp_2)
exp_2 <- as.data.frame(exp_2)
names(exp_2) <- as.matrix(exp_2[1,])
exp_2 <- exp_2[-1,]
group_2 <- c("C","C","C","C","C","C","SORc","SORc","SORc","SORc","SORf","SORf","SORf","SORf","WOc","WOc","WOc","WOc","WOf","WOf","WOf","WOf","WHEc","WHEc","WHEc","WHEc","WHEf","WHEf","WHEf","WHEf")
exp_2$group <- group_2

####add columns: pH and dummies####
pH_data <- read_excel("pH_data.xlsx")
exp_1$pH <- pH_data$pH[pH_data$exp=="Common set-up"]
exp_2$pH <- pH_data$pH[pH_data$exp=="Modified set-up"]
exp_1$exp2 <- 0
exp_1$exp3 <- 0
exp_2$exp2 <- 1
exp_2$exp3 <- 0
exp_1$BC <- c(0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1)
exp_2$BC <- c(0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
exp_1$fine <- c(0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1)
exp_2$fine <- c(0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1)
exp_1$coarse <- 0
exp_2$coarse <- c(0,0,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0)
exp_1$SOR <- c(0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0)
exp_2$SOR <- c(0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
exp_1$WO <- c(0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0)
exp_2$WO <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0)
exp_1$WHE <- c(0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1)
exp_2$WHE <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1)
#D# es gibt übrigens in R die rep() funktion, mit der man Zahlenreihen erzeugen kann...

####pre-tests####

#### plots ####
#prepare Tukey test 

##D## hab die funktion so erweitert, dass sie tukey und dunns test verwenden kann
generate_label <- function(test_res){
  if(grepl("tukey", deparse(substitute(test_res)), fixed = TRUE)){#if tukey test
    test.labels <- data.frame(multcompLetters(test_res[[1]][,4])['Letters'])
  }else if(grepl("dunn", deparse(substitute(test_res)), fixed = TRUE)){#if dunns test
    test.clean <- setNames(c(test_res[[2]][,4]), #create named num to match tukey
                           c(gsub(" ", "", test_res[[2]][,1], fixed=TRUE))) #remove whitespace from names
    test.labels <- data.frame(multcompLetters(test.clean)['Letters'])
  }else{}#do nothing (possibilty to add another test)
  test.labels$group=rownames(test.labels)
  test.labels=test.labels[order(test.labels$group) , ]
  names(test.labels)<-c("letter",'group')
  return(test.labels)}

#ph plot
pH_tukey <- TukeyHSD(aov(c(exp_1$pH,exp_2$pH) ~ c(exp_1$group,exp_2$group)), ordered = FALSE, conf.level = 0.95)
pH_label <-generate_label(pH_tukey)

pH_plot <- ggplot(pH_data,aes(x = group, y = pH))+
  theme_bw()+
  theme(text = element_text(size = 15))+
  ylim(cbind(6.4,7.5))+
  labs(caption="",x="",y="")+
  scale_x_discrete(limits=c("C","SORn","WOn","WHEn","SORc","SORf","WOc","WOf","WHEc","WHEf"))+
  geom_boxplot(color="#A2061C",fill="#DCD8E0")+
  geom_text(data = pH_label, aes(x = group, y = 6.4, label = letter), vjust=-.5,hjust=0, size=6.5) 

ggsave(filename="pH_plot.png",pH_plot, width = 15,height =8,dpi = 400)

#mean accumulation plots without modelling
means1 <- read_excel("1_means_full.xlsx")
means2 <- read_excel("2_means_full.xlsx")
#D# es gibt die Möglichkeit, Daten ins Langformat zu bringen, damit man nicht 10 x geom_point aufrufen muss...
#D# und wenn die x Achse gleich bleibt, kann man die einfach in ggplot() reinschreiben
plot_means1 <- ggplot()+
  theme_minimal()+
  theme(text = element_text(size = 15))+
  xlim(cbind(0,290))+
  ylim(cbind(0,60))+
  labs(caption="",x="",y="")+
  geom_point(aes(x=as.numeric(means1$day),y=as.numeric(means1$C_CH4)),shape= 1, size=3,stroke=1.5,color="black",fill="transparent")+
  geom_point(aes(x=as.numeric(means1$day),y=as.numeric(means1$SOR_CH4)),shape= 1, size=3,stroke=1.5,color="#A2061C",fill="transparent")+
  geom_point(aes(x=as.numeric(means1$day),y=as.numeric(means1$WO_CH4)),shape= 1, size=3,stroke=1.5,color="#A76E19",fill="transparent")+
  geom_point(aes(x=as.numeric(means1$day),y=as.numeric(means1$WHE_CH4)),shape= 1, size=3,stroke=1.5,color="#136699",fill="transparent")+
  geom_point(aes(x=as.numeric(means1$day),y=as.numeric(means1$C_CO2)),shape= 2, size=3,stroke=1.5,color="black",fill="transparent")+
  geom_point(aes(x=as.numeric(means1$day),y=as.numeric(means1$SOR_CO2)),shape= 2, size=3,stroke=1.5,color="#A2061C",fill="transparent")+
  geom_point(aes(x=as.numeric(means1$day),y=as.numeric(means1$WO_CO2)),shape= 2, size=3,stroke=1.5,color="#A76E19",fill="transparent")+
  geom_point(aes(x=as.numeric(means1$day),y=as.numeric(means1$WHE_CO2)),shape= 2, size=3,stroke=1.5,color="#136699",fill="transparent")
ggsave(filename="plot_means1.png",plot_means1, width = 15,height =6,dpi = 400)

  
plot_means2 <- ggplot()+
  theme_bw()+
  theme(text = element_text(size = 15))+
  xlim(cbind(0,185))+
  ylim(cbind(0,60))+
  labs(caption="",x="",y="")+
  geom_point(aes(x=as.numeric(means2$day),y=as.numeric(means2$C_CH4)),shape= 1, size=3,stroke=1.5,color="black",fill="transparent")+
  geom_point(aes(x=as.numeric(means2$day),y=as.numeric(means2$SORc_CH4)),shape= 1, size=3,stroke=1.5,color="#C76A77",fill="transparent")+
  geom_point(aes(x=as.numeric(means2$day),y=as.numeric(means2$SORf_CH4)),shape= 1, size=3,stroke=1.5,color="#A2061C",fill="transparent")+
  geom_point(aes(x=as.numeric(means2$day),y=as.numeric(means2$WOc_CH4)),shape= 1, size=3,stroke=1.5,color="#CAA875",fill="transparent")+
  geom_point(aes(x=as.numeric(means2$day),y=as.numeric(means2$WOf_CH4)),shape= 1, size=3,stroke=1.5,color="#A76E19",fill="transparent")+
  geom_point(aes(x=as.numeric(means2$day),y=as.numeric(means2$WHEc_CH4)),shape= 1, size=3,stroke=1.5,color="#71A3C2",fill="transparent")+
  geom_point(aes(x=as.numeric(means2$day),y=as.numeric(means2$WHEf_CH4)),shape= 1, size=3,stroke=1.5,color="#136699",fill="transparent")+
  geom_point(aes(x=as.numeric(means2$day),y=as.numeric(means2$C_CO2)),shape= 2, size=3,stroke=1.5,color="black",fill="transparent")+
  geom_point(aes(x=as.numeric(means2$day),y=as.numeric(means2$SORc_CO2)),shape= 2, size=3,stroke=1.5,color="#C76A77",fill="transparent")+
  geom_point(aes(x=as.numeric(means2$day),y=as.numeric(means2$SORf_CO2)),shape= 2, size=3,stroke=1.5,color="#A2061C",fill="transparent")+
  geom_point(aes(x=as.numeric(means2$day),y=as.numeric(means2$WOc_CO2)),shape= 2, size=3,stroke=1.5,color="#CAA875",fill="transparent")+
  geom_point(aes(x=as.numeric(means2$day),y=as.numeric(means2$WOf_CO2)),shape= 2, size=3,stroke=1.5,color="#A76E19",fill="transparent")+
  geom_point(aes(x=as.numeric(means2$day),y=as.numeric(means2$WHEc_CO2)),shape= 2, size=3,stroke=1.5,color="#71A3C2",fill="transparent")+
  geom_point(aes(x=as.numeric(means2$day),y=as.numeric(means2$WHEf_CO2)),shape= 2, size=3,stroke=1.5,color="#136699",fill="transparent")
ggsave(filename="plot_means2.png",plot_means2, width = 15,height =6,dpi = 400)

#density plots
Experiment <- c(as.factor(exp_1$exp),as.factor(exp_2$exp))
density_d <- ggplot() +
  geom_density(alpha=0.7,aes(x=c(as.numeric(exp_1$d_SE_rel),as.numeric(exp_2$d_SE_rel)),fill=Experiment))+
  geom_vline(aes(xintercept=median(as.numeric(exp_1$d_SE_rel))),linetype=2, linewidth=1,color="#136699")+
  geom_vline(aes(xintercept=median(as.numeric(exp_2$d_SE_rel))),linetype=2, linewidth=1,color="#A76E19")+
  geom_density(alpha=0.7,aes(x=c(as.numeric(exp_1$m_SE_rel),as.numeric(exp_2$m_SE_rel)),fill=Experiment))+
  geom_vline(aes(xintercept=median(as.numeric(exp_1$m_SE_rel))),linetype=2, linewidth=1,color="#136699")+
  geom_vline(aes(xintercept=median(as.numeric(exp_2$m_SE_rel))),linetype=2, linewidth=1,color="#A76E19")+
  
  scale_fill_manual(values=c("Common set-up"="#136699","Modified set-up"="#A76E19"))+
  theme_bw()+
  xlim(cbind(0,0.13))+
  scale_x_continuous()
  theme(text = element_text(size=20), legend.text = element_text(size=20))+
  labs(caption="",x="",y="")
ggsave(filename="density_d.png",density_d, width = 15,height =6,dpi = 400)

density_m <- ggplot() +
  geom_density(alpha=0.7,aes(x=c(as.numeric(exp_1$m_SE_rel),as.numeric(exp_2$m_SE_rel)),fill=Experiment))+
  geom_vline(aes(xintercept=median(as.numeric(exp_1$m_SE_rel))),linetype=2, linewidth=1,color="#136699")+
  geom_vline(aes(xintercept=median(as.numeric(exp_2$m_SE_rel))),linetype=2, linewidth=1,color="#A76E19")+
  scale_fill_manual(values=c("Common set-up"="#136699","Modified set-up"="#A76E19"))+
  theme_bw()+
  xlim(cbind(0,0.13))+
  theme(text = element_text(size=20), legend.text = element_text(size=20))+
  labs(caption="",x="",y="")
ggsave(filename="density_m.png",density_m, width = 15,height =6,dpi = 400)

median(c(as.numeric(exp_1$d_SE_rel),as.numeric(exp_2$d_SE_rel)))
median(c(as.numeric(exp_1$m_SE_rel),as.numeric(exp_2$m_SE_rel)))

#d
#D# pre-test: normality
shapiro.test(as.numeric(exp_1$d)) #D# p>0.05 -> data normally distributed
#D# pre-test: homogeneity of variance
leveneTest(as.numeric(exp_1$d) ~ exp_1$group) #D# p < 0.05 data has heterogeneous variances

#d_1_tukey <- TukeyHSD(aov(exp_1$d ~ exp_1$group), ordered = FALSE, conf.level = 0.95)
#D# use kruskal wallis instead
d_1_kw <- kruskal.test(as.numeric(exp_1$d) ~ exp_1$group)
#D# Dunns test as post-hoc
d_1_dunn <- dunnTest(as.numeric(exp_1$d) ~ exp_1$group)

####to-do D: change plot to use dunns test####
d_1_label <-generate_label(d_1_dunn)
d_1 <- ggplot(exp_1, aes(x = group, y = as.numeric(d)))+
    theme_bw()+
    theme(text = element_text(size = 15))+
    coord_flip()+
    ylim(cbind(0,250))+
    labs(caption="",x="",y="")+
    scale_x_discrete(limits=c("WHEn","WOn","SORn","C"))+
    geom_dotplot(binaxis = "y", stackdir = "center",stroke = 4,color="#24A177", stackratio = 0.5,fill="transparent",dotsize=0.5)+
    stat_summary(fun=mean, geom="point",shape=3,size=7,stroke=2,color="black")+
    geom_text(data = d_1_label, aes(x = group, y = 0, label = letter), vjust=-0.4,hjust=0, size=6.5)+
    geom_hline(aes(yintercept=mean(as.numeric(exp_1$d[exp_1$group=="C"]))), linetype=2, color="#F8A50E",linewidth=1)
ggsave(filename="d_1.png",d_1, width = 15,height =6,dpi = 400)

#D# pre-test: normality
shapiro.test(as.numeric(exp_2$d)) #D# p>0.05 -> data normally distributed
#D# pre-test: homogeneity of variance
leveneTest(as.numeric(exp_2$d) ~ exp_2$group) #D# p > 0.05 variance is homogenous

d_2_tukey <- TukeyHSD(aov(exp_2$d ~ exp_2$group), ordered = FALSE, conf.level = 0.95)
d_2_label <-generate_label(d_2_tukey)
d_2 <- ggplot(exp_2, aes(x = group, y = as.numeric(d)))+
    theme_bw()+
    theme(text = element_text(size = 15))+
    coord_flip()+
    ylim(cbind(0,250))+
    labs(caption="",x="",y="")+
    scale_x_discrete(limits=c("WHEf","WHEc","WOf","WOc","SORf","SORc","C"))+
    geom_dotplot(binaxis = "y", stackdir = "center",stroke = 4,color="#24A177", stackratio = 0.5,fill="transparent",dotsize=0.5)+
    stat_summary(fun=mean, geom="point",shape=3,size=7,stroke=2,color="black")+
    geom_text(data = d_2_label, aes(x = group, y = 0, label = letter), vjust=-0.4,hjust=0, size=6.5)+
    geom_hline(aes(yintercept=mean(as.numeric(exp_2$d[exp_2$group=="C"]))), linetype=2, color="#F8A50E",linewidth=1)
ggsave(filename="d_2.png",d_2, width = 15,height =6,dpi = 400)

#pmCO2
#D# pre-test: normality
shapiro.test(as.numeric(exp_1$pmCO2)) #D# p>0.05 -> data normally distributed
#D# pre-test: homogeneity of variance
leveneTest(as.numeric(exp_1$pmCO2) ~ exp_1$group) #D# p>0.05 variance is homogeneous

pmCO2_1_tukey <- TukeyHSD(aov(exp_1$pmCO2 ~ exp_1$group), ordered = FALSE, conf.level = 0.95)
pmCO2_1_label <-generate_label(pmCO2_1_tukey)
pmCO2_1 <- ggplot(exp_1, aes(x = group, y = as.numeric(pmCO2)))+
    theme_bw()+
    theme(text = element_text(size = 15))+
    coord_flip()+
    ylim(cbind(28,50))+
    labs(caption="",x="",y="")+
    scale_x_discrete(limits=c("WHEn","WOn","SORn","C"))+
    geom_dotplot(binaxis = "y", stackdir = "center",stroke = 4,color="#136699", stackratio = 0.5,fill="transparent",dotsize=0.5)+
    stat_summary(fun=mean, geom="point",shape=3,size=7,stroke=2,color="black")+
    geom_text(data = pmCO2_1_label, aes(x = group, y = 28, label = letter), vjust=-0.4,hjust=0, size=6.5)+
    geom_hline(aes(yintercept=mean(as.numeric(exp_1$pmCO2[exp_1$group=="C"]))), linetype=2, color="#F8A50E",linewidth=1)
ggsave(filename="pmCO2_1.png",pmCO2_1, width = 15,height =6,dpi = 400)

#D# pre-test: normality
shapiro.test(as.numeric(exp_2$pmCO2)) #D# p>0.05 -> data normally distributed
#D# pre-test: homogeneity of variance
leveneTest(as.numeric(exp_2$pmCO2) ~ exp_2$group) #D# variance is homogeneous

pmCO2_2_tukey <- TukeyHSD(aov(exp_2$pmCO2 ~ exp_2$group), ordered = FALSE, conf.level = 0.95)
pmCO2_2_label <-generate_label(pmCO2_2_tukey)
pmCO2_2 <- ggplot(exp_2, aes(x = group, y = as.numeric(pmCO2)))+
    theme_bw()+
    theme(text = element_text(size = 15))+
    coord_flip()+
    ylim(cbind(28,50))+
    labs(caption="",x="",y="")+
    scale_x_discrete(limits=c("WHEf","WHEc","WOf","WOc","SORf","SORc","C"))+
    geom_dotplot(binaxis = "y", stackdir = "center",stroke = 4,color="#136699", stackratio = 0.5,fill="transparent",dotsize=0.5)+
    stat_summary(fun=mean, geom="point",shape=3,size=7,stroke=2,color="black")+
    geom_text(data = pmCO2_2_label, aes(x = group, y = 28, label = letter), vjust=-0.4,hjust=0, size=6.5)+
    geom_hline(aes(yintercept=mean(as.numeric(exp_2$pmCO2[exp_2$group=="C"]))), linetype=2, color="#F8A50E",linewidth=1)
ggsave(filename="pmCO2_2.png",pmCO2_2, width = 15,height =6,dpi = 400)

#m
#D# pre-test: normality
shapiro.test(as.numeric(exp_1$m)) #D# p>0.05 -> data normally distributed
#D# pre-test: homogeneity of variance
leveneTest(as.numeric(exp_1$m) ~ exp_1$group) #D# p>0.05 -> variance is homogeneous


m_1_tukey <- TukeyHSD(aov(exp_1$m ~ exp_1$group), ordered = FALSE, conf.level = 0.95)
m_1_label <-generate_label(m_1_tukey)
m_1 <- ggplot(exp_1, aes(x = group, y = as.numeric(m)))+
    theme_bw()+
    theme(text = element_text(size = 15))+
    coord_flip()+
    ylim(cbind(0,0.11))+
    labs(caption="",x="",y="")+
    scale_x_discrete(limits=c("WHEn","WOn","SORn","C"))+
    geom_dotplot(binaxis = "y", stackdir = "center",stroke = 4,color="#A76E19", stackratio = 0.5,fill="transparent",dotsize=0.5)+
    stat_summary(fun=mean, geom="point",shape=3,size=7,stroke=2,color="black")+
    geom_text(data = m_1_label, aes(x = group, y = 0, label = letter), vjust=-0.4,hjust=0, size=6.5)+
    geom_hline(aes(yintercept=mean(as.numeric(exp_1$m[exp_1$group=="C"]))), linetype=2, color="#F8A50E",linewidth=1)
ggsave(filename="m_1.png",m_1, width = 15,height =6,dpi = 400)

#D# pre-test: normality
shapiro.test(as.numeric(exp_2$m)) #D# p < 0.05 -> data not normally distributed
#D# pre-test: homogeneity of variance
leveneTest(as.numeric(exp_2$m) ~ exp_2$group) #D# p>0.05 -> variance is homogeneous

#D# use kruskal wallis instead
m_2_kw <- kruskal.test(exp_2$m ~ exp_2$group)
#D# Dunns test as post-hoc
m_2_dunn <- dunnTest(as.numeric(exp_2$m) ~ exp_2$group)

#m_2_tukey <- TukeyHSD(aov(exp_2$m ~ exp_2$group), ordered = FALSE, conf.level = 0.95)
m_2_label <-generate_label(m_2_tukey)
m_2 <- ggplot(exp_2, aes(x = group, y = as.numeric(m)))+
    theme_bw()+
    theme(text = element_text(size = 15))+
    coord_flip()+
    ylim(cbind(0,0.11))+
    labs(caption="",x="",y="")+
    scale_x_discrete(limits=c("WHEf","WHEc","WOf","WOc","SORf","SORc","C"))+
    geom_dotplot(binaxis = "y", stackdir = "center",stroke = 4,color="#A76E19", stackratio = 0.5,fill="transparent",dotsize=0.5)+
    stat_summary(fun=mean, geom="point",shape=3,size=7,stroke=2,color="black")+
    geom_text(data = m_2_label, aes(x = group, y = 0, label = letter), vjust=-0.4,hjust=0, size=6.5)+
    geom_hline(aes(yintercept=mean(as.numeric(exp_2$m[exp_2$group=="C"]))), linetype=2, color="#F8A50E",linewidth=1)
ggsave(filename="m_2.png",m_2, width = 15,height =6,dpi = 400)

####OLS: biochar focus####
#response variables
r_d <-cbind(c(as.numeric(exp_1$d),as.numeric(exp_2$d)))
r_m <-cbind(c(as.numeric(exp_1$m),as.numeric(exp_2$m)))
r_pmCO2 <-cbind(c(as.numeric(exp_1$pmCO2),as.numeric(exp_2$pmCO2)))
#explanatory variables
e_pH <-c(as.numeric(exp_1$pH),as.numeric(exp_2$pH))
e_pmCO2 <-c(as.numeric(exp_1$pmCO2),as.numeric(exp_2$pmCO2))
e_b <-c(as.numeric(exp_1$b),as.numeric(exp_2$b))
e_a <-c(as.numeric(exp_1$a),as.numeric(exp_2$a))
e_k <-c(as.numeric(exp_1$k),as.numeric(exp_2$k))
e_a_k <-c(as.numeric(exp_1$a_k),as.numeric(exp_2$a_k))
e_d <-c(as.numeric(exp_1$d),as.numeric(exp_2$d))
e_omCO2r <-c(as.numeric(exp_1$omCO2r),as.numeric(exp_2$omCO2r))
e_nmCO2 <-c(as.numeric(exp_1$nmCO2),as.numeric(exp_2$nmCO2))
#dummies
e_exp2 <-c(as.numeric(exp_1$exp2),as.numeric(exp_2$exp2))

e_BC <- c(as.numeric(exp_1$BC),as.numeric(exp_2$BC))
e_fine <- c(as.numeric(exp_1$fine),as.numeric(exp_2$fine))
e_coarse <- c(as.numeric(exp_1$coarse),as.numeric(exp_2$coarse))
e_fine <- c(as.numeric(exp_1$coarse),as.numeric(exp_2$coarse))
e_SOR <- c(as.numeric(exp_1$SOR),as.numeric(exp_2$SOR))
e_WO <- c(as.numeric(exp_1$WO),as.numeric(exp_2$WO))
e_WHE <- c(as.numeric(exp_1$WHE),as.numeric(exp_2$WHE))
#for diagnostic plots:
par(mfrow = c(2, 2))
#interpretation of values
#only r is log
#(exp(estimate)-1)*100
#both r and e are log ... 10 percent increase mean...
#(1.10^(estimate)-1)*100

#####d#####
#1: Does biochar in general have an effect?
biochar <- data.frame("r_d" = r_d, "e_BC"=e_BC, "e_exp2" = e_exp2)
ggplot(biochar, aes(y = r_d, x = e_BC+e_exp2)) + 
  geom_point() + theme_bw()+
  stat_smooth(method = "lm", col = "red")

d_BC <- lm(r_d~ e_BC+e_exp2)

#glm(r_d~ e_BC+e_exp2, family=binomial)

summary(d_BC)
plot(d_BC)

#D# pre-test: normality
shapiro.test(resid(d_BC)) #D# p>0.05 -> normally distributed
#heteroscedasicity from plot
#could try to transform the data
#could use weighted least squares linear regression 
#could use non-parametric regression

#2: Does one of the feedstock have an effect?
feedstock <- data.frame("r_d" = r_d, "e_SOR"=e_SOR, "e_WHE"=e_WHE, "e_exp2" = e_exp2)
ggplot(feedstock, aes(y = r_d, x = e_SOR+e_WHE+e_exp2)) + 
  geom_point() + theme_bw()+
  stat_smooth(method = "lm", col = "red")

d_feedstock <- lm(r_d~ e_SOR+e_WO+e_WHE+e_exp2)
summary(d_feedstock)
plot(d_feedstock)

#D# pre-test: normality
shapiro.test(resid(d_feedstock)) #D# p>0.05 -> normally distributed

#3: Does biochar texture have an effect among BC treatments of exp 2&3?
biochar_texture <- data.frame("r_d_onlytexture" = r_d_onlytexture, "e_coarse_onlytexture"=e_coarse_onlytexture)
ggplot(biochar_texture, aes(y= r_d_onlytexture, x=e_coarse_onlytexture)) + 
  geom_point() + theme_bw()+
  stat_smooth(method = "lm", col = "red")

r_d_onlytexture <-cbind(c(as.numeric(exp_2$d[exp_2$BC==1])))
e_coarse_onlytexture <- c(as.numeric(exp_2$coarse[exp_2$BC==1]))

d_texture <- lm(r_d_onlytexture~ e_coarse_onlytexture)
summary(d_texture)
plot(d_texture)

#as percentage:
d_texture_percent <- lm(log(r_d_onlytexture)~ e_coarse_onlytexture)
summary(d_texture_percent)
(exp(-0.16920)-1)*100

#D# pre-test: normality
shapiro.test(resid(d_texture)) #D# p<0.05 -> not normally distributed

#Can a high d be described by high pH and low k? Yes!
high_d <- data.frame("log_r_d" = log(r_d), "log_e_k"= log(e_k),
                              "e_pH"=e_pH)
ggplot(high_d, aes(y= log_r_d, x=log_e_k+e_pH)) + 
  geom_point() + theme_bw()+
  stat_smooth(method = "lm", col = "red")

d_pH_k <- lm(log(r_d)~log(e_k)+e_pH)
summary(d_pH_k)
plot(d_pH_k)

#D# pre-test: normality
shapiro.test(resid(d_pH_k)) #D# p<0.05 -> not normally distributed

#####pmCO2#####
#1: Does biochar in general have an effect?
biochar_ge <- data.frame("r_pmCO2" = r_pmCO2, "e_BC"= e_BC,
                     "e_exp2"=e_exp2)
ggplot(biochar_ge, aes(y= r_pmCO2, x= e_BC+e_exp2)) + 
  geom_point() + theme_bw()+
  stat_smooth(method = "lm", col = "red")

pmCO2_BC <- lm(r_pmCO2~ e_BC+e_exp2)
summary(pmCO2_BC)
plot(pmCO2_BC)

#as percentage:
pmCO2_BC_percent <- lm(log(r_pmCO2)~ e_BC+e_exp2)
summary(pmCO2_BC_percent)
(exp(0.06818)-1)*100

#D# pre-test: normality
shapiro.test(resid(pmCO2_BC)) #D# p<0.05 -> normally distributed

#2: Does one of the feedstock have an effect?
feedstock_one <- data.frame("r_pmCO2" = r_pmCO2, "e_SOR"= e_SOR,
                                       "e_WO"=e_WO, "e_WHE"=e_WHE, 
                            "e_exp2" = e_exp2)
ggplot(feedstock_one, aes(y= r_pmCO2, x= e_SOR+e_WO+e_WHE+e_exp2)) + 
  geom_point() + theme_bw()+
  stat_smooth(method = "lm", col = "red")

pmCO2_feedstock <- lm(r_pmCO2~ e_SOR+e_WO+e_WHE+e_exp2)
summary(pmCO2_feedstock)
plot(pmCO2_feedstock)
#D# pre-test: normality
shapiro.test(resid(pmCO2_feedstock)) #D# p >0.05 -> normally distributed

#as percentage:
pmCO2_feedstock_percent <- lm(log(r_pmCO2)~ e_SOR+e_WO+e_WHE+e_exp2)
summary(pmCO2_feedstock_percent)
(exp(0.12412)-1)*100

#3: Does biochar texture have an effect among BC treatments of exp 2&3?
biochar_texture <- data.frame("r_pmCO2_onlytexture" = r_pmCO2_onlytexture, "e_coarse_onlytexture"=e_coarse_onlytexture)
ggplot(biochar_texture, aes(y= r_pmCO2_onlytexture, x=e_coarse_onlytexture)) + 
  geom_point() + theme_bw()+
  stat_smooth(method = "lm", col = "red")

r_pmCO2_onlytexture <-cbind(c(as.numeric(exp_2$pmCO2[exp_2$BC==1])))
e_coarse_onlytexture <- c(as.numeric(exp_2$coarse[exp_2$BC==1]))

pmCO2_texture <- lm(r_pmCO2_onlytexture~ e_coarse_onlytexture)
summary(pmCO2_texture)
plot(pmCO2_texture)

#as percentage:
pmCO2_texture_percent <- lm(log(r_pmCO2_onlytexture)~ e_coarse_onlytexture)
summary(pmCO2_texture_percent)
(exp(-0.07579)-1)*100

#D# pre-test: normality
shapiro.test(resid(pmCO2_texture)) #D# p >0.05 -> normally distributed

#####m#####
#1: Does biochar in general have an effect?
biochar_ge <- data.frame("r_m" = r_m, "e_BC"= e_BC,
                         "e_exp2"=e_exp2)
ggplot(biochar_ge, aes(y= r_m, x= e_BC+e_exp2)) + 
  geom_point() + theme_bw()+
  stat_smooth(method = "lm", col = "red")

m_BC <- lm(r_m~ e_BC+e_exp2)
summary(m_BC)
plot(m_BC)

#D# pre-test: normality
shapiro.test(resid(m_BC)) #D# p <0.05 -> not normally distributed

#2: Does one of the feedstock have an effect?
feedstock_one <- data.frame("r_m" = r_m, "e_SOR"= e_SOR,
                            "e_WO"=e_WO, "e_WHE"=e_WHE, 
                            "e_exp2" = e_exp2)
ggplot(feedstock_one, aes(y= r_m, x= e_SOR+e_WO+e_WHE+e_exp2)) + 
  geom_point() + theme_bw()+
  stat_smooth(method = "lm", col = "red")

m_feedstock <- lm(r_m~ e_SOR+e_WO+e_WHE+e_exp2)
summary(m_feedstock)
plot(m_feedstock)

#D# pre-test: normality
shapiro.test(resid(m_feedstock)) #D# p <0.05 -> not normally distributed

#3: Does biochar texture have an effect among BC treatments of exp 2&3?
biochar_texture <- data.frame("r_m_onlytexture" = r_m_onlytexture, "e_coarse_onlytexture"=e_coarse_onlytexture)
ggplot(biochar_texture, aes(y= r_m_onlytexture, x=e_coarse_onlytexture)) + 
  geom_point() + theme_bw()+
  stat_smooth(method = "lm", col = "red")

r_m_onlytexture <-cbind(c(as.numeric(exp_2$m[exp_2$BC==1])))
e_coarse_onlytexture <- c(as.numeric(exp_2$coarse[exp_2$BC==1]))

m_texture <- lm(r_m_onlytexture~ e_coarse_onlytexture)
summary(m_texture)
plot(m_texture)

#as percentage:
m_texture_percent <- lm(log(r_m_onlytexture)~ e_coarse_onlytexture)
summary(m_texture_percent)
(exp(0.15883)-1)*100

#D# pre-test: normality
shapiro.test(resid(m_texture)) #D# p <0.05 -> not normally distributed

#Can the low m be described by high non-methanogenic respiration and low omCO2r? Yes!
low_m <- data.frame("r_m_log" = log(r_m), "e_omCO2r_log"= log(e_omCO2r),
                                   "e_nmCO2_log"=log(e_nmCO2), "e_pH"= e_pH)
ggplot(low_m, aes(y= r_m_log, x= e_omCO2r_log+e_nmCO2_log+e_pH)) + 
  geom_point() + theme_bw()+
  stat_smooth(method = "lm", col = "red")

m_omCO2r_nmCO2 <- lm(log(r_m)~ log(e_omCO2r)+log(e_nmCO2)+e_pH)
summary(m_omCO2r_nmCO2)
plot(m_omCO2r_nmCO2)

#D# pre-test: normality
shapiro.test(resid(m_omCO2r_nmCO2)) #D# p >0.05 -> normally distributed

####OLS: mechanistic model with higher data quality####
#set limit for data quality
#SE_limit <- 0.05
#response variable
#r_d_l <-cbind(c(as.numeric(exp_1$d[exp_1$m_SE_rel<SE_limit]),as.numeric(exp_2$d[exp_2$m_SE_rel<SE_limit]),as.numeric(exp_3$d[exp_3$m_SE_rel<SE_limit])))
#r_m_l <-cbind(c(as.numeric(exp_1$m[exp_1$m_SE_rel<SE_limit]),as.numeric(exp_2$m[exp_2$m_SE_rel<SE_limit]),as.numeric(exp_3$m[exp_3$m_SE_rel<SE_limit])))
#explanatory variables
#e_pH_l <-c(as.numeric(exp_1$pH[exp_1$m_SE_rel<SE_limit]),as.numeric(exp_2$pH[exp_2$m_SE_rel<SE_limit]),as.numeric(exp_3$pH[exp_3$m_SE_rel<SE_limit]))
#e_p_l <-c(as.numeric(exp_1$omCO2[exp_1$m_SE_rel<SE_limit]),as.numeric(exp_2$pmCO2[exp_2$m_SE_rel<SE_limit]),as.numeric(exp_3$pmCO2[exp_3$m_SE_rel<SE_limit]))
#e_b_l <-c(as.numeric(exp_1$b[exp_1$m_SE_rel<SE_limit]),as.numeric(exp_2$b[exp_2$m_SE_rel<SE_limit]),as.numeric(exp_3$b[exp_3$m_SE_rel<SE_limit]))
#e_a_l <-c(as.numeric(exp_1$a[exp_1$m_SE_rel<SE_limit]),as.numeric(exp_2$a[exp_2$m_SE_rel<SE_limit]),as.numeric(exp_3$a[exp_3$m_SE_rel<SE_limit]))
#e_k_l <-c(as.numeric(exp_1$k[exp_1$m_SE_rel<SE_limit]),as.numeric(exp_2$k[exp_2$m_SE_rel<SE_limit]),as.numeric(exp_3$k[exp_3$m_SE_rel<SE_limit]))
#e_a_k_l <-c(as.numeric(exp_1$a_k[exp_1$m_SE_rel<SE_limit]),as.numeric(exp_2$a_k[exp_2$m_SE_rel<SE_limit]),as.numeric(exp_3$a_k[exp_3$m_SE_rel<SE_limit]))
#e_d_l <-c(as.numeric(exp_1$d[exp_1$m_SE_rel<SE_limit]),as.numeric(exp_2$d[exp_2$m_SE_rel<SE_limit]),as.numeric(exp_3$d[exp_3$m_SE_rel<SE_limit]))
#e_c_l <-c(as.numeric(exp_1$omCO2r[exp_1$m_SE_rel<SE_limit]),as.numeric(exp_2$omCO2r[exp_2$m_SE_rel<SE_limit]),as.numeric(exp_3$omCO2r[exp_3$m_SE_rel<SE_limit]))

#####d#####
#option 1
#summary(lm(log(r_d_l)~log(e_a_k_l)+log(e_k_l)+e_pH_l))
#option 2
#summary(lm(log(r_d_l)~log(e_b_l)+log(e_a_k_l)+e_pH_l))

#####m#####
#summary(lm(log(r_m_l)~log(e_c_l)+log(e_k_l)+log(e_d_l)))
#summary(lm(log(r_m_l)~log(e_a_l)+log(e_p_l)+log(e_d_l)))
#good fit, but including b (possibly endogenous):
#summary(lm(log(r_m_l)~log(e_a_l)+log(e_b_l)+log(e_k_l)))
#summary(lm(log(r_m_l)~log(e_a_k_l)+log(e_b_l)+log(e_d_l)))
#summary(lm(log(r_m_l)~log(e_a_k_l)+log(e_p_l)+log(e_b_l)))
#summary(lm(log(r_m_l)~log(e_b_l)+log(e_c_l)+log(e_d_l)))
#summary(lm(log(r_m_l)~log(e_b_l)+log(e_p_l)+log(e_k_l)))
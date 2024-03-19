
master <- read.csv("C:/Users/Katrina/Box/Duke/Kat-PhD/projects/fall-winter-2020-exp/Chapter.1/00_raw_data/Master.csv", na.strings="NA",stringsAsFactors = FALSE)
head(master)

library(ggplot2)
library(gridExtra)
library(grid)
#### need to run resp_ana before using combined data 
colnames(respiration.master) <- c("Jar.number", "Respiration")
head(respiration.master)

master <- merge(master, respiration.master, by = "Jar.number")

master$Temp
master$Nutrient


master$Nutrient <- as.factor(master$Nutrient)
master$Protist <- as.factor(master$Protist)
master$Temp <- as.factor(master$Temp)

master$Nutrient <- factor(master$Nutrient, levels = c("Half", "Full"))
head(master)
########################################## 
##---------------------------------------------------------------------------------------------------------------------
## The code below does preliminary data analysis and model selection for OD600
##-------------------------------------------------------------------------------------------------------------------------------
library(MuMIn)
full.model.1 <- lm(Spec~Temp*Protist*Nutrient, na.action = "na.fail", data= master)
full.model.log <- lm(log(Spec+1)~Temp*Protist*Nutrient, na.action = "na.fail", data= master)
summary(full.model.1)
summary(full.model.log)

par(mfrow=c(2,2))
plot(full.model.log)

full.model.1 <- lm(Spec~Temp*Protist*Nutrient, na.action = "na.fail", data= master)
dd1 <- dredge(full.model.1)
dd1
plot_dd1 <- sw(dd1)
barplot(t(plot_dd1))


model.sel(dd1, rank = AIC, rank.args = alist(k = log(nobs(x))))

#### indi
temp.mod <- summary(lm(Spec~Temp, data = master))
cf <- coef(temp.mod)
Intercept <- cf[1]
Slope <- cf[2]
temp.mod 

boxplot((Spec)~Temp, data = master)

nut.mod <- summary(lm(Spec~Nutrient, data = master))
nut.mod

pro.mod <- summary(lm(Spec~Protist, data = master))
pro.mod

summary(lm(Spec~Temp+Protist+Nutrient, data= master)) #?????????????????????????

#### interactions

full.model.1 <- lm(Spec~Temp*Protist*Nutrient, na.action = "na.fail", data= master)
summary(full.model.1)


boxplot((Spec)~Temp*Nutrient*Protist, data = master)
summary(lm(Spec~Temp*Nutrient*Protist,na.action = "na.fail", data= master))  #use this


library("tidyverse")
master_new<-master %>%
  group_by(Nutrient, Temp, Protist) %>%
  summarize(mean_spec=mean(Spec), sd_spec=sd(Spec)) %>%
  mutate(lower=mean_spec-sd_spec, upper=mean_spec+sd_spec) %>%
  ungroup()

###################################
############################fig 1


finally1 <- ggplot(data=master_new, aes(Temp, mean_spec, shape=Protist, colour=Protist)) +
  geom_errorbar(aes(ymin=lower, ymax=upper), linewidth=1, width=0.5, alpha=0.5,  position = "dodge") +
  geom_point(size=3) +
  ggtitle("D: Low               E: High")+
  geom_line(aes(group=Protist), linetype=2) +
  scale_color_manual(values=c("deeppink", "purple"))+
  facet_wrap(~Nutrient, scale="free_y") +
  coord_cartesian(ylim = c(0, 0.07))+
  theme_bw()
finally.od <- finally1 + labs(x="Temperature (\u00B0C)", y="Microbial Biomass (OD600)") + theme(text = element_text(size = 25)) 
finally.od <- finally.od + theme(legend.position = c(0.08, 0.80),legend.background = element_rect(fill = "white", color = "black"))+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(-1.4, "mm"))
finally.od.1 <- finally.od + theme(strip.background = element_blank(), axis.title.y = element_blank(), strip.text.x = element_blank(), axis.text.y=element_blank(),   axis.text.x = element_text(size = 20),  axis.title.x = element_blank(),  plot.title = element_text(size = 25))  
finally.od.1


master$Temp <- as.factor(master$Temp)

t<- master %>%
  group_by(Temp) %>%
  mutate(Spec.show = as.numeric(  # so ggplot doesn't complain about alpha being discrete
    between(Spec, 
            quantile(Spec)[2] - 1.5*IQR(Spec),
            quantile(Spec)[4] + 1.5*IQR(Spec)))) %>% 
  ggplot(aes(Temp, Spec)) + 
  stat_boxplot(geom = "errorbar",
               width = 0.15) + 
  geom_boxplot(outlier.shape = NA, aes(fill=Temp), fatten =3,  show.legend=FALSE) + 
  scale_fill_manual(values=c("dodgerblue", "tomato2"))+
  geom_jitter(aes(alpha=Spec.show), show.legend=FALSE, size=3, width = 0.35, shape =1, stroke =1.05) +
  scale_alpha_continuous(range = c(0, 1)) +
  labs(x="Temperature (\u00B0C)", y="Microbial Biomass (OD600)") +
  ggtitle("A")+
  theme_bw() +
  theme( axis.title.x = element_blank())+
  theme(
    # Hide panel borders and remove grid lines
    # Hide panel borders and remove grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(-1.4, "mm"),
    axis.title.x = element_blank()) +
  theme(legend.background = element_rect(fill = "white", color = "black"))+
  theme(text = element_text(size = 25)) +
  coord_cartesian(ylim = c(0, 0.07)) %>% {.}
t.1 <- t
t.1 


n<- master %>%
  group_by(Nutrient) %>%
  mutate(Spec.show = as.numeric(  # so ggplot doesn't complain about alpha being discrete
    between(Spec, 
            quantile(Spec)[2] - 1.5*IQR(Spec),
            quantile(Spec)[4] + 1.5*IQR(Spec)))) %>% 
  ggplot(aes(Nutrient, Spec)) + 
  stat_boxplot(geom = "errorbar",
               width = 0.15) + 
  geom_boxplot(outlier.shape = NA, aes(fill= Nutrient), fatten=3,  show.legend=FALSE) + 
  geom_jitter(aes(alpha=Spec.show), show.legend=FALSE, size=3, width = 0.35, shape =1, stroke =1.05) +
  scale_alpha_continuous(range = c(0, 1)) +
  scale_fill_manual(values=c("tan4", "darkolivegreen3"))+
  labs(x="Nutrient", y="Microbial Biomass (OD600)") +
  ggtitle("B")+
  theme_bw() +
  theme(axis.title.y = element_blank(),axis.text.y=element_blank(), axis.title.x = element_blank())+
  theme(
    # Hide panel borders and remove grid lines
    panel.grid.major = element_blank(),
    axis.ticks.length = unit(-1.4, "mm"),
    panel.grid.minor = element_blank())+
  theme(text = element_text(size = 25)) +
  theme(legend.background = element_rect(fill = "white", color = "black"))+
  coord_cartesian(ylim = c(0, 0.07)) %>% {.}
n.1 <- n 
n.1 


p<- master %>%
  group_by(Protist) %>%
  mutate(Spec.show = as.numeric(  # so ggplot doesn't complain about alpha being discrete
    between(Spec, 
            quantile(Spec)[2] - 1.5*IQR(Spec),
            quantile(Spec)[4] + 1.5*IQR(Spec)))) %>% 
  ggplot(aes(Protist, Spec)) + 
  stat_boxplot(geom = "errorbar",
               width = 0.15) + 
  geom_boxplot(outlier.shape = NA, fatten =3, aes(fill = Protist), show.legend=FALSE) + 
  geom_jitter(aes(alpha=Spec.show), show.legend=FALSE, size=3, width = 0.35, shape =1, stroke =1.05) +
  scale_alpha_continuous(range = c(0, 1)) +
  scale_fill_manual(values=c("deeppink", "purple"))+
  labs(x="Protist", y="Microbial Biomass (OD600)") +
  ggtitle("C")+
  theme_bw() +
  theme( axis.title.x = element_blank())+
  theme(
    # Hide panel borders and remove grid lines
    # Hide panel borders and remove grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(-1.4, "mm"),
    axis.title.x = element_blank())+
  theme(text = element_text(size =25 )) +
  coord_cartesian(ylim = c(0, 0.07)) %>% {.}
p.1 <- p 
p.1 

library(cowplot)
figure.1.pdf <- plot_grid(t.1, n.1, p.1, finally.od.1, align = "h")
figure.1.pdf
#pdf(file= "figure.1.pdf", width=16, height=16, useDingbats=FALSE )
#plot(figure.1.pdf)
dev.off()


##-----------------------------------------------------------------------------------------------------------------------
## The code below does preliminary data analysis and model selection for Respiration
##-------------------------------------------------------------------------------------------------------------------------

library(MuMIn)
par(mfrow=c(2,2))

full.model <- lm(Respiration~Temp*Protist*Nutrient, na.action = "na.fail", data= master)
full.model.log <- lm(log(Respiration+1)~Temp*Protist*Nutrient, na.action = "na.fail", data= master)
summary(full.model)
summary(full.model.log)
dd <- dredge(full.model)
plot_dd <- sw(dd)
barplot(t(plot_dd))

plot(full.model)

model.sel(dd, rank = AIC, rank.args = alist(k = log(nobs(x))))

fmList <- get.models(dd, 1:7)
summary(model.avg(fmList))

#### indi
model.temp <- summary(lm(Respiration~Temp, data = master))
model.temp

model.pro <-summary(lm(Respiration~Protist, data = master))
model.pro

model.nut <-summary(lm(Respiration~Nutrient, data = master))
model.nut

model.allwithone <- lm(Respiration~Temp+Protist+Nutrient, data= master)
summary(model.allwithone)



#############################################################
plot((Respiration)~Temp*Nutrient*Protist, data = master) #report these

interactive <-summary(lm(Respiration)~Temp*Nutrient*Protist, data = master) #report these for two way
interactive

#########################Figure 2 

library("tidyverse")
master_new<-master %>%
  group_by(Nutrient, Temp, Protist) %>%
  summarize(mean_Respiration=mean(Respiration), sd_Respiration=sd(Respiration)) %>%
  mutate(lower=mean_Respiration-sd_Respiration, upper=mean_Respiration+sd_Respiration) %>%
  ungroup()


finally1 <- ggplot(data=master_new, aes(Temp, mean_Respiration, shape=Protist, colour=Protist)) +
  geom_errorbar(aes(ymin=lower, ymax=upper), linewidth=1, width=0.5, alpha=0.5,  position = "dodge") +
  geom_point(size=3) +
  ggtitle("D: Low               E: High")+
  geom_line(aes(group=Protist), linetype=2) +
  scale_color_manual(values=c("deeppink", "purple"))+
  facet_wrap(~Nutrient, scale="free_y") +
  coord_cartesian(ylim = c(0, 0.0012))+
  theme_bw()
finally.resp <- finally1 + labs(x="Temperature (\u00B0C)", y="Respiration (umol O[2]/min)") + theme(text = element_text(size = 25)) 
finally.resp <- finally.resp + theme(legend.position = c(0.08, 0.80),axis.text.y=element_blank(),legend.background = element_rect(fill = "white", color = "black"))+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(-1.4, "mm"))
finally.resp.1 <- finally.resp + theme(strip.background = element_blank(),axis.text.y=element_blank(), strip.text.x = element_blank(),   axis.text.x = element_text(size = 20),  axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_text(size = 25))  
finally.resp.1
finally.resp.2  <- finally.resp.1 
finally.resp.2 



master$Temp <- as.factor(master$Temp)

t.resp <- master %>%
  group_by(Temp) %>%
  mutate(Respiration.show = as.numeric(  # so ggplot doesn't complain about alpha being discrete
    between(Respiration, 
            quantile(Respiration)[2] - 1.5*IQR(Respiration),
            quantile(Respiration)[4] + 1.5*IQR(Respiration)))) %>% 
  ggplot(aes(Temp, Respiration)) + 
  stat_boxplot(geom = "errorbar",
               width = 0.15) +
  geom_boxplot(outlier.shape = NA, aes(fill=Temp), fatten =3, show.legend=FALSE) + 
  geom_jitter(aes(alpha=Respiration.show), show.legend=FALSE, shape = 1,size=3, width = 0.35,stroke =1.05) +
  scale_alpha_continuous(range = c(0, 1)) +
  scale_fill_manual(values=c("dodgerblue", "tomato2"))+
  labs(x="Temperature (\u00B0C)", y="Respiration (umol O[2]/min)") +
  ylab(expression(Respiration~(umol~o[2]/min)))+
  ggtitle("A ")+
  theme_bw() +
  theme(
    # Hide panel borders and remove grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(-1.4, "mm"),
    axis.title.x = element_blank())+
  theme(text = element_text(size = 25)) +
  coord_cartesian(ylim = c(0, 0.0012)) %>% {.}
t.resp1 <- t.resp
t.resp1 ##p value needs to be fixed


Protist.resp <- master %>%
  group_by(Protist) %>%
  mutate(Respiration.show = as.numeric(  # so ggplot doesn't complain about alpha being discrete
    between(Respiration, 
            quantile(Respiration)[2] - 1.5*IQR(Respiration),
            quantile(Respiration)[4] + 1.5*IQR(Respiration)))) %>% 
  ggplot(aes(Protist, Respiration)) + 
  stat_boxplot(geom = "errorbar",
               width = 0.15) +
  geom_boxplot(outlier.shape = NA, aes(fill = Protist), show.legend = FALSE, fatten = 3) + 
  geom_jitter(aes(alpha=Respiration.show), show.legend=FALSE, size=3, shape =1, stroke = 1.05, width = 0.35,) +
  scale_alpha_continuous(range = c(0, 1)) +
  scale_fill_manual(values=c("deeppink", "purple"))+
  labs(x="Protist", y="Respiration (umol O[2]/min)") +
  ylab(expression(Respiration~(umol~o[2]/min)))+
  ggtitle("C")+
  theme_bw() +
  theme(
    # Hide panel borders and remove grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(-1.4, "mm"),
    axis.title.x = element_blank())+
  theme(text = element_text(size = 25)) +
  coord_cartesian(ylim = c(0, 0.0012)) %>% {.}
p.resp1 <- Protist.resp  
p.resp1 ##p value needs to be fixed

Nutrient.resp <- master %>%
  group_by(Nutrient) %>%
  mutate(Respiration.show = as.numeric(  # so ggplot doesn't complain about alpha being discrete
    between(Respiration, 
            quantile(Respiration)[2] - 1.5*IQR(Respiration),
            quantile(Respiration)[4] + 1.5*IQR(Respiration)))) %>% 
  ggplot(aes(Nutrient, Respiration)) + 
  stat_boxplot(geom = "errorbar",
               width = 0.15) + 
  geom_boxplot(outlier.shape = NA, aes(fill =Nutrient), fatten = 3, show.legend=FALSE ) + 
  geom_jitter(aes(alpha=Respiration.show), show.legend=FALSE, size=3, shape =1, stroke= 1.05, width = 0.35,) +
  scale_alpha_continuous(range = c(0, 1)) +
  scale_fill_manual(values=c("tan4", "darkolivegreen3"))+
  labs(x="Nutrient", y="Respiration (umol O[2]/min)") +
  ylab(expression(Respiration~(umol~o[2]/min)))+
  ggtitle("B")+
  theme_bw() +
  theme(axis.title.y = element_blank(),axis.text.y=element_blank(), axis.title.x = element_blank())+
  theme(
    # Hide panel borders and remove grid lines
    panel.grid.major = element_blank(),
    axis.ticks.length = unit(-1.4, "mm"),
    panel.grid.minor = element_blank())+
  theme(text = element_text(size = 25)) +
  coord_cartesian(ylim = c(0, 0.0012)) %>% {.}
n.resp1 <- Nutrient.resp
n.resp1 ##p value needs to be fixed

grid.arrange(arrangeGrob(t.resp1,n.resp1,p.resp1, finally.resp.2, ncol=2), widths=1)

library(cowplot)
figure.2 <- plot_grid(t.resp1, n.resp1, p.resp1, finally.resp.2, align = "h")



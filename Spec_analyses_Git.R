
master <- read.csv("C:/Users/Katrina/Box/Duke/Kat-PhD/projects/fall-winter-2020-exp/Chapter.1/00_raw_data/Master.csv", na.strings="NA",stringsAsFactors = FALSE)
master <- read.csv("C:/Users/18568/Box/Duke/Kat-PhD/projects/fall-winter-2020-exp/Chapter.1/00_raw_data/Master.csv", na.strings="NA",stringsAsFactors = FALSE)
master <- read.csv("C:/Users/kd239/Box/Duke/Kat-PhD/projects/fall-winter-2020-exp/Chapter.1/00_raw_data/Master.csv", na.strings="NA",stringsAsFactors = FALSE)

head(master)

library(ggplot2)
library(gridExtra)
library(grid)
#### need to run resp_ana and sem before using combined data 
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

# Calculate Z-scores
sem.new$Spec1_z <- scale(sem.new$Spec1)
sem.new$biomass.prot1_z <- scale(sem.new$biomass.prot1)

lm_model <- lm(Spec1_z ~ biomass.prot1_z, data = sem.new)

# Summary of the model
summary(lm_model)

sem.new$lm_pred <- predict(lm_model)

# Plot with both 1:1 line and regression prediction
ggplot(sem.new, aes(x = biomass.prot1_z, y = Spec1_z)) +
  geom_point(alpha = 0.7) +
  geom_line(aes(y = lm_pred), color = "blue", linetype = "solid", linewidth = 1.2) + # Regression line
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "red") + # 1:1 line
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  labs(
    x = "Biomass.prot1 Z-score",
    y = "Spec1 Z-score",
    title = "Linear Model Prediction vs. 1:1 Line"
  ) +
  theme_minimal()



# Identify potential outliers using Cook's distance
cooksd <- cooks.distance(lm_model)

# Plot Cook's distance to visualize potential outliers
plot(cooksd, pch="*", cex=2, main="Cook's Distance for Outlier Detection")
abline(h = 4/length(cooksd), col="red") # Cutoff threshold

# Remove outlier (e.g., if Cook's distance is above threshold)
outlier_threshold <- 4/length(cooksd)
sem.new_clean <- sem.new[cooksd < outlier_threshold, ]

# Rerun the model without the outlier
lm_model_clean <- lm(Spec1_z ~ biomass.prot1_z, data = sem.new_clean)

summary(lm_model_clean)

# Plot again without the outlier
sem.new_clean$lm_pred_clean <- predict(lm_model_clean)

ggplot(sem.new_clean, aes(x = biomass.prot1_z, y = Spec1_z)) +
  geom_point(alpha = 0.7) +
  geom_line(aes(y = lm_pred_clean), color = "blue", linetype = "solid", linewidth = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  labs(
    x = "Biomass.prot1 Z-score",
    y = "Spec1 Z-score",
    title = "Linear Model Without Outliers"
  ) +
  theme_minimal()





###############raw fig 1 ana
# Fit model using raw (non-standardized) values
lm_raw <- lm(Spec ~ biomass.pro1, data = sem.new)

plot((Spec)~biomass.pro1, data = sem.new)
abline(lm_raw)
# Show summary
summary(lm_raw)

master <- read.csv("C:/Users/kd239/Box/Duke/Kat-PhD/projects/fall-winter-2020-exp/Chapter.1/00_raw_data/Master.csv", na.strings="NA",stringsAsFactors = FALSE)
master <- read.csv("C:/Users/18568/Box/Duke/Kat-PhD/projects/fall-winter-2020-exp/Chapter.1/00_raw_data/Master.csv", na.strings="NA",stringsAsFactors = FALSE)



prot.biomassdata <- sem.new %>% select(Jar.number, biomass.pro1)

master <- merge(master, prot.biomassdata, by = "Jar.number", all.x = TRUE)
master[is.na(master)] <- 0
master$biomass.pro1[master$Jar.number == 34] <- NA

class(master$biomass.pro1)


lm_master <- lm(Spec ~ biomass.pro1, data = master, na.action = na.exclude)
plot(Spec ~ biomass.pro1, data = master)
abline(lm_master)
summary(lm_master)

master$OD600_controlled <- resid(lm_master)

# Step 2: Fit full model on unfiltered data
full.model.1 <- lm(OD600_controlled ~ Temp * Protist * Nutrient,
                   data = master,
                   na.action = na.exclude)
summary(full.model.1)

full.model.pro <- lm(OD600_controlled ~ Protist, data = master, na.action = na.exclude)
summary(full.model.pro)

full.model.temp <- lm(OD600_controlled ~ Temp, data = master, na.action = na.exclude)
summary(full.model.temp)

full.model.nut <- lm(OD600_controlled ~ Nutrient, data = master, na.action = na.exclude)
summary(full.model.nut)

# Ensure factors are set correctly
master$Nutrient <- factor(master$Nutrient, levels = c("Half", "Full"))
master$Temp <- as.factor(master$Temp)


library(tidyverse)
library(cowplot)
# Step 2: Summarized data (mean ± SD) by group (if needed for future error bars)
master_new <- master %>%
  group_by(Nutrient, Temp, Protist) %>%
  summarize(
    mean_spec = mean(OD600_controlled, na.rm = TRUE),
    sd_spec = sd(OD600_controlled, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    lower = mean_spec - sd_spec,
    upper = mean_spec + sd_spec
  )

master$Nutrient <- recode_factor(master$Nutrient,
                                 "Half" = "Low",
                                 "Full" = "High")
# Step 3: Faceted Temp x Protist Plot by Nutrient

# Create ggplot with the second graph's data and first graph's styling
finally1 <- master %>%
  filter(!is.na(OD600_controlled)) %>%
  ggplot(aes(x = Temp, y = OD600_controlled, colour = Protist)) +
  scale_color_manual(values = c("deeppink", "purple")) +
  stat_boxplot(geom = 'errorbar') + 
  geom_boxplot(outlier.shape = NA, outlier.colour = NULL) +  # Match first graph style (hide outliers)
  stat_summary(
    fun = mean,
    geom = "crossbar",
    aes(group = Protist),
    width = 0.5,
    color = "black",
    fatten = 2,
    position = position_dodge(width = 0.75)
  ) +
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.5,
      jitter.height = 0
    )
  ) +
  coord_cartesian(ylim = c(-0.02, 0.06))  +
  ggtitle("B: Low Nutrient               C: High Nutrient") +
  facet_wrap(~Nutrient, scale = "free_y") +
  theme_bw()

# Styling as in the original first graph
finally.od <- finally1 +
  labs(x = "Temperature (\u00B0C)", y = "Prokaryotic Biomass (OD600-residuals)") +
  theme(text = element_text(size = 25)) +
  theme(
    legend.position = c(0.08, 0.80),
    legend.background = element_rect(fill = "white", color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(-1.4, "mm")
  )

# Additional cleanup (optional, as in first graph’s final version)
finally.od.1 <- finally.od +
  theme(
    strip.background = element_blank(),
    axis.title.y = element_blank(),
    strip.text.x = element_blank(),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    plot.title = element_text(size = 25)
  )

# Display plot
finally.od.1
t.1 <- master %>%
  filter(!is.na(OD600_controlled)) %>%
  group_by(Temp) %>%
  mutate(Spec.show = as.numeric(OD600_controlled)) %>%
  ggplot(aes(Temp, OD600_controlled)) + 
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot(outlier.shape = NA, aes(fill = Temp), fatten = 3, show.legend = FALSE) + 
  scale_fill_manual(values = c("dodgerblue", "tomato2")) +
  geom_jitter(aes(alpha = 1), show.legend = FALSE, size = 3, width = 0.35, shape = 1, stroke = 1.05) +
  scale_alpha_continuous(range = c(0, 1)) +
  labs(x = "Temperature (\u00B0C)", y = "Prokaryotic Biomass (OD600-residuals)") +
  ggtitle("A: Temperature") +
  theme_bw() +
  theme(
    text = element_text(size = 25),
    axis.title.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(-1.4, "mm"),
    legend.background = element_rect(fill = "white", color = "black")
  ) +
  coord_cartesian(ylim = c(-0.02, 0.06))

n.1 <- master %>%
  filter(!is.na(OD600_controlled)) %>%
  group_by(Nutrient) %>%
  mutate(Spec.show = as.numeric(OD600_controlled)) %>%
  ggplot(aes(Nutrient, OD600_controlled)) + 
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot(outlier.shape = NA, aes(fill = Nutrient), fatten = 3, show.legend = FALSE) + 
  scale_fill_manual(values = c("tan4", "darkolivegreen3")) +
  geom_jitter(aes(alpha = 1), size = 3, width = 0.35, shape = 1, stroke = 1.05, show.legend = FALSE) +
  scale_alpha_continuous(range = c(0, 1)) +
  labs(x = "Nutrient", y = "Prokaryotic Biomass (OD600-residuals)") +
  ggtitle("B: Nutrient") +
  coord_cartesian(ylim = c(-0.02, 0.06)) +
  theme_bw() +
  theme(
    text = element_text(size = 25),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(-1.4, "mm"),
    legend.background = element_rect(fill = "white", color = "black")
  )


# Protist Panel


p.1 <- master %>%
  filter(!is.na(OD600_controlled)) %>%
  group_by(Protist) %>%
  mutate(Spec.show = as.numeric(OD600_controlled)) %>%
  ggplot(aes(Protist, OD600_controlled)) + 
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot(outlier.shape = NA, aes(fill = Protist), fatten = 3, show.legend = FALSE) + 
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, fatten = 2, color = "black") +
  geom_jitter(aes(alpha = 1), size = 3, width = 0.35, shape = 1, stroke = 1.05, show.legend = FALSE) +
  scale_alpha_continuous(range = c(0, 1)) +
  scale_fill_manual(values = c("deeppink", "purple")) +
  labs(x = "Protist", y = "Prokaryotic Biomass (OD600-residuals)") +
  ggtitle("A: Ciliate Presence") +
  coord_cartesian(ylim = c(-0.02, 0.06)) +
  theme_bw() +
  theme(
    text = element_text(size = 25),
    axis.title.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(-1.4, "mm"),
    legend.background = element_rect(fill = "white", color = "black")
  )

# Export combined figure (Protist panel + Fig 1 combo)
library(cowplot)

pdf(file = "C:/Users/kd239/Box/Duke/Fig1.pdf", width = 20, height = 10)

figure.1.pdf <- plot_grid(p.1, finally.od.1, align = "vh", axis = "bt")
figure.1.pdf
dev.off()


app.temp.nut.od.pdf <- plot_grid(t.1, n.1, align = "vh")
app.temp.nut.od.pdf

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

interactive <- lm(Respiration ~ Temp * Nutrient * Protist, data = master)
summary (interactive)

#########################Figure 2 

library("tidyverse")
master_new<-master %>%
  group_by(Nutrient, Temp, Protist) %>%
  summarize(mean_Respiration=mean(Respiration), sd_Respiration=sd(Respiration)) %>%
  mutate(lower=mean_Respiration-sd_Respiration, upper=mean_Respiration+sd_Respiration) %>%
  ungroup()


finally1 <- ggplot(data=master, aes(Temp, Respiration, colour=Protist))  +
  scale_color_manual(values=c("deeppink", "purple"))+
  stat_boxplot(geom ='errorbar') + 
  geom_boxplot(outlier.shape = NA, outlier.colour = "white")+
  geom_point(position=position_jitterdodge(jitter.width = 0.5, jitter.height = 0))+
  ggtitle("B: Low Nutrient               C: High Nutrient")+
  facet_wrap(~Nutrient, scale="free_y") +
  coord_cartesian(ylim = c(0, 0.0020))+
  theme_bw()
finally.resp <- finally1 + labs(x="Temperature (\u00B0C)", y="Total Community Respiration (umol O[2]/min)") + theme(text = element_text(size = 25)) 
finally.resp <- finally.resp + theme(legend.position = c(0.08, 0.80),axis.text.y=element_blank(),legend.background = element_rect(fill = "white", color = "black"))+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(-1.4, "mm"))
finally.resp.1 <- finally.resp + theme(strip.background = element_blank(),axis.text.y=element_blank(), strip.text.x = element_blank(), axis.title.y = element_blank(),  axis.text.x = element_text(size = 20),  axis.title.x = element_blank(), plot.title = element_text(size = 25))  
finally.resp.1
finally.resp.2  <- finally.resp.1 
finally.resp.2 


master$Temp <- as.factor(master$Temp)




t.resp <- master %>%
  group_by(Temp) %>%
  mutate(Respiration.show = as.numeric(Respiration)) %>%  # so ggplot doesn't complain about alpha being discrete
    #between(Respiration, 
     #       quantile(Respiration)[2] - 1.5*IQR(Respiration),
      #      quantile(Respiration)[4] + 1.5*IQR(Respiration)))) %>% 
  ggplot(aes(Temp, Respiration)) + 
  stat_boxplot(geom = "errorbar",
               width = 0.15) +
  geom_boxplot(outlier.shape = NA, aes(fill=Temp), fatten =3, show.legend=FALSE) + 
  geom_jitter(aes(alpha=1), show.legend=FALSE, shape = 1,size=3, width = 0.35,stroke =1.05) +
  scale_alpha_continuous(range = c(0, 1)) +
  scale_fill_manual(values=c("dodgerblue", "tomato2"))+
  labs(x="Temperature (\u00B0C)", y="Total Community Respiration (umol O[2]/min)") +
  ylab(expression(Total~Community~Respiration~(umol~o[2]/min)))+
  ggtitle("A: Temperature ")+
  theme_bw() +
  theme(
    # Hide panel borders and remove grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(-1.4, "mm"),
    axis.title.x = element_blank())+
  theme(text = element_text(size = 25)) +
  coord_cartesian(ylim = c(0, 0.0020)) %>% {.}
t.resp1 <- t.resp
t.resp1 ##p value needs to be fixed


Protist.resp <- master %>%
  group_by(Protist) %>%
  mutate(Respiration.show = as.numeric(Respiration))  %>%  # so ggplot doesn't complain about alpha being discrete
   # between(Respiration, 
    #        quantile(Respiration)[2] - 1.5*IQR(Respiration),
     #       quantile(Respiration)[4] + 1.5*IQR(Respiration)))) %>% 
  ggplot(aes(Protist, Respiration)) + 
  stat_boxplot(geom = "errorbar",
               width = 0.15) +
  geom_boxplot(outlier.shape = NA, aes(fill = Protist), show.legend = FALSE, fatten = 3) + 
  geom_jitter(aes(alpha=1), show.legend=FALSE, size=3, shape =1, stroke = 1.05, width = 0.35,) +
  scale_alpha_continuous(range = c(0, 1)) +
  scale_fill_manual(values=c("deeppink", "purple"))+
  labs(x="Protist", y="Respiration (umol O[2]/min)") +
  ylab(expression(Total~Community~Respiration~(umol~o[2]/min)))+
  ggtitle("A: Ciliate Presence")+
  theme_bw() +
  theme(
    # Hide panel borders and remove grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(-1.4, "mm"),
    axis.title.x = element_blank())+
  theme(text = element_text(size = 25)) +
  coord_cartesian(ylim = c(0, 0.0020)) %>% {.}
p.resp1 <- Protist.resp  
p.resp1 ##p value needs to be fixed

Nutrient.resp <- master %>%
  group_by(Nutrient) %>%
  mutate(Respiration.show = as.numeric(Respiration)) %>%  # so ggplot doesn't complain about alpha being discrete
 #   between(Respiration, 
  #          quantile(Respiration)[2] - 1.5*IQR(Respiration),
   #         quantile(Respiration)[4] + 1.5*IQR(Respiration)))) %>% 
  ggplot(aes(Nutrient, Respiration)) + 
  stat_boxplot(geom = "errorbar",
               width = 0.15) + 
  geom_boxplot(outlier.shape = NA, aes(fill =Nutrient), fatten = 3, show.legend=FALSE ) + 
  geom_jitter(aes(alpha=1), show.legend=FALSE, size=3, shape =1, stroke= 1.05, width = 0.35,) +
  scale_alpha_continuous(range = c(0, 1)) +
  scale_fill_manual(values=c("tan4", "darkolivegreen3"))+
  labs(x="Nutrient", y="Respiration (umol O[2]/min)") +
  ylab(expression(Respiration~(umol~o[2]/min)))+
  ggtitle("B: Nutrient")+
  theme_bw() +
  theme(axis.title.y = element_blank(),axis.text.y=element_blank(), axis.title.x = element_blank())+
  theme(
    # Hide panel borders and remove grid lines
    panel.grid.major = element_blank(),
    axis.ticks.length = unit(-1.4, "mm"),
    panel.grid.minor = element_blank())+
  theme(text = element_text(size = 25)) +
  coord_cartesian(ylim = c(0, 0.0020)) %>% {.}
n.resp1 <- Nutrient.resp
n.resp1 ##p value needs to be fixed

grid.arrange(arrangeGrob(t.resp1,n.resp1,p.resp1, finally.resp.2, ncol=2), widths=1)



library(cowplot)

pdf(file = "C:/Users/Katrina/Box/Duke/fig2.pdf",   # The directory you want to save the file in
    width = 20, # The width of the plot in inches
    height = 10) # The height of the plot in inches

figure.2.pdf<- plot_grid(p.resp1, finally.resp.2, align = "vh", axis = "bt")
figure.2.pdf
#pdf(file= "figure.2.pdf", width=16, height=16, useDingbats=FALSE )
#plot(figure.2.pdf)
dev.off()


pdf(file = "C:/Users/Katrina/Box/Duke/resp.pdf",   # The directory you want to save the file in
    width = 20, # The width of the plot in inches
    height = 10) # The height of the plot in inches


appendix.nut.temp.resp.pdf <- plot_grid(t.resp1, n.resp1, align = "vh")
appendix.nut.temp.resp.pdf

dev.off()


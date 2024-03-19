rm(list=ls())

# Load packages, set working directory

library("readxl")
library("stringr")
library("dplyr")

## DAY 1 OF RESP DATA
setwd("C:/Users/Katrina/Box/Duke/Kat-PhD/projects")
setwd("C:/Users/Katrina/Box/Duke/Kat-PhD/projects/fall-winter-2020-exp/Chapter.1/00_raw_data/resp_data/12.14.2020")


##--------------------------------------------------------------
# Load and prep respirometry RAW data (found in https://github.com/kmd304/Predation-Mediates-TempXNut/tree/main/Raw%20data/Resp_data)
##--------------------------------------------------------------

## These loops load and prep the data, and produce Figs S1, S2, s3, and s4 of the Appendix.
# Data needs to be separated in the three folders 12.14.2020, 12.15.2020, and 12.16.2020 for this code to run

dir()

names <- dir()
respirometry <- rep(0,(length(names)-12))
##

non_control <- grep(pattern="[a-z]\\.([1-9]|1[0-2])\\.[[:digit:]].",names)
control <- setdiff(seq(1,length(names),1),non_control)
respiration <- rep(c(0,0,0,0),12)
jar <- rep(c(0,0,0,0),12)

par(mar = c(1,1,1,1))
par(mfrow=c(8,6))

for(i in 1:12){
  # The next three lines of code run through all 1, where 1 is an entire run in day 1. 
  # An entire run encompasses three experimental jrs and one control.
  data1 <- read_excel(names[(4*i-1)], sheet=6, col_names=TRUE) 
  data2 <- read_excel(names[(4*i-2)], sheet=6, col_names=TRUE)
  data3 <- read_excel(names[(4*i-3)], sheet=6, col_names=TRUE)
  
  ##### Extracted number of jar here
  test1<-sub(".xlsx","",sub("control.","",names[(4*i-1)])) ###  first get rid of control and .xlsx
  test1<- regmatches(test1,regexpr("\\.[^\\.]*$",test1))     ### gets rid of number. to find the jar name 
  jar[(4*i-1)] <- sub("\\.","",test1)                           ### gets rid of .
  
  test2<-sub(".xlsx","",sub("control.","",names[(4*i-2)])) ###  first get rid of control and .xlsx
  test2<- regmatches(test2,regexpr("\\.[^\\.]*$",test2))     ### gets rid of number. to find the jar name 
  jar[(4*i-2)] <- sub("\\.","",test2)                           ### gets rid of .
  
  
  test3<-sub(".xlsx","",sub("control.","",names[(4*i-3)])) ###  first get rid of control and .xlsx
  test3<- regmatches(test3,regexpr("\\.[^\\.]*$",test3))     ### gets rid of number. to find the jar name 
  jar[(4*i-3)] <- sub("\\.","",test3)                           ### gets rid of .
  
 # This is where we store the control file associated with the three jars saved as data1-data3
  data_control1 <- read_excel(names[4*i], sheet=6, col_names=TRUE)
 
  # Prep data by removing first 10 mins of run and using controls to calculate slope
  data1 <- data1 %>%
    # Slopes need to be calculated relative to control and air partial O2 pressure following Tomlinson et al 2018
    mutate(controlled_02 = data1[,7]*0.21*(1/data_control1[,7]), Time=data1[,3]*1) %>%
    filter(Time > 10)
  
  data2 <- data2 %>%
    # Slopes need to be calculated relative to control and air partial O2 pressure following Tomlinson et al 2018
    mutate(controlled_02 = data2[,7]*0.21*(1/data_control1[,7]), Time=data2[,3]*1) %>%
    filter(Time > 10)
  
  data3 <- data3 %>%
    # Slopes need to be calculated relative to control and air partial O2 pressure following Tomlinson et al 2018
    mutate(controlled_02 = data3[,7]*0.21*(1/data_control1[,7]), Time=data3[,3]*1) %>%
    filter(Time > 10)
  
  respiration[(4*i-1)] <- abs(summary(lm(t(t(data1$controlled_02))~t(t(data1$Time))))$coefficients[2])
  respiration[(4*i-2)] <- abs(summary(lm(t(t(data2$controlled_02))~t(t(data2$Time))))$coefficients[2])
  respiration[(4*i-3)] <- abs(summary(lm(t(t(data3$controlled_02))~t(t(data3$Time))))$coefficients[2])
  
  plot(t(t(data1$controlled_02))~t(t(data1$Time)), ylab="Oxygen concentraion (micro mol/L)", xlab="Time (minutes)")
  abline(lm(t(t(data1$controlled_02))~t(t(data1$Time))), col="red", lwd=2)	
  
 plot(t(t(data2$controlled_02))~t(t(data2$Time)), ylab="Oxygen concentraion (micro mol/L)", xlab="Time (minutes)")
 abline(lm(t(t(data2$controlled_02))~t(t(data2$Time))), col="red", lwd=2)	
  
  plot(t(t(data3$controlled_02))~t(t(data3$Time)), ylab="Oxygen concentraion (micro mol/L)", xlab="Time (minutes)")
  abline(lm(t(t(data3$controlled_02))~t(t(data3$Time))), col="red", lwd=2)
}  

jar <- jar[-which(jar == "0")]

################################################### 
respiration <- respiration[-which(respiration == 0)]   # Replace 0 with NA

new.data <- data.frame(jar = jar ,respiration = respiration)
                        


####---------------------------------------------------------------------------------------
## DAY 2 OF RESP DATA
setwd("C:/Users/Katrina/Box/Duke/Kat-PhD/projects/fall-winter-2020-exp/Chapter.1/00_raw_data/resp_data/12.15.2020")

names2 <- dir()

non_control2 <- grep(pattern="[a-z]\\.([1-9]|1[0-2])\\.[[:digit:]].",names2)
control2 <- setdiff(seq(1,length(names2),1),non_control2)
respiration2 <- rep(c(0,0,0,0),15)
jar2 <- rep(c(0,0,0,0),15)

par(mar = c(1,1,1,1))
par(mfrow=c(8,6))

for(i in 1:15){
  data1 <- read_excel(names2[(4*i-1)], sheet=6, col_names=TRUE) 
  data2 <- read_excel(names2[(4*i-2)], sheet=6, col_names=TRUE)
  data3 <- read_excel(names2[(4*i-3)], sheet=6, col_names=TRUE)
  
  ##### Extracted number of jar here
  test1<-sub(".xlsx","",sub("control.","",names2[(4*i-1)])) ###  first get rid of control and .xlsx
  test1<- regmatches(test1,regexpr("\\.[^\\.]*$",test1))     ### gets rid of number. to find the jar name 
  jar2[(4*i-1)] <- sub("\\.","",test1)                           ### gets rid of .
  
  test2<-sub(".xlsx","",sub("control.","",names2[(4*i-2)])) ###  first get rid of control and .xlsx
  test2<- regmatches(test2,regexpr("\\.[^\\.]*$",test2))     ### gets rid of number. to find the jar name 
  jar2[(4*i-2)] <- sub("\\.","",test2)                           ### gets rid of .
  
  
  test3<-sub(".xlsx","",sub("control.","",names2[(4*i-3)])) ###  first get rid of control and .xlsx
  test3<- regmatches(test3,regexpr("\\.[^\\.]*$",test3))     ### gets rid of number. to find the jar name 
  jar2[(4*i-3)] <- sub("\\.","",test3)                           ### gets rid of .
  
  
  
  # This is where we store the control file associated with the three jars saved as data1-data3
  data_control1 <- read_excel(names2[4*i], sheet=6, col_names=TRUE)
  
  # Prep data by removing first 10 mins of run and using controls to calculate slope
  data1 <- data1 %>%
    # Slopes need to be calculated relative to control and air partial O2 pressure following Tomlinson et al 2018
    mutate(controlled_02 = data1[,7]*0.21*(1/data_control1[,7]), Time=data1[,3]*1) %>%
    filter(Time > 10)
  
  data2 <- data2 %>%
    # Slopes need to be calculated relative to control and air partial O2 pressure following Tomlinson et al 2018
    mutate(controlled_02 = data2[,7]*0.21*(1/data_control1[,7]), Time=data2[,3]*1) %>%
    filter(Time > 10)
  
  data3 <- data3 %>%
    # Slopes need to be calculated relative to control and air partial O2 pressure following Tomlinson et al 2018
    mutate(controlled_02 = data3[,7]*0.21*(1/data_control1[,7]), Time=data3[,3]*1) %>%
    filter(Time > 10)
  
  respiration2[(4*i-1)] <- abs(summary(lm(t(t(data1$controlled_02))~t(t(data1$Time))))$coefficients[2])
  respiration2[(4*i-2)] <- abs(summary(lm(t(t(data2$controlled_02))~t(t(data2$Time))))$coefficients[2])
  respiration2[(4*i-3)] <- abs(summary(lm(t(t(data3$controlled_02))~t(t(data3$Time))))$coefficients[2])
  
  plot(t(t(data1$controlled_02))~t(t(data1$Time)))
  abline(lm(t(t(data1$controlled_02))~t(t(data1$Time))), col="red", lwd=2)	
  
  plot(t(t(data2$controlled_02))~t(t(data2$Time)), ylab="", xlab="")
  abline(lm(t(t(data2$controlled_02))~t(t(data2$Time))), col="red", lwd=2)	
  
  plot(t(t(data3$controlled_02))~t(t(data3$Time)), ylab="", xlab="")
  abline(lm(t(t(data3$controlled_02))~t(t(data3$Time))), col="red", lwd=2)
}  


jar2 <- jar2[-which(jar2 == "0")]

################################################### 

respiration2 <- respiration2[-which(respiration2 == 0)]   # Replace 0 with NA
new.data2<- data.frame(jar = jar2 ,respiration = respiration2)

##---------------------------------------------------------------------------------------
## DAY 3 OF RESP DATA

setwd("C:/Users/Katrina/Box/Duke/Kat-PhD/projects/fall-winter-2020-exp/Chapter.1/00_raw_data/resp_data/12.16.2020")

names3 <- dir()


non_control3 <- grep(pattern="[a-z]\\.([1-9]|1[0-2])\\.[[:digit:]].",names3)
control3 <- setdiff(seq(1,length(names3),1),non_control3)
respiration3 <- rep(c(0,0,0,0),13)
jar3 <- rep(c(0,0,0,0),13)

par(mar = c(1,1,1,1))
par(mfrow=c(8,6))

for(i in 1:13){
  data1 <- read_excel(names3[(4*i-1)], sheet=6, col_names=TRUE) 
  data2 <- read_excel(names3[(4*i-2)], sheet=6, col_names=TRUE)
  data3 <- read_excel(names3[(4*i-3)], sheet=6, col_names=TRUE)
  
  
  ##### Extracted number of jar here
  test1<-sub(".xlsx","",sub("control.","",names3[(4*i-1)])) ###  first get rid of control and .xlsx
  test1<- regmatches(test1,regexpr("\\.[^\\.]*$",test1))     ### gets rid of number. to find the jar name 
  jar3[(4*i-1)] <- sub("\\.","",test1)                           ### gets rid of .
  
  test2<-sub(".xlsx","",sub("control.","",names3[(4*i-2)])) ###  first get rid of control and .xlsx
  test2<- regmatches(test2,regexpr("\\.[^\\.]*$",test2))     ### gets rid of number. to find the jar name 
  jar3[(4*i-2)] <- sub("\\.","",test2)                           ### gets rid of .
  
  
  test3<-sub(".xlsx","",sub("control.","",names3[(4*i-3)])) ###  first get rid of control and .xlsx
  test3<- regmatches(test3,regexpr("\\.[^\\.]*$",test3))     ### gets rid of number. to find the jar name 
  jar3[(4*i-3)] <- sub("\\.","",test3)                           ### gets rid of .
  
  # This is where we store the control file associated with the three jars saved as data1-data3
  data_control1 <- read_excel(names3[4*i], sheet=6, col_names=TRUE)
  
  # Prep data by removing first 10 mins of run and using controls to calculate slope
  data1 <- data1 %>%
    # Slopes need to be calculated relative to control and air partial O2 pressure following Tomlinson et al 2018
    mutate(controlled_02 = data1[,7]*0.21*(1/data_control1[,7]), Time=data1[,3]*1) %>%
    filter(Time > 10)
  
  data2 <- data2 %>%
    # Slopes need to be calculated relative to control and air partial O2 pressure following Tomlinson et al 2018
    mutate(controlled_02 = data2[,7]*0.21*(1/data_control1[,7]), Time=data2[,3]*1) %>%
    filter(Time > 10)
  
  data3 <- data3 %>%
    # Slopes need to be calculated relative to control and air partial O2 pressure following Tomlinson et al 2018
    mutate(controlled_02 = data3[,7]*0.21*(1/data_control1[,7]), Time=data3[,3]*1) %>%
    filter(Time > 10)
  
  respiration3[(4*i-1)] <- abs(summary(lm(t(t(data1$controlled_02))~t(t(data1$Time))))$coefficients[2])
  respiration3[(4*i-2)] <- abs(summary(lm(t(t(data2$controlled_02))~t(t(data2$Time))))$coefficients[2])
  respiration3[(4*i-3)] <- abs(summary(lm(t(t(data3$controlled_02))~t(t(data3$Time))))$coefficients[2])
  
  plot(t(t(data1$controlled_02))~t(t(data1$Time)))
    abline(lm(t(t(data1$controlled_02))~t(t(data1$Time))), col="red", lwd=2)	
  
   plot(t(t(data2$controlled_02))~t(t(data2$Time)), ylab="", xlab="")
  abline(lm(t(t(data2$controlled_02))~t(t(data2$Time))), col="red", lwd=2)	
  
    plot(t(t(data3$controlled_02))~t(t(data3$Time)), ylab="", xlab="")
   abline(lm(t(t(data3$controlled_02))~t(t(data3$Time))), col="red", lwd=2)
}  


jar3 <- jar3[-which(jar3 == "0")]

################################################### 
respiration3 <- respiration3[-which(respiration3 == 0)]   # Replace 0 with NA
new.data3<- data.frame(jar = jar3 ,respiration = respiration3)
respiration.master <- rbind(new.data, new.data2, new.data3)


dev.off()
# Data curating

setwd("C:/Users/Katrina/Box/Duke/Kat-PhD/projects/fall-winter-2020-exp/Chapter.1/00_raw_data/resp_data/12.14.2020")
curation_1 <- read_excel('control.1.72.xlsx', sheet=6, col_names=TRUE)
control_1 <- read_excel('control.1.xlsx', sheet=6, col_names=TRUE)
curation_1 <- curation_1 %>%
  # Slopes need to be calculated relative to control and air partial O2 pressure following Tomlinson et al 2018
  mutate(controlled_02 = curation_1[,7]*0.21*(1/control_1[,7]), Time=curation_1[,3]*1) %>%
  filter(Time > 15)
curation_1_resp <- abs(summary(lm(t(t(curation_1$controlled_02))~t(t(curation_1$Time))))$coefficients[2])

par(mar = c(1,1,1,1))
par(mfrow=c(3,3))

plot(t(t(curation_1$controlled_02))~t(t(curation_1$Time)), ylab="", xlab="")
abline(lm(t(t(curation_1$controlled_02))~t(t(curation_1$Time))), col="red", lwd=2)


respiration.master$respiration[3] <- curation_1_resp # Load data
names

respiration.master$respiration[which(respiration.master$jar==72)] <- curation_1_resp

# Data curating

setwd("C:/Users/Katrina/Box/Duke/Kat-PhD/projects/fall-winter-2020-exp/Chapter.1/00_raw_data/resp_data/12.15.2020")
curation_2 <- read_excel('control.13.92.xlsx', sheet=6, col_names=TRUE)
control_2 <- read_excel('control.13.xlsx', sheet=6, col_names=TRUE)
curation_2 <- curation_2 %>%
  # Slopes need to be calculated relative to control and air partial O2 pressure following Tomlinson et al 2018
  mutate(controlled_02 = curation_2[,7]*0.21*(1/control_2[,7]), Time=curation_2[,3]*1) %>%
  filter(Time > 14)
curation_2_resp <- abs(summary(lm(t(t(curation_2$controlled_02))~t(t(curation_2$Time))))$coefficients[2])

respiration.master$respiration[which(respiration.master$jar==92)] <- curation_2_resp

names2

plot(t(t(curation_2$controlled_02))~t(t(curation_2$Time)), ylab="", xlab="")
abline(lm(t(t(curation_2$controlled_02))~t(t(curation_2$Time))), col="red", lwd=2)


# Data curating

  setwd("C:/Users/Katrina/Box/Duke/Kat-PhD/projects/fall-winter-2020-exp/Chapter.1/00_raw_data/resp_data/12.16.2020")
curation_3 <- read_excel('control.32.28.xlsx', sheet=6, col_names=TRUE)
control_3 <- read_excel('control.32.xlsx', sheet=6, col_names=TRUE)
curation_3 <- curation_3 %>%
  # Slopes need to be calculated relative to control and air partial O2 pressure following Tomlinson et al 2018
  mutate(controlled_02 = curation_3[,7]*0.21*(1/control_3[,7]), Time=curation_3[,3]*1) %>%
  filter(Time > 14)
curation_3_resp <- abs(summary(lm(t(t(curation_3$controlled_02))~t(t(curation_3$Time))))$coefficients[2])

respiration.master$respiration[which(respiration.master$jar==28)]  <- curation_3_resp

names3

plot(t(t(curation_3$controlled_02))~t(t(curation_3$Time)), ylab="", xlab="")
abline(lm(t(t(curation_3$controlled_02))~t(t(curation_3$Time))), col="red", lwd=2)


# Data curating
setwd("C:/Users/Katrina/Box/Duke/Kat-PhD/projects/fall-winter-2020-exp/Chapter.1/00_raw_data/resp_data/12.16.2020")
curation_4 <- read_excel('control.34.22.xlsx', sheet=6, col_names=TRUE)
control_4 <- read_excel('control.34.xlsx', sheet=6, col_names=TRUE)
curation_4 <- curation_4 %>%
  # Slopes need to be calculated relative to control and air partial O2 pressure following Tomlinson et al 2018
  mutate(controlled_02 = curation_4[,7]*0.21*(1/control_4[,7]), Time=curation_4[,3]*1) %>%
  filter(Time > 14)
curation_4_resp <- abs(summary(lm(t(t(curation_4$controlled_02))~t(t(curation_4$Time))))$coefficients[2])

respiration.master$respiration[which(respiration.master$jar==22)]  <- curation_4_resp

plot(t(t(curation_4$controlled_02))~t(t(curation_4$Time)), ylab="", xlab="")
abline(lm(t(t(curation_4$controlled_02))~t(t(curation_4$Time))), col="red", lwd=2)



# Data curating
setwd("C:/Users/Katrina/Box/Duke/Kat-PhD/projects/fall-winter-2020-exp/Chapter.1/00_raw_data/resp_data/12.16.2020")
curation_5 <- read_excel('control.34.25.xlsx', sheet=6, col_names=TRUE)
control_5 <- read_excel('control.34.xlsx', sheet=6, col_names=TRUE)
curation_5 <- curation_5 %>%
  # Slopes need to be calculated relative to control and air partial O2 pressure following Tomlinson et al 2018
  mutate(controlled_02 = curation_5[,7]*0.21*(1/control_5[,7]), Time=curation_5[,3]*1) %>%
  filter(Time > 16)
curation_5_resp <- abs(summary(lm(t(t(curation_5$controlled_02))~t(t(curation_5$Time))))$coefficients[2])

respiration.master$respiration[which(respiration.master$jar==25)] <- curation_5_resp

plot(t(t(curation_5$controlled_02))~t(t(curation_5$Time)), ylab="", xlab="")
abline(lm(t(t(curation_5$controlled_02))~t(t(curation_5$Time))), col="red", lwd=2)




# Data curating
setwd("C:/Users/Katrina/Box/Duke/Kat-PhD/projects/fall-winter-2020-exp/Chapter.1/00_raw_data/resp_data/12.16.2020")
curation_6 <- read_excel('control.35.42.xlsx', sheet=6, col_names=TRUE)
control_6 <- read_excel('control.35.xlsx', sheet=6, col_names=TRUE)
curation_6 <- curation_6 %>%
  # Slopes need to be calculated relative to control and air partial O2 pressure following Tomlinson et al 2018
  mutate(controlled_02 = curation_6[,7]*0.21*(1/control_6[,7]), Time=curation_6[,3]*1) %>%
  filter(Time > 16)
curation_6_resp <- abs(summary(lm(t(t(curation_6$controlled_02))~t(t(curation_6$Time))))$coefficients[2])

respiration.master$respiration[which(respiration.master$jar==42)] <-  curation_6_resp

plot(t(t(curation_6$controlled_02))~t(t(curation_6$Time)), ylab="", xlab="")
abline(lm(t(t(curation_6$controlled_02))~t(t(curation_6$Time))), col="red", lwd=2)




# Data curating
setwd("C:/Users/Katrina/Box/Duke/Kat-PhD/projects/fall-winter-2020-exp/Chapter.1/00_raw_data/resp_data/12.16.2020")
curation_7 <- read_excel('control.37.44.xlsx', sheet=6, col_names=TRUE)
control_7 <- read_excel('control.37.xlsx', sheet=6, col_names=TRUE)
curation_7 <- curation_7 %>%
  # Slopes need to be calculated relative to control and air partial O2 pressure following Tomlinson et al 2018
  mutate(controlled_02 = curation_7[,7]*0.21*(1/control_7[,7]), Time=curation_7[,3]*1) %>%
  filter(Time > 16)
curation_7_resp <- abs(summary(lm(t(t(curation_7$controlled_02))~t(t(curation_7$Time))))$coefficients[2])

respiration.master$respiration[which(respiration.master$jar==44)]  <- curation_7_resp

plot(t(t(curation_7$controlled_02))~t(t(curation_7$Time)), ylab="", xlab="")
abline(lm(t(t(curation_7$controlled_02))~t(t(curation_7$Time))), col="red", lwd=2)





# Data curating
setwd("C:/Users/Katrina/Box/Duke/Kat-PhD/projects/fall-winter-2020-exp/Chapter.1/00_raw_data/resp_data/12.16.2020")
curation_8 <- read_excel('control.34.15.xlsx', sheet=6, col_names=TRUE)
control_8 <- read_excel('control.34.xlsx', sheet=6, col_names=TRUE)
curation_8 <- curation_8 %>%
  # Slopes need to be calculated relative to control and air partial O2 pressure following Tomlinson et al 2018
  mutate(controlled_02 = curation_8[,7]*0.21*(1/control_8[,7]), Time=curation_8[,3]*1) %>%
  filter(Time > 16)
curation_8_resp <- abs(summary(lm(t(t(curation_8$controlled_02))~t(t(curation_8$Time))))$coefficients[2])

respiration.master$respiration[which(respiration.master$jar==15)]  <- curation_8_resp  ###line in master 

plot(t(t(curation_8$controlled_02))~t(t(curation_8$Time)), ylab="", xlab="")
abline(lm(t(t(curation_8$controlled_02))~t(t(curation_8$Time))), col="red", lwd=2)


# Data curating
setwd("C:/Users/Katrina/Box/Duke/Kat-PhD/projects/fall-winter-2020-exp/Chapter.1/00_raw_data/resp_data/12.16.2020")
curation_9 <- read_excel('control.35.33.xlsx', sheet=6, col_names=TRUE)
control_9 <- read_excel('control.35.xlsx', sheet=6, col_names=TRUE)
curation_9 <- curation_9 %>%
  # Slopes need to be calculated relative to control and air partial O2 pressure following Tomlinson et al 2018
  mutate(controlled_02 = curation_9[,7]*0.21*(1/control_9[,7]), Time=curation_9[,3]*1) %>%
  filter(Time > 16)
curation_9_resp <- abs(summary(lm(t(t(curation_9$controlled_02))~t(t(curation_9$Time))))$coefficients[2])

respiration.master$respiration[which(respiration.master$jar==33)]  <- curation_9_resp  ###line in master 

plot(t(t(curation_9$controlled_02))~t(t(curation_9$Time)), ylab="", xlab="")
abline(lm(t(t(curation_9$controlled_02))~t(t(curation_9$Time))), col="red", lwd=2)


### THE END

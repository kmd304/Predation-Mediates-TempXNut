library(piecewiseSEM)
# https://murphymv.github.io/semEff/articles/semEff.html
# https://jslefche.github.io/sem_book/composite-variables.html#what-is-a-composite-variable
#https://www.youtube.com/watch?v=rUUzF1_yaXg

sem.new <- read.csv("C:/Users/Katrina/Box/Duke/Kat-PhD/projects/fall-winter-2020-exp/Chapter.1/sem.new.csv", na.strings="NA",stringsAsFactors = FALSE)
head(sem.new)

sem.new$Temp1  <- as.numeric(sem.new$Temp)
sem.new$Spec1 <- log(sem.new$Spec+1)

head(sem.new)
sem.new <-na.omit(sem.new)

resp.t <- lm(Respiration ~Temp1*Nutrient + shannon.pro + shannon.bac +biomass.pro + Spec1, na.action = "na.fail", data= sem.new)
full.model.1 <- lm(Spec1~Temp1*Nutrient +shannon.pro + shannon.bac + biomass.pro, na.action = "na.fail", data= sem.new)
a <- lm(shannon.pro~Temp1*Nutrient + biomass.pro, na.action = "na.fail", data= sem.new)
b <- lm(shannon.bac~Temp1*Nutrient, na.action = "na.fail", data= sem.new)
c <- lm(biomass.pro~Temp1*Nutrient, na.action = "na.fail", data= sem.new)

model<- psem(resp.t, full.model.1, a, b, c) 
summary(model)

dSep(model, interactions = TRUE, .progressBar = TRUE)
AIC(model, AIC.type = "dsep")
fisherC(model, interactions = TRUE) ##a low p value is bad here so it is telling us that our model does match the assoctions that are impiled by our data
coefs(model)

plot(model)

summary(model)


piecewiseSEM:::plot.psem(
  piecewiseSEM::as.psem(model),
  node_attrs = data.frame(shape = "rectangle", color = "black", fillcolor = "grey"), ns_dashed = F,
  layout = "tree"
)

piecewiseSEM:::plot.psem(
  piecewiseSEM::as.psem(model),
  node_attrs = data.frame(shape = "rectangle", color = "black", fillcolor = "white"), ns_dashed = T, show = "std", add_edge_label_spaces = TRUE,
  layout = "tree",   edge_attrs = data.frame(style = "solid", color = "black")
)


library(semEff)
system.time(
  keeley.sem.boot <- bootEff(model, R = 500, seed = 13, parallel = "no"))
(keeley.sem.eff <- semEff(keeley.sem.boot))
summary(keeley.sem.eff, response = "biomass.pro")
summary(keeley.sem.eff, response = "Respiration")
summary(keeley.sem.eff, response = "shannon.bac")
summary(keeley.sem.eff, response = "shannon.pro")
summary(keeley.sem.eff, response = "Spec1")





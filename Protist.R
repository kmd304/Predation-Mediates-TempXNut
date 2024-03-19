rm(list=ls()) # remove everything currently held in the R memory
library(tidyverse)
library(vegan)
library(viridis)
library(ggpubr)
library(cowplot)
library(gridExtra)
library(grid)
library(ggplot2)
library(ggfortify)
library(factoextra)
library(corrplot)


trying <- read.csv("C:/Users/Katrina/Box/Duke/Kat-PhD/projects/fall-winter-2020-exp/Chapter.1/00_raw_data/trying.csv", na.strings="NA",stringsAsFactors = FALSE)
trying <-na.omit(trying)

# Here we transform variables nutrient, protist, and temp to factor:
trying$Nutrient1 <- as.factor(trying$Nutrient)
trying$Protist1 <- as.factor(trying$Protist)
trying$Temp1 <- as.factor(trying$Temp)


#############################################
#shannon diversity
H <- diversity(trying[,7:17])

#Observed Richness
richness <- specnumber(trying[,7:17])

# Pielou's Evenness
evenness <- H/log(trying[,7:17])

trying <- cbind( trying, cbind(shannon = H, richness = richness, pielou = evenness, trying)[,1:2])

###################################################################################

trying$Temp1 <- as.factor(trying$Temp)
trying$Nutrient1 <- as.factor(trying$Nutrient)


##################################################################
#this is richness for only protist jars 

hist(trying$richness)

mod<- lm(richness ~ Temp1*Nutrient1, data=trying)
summary(mod)
anova(mod)

#####################################
#this is shannon for only protist jars 

hist(trying$shannon)

mod2<- lm(shannon ~ Temp1*Nutrient1, data=trying)
summary(mod2)
anova_result <- aov(lm(shannon ~ Temp1*Nutrient1, data=trying))
anova_result


mod2<- lm(shannon ~ Temp1+Nutrient1, data=trying)
summary(mod2)
anova(mod2)

mod2<- lm(shannon ~ Temp1, data=trying)
summary(mod2)
anova(mod2)

mod2<- lm(shannon ~ Nutrient1, data=trying)
summary(mod2)
anova(mod2)


##############################################
##################################################################

trying.changed <- gather(trying, key = "species", value = "count",
         Blepharisma.sp, Colpidium.sp, Glaucoma.sp, Paramecium.bursaria, Tillina.magna, Colpoda.steinii, Strombidium.sp, Halteria.grandinella, Paramecium.aurelia, Cyclidium.glaucoma, Tetrahymena.pyriformis)
dim(trying.changed)

trying.changed$log.count <- log(trying.changed$count + 1)

trying.changed <-na.omit(trying.changed)

is.numeric(trying.changed$log.count)

ggplot(trying.changed, aes(x = Nutrient, y = shannon, group =Nutrient)) + 
  geom_boxplot()

ggplot(trying.changed, aes(x = Temp1, y = richness, group =Temp1)) + 
  geom_boxplot()


#######################################################################
trying.changed.new<-trying.changed %>%
  group_by(Nutrient1, Temp1, Protist1) %>%
  ungroup()

trying.changed.new$Nutrient1 <- factor(trying.changed.new$Nutrient1, levels = c("Half", "Full"))
trying$Nutrient1 <- factor(trying$Nutrient1, levels = c("Half", "Full"))


p2 <- ggplot(trying, aes(x=Nutrient1, y=shannon, group= Nutrient1,fill = Nutrient1)) +
  stat_boxplot(geom = "errorbar",
               width = 0.15) + 
  geom_boxplot(outlier.shape = NA, fatten =3, show.legend=FALSE)+
  geom_jitter( show.legend=FALSE, size=3, shape =1, stroke= 1.05, width = 0.35,) +
  scale_fill_manual(values = c("tan4", "darkolivegreen3")) +
  coord_cartesian(ylim = c(0, 1.5)) +
  theme_bw()+
  labs( y="Shannon Diversity") + theme(text = element_text(size = 25)) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(-1.4, "mm"))+
  theme(legend.key = element_blank(), strip.background = element_rect(colour="black", fill="white") ) +
  theme(legend.position = c(0.1, 0.90),legend.background = element_rect(fill = "white", color = "black"))

p2

######################################################################################################

trying.changed$Nutrient1 <- factor(trying.changed$Nutrient1, levels = c("Half", "Full"))

library(ggh4x)
every.jar<- subset(trying.changed, log.count>0)###remove NAs from column by subsetting values above zero
every.jar.1 <-ggplot(trying.changed, aes(fill=species, y=log.count, x= reorder(Jar.number, -log.count))) +  xlab("Jar") + ylab("Log Density") +
  facet_nested(.~Temp1+Nutrient1, scale="free") +
  geom_bar(position="stack", stat="identity", width = 1) +
  scale_y_continuous(expand=c(0,0),limits=c(0,14))

every.jar.2 <- every.jar.1 + labs(fill = "Species") +  theme_bw() + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(), axis.ticks.length = unit(-1.4, "mm"), axis.title.x = element_blank(),
                                               axis.text.x=element_blank(),
                                                                         axis.ticks.x=element_blank())  + theme(text = element_text(size = 25)) + theme(legend.key = element_blank(), strip.background = element_rect(colour="black", fill="white") ) 

every.jar.2

#####################################################################################
##################################### PCA

trying.21 <- trying[1:18]
trying.21 <-na.omit(trying.21)



trying.21$Nutrient <- as.factor(trying.21$Nutrient)
trying.21$Protist <- as.factor(trying.21$Protist)
trying.21$Temp  <- as.factor(trying.21$Temp)
trying.21$treatment  <- as.factor(trying.21$treatment)



###Preselecting variables########

all <-dplyr::select(trying.21, Jar.number, Nutrient, Temp, Protist, treatment, Spec, Respiration, Blepharisma.sp, Colpidium.sp,
                    Glaucoma.sp, Paramecium.bursaria , Tillina.magna, Colpoda.steinii , Strombidium.sp, Halteria.grandinella , Paramecium.aurelia, Cyclidium.glaucoma ,
                    Tetrahymena.pyriformis )


###Omitting NAs in data
all <-na.omit(all)


all.1<- all[, c(8:16)]+1
final_data <- cbind(all.1, all[,c(2:3)])


res.pca <- prcomp(log(all.1), center = TRUE, scale = TRUE) #got rid of everything but protists
summary(res.pca)

fviz_eig(res.pca)

eig.val <- get_eigenvalue(res.pca)
eig.val


fviz_eig(res.pca, choice = "eigenvalue", 
         addlabels=TRUE)

fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

fviz_pca_biplot(res.pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables/species color
                col.ind = "#696969"  # Individuals/each jar color
)

autoplot(res.pca, loadings = TRUE, loadings.label = TRUE)

autoplot(res.pca, loadings = TRUE, loadings.label = TRUE,
         data = all, colour = 'treatment')


res.var <- get_pca_var(res.pca) ###species
res.var$coord          # Coordinates
res.var$contrib[, 1:4]        # Contributions to the PCs
res.var$cos2           # Quality of representation 



# Results for individuals
res.ind <- get_pca_ind(res.pca)
res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation 

corrplot(res.var$cos2, is.corr=FALSE)

corrplot(res.var$contrib, is.corr=FALSE)

pca_points <- 
  # first convert the pca results to a tibble
  as_tibble(res.pca$x) %>% 
  bind_cols(all)

head(pca_points)


str(res.pca)
res.pca$x
varspec.PCAs <- cbind(all, res.pca$x[,1:3])

basic_plot <- 
  ggplot(pca_points, aes(x = PC1, y = PC2)) +
  geom_point(aes(colour = treatment), size= 3) +
  theme_light()

basic_plot


pca_hull <- 
  pca_points %>% 
  group_by(treatment) %>% 
  slice(chull(PC1, PC2))


chull_plot <- 
  basic_plot +
  geom_polygon(data = pca_hull,
               aes(fill = treatment, colour = treatment), linewidth=1, fill =NA,
               alpha = 0.3,show.legend = FALSE) 



chull_plot1 <- chull_plot +  theme(panel.border = element_blank(),
                                   panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(),
                                   axis.line = element_line(linewidth = 0.5, linetype = "solid",
                                                            colour = "black"))

chull_plot1 


pca_load <-  as_tibble(res.pca$rotation, rownames = 'variable') 
# we can rename the variables so they look nicer on the figure


head(pca_load)

library(ggrepel)
p <- chull_plot1 +
  geom_segment(data = pca_load, 
               aes(x = 0, y = 0, 
                   xend = (PC1*5),
                   yend = (PC2*5)),
               arrow = arrow(length = unit(1/2, 'picas')), linewidth = 0.5) + theme_bw()+ theme(panel.grid.minor = element_blank(),
                                                                                                panel.grid.major = element_blank(), axis.ticks.length = unit(-1.4, "mm")) +  
  theme(legend.position = c(0.8, 0.75)) + theme(text = element_text(size = 25)) +
  annotate("text", x = (pca_load$PC1*5), y = (pca_load$PC2*5),
           label = pca_load$variable)

p



figure.4.pdf<-ggarrange(every.jar.2,                                                 # First row with scatter plot
                        ggarrange(p2, p, ncol = 2, labels = c("B", "C")), # Second row with box and dot plots
                        nrow = 2, 
                        labels = "A"                                        # Labels of the scatter plot
) 


pdf(file= "figure.4.pdf", width=16, height=16, useDingbats=FALSE )
plot(figure.4.pdf)
dev.off()


pc.jars <- scores(res.pca, display = "sites")
pc.species <- scores(res.pca, display = "species")


d.permanova <- adonis2(pc.jars ~ Nutrient*Temp, method = "euclidean", data= all, permutations = 999)
d.permanova

#####################################################################################
##################################### INDI SPECIES

pc1_mod <- lm(Blepharisma.sp ~ Temp * Nutrient, data = trying.21)
summary(pc1_mod)
anova(pc1_mod)
boxplot((Blepharisma.sp)~Nutrient*Temp, data = trying.21)


pc1_mod <- lm(Colpidium.sp ~ Temp*Nutrient, data = trying.21)
summary(pc1_mod)
anova(pc1_mod)

pc1_mod <- lm(Glaucoma.sp ~ Temp*Nutrient, data = trying.20)
summary(pc1_mod)
anova(pc1_mod)



pc1_mod <- lm(Paramecium.aurelia ~ Temp, data = trying.20)
summary(pc1_mod)
anova(pc1_mod)


pc1_mod <- lm(Paramecium.aurelia ~ Temp*Nutrient, data = trying.20)
summary(pc1_mod)
anova(pc1_mod)


pc1_mod <- lm(Paramecium.bursaria ~ Temp*Nutrient, data = trying.20)
summary(pc1_mod)
anova(pc1_mod)


pc1_mod1 <- lm(Tillina.magna ~ Temp*Nutrient, data = trying.20)
summary(pc1_mod1)
anova(pc1_mod1)


pc1_mod <- lm(Colpoda.steinii ~ Temp*Nutrient, data = trying.20)
summary(pc1_mod)
anova(pc1_mod)


pc1_mod <- lm(Strombidium.sp ~ Temp*Nutrient, data = trying.20)
summary(pc1_mod)
anova(pc1_mod)


pc1_mod <- lm(Halteria.grandinella ~ Temp*Nutrient, data = trying.20)
summary(pc1_mod)
anova(pc1_mod)




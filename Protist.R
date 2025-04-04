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
library(dplyr)

biomass.pro <- read.csv("C:/Users/18568/Box/Duke/Kat-PhD/projects/fall-winter-2020-exp/Chapter.1/00_raw_data/master.pro.csv", na.strings="NA",stringsAsFactors = FALSE)
biomass.pro <- read.csv("C:/Users/kd239/Box/Duke/Kat-PhD/projects/fall-winter-2020-exp/Chapter.1/00_raw_data/master.pro.csv", na.strings="NA",stringsAsFactors = FALSE)
biomass.pro <- read.csv("C:/Users/Katrina/Box/Duke/Kat-PhD/projects/fall-winter-2020-exp/Chapter.1/00_raw_data/master.pro.csv", na.strings="NA",stringsAsFactors = FALSE)

trying <- read.csv("C:/Users/kd239/Box/Duke/Kat-PhD/projects/fall-winter-2020-exp/Chapter.1/00_raw_data/protist.csv", na.strings="NA",stringsAsFactors = FALSE)
trying <- read.csv("C:/Users/Katrina/Box/Duke/Kat-PhD/projects/fall-winter-2020-exp/Chapter.1/00_raw_data/protist.csv", na.strings="NA",stringsAsFactors = FALSE)
trying <- read.csv("C:/Users/18568/Box/Duke/Kat-PhD/projects/fall-winter-2020-exp/Chapter.1/00_raw_data/protist.csv", na.strings="NA",stringsAsFactors = FALSE)
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


J <- H/log(9)

trying <- cbind( trying, cbind(shannon = H, richness = richness, eve= J, trying)[,1:3])

###################################################################################

trying$Temp1 <- as.factor(trying$Temp)
trying$Nutrient1 <- as.factor(trying$Nutrient)


##################################################################
#this is richness for only protist jars 

hist(trying$richness)

mod<- lm(richness ~ Temp1*Nutrient1, data=trying)
summary(mod)
anova(mod)

modadd<- lm(richness ~ Temp1+Nutrient1, data=trying)
summary(mod)
anova(mod)

hist(trying$eve)
even<- lm(eve ~ Temp1*Nutrient1, data=trying)
summary(even)

#####################################
#this is shannon for only protist jars 


mod3<- lm(shannon ~ Temp1*Nutrient1, data=trying)
summary(mod3)

##############################################
##################################################################

trying.changed <- gather(trying, key = "species", value = "count",
         Blepharisma.sp, Colpidium.sp, Glaucoma.sp, Paramecium.bursaria, Tillina.magna, Colpoda.steinii, Strombidium.sp, Halteria.grandinella, Paramecium.aurelia, Cyclidium.glaucoma, Tetrahymena.pyriformis)
dim(trying.changed)

trying.changed$log.count <- log(trying.changed$count + 1)

trying.changed <-na.omit(trying.changed)

is.numeric(trying.changed$log.count)


final.abundances <- trying.changed %>%
  group_by(Jar.number, Nutrient, Temp, Respiration,Spec) %>%
  summarise(abundance = sum(count)) %>%
  ungroup() 
  
final.abundances$log.abundance <- log(final.abundances$abundance)

final.abundances$Nutrient <- as.factor(final.abundances$Nutrient)
final.abundances$Temp <- as.factor(final.abundances$Temp)

final.abundances$Nutrient <- factor(final.abundances$Nutrient, levels = c("Half", "Full"))
head(final.abundances)

# Add the new column
final.abundances$per.ml <- (1 / 0.31) * final.abundances$abundance

final.abundances


##better model
nl <- lm(log.abundance ~ Nutrient*Temp, data = final.abundances[-which(final.abundances$abundance > 200),])
summary(nl)

ggplot(final.abundances, aes(x=Nutrient, y= log.abundance, group = Nutrient))+
  geom_boxplot()

ggplot(final.abundances, aes(x=Temp, y= log.abundance, group = Temp))+
  geom_boxplot()

# Summarize total counts by species
species_abundance <- trying.changed %>%
  group_by(species) %>%
  summarise(total_count = sum(count)) %>%
  arrange(desc(total_count))

# View the most abundant species
head(species_abundance)

species_treatment_abundance <- trying.changed %>%
  group_by(Nutrient, Temp, species) %>%
  summarise(total_count = sum(count), .groups = "drop")

# Find the most abundant species per treatment
most_abundant_per_treatment <- species_treatment_abundance %>%
  group_by(Nutrient, Temp) %>%
  filter(total_count == max(total_count)) %>%
  arrange(Nutrient, Temp)

# View the results
most_abundant_per_treatment

ggplot(most_abundant_per_treatment, aes(x = interaction(Nutrient, Temp), y = total_count, fill = species)) +
  geom_bar(stat = "identity") +
  labs(x = "Treatment (Nutrient x Temp)", y = "Total Count", title = "Most Abundant Species per Treatment") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Summarize counts of each species by treatment
species_by_treatment <- trying.changed %>%
  group_by(Nutrient, Temp, species) %>%
  summarise(total_count = sum(count), .groups = "drop")

# View the summarized data
head(species_by_treatment)
print(species_by_treatment)

# Convert Nutrient and Temp to factors with appropriate levels
species_by_treatment$Nutrient <- factor(species_by_treatment$Nutrient, levels = c("Half", "Full"))

species_by_treatment$Temp <- factor(species_by_treatment$Temp, levels = c("22", "25"))

levels(species_by_treatment$Nutrient)
levels(species_by_treatment$Temp)

# Filter for each treatment and plot
p1 <- ggplot(species_by_treatment %>% filter(Nutrient == "Half", Temp == "22"),
             aes(x = species, y = total_count, fill = species)) +
  geom_bar(stat = "identity") +
  labs(title = "Half Nutrient, Low Temp", x = "Species", y = "Total Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2 <- ggplot(species_by_treatment %>% filter(Nutrient == "Half", Temp == "25"),
             aes(x = species, y = total_count, fill = species)) +
  geom_bar(stat = "identity") +
  labs(title = "Half Nutrient, High Temp", x = "Species", y = "Total Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p3 <- ggplot(species_by_treatment %>% filter(Nutrient == "Full", Temp == "22"),
             aes(x = species, y = total_count, fill = species)) +
  geom_bar(stat = "identity") +
  labs(title = "Full Nutrient, Low Temp", x = "Species", y = "Total Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p4 <- ggplot(species_by_treatment %>% filter(Nutrient == "Full", Temp == "25"),
             aes(x = species, y = total_count, fill = species)) +
  geom_bar(stat = "identity") +
  labs(title = "Full Nutrient, High Temp", x = "Species", y = "Total Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Display plots
library(gridExtra)
grid.arrange(p1, p2, p3, p4, ncol = 2)


finallypro <- ggplot(data=final.abundances, aes(Temp, log.abundance, group = Temp, fill = Temp)) +
  scale_fill_manual(values=c("dodgerblue", "tomato2"))+
  stat_boxplot(geom ='errorbar') + 
  geom_boxplot(outlier.shape = NA, outlier.colour = NULL)+
  ggtitle("B: Low Nutrient               C: High Nutrient") +
  facet_wrap(~Nutrient) +
  coord_cartesian(ylim = c(1.5, 5))+
  theme_bw()

finallypro <- finallypro + labs(x="Temperature (\u00B0C)", y="Log Abundance of Ciliates") + theme(text = element_text(size = 25)) 

finallypro <- finallypro + theme(legend.background = element_rect(fill = "white", color = "black"))+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(-1.4, "mm"))


finallypro <- finallypro + theme(strip.background = element_blank(), strip.text.x = element_blank(),    axis.text.x = element_text(size = 20),  axis.title.x = element_blank(),  plot.title = element_text(size = 25))  
finallypro



finallypronut <- ggplot(data=final.abundances, aes(Nutrient, log.abundance, group = Nutrient, fill = Nutrient)) +
  scale_fill_manual(values=c("tan4", "darkolivegreen3"))+
  stat_boxplot(geom ='errorbar') + 
  geom_boxplot(outlier.shape = NA, outlier.colour = NULL)+
  geom_jitter( show.legend=FALSE, size=3, shape =1, stroke= 1.05, width = 0.35,) +
  facet_wrap(~Temp) +
  coord_cartesian(ylim = c(1.5, 5))+
  theme_bw()

finallypronut <- finallypronut + labs(x="Nutrient", y="Log Abundance of Ciliates") + theme(text = element_text(size = 25)) 

finallypronut <- finallypronut + theme(legend.background = element_rect(fill = "white", color = "black"))+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(-1.4, "mm"))
finallypronut <- finallypronut +  theme_bw() +  theme(legend.position = "none", axis.text.x=element_blank()) + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(), axis.ticks.length = unit(-1.4, "mm"), axis.title.x = element_blank(),
                                                     axis.text.x=element_blank(),
                                                     axis.ticks.x=element_blank())  + theme(text = element_text(size = 25)) + theme(legend.key = element_blank(), strip.background = element_rect(colour="black", fill="white") ) 

  

#######################################################################
trying.changed.new<-trying.changed %>%
  group_by(Nutrient1, Temp1, Protist1) %>%
  ungroup()

trying.changed.new$Nutrient1 <- factor(trying.changed.new$Nutrient1, levels = c("Half", "Full"))
trying$Nutrient1 <- factor(trying$Nutrient1, levels = c("Half", "Full"))


p3 <- ggplot(trying, aes(x=Nutrient1, y=shannon, group= Nutrient1,fill = Nutrient1)) +
  stat_boxplot(geom = "errorbar",
               width = 0.15) + 
  geom_boxplot(outlier.shape = NA, fatten =3, show.legend=FALSE)+
  geom_jitter( show.legend=FALSE, size=3, shape =1, stroke= 1.05, width = 0.35,) +
  scale_fill_manual(values = c("tan4", "darkolivegreen3")) +
  coord_cartesian(ylim = c(0, 1.5)) +
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(-1.4, "mm"))+
  theme(legend.key = element_blank(), strip.background = element_rect(colour="black", fill="white") ) +
  theme(legend.position.inside = c(0.1, 0.90),axis.title.x=element_blank(),legend.background = element_rect(fill = "white", color = "black"))

p3


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



##################################
# Step 1: Calculate average dimensions and total count for each class within each jar
biomass_class_jar <- biomass.pro %>%
  group_by(Jar.number, Class) %>%
  summarise(
    avg_length = mean(Length, na.rm = TRUE),
    avg_width = mean(Width, na.rm = TRUE),
    count = n()  # This will count the number of rows (particles) for each class in each jar
  ) %>%
  ungroup()

# Step 2: Calculate ellipsoid volume using averaged dimensions and multiply by count
biomass_class_jar <- biomass_class_jar %>%
  mutate(
    # Calculate ellipsoid volume based on average dimensions in µm³
    avg_ellipsoid_volume_um3 = (4/3) * pi * (avg_length / 2) * (avg_width / 2)^2,
    
    # Convert volume from µm³ to cm³ (1 µm³ = 1e-12 cm³)
    avg_ellipsoid_volume_cm3 = avg_ellipsoid_volume_um3 * 1e-12,
    
    # Multiply the volume by the count to get total volume per class within each jar
    total_volume_cm3_per_class_jar = avg_ellipsoid_volume_cm3 * count
  )

# Step 3: Sum the total count and biomass per jar
total_biomass_and_count_per_jar <- biomass_class_jar %>%
  group_by(Jar.number) %>%
  summarise(
    total_biomass_cm3 = sum(total_volume_cm3_per_class_jar, na.rm = TRUE),  # Sum total biomass
    total_count = sum(count, na.rm = TRUE)                                  # Sum total count of particles
  ) %>%
  ungroup()

# Since density is 1g/cm³, total biomass in grams is the same as total volume in cm³
total_biomass_and_count_per_jar$total_biomass_g = total_biomass_and_count_per_jar$total_biomass_cm3

# View the result
print(total_biomass_and_count_per_jar)


# Ensure both datasets have the Jar.number column for merging
total_biomass_and_count_per_jar <- total_biomass_and_count_per_jar %>%
  right_join(trying %>% select(Jar.number, Nutrient, Temp), by = "Jar.number")

# Check if the merge was successful
head(total_biomass_and_count_per_jar)

total_biomass_and_count_per_jar$Temp <- as.factor(total_biomass_and_count_per_jar$Temp)
total_biomass_and_count_per_jar$Nutrient <- as.factor(total_biomass_and_count_per_jar$Nutrient)

# Linear model to analyze the effect of Temp and Nutrient on total biomass
biomass_model <- lm(total_biomass_g ~ Temp * Nutrient, data = total_biomass_and_count_per_jar)

# Summary of the model
summary(biomass_model)

# ANOVA to check for significance
anova(biomass_model)


# Ensure "Low" comes first and "High" comes second
total_biomass_and_count_per_jar$Nutrient <- factor(total_biomass_and_count_per_jar$Nutrient, 
                                                   levels = c("Half", "Full"),  # Reversing order
                                                   labels = c("Low", "High"))  # Matching labels

# Generate the plot with the correct order
ggplot(total_biomass_and_count_per_jar, aes(x = Nutrient, y = total_biomass_g)) +
  geom_boxplot() +
  facet_wrap(~ Temp) +
  labs(title = "Total Ciliate Biomass Across Nutrient Levels by Temperature",
       x = "Nutrient Treatment", y = "Total Biomass (g)") +
  theme_bw()


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


pca1 <- res.pca$x[,1]


pca2 <- res.pca$x[,2]

pca3 <- res.pca$x[,3]

res.var <- get_pca_var(res.pca) ###species
res.var$coord          # Coordinates
res.var$contrib[, 1:4]        # Contributions to the PCs
res.var$cos2           # Quality of representation 



# Results for individuals- jar/samples
res.ind <- get_pca_ind(res.pca)
res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation 

corrplot(res.var$cos2, is.corr=FALSE)

corrplot(res.var$contrib[, 1:4], is.corr=FALSE)

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
  theme(legend.position.inside  = c(0.8, 0.75)) + theme(text = element_text(size = 25)) +
  annotate("text", x = (pca_load$PC1*5), y = (pca_load$PC2*5),
           label = pca_load$variable)

p


library(grid)
# Move to a new page
grid.newpage()
# Create layout : nrow = 3, ncol = 2
pushViewport(viewport(layout = grid.layout(nrow = 3, ncol = 2)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 

print(finallypronut , vp = define_region(row = 1, col = 1:2))
print(every.jar.2, vp = define_region(row = 2, col = 1:2))
print(p3, vp = define_region(row = 3, col = 1))
print(p, vp = define_region(row = 3, col = 2))

dev.off()


pc.jars <- scores(res.pca, display = "sites")
pc.jars <- res.pca$x  # Coordinates of individuals (samples or jars)

pc.species <- scores(res.pca, display = "species")
pc.species <- res.pca$rotation  # Coordinates of variables (species)





d.permanova <- adonis2(pc.jars ~ Nutrient*Temp, method = "euclidean", data= all, permutations = 999)
d.permanova

d.permanova1 <- adonis2(pca1 ~ Nutrient*Temp, method = "euclidean", data= all, permutations = 999)
d.permanova1

summary(lm(pca1 ~ Nutrient*Temp, data= all))

d.permanova2 <- adonis2(pca2 ~ Nutrient*Temp, method = "euclidean", data= all, permutations = 999)
d.permanova2

summary(lm(pca2 ~ Nutrient*Temp, data= all))

d.permanova3 <- adonis2(pca3 ~ Nutrient*Temp, method = "euclidean", data= all, permutations = 999)
d.permanova3

summary(lm(pca3 ~ Nutrient*Temp, data= all))

# Calculate Bray-Curtis dissimilarity
bray_curtis <- vegdist(all.1, method = "bray")

# Convert the distance object to a matrix for easier manipulation
bray_curtis_matrix <- as.matrix(bray_curtis)

mean_dissimilarity_per_sample <- apply(bray_curtis_matrix, 1, mean)

# Combine the mean dissimilarity with the original data (optional)
final_data_with_dissimilarity <- cbind(all, MeanDissimilarity = mean_dissimilarity_per_sample)

# View the final data with dissimilarity values
head(final_data_with_dissimilarity)

library(writexl)
write_xlsx(final_data_with_dissimilarity, path = "sfinal_data_with_dissimilarity.xlsx")



#####################################################################################
##################################### INDI SPECIES

pc1_mod <- lm(Blepharisma.sp ~ Temp * Nutrient *Respiration, data = trying.21)
summary(pc1_mod)
anova(pc1_mod)
boxplot((Blepharisma.sp)~Nutrient*Temp, data = trying.21)


pc1_mod <- lm(Colpidium.sp ~ Temp*Nutrient*Respiration, data = trying.21)
summary(pc1_mod)
anova(pc1_mod)

pc1_mod <- lm(Glaucoma.sp ~ Temp*Nutrient*Respiration, data = trying.21)
summary(pc1_mod)
anova(pc1_mod)


pc1_mod <- lm(Paramecium.aurelia ~ Temp*Nutrient*Respiration, data = trying.21)
summary(pc1_mod)
anova(pc1_mod)


pc1_mod <- lm(Paramecium.bursaria ~ Temp*Nutrient*Respiration, data = trying.21)
summary(pc1_mod)
anova(pc1_mod)


pc1_mod1 <- lm(Tillina.magna ~ Temp*Nutrient*Respiration, data = trying.21)
summary(pc1_mod1)
anova(pc1_mod1)


pc1_mod <- lm(Colpoda.steinii ~ Temp*Nutrient*Respiration, data = trying.21)
summary(pc1_mod)
anova(pc1_mod)


pc1_mod <- lm(Strombidium.sp ~ Temp*Nutrient*Respiration, data = trying.21)
summary(pc1_mod)
anova(pc1_mod)


pc1_mod <- lm(Halteria.grandinella ~ Temp*Nutrient*Respiration, data = trying.21)
summary(pc1_mod)
anova(pc1_mod)


####################################### on function
library(car)
something2 <- lm(cbind(shannon,Respiration)~Temp1*Nutrient1, data=trying)
summary(something2)
Anova(something2)
vcov(something2)


plot(pca1~Respiration, data = all)
plot(pca2~Respiration, data = all)


something1 <- lm(cbind(log.abundance,Respiration)~Temp*Nutrient, data = final.abundances[-which(final.abundances$abundance > 200),])
summary(something1)
Anova(something1)

cor.test(final.abundances[-which(final.abundances$abundance > 200),]$Respiration,final.abundances[-which(final.abundances$abundance > 200),]$log.abundance)



cor.test(final.abundances[-which(final.abundances$abundance > 200),]$Spec,
         final.abundances[-which(final.abundances$abundance > 200),]$log.abundance)

lm_spec_abundance <- lm(Spec ~ log.abundance, data = final.abundances[-which(final.abundances$abundance > 200),])
summary(lm_spec_abundance)


ggplot(final.abundances[-which(final.abundances$abundance > 200),], aes(x = log.abundance, y = Spec)) +
  geom_point() +
  geom_smooth(method = "lm", col = "blue") +
  labs(title = "Relationship between Biomass (OD600) and Abundance",
       x = "Log Abundance",
       y = "Biomass (OD600)")


ng <-lm(cbind(pca1,pca2,pca3,Respiration)~Nutrient*Temp, data = all)
summary(ng)
Anova(ng)

all <-cbind(all, pc.jars)

cor(all[, c("PC1", "PC2", "PC3", "Respiration")], use = "complete.obs")

vars <- all[, c("PC1", "PC2", "PC3", "Respiration")]
cor.test(vars$PC1, vars$Respiration)
cor.test(vars$PC2, vars$Respiration)
cor.test(vars$PC3, vars$Respiration)



cor.test(trying$Respiration,trying$shannon)
summary(lm(Respiration~shannon, data = trying))






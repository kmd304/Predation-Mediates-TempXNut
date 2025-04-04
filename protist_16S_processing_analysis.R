rm(list=ls())
#import QIIME2 .qza files and create phyloseq object for 16S
library(qiime2R)
library(tidyverse)
library(phyloseq)
library(ggplot2)
library(ggh4x)
library(vegan)
library(car)
library(dplyr)     

setwd("C:/Users/Katrina/Box/Duke/Kat-PhD/projects/fall-winter-2020-exp/code from alyssa")
setwd("C:/Users/18568/Box/Duke/Kat-PhD/projects/fall-winter-2020-exp/code from alyssa")
setwd("C:/Users/kd239/Box/Duke/Kat-PhD/projects/fall-winter-2020-exp/code from alyssa")


#import QIIME2 .qza files and create phyloseq object for 16S

mapping<-read.delim("Book2.txt")
SVs16S<-read_qza("16S-dada2table_all.qza")
taxonomy16S<-read_qza("16S-taxonomy-all.qza")
tax16S<-taxonomy16S$data %>% as.tibble() %>%mutate(Taxon=gsub("[a-z]__", "", Taxon)) %>% separate(Taxon, sep=";", c("Kingdom","Phylum","Class","Order","Family","Genus","Species")) ###remove SILVA notation and change column names to classifications

#generate phyloseq object
physeq16S<-phyloseq(
  otu_table(SVs16S$data, taxa_are_rows = T), 
  tax_table(as.data.frame(tax16S) %>% 
              select(-Confidence) %>% 
              column_to_rownames("Feature.ID") %>% 
              as.matrix()), #moving the taxonomy to the way phyloseq wants it
  sample_data(mapping %>% 
                as.data.frame() %>% 
                column_to_rownames("SampleID"))) %>%
  subset_taxa(Family!= "Mitochondria" | is.na(Family)) %>%
  subset_taxa(Class!="Chloroplast" | is.na(Class))%>%
  subset_taxa(Order!="Chloroplast" | is.na(Order))

physeq16S<-subset_samples(physeq16S, analyze=="keep")%>%
  rarefy_even_depth()

physeq16S@sam_data$Nutrient <- factor(physeq16S@sam_data$Nutrient, levels = c("Half", "Full"))


#summariza taxa at phylum level for all samples
phyla16S<-physeq16S %>%
  tax_glom(taxrank = "Phylum")%>%
  transform_sample_counts(function(x){x/sum(x)})%>%
  psmelt()%>%
  filter(Abundance>0.02)%>%
  arrange(Phylum)
newp_16S<-aggregate(Abundance~Sample,data=phyla16S,sum)
newp_16S$Abundance<-1-newp_16S$Abundance
newp_16S$OTU<-"Other"
newp_16S$Phylum<-"Other"
newp_16S$Kingdom<-"Other"
newp_16S2<-merge(newp_16S,mapping,by.x=c("Sample"),by.y=c("SampleID"))
phyla16S2<-merge(phyla16S,newp_16S2,all=TRUE)



#plot taxonomy at phylum level
plot1 <- ggplot(phyla16S2,aes(x=Sample, y=Abundance, fill=Phylum)) +
  geom_bar(stat="identity", position="stack", color="black")+ theme_bw() +
  theme(axis.text=element_text(size=8),
        panel.grid.major = element_blank(),
        axis.title=element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        legend.text=element_text(size=8),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(color="black",size=8, face="bold")) +
  facet_nested(.~Protist+Temp+Nutrient,scales="free")+
  scale_fill_manual(values=c("#ebac23","#00c6f8","darkgreen","#dab1da", "#b24502","red","lightpink","palegreen", "royalblue3", "darkorange","#00a76c", "black", "brown", "purple"), name="Bacterial Phylum")+
  xlab("Samples") +
  ylab("Relative Abundance Phylum >2%") +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(color="black",size=8, face="bold"))+
  guides(fill=guide_legend(keywidth = 0.5,keyheight = 0.5))+ylab("Relative Frequency")

plot1

phylum_summary <- phyla16S2 %>%
  group_by(Phylum) %>%
  summarize(TotalAbundance = sum(Abundance, na.rm = TRUE)) %>%
  arrange(desc(TotalAbundance))

# View the most abundant phylum
head(phylum_summary)

# calculate alpha diversity
library(data.table)
set.seed(092)
diversity16S<-estimate_richness(physeq16S, measures=c("Observed","Shannon"))
diversity16S<-setDT(diversity16S,keep.rownames = "SampleID")[]
diversity16S$SampleID<-str_replace_all(diversity16S$SampleID,"\\.","-")

alpha16S<-diversity16S%>%
  left_join(mapping)%>%
  filter(Shannon!=0)


# ANOVA on alpha
aov.16S.observed<-aov(Observed~Temp*Nutrient*Protist, data=alpha16S) %>%
  Anova( type = 3)


aov.16S.shannon<-aov(Shannon~Temp*Nutrient*Protist, data=alpha16S) %>%
  Anova( type = 3)




alpha16S$Nutrient <- factor(alpha16S$Nutrient, levels = c("Half", "Full"))


# plot alpha diversity
alpha <-ggplot(alpha16S, 
       aes(x=Protist,
           y=Shannon,
           group=Protist))+
  stat_boxplot(geom = 'errorbar')+
  geom_boxplot(aes(fill=Protist),outlier.shape=NA) +
  facet_nested(.~Temp+Nutrient,scales="free") +
  theme_bw()+
  theme(axis.text=element_text(size=8),
        axis.title.x=element_blank(),
        panel.grid.major = element_blank(),
        axis.title=element_text(size=8,face="bold"),
        axis.text.x = element_blank(),
        legend.text=element_text(size=8),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(color="black",size=8, face="bold")) +
  scale_fill_manual(values=c("deeppink","purple"), name="Protist") +
  guides(fill=guide_legend(keywidth = 0.5,keyheight = 0.5))


# calculate distance matrix
GPr16S  = transform_sample_counts(physeq16S, function(x) x / sum(x) )
physeq16S2 <- filter_taxa(GPr16S, function(x) sum(x) > 1e-5, TRUE)
ord16S_bc<-ordinate(physeq16S2, "NMDS", "bray")



# permanova on betadiv
bc16S<-phyloseq::distance(physeq16S2, method="bray")
sampled16S<-data.frame(sample_data(physeq16S2))
BC16S.adns<-adonis2(bc16S~Temp*Nutrient*Protist, data=sampled16S)
BC16S.adns


bray_curtis_matrix_bac <- as.matrix(bc16S)
print(bray_curtis_matrix_bac)
# Calculate mean dissimilarity per sample
mean_dissimilarity_per_sample.bac <- apply(bray_curtis_matrix_bac, 1, mean)

# Combine the mean dissimilarity with the original data (optional)
final_data_with_dissimilarity.bac <- cbind(sampled16S, MeanDissimilarity = mean_dissimilarity_per_sample.bac)

# View the final data with dissimilarity values
head(final_data_with_dissimilarity.bac)
library(dplyr)

# Rename the column
final_data_with_dissimilarity.bac <- final_data_with_dissimilarity.bac %>%
  rename(Jar.number = samp)

final_data_with_dissimilarity.bac <- final_data_with_dissimilarity.bac %>%
  rename(MeanDissimilaritybacteria = MeanDissimilarity)

# View the updated data frame
head(final_data_with_dissimilarity.bac)

# Remove the leading 'S' from Jar.number
final_data_with_dissimilarity.bac$Jar.number <- gsub("^S", "", final_data_with_dissimilarity.bac$Jar.number)

# View the updated data frame
head(final_data_with_dissimilarity.bac)

final_data_with_dissimilarity.bac$Jar.number <- as.factor(final_data_with_dissimilarity.bac$Jar.number)


library(writexl)
write_xlsx(final_data_with_dissimilarity.bac, path = "sfinal_data_with_dissimilarity.bac.xlsx")

#ng <-lm(cbind(bc16S,Respiration)~Nutrient*Temp, data = sampled16S)
#summary(ng)
#Anova(ng)

#plot betadiversity
p_nmds_bc_physeq16S<- plot_ordination(physeq16S2, 
                                      ord16S_bc, 
                                      type ="sample", 
                                      color = "Protist", 
                                      shape="Nutrient_Temp") +
  theme_bw()+ 
  labs(title="16S-Bray Curtis NMDS") +
  theme(plot.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_shape_manual(values=c(0, 15, 2, 17)) +
  scale_color_manual(values=c("deeppink","purple"), name="Treatment") +
  geom_point(size = 4)+
  ylim(-1, 1)

p_nmds_bc_physeq16S$layers <- p_nmds_bc_physeq16S$layers[-1]

p_nmds_bc_physeq16S + xlim(-1,1) + ylim(-1,1)

library(grid)
# Move to a new page
grid.newpage()
# Create layout : nrow = 3, ncol = 2
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 2)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 

print(plot1 , vp = define_region(row = 1, col = 1:2))
print(alpha, vp = define_region(row = 2, col = 1))
print(p_nmds_bc_physeq16S, vp = define_region(row = 2, col = 2))

dev.off()



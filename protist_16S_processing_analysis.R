
#import QIIME2 .qza files and create phyloseq object for 16S

mapping<-read.delim("mapping_16S_all.txt")
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
ggplot(phyla16S,aes(x=Sample, y=Abundance, fill=Phylum)) +
  geom_bar(stat="identity",color="black")+ theme_bw() +
  theme(axis.text=element_text(size=8),
        panel.grid.major = element_blank(),
        axis.title=element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        legend.text=element_text(size=8),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(color="black",size=8, face="bold")) +
  facet_nested(.~Protist+Temp+Nutrient,scales="free")+
  scale_fill_manual(values=c("#ebac23","deeppink","#008cf9","#006e00","#00bbad","#d163e6","#b24502","#ff9287","#00c6f8","purple","#878500","#00a76c"), name="Bacterial Phylum")+
  xlab("Samples") +
  ylab("Relative Abundance Phylum >2%") +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(color="black",size=8, face="bold"))+
  guides(fill=guide_legend(keywidth = 0.5,keyheight = 0.5))+ylab("Relative Frequency")


# calculate alpha diversity
library(data.table)
set.seed(092)
diversity16S<-estimate_richness(physeq16S, measures=c("Observed","Shannon"))
diversity16S<-setDT(diversity16S,keep.rownames = "SampleID")[]
alpha16S<-diversity16S%>%
  left_join(mapping)%>%
  filter(Shannon!=0)


# ANOVA on alpha
aov.16S.observed<-aov(Observed~Temp*Nutrient*Protist, data=alpha16S) %>%
Anova( type = 3)

aov.16S.shannon<-aov(Shannon~Temp*Nutrient*Protist, data=alpha16S) %>%
  Anova( type = 3)


# plot alpha diversity
ggplot(alpha16S, 
       aes(x=Protist,
           y=Shannon,
          group=Protist))+
  stat_boxplot(geom = 'errorbar')+
  geom_boxplot(aes(fill=Protist),outlier.shape=NA) +
  facet_nested(gene~Temp+Nutrient,scales="free") +
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
sampled16S<-data.frame(sample_data(physeq16S2))
BC16S.adns<-adonis2(bc16S~Temp*Nutrient*Protist, data=sampled16S)

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
  scale_color_manual(values=c("purple","deeppink"), name="Treatment") +
  geom_point(size = 4)


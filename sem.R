library(lavaan)

# https://murphymv.github.io/semEff/articles/semEff.html
# https://jslefche.github.io/sem_book/composite-variables.html#what-is-a-composite-variable
#https://www.youtube.com/watch?v=rUUzF1_yaXg

#run protist.R and protist_16S_processing_analysis.R first
sem.new <- read.csv("C:/Users/Katrina/Box/Duke/Kat-PhD/projects/fall-winter-2020-exp/Chapter.1/sem.new1.csv", na.strings="NA",stringsAsFactors = FALSE)
sem.new <- read.csv("C:/Users/18568/Box/Duke/Kat-PhD/projects/fall-winter-2020-exp/Chapter.1/sem.new1.csv", na.strings="NA",stringsAsFactors = FALSE)
sem.new <- read.csv("C:/Users/kd239/Box/Duke/Kat-PhD/projects/fall-winter-2020-exp/Chapter.1/sem.new1.csv", na.strings="NA",stringsAsFactors = FALSE)
head(sem.new)


sem.new <- sem.new %>%
  left_join(final_data_with_dissimilarity %>%
              select(Jar.number, MeanDissimilarity), by = "Jar.number")

sem.new$Spec1 <- log(sem.new$Spec+1)
sem.new$Temp1 <- as.numeric(sem.new$Temp)


head(sem.new)
sem.new <-na.omit(sem.new)
# Convert Jar.number to a factor
sem.new$Jar.number <- as.factor(sem.new$Jar.number)

sem.new <- sem.new %>%
  left_join(final_data_with_dissimilarity.bac %>%
              select(Jar.number, MeanDissimilaritybacteria), by = "Jar.number")

head(sem.new)


sem.forpro <- sem.new %>%
  select(-Spec, -Spec1, -Temp1)

head(sem.forpro)

write.csv(sem.forpro, "sem.forpro.csv", row.names = FALSE)


sem.new$biomass.pro1 <- total_biomass_and_count_per_jar$total_biomass_g

head(sem.new)
###################################################################
sem.new$Temperature <- scale(sem.new$Temp1)
sem.new$Nutrient <- scale(sem.new$Nutrient)
sem.new$Respiration <- scale(sem.new$Respiration)

sem.new$Ciliate_community_composition <- scale(sem.new$MeanDissimilarity)
sem.new$Prokaryotic_community_composition <- scale(sem.new$MeanDissimilaritybacteria)

sem.new$Ciliate_diversity <- scale(sem.new$shannon.pro)
sem.new$Prokaryote_diversity <- scale(sem.new$shannon.bac)

sem.new$biomass_prokaryote <- scale(sem.new$Spec1)
sem.new$biomass_ciliate <- scale(sem.new$biomass.pro1)

# Create composite biomass variable
sem.new$total_biomass <- sem.new$biomass_ciliate + sem.new$biomass_prokaryote

# Define model
model.1 <- '
  Ciliate_diversity ~ Temperature + Nutrient + Temperature : Nutrient 
  Prokaryote_diversity ~ Temperature + Nutrient + Temperature : Nutrient 
  Ciliate_community_composition ~ Temperature + Nutrient + Ciliate_diversity + Temperature : Nutrient + Prokaryotic_community_composition 
  Prokaryotic_community_composition ~ Temperature + Nutrient + Prokaryote_diversity + Temperature : Nutrient + Ciliate_community_composition  
  total_biomass ~ Temperature + Nutrient + Temperature : Nutrient + Prokaryote_diversity + Ciliate_diversity + Ciliate_community_composition + Prokaryotic_community_composition
  Respiration ~ Temperature + Nutrient + Ciliate_diversity + Prokaryote_diversity + total_biomass + Ciliate_community_composition + Prokaryotic_community_composition + Temperature : Nutrient
'

# Fit model
fit.1 <- sem(model.1, data = sem.new)

# Show summary with fit measures
summary(fit.1, fit.measures = TRUE, standardized = TRUE)

model.2 <- '
  # Regressions
  Ciliate_diversity ~ Temperature + Nutrient + Temperature : Nutrient 
  Prokaryote_diversity ~ Temperature + Nutrient + Temperature : Nutrient 
  Ciliate_community_composition ~ Temperature + Nutrient + Ciliate_diversity + Temperature : Nutrient + Prokaryotic_community_composition 
  Prokaryotic_community_composition ~ Temperature + Nutrient + Prokaryote_diversity + Temperature : Nutrient + Ciliate_community_composition  
  total_biomass ~ Temperature + Nutrient + Temperature : Nutrient + Prokaryote_diversity + Ciliate_diversity + Ciliate_community_composition + Prokaryotic_community_composition
  Respiration ~ Temperature + Nutrient + Ciliate_diversity + total_biomass + Ciliate_community_composition + Prokaryotic_community_composition
  
  # Covariance between Ciliate_diversity and Ciliate_community_composition
  # Ciliate_diversity ~~ Ciliate_community_composition
'


# Fit the model
fit.2 <- sem(model.2, data = sem.new)

# Summary of the model with fit measures
summary(fit.2, fit.measures = TRUE, standardized = TRUE)

model.3 <- '
  # Regressions
  Ciliate_diversity ~ Temperature + Nutrient + Temperature : Nutrient 
  Prokaryote_diversity ~ Temperature + Temperature : Nutrient 
  Ciliate_community_composition ~ Temperature + Nutrient + Ciliate_diversity + Temperature : Nutrient + Prokaryotic_community_composition 
  Prokaryotic_community_composition ~ Temperature + Nutrient + Prokaryote_diversity + Temperature : Nutrient + Ciliate_community_composition  
  total_biomass ~ Temperature + Nutrient + Temperature : Nutrient + Prokaryote_diversity + Ciliate_diversity + Ciliate_community_composition + Prokaryotic_community_composition
  Respiration ~ Temperature + Nutrient + Ciliate_diversity + total_biomass + Ciliate_community_composition + Prokaryotic_community_composition
  
  # Covariance between Ciliate_diversity and Ciliate_community_composition
  # Ciliate_diversity ~~ Ciliate_community_composition
'


# Fit the model
fit.3 <- sem(model.3, data = sem.new)

# Summary of the model with fit measures
summary(fit.3, fit.measures = TRUE, standardized = TRUE)

model.4 <- '
  Ciliate_diversity ~ Temperature + Nutrient + Temperature : Nutrient 
  Prokaryote_diversity ~ Temperature + Temperature : Nutrient 
  Ciliate_community_composition ~ Temperature + Nutrient + Ciliate_diversity + Temperature : Nutrient + Prokaryotic_community_composition 
  Prokaryotic_community_composition ~ Temperature + Nutrient + Prokaryote_diversity + Temperature : Nutrient + Ciliate_community_composition  
  total_biomass ~ Temperature + Nutrient + Temperature : Nutrient + Prokaryote_diversity + Ciliate_diversity + Ciliate_community_composition + Prokaryotic_community_composition
  Respiration ~ Nutrient + Ciliate_diversity + total_biomass
'

# Fit the model
fit.4 <- sem(model.4, data = sem.new)

# Summary of the model with fit measures
summary(fit.4, fit.measures = TRUE, standardized = TRUE)

model.5 <- '
  Ciliate_diversity ~ Temperature + Nutrient + Temperature : Nutrient 
  Prokaryote_diversity ~ Temperature + Temperature : Nutrient 
  Ciliate_community_composition ~ Temperature + Nutrient + Ciliate_diversity + Temperature : Nutrient + Prokaryotic_community_composition 
  Prokaryotic_community_composition ~ Temperature + Prokaryote_diversity + Temperature : Nutrient + Ciliate_community_composition  
  total_biomass ~ Temperature + Nutrient + Temperature : Nutrient + Prokaryote_diversity + Ciliate_community_composition + Prokaryotic_community_composition
  Respiration ~ Nutrient + Ciliate_diversity + total_biomass
'

# Fit the model
fit.5 <- sem(model.5, data = sem.new)

# Summary of the model with fit measures
summary(fit.5, fit.measures = TRUE, standardized = TRUE)

model.6 <- '
  Ciliate_diversity ~ Temperature + Nutrient + Temperature : Nutrient 
  Prokaryote_diversity ~ Temperature + Temperature : Nutrient 
  Ciliate_community_composition ~ Temperature + Nutrient + Ciliate_diversity + Temperature : Nutrient + Prokaryotic_community_composition 
  Prokaryotic_community_composition ~ Prokaryote_diversity + Temperature : Nutrient + Ciliate_community_composition  
  total_biomass ~ Temperature + Nutrient + Temperature : Nutrient + Prokaryote_diversity + Ciliate_diversity + Ciliate_community_composition + Prokaryotic_community_composition
  Respiration ~ Temperature + Nutrient + Ciliate_diversity + total_biomass + Ciliate_community_composition + Prokaryotic_community_composition
'

# Fit the model
fit.6 <- sem(model.6, data = sem.new)

# Summary of the model with fit measures
summary(fit.6, fit.measures = TRUE, standardized = TRUE)


model.7 <- '
  Ciliate_diversity ~ Temperature + Nutrient + Temperature : Nutrient 
  Prokaryote_diversity ~ Temperature + Temperature : Nutrient 
  Ciliate_community_composition ~ Temperature + Nutrient + Ciliate_diversity + Temperature : Nutrient + Prokaryotic_community_composition 
  Prokaryotic_community_composition ~ Temperature + Nutrient + Prokaryote_diversity + Temperature : Nutrient + Ciliate_community_composition  
  biomass_ciliate ~ Temperature + Nutrient + Ciliate_community_composition + Ciliate_diversity + Temperature : Nutrient 
  biomass_prokaryote ~ Temperature + Nutrient + Ciliate_diversity + Prokaryote_diversity + biomass_ciliate + Temperature : Nutrient + Prokaryotic_community_composition
  Respiration ~ Nutrient + Ciliate_diversity + biomass_prokaryote + biomass_ciliate
'
# Fit the model
fit.7 <- sem(model.7, data = sem.new)


model.8 <- '
  Ciliate_diversity ~ Temperature + Nutrient + Temperature : Nutrient 
  Prokaryote_diversity ~ Temperature + Temperature : Nutrient 
  Ciliate_community_composition ~ Temperature + Nutrient + Ciliate_diversity + Temperature : Nutrient + Prokaryotic_community_composition 
  Prokaryotic_community_composition ~ Temperature + Nutrient + Prokaryote_diversity + Temperature : Nutrient + Ciliate_community_composition  
  biomass_ciliate ~ Temperature + Nutrient + Ciliate_community_composition + Temperature : Nutrient 
  biomass_prokaryote ~ Temperature + Nutrient + biomass_ciliate + Temperature : Nutrient 
  Respiration ~ Nutrient + Ciliate_diversity + biomass_prokaryote + biomass_ciliate
'

# Fit the model
fit.8 <- sem(model.8, data = sem.new)
# Summary of the model with fit measures
summary(fit.8, fit.measures = TRUE, standardized = TRUE)

model.9 <- '
  Ciliate_diversity ~ a1 * Temperature + a2 * Nutrient
  Prokaryote_diversity ~ b1 * Temperature + b2 * Temperature:Nutrient
  Ciliate_community_composition ~ c1 * Temperature + c2 * Nutrient + c3 * Ciliate_diversity + c4 * Prokaryotic_community_composition
  Prokaryotic_community_composition ~ d1 * Prokaryote_diversity + d2 * Ciliate_community_composition  
  total_biomass ~ e1 * Temperature + e2 * Nutrient + e3 * Temperature:Nutrient + e4 * Prokaryote_diversity + e5 * Ciliate_community_composition + e6 * Prokaryotic_community_composition
  Respiration ~ f1 * Nutrient + f2 * Ciliate_diversity + f3 * total_biomass

  indirect_Ciliate_diversity_to_Ciliate_community_composition_via_total_biomass := e5 * f3
  indirect_Ciliate_diversity_to_Respiration_via_total_biomass := e5 * f3
  indirect_Ciliate_diversity_to_Respiration_via_Ciliate_community_composition := c3 * e5 * f3
  indirect_Prokaryote_diversity_to_Ciliate_community_composition_via_Prokaryotic_community_composition := d1 * c4
  indirect_Prokaryote_diversity_to_Respiration_via_total_biomass := e4 * f3
  indirect_Temperature_to_Respiration_via_Ciliate_diversity := a1 * f2
  indirect_Temperature_to_Respiration_via_total_biomass := e1 * f3
  indirect_Nutrient_to_Respiration_via_Ciliate_diversity := a2 * f2
  indirect_Nutrient_to_Respiration_via_total_biomass := e2 * f3
  indirect_Temperature_Nutrient_to_Respiration_via_Prokaryote_diversity := b2 * f3
  indirect_Temperature_Nutrient_to_Respiration_via_total_biomass := e3 * f3
  indirect_Ciliate_community_composition_to_Respiration_via_total_biomass := e5 * f3
  indirect_Prokaryotic_community_composition_to_Respiration_via_total_biomass := e6 * f3
'

# Fit the model
fit.9 <- sem(model.9, data = sem.new)

# Summary of the model with fit measures
summary(fit.9, fit.measures = TRUE, standardized = TRUE)


# Extract parameter estimates, including defined indirect effects
indirect_effects1 <- parameterEstimates(fit.9, standardized = TRUE) %>%
  subset(op == ":=") # Only defined indirect effects are extracted

# Display indirect effects
print(indirect_effects1)


model.10 <- '
  Ciliate_diversity ~ a1 * Temperature + a2 * Nutrient + a3 * Temperature:Nutrient
  Prokaryote_diversity ~ b1 * Temperature + b2 * Temperature:Nutrient
  Ciliate_community_composition ~ c1 * Temperature + c2 * Nutrient + c3 * Ciliate_diversity + c4 * Temperature:Nutrient + c5 * Prokaryotic_community_composition
  Prokaryotic_community_composition ~ d1 * Prokaryote_diversity + d2 * Ciliate_community_composition  
  total_biomass ~ e1 * Temperature + e2 * Nutrient + e3 * Temperature:Nutrient + e4 * Prokaryote_diversity + e5 * Ciliate_community_composition + e6 * Prokaryotic_community_composition + e7 * Ciliate_diversity
  Respiration ~ f1 * Nutrient + f2 * Ciliate_diversity + f3 * total_biomass + f4 * Prokaryote_diversity

  indirect_Ciliate_diversity_to_Ciliate_community_composition_via_total_biomass := e7 * e5
  indirect_Ciliate_diversity_to_Respiration_via_total_biomass := e7 * f3
  indirect_Ciliate_diversity_to_Respiration_via_Ciliate_community_composition := c3 * e5 * f3
  indirect_Prokaryote_diversity_to_Ciliate_community_composition_via_Prokaryotic_community_composition := d1 * c5
  indirect_Prokaryote_diversity_to_Respiration_via_total_biomass := e4 * f3
  indirect_Temperature_to_Respiration_via_Ciliate_diversity := a1 * f2
  indirect_Temperature_to_Respiration_via_total_biomass := e1 * f3
  indirect_Nutrient_to_Respiration_via_Ciliate_diversity := a2 * f2
  indirect_Nutrient_to_Respiration_via_total_biomass := e2 * f3
  indirect_Temperature_Nutrient_to_Respiration_via_Prokaryote_diversity := b2 * f4
  indirect_Temperature_Nutrient_to_Respiration_via_total_biomass := e3 * f3
  indirect_Ciliate_community_composition_to_Respiration_via_total_biomass := e5 * f3
  indirect_Prokaryotic_community_composition_to_Respiration_via_total_biomass := e6 * f3
'

# Fit the model
fit.10 <- sem(model.10, data = sem.new)

# Summary of the model with fit measures
summary(fit.10, fit.measures = TRUE, standardized = TRUE)

# Extract parameter estimates, including defined indirect effects
indirect_effects <- parameterEstimates(fit.10, standardized = TRUE) %>%
  subset(op == ":=") # Only defined indirect effects are extracted

# Display indirect effects
print(indirect_effects)



# Extract fit indices for each model
models <- list(fit.1, fit.2, fit.3, fit.4, fit.5, fit.6, fit.7, fit.8, fit.9, fit.10)
model_names <- c("fit.1", "fit.2", "fit.3", "fit.4", "fit.5", "fit.6", "fit.7", "fit.8", "fit.9", "fit.10")

fit_measures <- lapply(models, function(model) {
  fitMeasures(model, c("chisq", "aic", "bic", "cfi", "tli", "rmsea", "srmr", "pvalue"))
})

# Convert to data frame
fit_comparison <- data.frame(
  Model = model_names,
  Chi_Square = sapply(fit_measures, function(x) x["chisq"]),
  AIC = sapply(fit_measures, function(x) x["aic"]),
  BIC = sapply(fit_measures, function(x) x["bic"]),
  CFI = sapply(fit_measures, function(x) x["cfi"]),
  TLI = sapply(fit_measures, function(x) x["tli"]),
  RMSEA = sapply(fit_measures, function(x) x["rmsea"]),
  SRMR = sapply(fit_measures, function(x) x["srmr"]),
  P_Value = sapply(fit_measures, function(x) x["pvalue"])
)

# Display the comparison table
print(fit_comparison)

fit_comparison_sorted <- fit_comparison[order(-fit_comparison$AIC), ]

# Print the sorted dataframe
print(fit_comparison_sorted)


# Load the package
library(writexl)

# Save the dataframe to an Excel file
write_xlsx(fit_comparison_sorted, "fit_comparison.xlsx")

# The file will be saved in your working directory
getwd()



library(lavaanPlot)
lavaanPlot(
  model = fit.9,
  node_options = list(shape = "box", fontname = "Helvetica"),
  edge_options = list(color = "grey"),
  coefs = FALSE
)

# Plot with significant paths and labels
lavaanPlot(
  model = fit.9,
  node_options = list(shape = "box", fontname = "Helvetica"),
  edge_options = list(color = "grey"),
  coefs = TRUE,
  sig = .05
)

lavaanPlot(model = fit.9, coefs = TRUE, stand = TRUE, sig = 0.05) #standardized regression paths, showing only paths with p<= .05



###############################################


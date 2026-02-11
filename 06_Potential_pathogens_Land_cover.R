##################################################################################################
# Code to perform analysis for                                                               #####
# Shapiro et al.                                                                             #####
# Gut microbiomes of the little free-tailed bat (Mops pumilus) show high                     #####
# prevalence of potential bacterial pathogens but limited responses to land cover            #####
#                                                                                            #####
#                                                                                            #####
# Script 06 : Spatial models to determine the association between the                        #####
#             presence / absence of each potential pathogen at the roost level               #####
#             and land-cover metrics                                                         #####
#                                                                                            #####
#                                                                                            #####
# Script tested for R version 4.1.2                                                          #####
#                                                                                            #####
##################################################################################################


##################################################################################################
#### Load libraries
##################################################################################################
library(tidyverse)
library(lme4)
library(MuMIn)
##################################################################################################


##################################################################################################
#### --- Data --- ####
##################################################################################################
# Data with land cover metrics : 
land.cover.roosts <- read.csv("land.cover.roosts.csv", header=TRUE)

# Data files created in previous scripts 
rich.metadat.geo.roost # Created in Script04
bart.asvs # Created in Script05
rick.asvs # Created in Script05
camp.asvs # Created in Script05
salm.asvs # Created in Script05
mycop.asvs # Created in Script05
##################################################################################################


##################################################################################################
#### --- Data preparation --- ####
##################################################################################################
# Filter the full dataframe with land cover data to keep only the potential pathogens 
# Then create a data frame for each pathogen Genus

# Bartonella:
bart.prev.geo.roost <- rich.metadat.geo.roost %>%
  ungroup() %>%
  dplyr::select(Bat_ID:Roost, PropSav120:PropSug2000,
                matches(bart.asvs$asv)) %>%
  mutate(bart.pres = across(starts_with("seq")) %>% rowSums) %>%
  mutate(bart.pres = ifelse(bart.pres > 0, 1, 0)) %>%
  group_by(Site) %>%
  mutate(bart.prev = sum(bart.pres)) %>%
  mutate(bart.binom = ifelse(bart.prev > 0, 1, 0)) %>%
  dplyr::select(Site, Roost, Roost.Size,bart.prev, bart.binom, PropSav120:PropSug2000) %>%
  distinct()

# Rickettsia:
rick.prev.geo.roost <- rich.metadat.geo.roost %>%
  ungroup() %>%
  dplyr::select(Bat_ID:Roost, PropSav120:PropSug2000,
                matches(rick.asvs$asv)) %>%
  mutate(rick.pres = across(starts_with("seq")) %>% rowSums) %>%
  mutate(rick.pres = ifelse(rick.pres > 0, 1, 0)) %>%
  group_by(Site) %>%
  mutate(rick.prev = sum(rick.pres)) %>%
  mutate(rick.binom = ifelse(rick.prev > 0, 1, 0)) %>%
  dplyr::select(Site, Roost, Roost.Size,rick.prev, rick.binom,PropSav120:PropSug2000) %>%
  distinct()  

# Salmonella:
salm.prev.geo.roost <- rich.metadat.geo.roost %>%
  ungroup() %>%
  dplyr::select(Bat_ID:Roost, PropSav120:PropSug2000,
                matches(salm.asvs$asv)) %>%
  mutate(salm.pres = across(starts_with("seq")) %>% rowSums) %>%
  mutate(salm.pres = ifelse(salm.pres > 0, 1, 0)) %>%
  group_by(Site) %>%
  mutate(salm.prev = sum(salm.pres)) %>%
  mutate(salm.binom = ifelse(salm.prev > 0, 1, 0)) %>%
  dplyr::select(Site, Roost, Roost.Size,salm.prev,salm.binom, PropSav120:PropSug2000) %>%
  distinct()  

# Campylobacter:
camp.prev.geo.roost <- rich.metadat.geo.roost %>%
  ungroup() %>%
  dplyr::select(Bat_ID:Roost, PropSav120:PropSug2000,
                matches(camp.asvs$asv)) %>%
  mutate(camp.pres = across(starts_with("seq")) %>% rowSums) %>%
  mutate(camp.pres = ifelse(camp.pres > 0, 1, 0)) %>%
  group_by(Site) %>%
  mutate(camp.prev = sum(camp.pres)) %>%
  mutate(camp.binom = ifelse(camp.prev > 0, 1, 0)) %>%
  dplyr::select(Site, Roost, Roost.Size,camp.prev,camp.binom, PropSav120:PropSug2000) %>%
  distinct()  

# Mycoplasma:
mycop.prev.geo.roost <- rich.metadat.geo.roost %>%
  ungroup() %>%
  dplyr::select(Bat_ID:Roost, PropSav120:PropSug2000,
                matches(mycop.asvs$asv)) %>%
  mutate(mycop.pres = across(starts_with("seq")) %>% rowSums) %>%
  mutate(mycop.pres = ifelse(mycop.pres > 0, 1, 0)) %>%
  group_by(Site) %>%
  mutate(mycop.prev = sum(mycop.pres)) %>%
  mutate(mycop.binom = ifelse(mycop.prev > 0, 1, 0)) %>%
  dplyr::select(Site, Roost, Roost.Size,mycop.prev,mycop.binom, PropSav120:PropSug2000) %>%
  distinct() 
##################################################################################################


##################################################################################################
#### --- Modeling and model selection --- ####
##################################################################################################

# For each potential pathogen, create the list of candidate models, 
# run them, and perform model selection

# Bartonella:
# List of candidate models:
bart.glms.bin <- list(
  propsav120 = glm(bart.binom ~ scale(PropSav120), data=bart.prev.geo.roost, family=binomial),
  propsug120 = glm(bart.binom ~ scale(PropSug120), data=bart.prev.geo.roost, family=binomial),
  proprur120 = glm(bart.binom ~ scale(PropRur120), data=bart.prev.geo.roost, family=binomial),
  edge120 = glm(bart.binom ~ scale(SavEdge120), data=bart.prev.geo.roost, family=binomial),
  propsav2000 = glm(bart.binom ~ scale(PropSav2000), data=bart.prev.geo.roost, family=binomial),
  propsug2000 = glm(bart.binom ~ scale(PropSug2000), data=bart.prev.geo.roost, family=binomial),
  proprur2000 = glm(bart.binom ~ scale(PropRur2000), data=bart.prev.geo.roost, family=binomial),
  propwat2000 = glm(bart.binom ~ scale(PropWat2000), data=bart.prev.geo.roost, family=binomial),
  edge2000 = glm(bart.binom ~ scale(SavEdge2000), data=bart.prev.geo.roost, family=binomial),
  split2000 = glm(bart.binom ~ scale(SavSplit2000), data=bart.prev.geo.roost, family=binomial),
  null = glm(bart.binom ~ 1, data=bart.prev.geo.roost, family=binomial)
  
)

# Model selection:
bart.modsel.bin <- model.sel(bart.glms.bin)

#Summary of top model:
summary(bart.glms.bin$split2000)

# Confidence intervals of beta coefficients in the top model:
confint(bart.glms.bin$split2000)

# Repeat for competing models:
summary(bart.glms.bin$split2000)
confint(bart.glms.bin$split2000)

summary(bart.glms.bin$propsav120)
confint(bart.glms.bin$propsav120)

summary(bart.glms.bin$propsav2000)
confint(bart.glms.bin$propsav2000)


# Rickettsia
# List of candidate models:
rick.glms.bin <- list(
  propsav120 = glm(rick.binom ~ scale(PropSav120), data=rick.prev.geo.roost, family=binomial),
  propsug120 = glm(rick.binom ~ scale(PropSug120), data=rick.prev.geo.roost, family=binomial),
  proprur120 = glm(rick.binom ~ scale(PropRur120), data=rick.prev.geo.roost, family=binomial),
  edge120 = glm(rick.binom ~ scale(SavEdge120), data=rick.prev.geo.roost, family=binomial),
  propsav2000 = glm(rick.binom ~ scale(PropSav2000), data=rick.prev.geo.roost, family=binomial),
  propsug2000 = glm(rick.binom ~ scale(PropSug2000), data=rick.prev.geo.roost, family=binomial),
  proprur2000 = glm(rick.binom ~ scale(PropRur2000), data=rick.prev.geo.roost, family=binomial),
  propwat2000 = glm(rick.binom ~ scale(PropWat2000), data=rick.prev.geo.roost, family=binomial),
  edge2000 = glm(rick.binom ~ scale(SavEdge2000), data=rick.prev.geo.roost, family=binomial),
  split2000 = glm(rick.binom ~ scale(SavSplit2000), data=rick.prev.geo.roost, family=binomial),
  null = glm(rick.binom ~ 1, dat=rick.prev.geo.roost, family=binomial)
  
)

# Model selection:
rick.modsel.bin <- model.sel(rick.glms.bin)
# Null is best


# Salmonella
# List of candidate models:
salm.glms.bin <- list(
  propsav120 = glm(salm.binom ~ scale(PropSav120), data=salm.prev.geo.roost, family=binomial),
  propsug120 = glm(salm.binom ~ scale(PropSug120), data=salm.prev.geo.roost, family=binomial),
  proprur120 = glm(salm.binom ~ scale(PropRur120), data=salm.prev.geo.roost, family=binomial),
  edge120 = glm(salm.binom ~ scale(SavEdge120), data=salm.prev.geo.roost, family=binomial),
  propsav2000 = glm(salm.binom ~ scale(PropSav2000), data=salm.prev.geo.roost, family=binomial),
  propsug2000 = glm(salm.binom ~ scale(PropSug2000), data=salm.prev.geo.roost, family=binomial),
  proprur2000 = glm(salm.binom ~ scale(PropRur2000), data=salm.prev.geo.roost, family=binomial),
  propwat2000 = glm(salm.binom ~ scale(PropWat2000), data=salm.prev.geo.roost, family=binomial),
  edge2000 = glm(salm.binom ~ scale(SavEdge2000), data=salm.prev.geo.roost, family=binomial),
  split2000 = glm(salm.binom ~ scale(SavSplit2000), data=salm.prev.geo.roost, family=binomial),
  null = glm(salm.binom ~ 1, dat=salm.prev.geo.roost, family=binomial)
  
)

# Model selection
salm.modsel.bin <- model.sel(salm.glms.bin)
salm.modsel.bin
# Null is best


# Campylobacter
# List of candidate models:
camp.glms.bin <- list(
  propsav120 = glm(camp.binom ~ scale(PropSav120), data=camp.prev.geo.roost, family=binomial),
  propsug120 = glm(camp.binom ~ scale(PropSug120), data=camp.prev.geo.roost, family=binomial),
  proprur120 = glm(camp.binom ~ scale(PropRur120), data=camp.prev.geo.roost, family=binomial),
  edge120 = glm(camp.binom ~ scale(SavEdge120), data=camp.prev.geo.roost, family=binomial),
  propsav2000 = glm(camp.binom ~ scale(PropSav2000), data=camp.prev.geo.roost, family=binomial),
  propsug2000 = glm(camp.binom ~ scale(PropSug2000), data=camp.prev.geo.roost, family=binomial),
  proprur2000 = glm(camp.binom ~ scale(PropRur2000), data=camp.prev.geo.roost, family=binomial),
  propwat2000 = glm(camp.binom ~ scale(PropWat2000), data=camp.prev.geo.roost, family=binomial),
  edge2000 = glm(camp.binom ~ scale(SavEdge2000), data=camp.prev.geo.roost, family=binomial),
  split2000 = glm(camp.binom ~ scale(SavSplit2000), data=camp.prev.geo.roost, family=binomial),
  null = glm(camp.binom ~ 1, data=camp.prev.geo.roost, family=binomial)
  
)

# Model selection:
camp.modsel.bin <- model.sel(camp.glms.bin)
# Best is null


# Mycoplasma
# List of candidate models:
mycop.glms.bin <- list(
  propsav120 = glm(mycop.binom ~ scale(PropSav120), data=mycop.prev.geo.roost, family=binomial),
  propsug120 = glm(mycop.binom ~ scale(PropSug120), data=mycop.prev.geo.roost, family=binomial),
  proprur120 = glm(mycop.binom ~ scale(PropRur120), data=mycop.prev.geo.roost, family=binomial),
  edge120 = glm(mycop.binom ~ scale(SavEdge120), data=mycop.prev.geo.roost, family=binomial),
  propsav2000 = glm(mycop.binom ~ scale(PropSav2000), data=mycop.prev.geo.roost, family=binomial),
  propsug2000 = glm(mycop.binom ~ scale(PropSug2000), data=mycop.prev.geo.roost, family=binomial),
  proprur2000 = glm(mycop.binom ~ scale(PropRur2000), data=mycop.prev.geo.roost, family=binomial),
  propwat2000 = glm(mycop.binom ~ scale(PropWat2000), data=mycop.prev.geo.roost, family=binomial),
  edge2000 = glm(mycop.binom ~ scale(SavEdge2000), data=mycop.prev.geo.roost, family=binomial),
  split2000 = glm(mycop.binom ~ scale(SavSplit2000), data=mycop.prev.geo.roost, family=binomial),
  null = glm(mycop.binom ~ 1, data=mycop.prev.geo.roost, family=binomial)
  
)

# Model selection:
mycop.modsel.bin <- model.sel(mycop.glms.bin)
# Null is best
##################################################################################################


##################################################################################################
#### --- End of script, proceed to Script 06 --- ####
##################################################################################################
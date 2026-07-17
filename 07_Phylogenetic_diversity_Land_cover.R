##################################################################################################
# Code to perform analysis for                                                               #####
# Shapiro et al.                                                                             #####
# Gut microbiomes of the little free-tailed bat (Mops pumilus) show high                     #####
# prevalence of potential bacterial pathogens but limited responses to land cover            #####
#                                                                                            #####
#                                                                                            #####
# Script 07 : Spatial models to determine the association between phylogenetic diversity     #####
#             and land-cover metrics                                                         #####
#                                                                                            #####
#                                                                                            #####
# Script tested for R version 4.1.2                                                          #####
#                                                                                            #####
##################################################################################################


##################################################################################################
#### --- Load libraries --- ####
##################################################################################################
library(tidyverse)
library(lme4)
library(MuMIn)
library(lmerTest)
##################################################################################################

##################################################################################################
#### --- Data ---####
##################################################################################################
# Data created in other scripts:
rich.metadat # Created in Script 03

phylo_div # Created in Script 05

# Join phylogenetic diversity data to landsdcape metric data
phylometadat.roost <- phylo_div %>%
  rownames_to_column(var="Bat_ID") %>%
  mutate(Bat_ID = as.integer(Bat_ID)) %>%
  left_join(rich.metadat.geo.roost)


##################################################################################################
#### --- Modeling and model selection --- ####
##################################################################################################
# Create the list of candidate models:
bact.phylo.div <- list(
  propsug2000 = lmer(PD ~ scale(PropSug2000)+ (1|Site), data=phylometadat.roost),
  proprur2000 = lmer(PD ~ scale(PropRur2000)+ (1|Site), data=phylometadat.roost),
  propwat2000 = lmer(PD ~ scale(PropWat2000)+ (1|Site), data=phylometadat.roost),
  edge2000 = lmer(PD ~ scale(SavEdge2000)+ (1|Site), data=phylometadat.roost),
  split2000 = lmer(PD ~ scale(SavSplit2000)+ (1|Site), data=phylometadat.roost),
  propsav120 = lmer(PD ~ scale(PropSav120)+ (1|Site), data=phylometadat.roost),
  propsug120 = lmer(PD ~ scale(PropSug120)+ (1|Site), data=phylometadat.roost),
  proprur120 = lmer(PD ~ scale(PropRur120)+ (1|Site), data=phylometadat.roost),
  null = lmer(PD ~ 1 + (1|Site), data=phylometadat.roost)
  
)

# Run models and model selection :
bact.phylo.modsel <- model.sel(bact.phylo.div)

# Best model: 
summary(bact.phylo.div$propwat2000)

lmm <- lmer(PD ~ scale(PropWat2000)+ (1|Site), data=phylometadat.roost)
summary(lmm)
anova(lmm)

# Confidence interval of coefficient for proportion water
confint(bact.phylo.div$propwat2000)
# Not significant, no further analysis
##################################################################################################
#### --- End of script, proceed to Script 08 --- ####
##################################################################################################
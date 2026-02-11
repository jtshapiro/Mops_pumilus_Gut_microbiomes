##################################################################################################
# Code to perform analysis for                                                               #####
# Shapiro et al.                                                                             #####
# Gut microbiomes of the little free-tailed bat (Mops pumilus) show high                     #####
# prevalence of potential bacterial pathogens but limited responses to land cover            #####
#                                                                                            #####
#                                                                                            #####
# Script 04 : Spatial models to determine the association between ASV richness and           #####
#            land-cover metrics                                                              #####
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
##################################################################################################

##################################################################################################
#### --- Data ---####
##################################################################################################
# Data with land cover metrics : 
land.cover.roosts <- read.csv("land.cover.roosts.csv", header=TRUE)

# Data created in other scripts:
rich.metadat # Created in Script 03

# Combine the two dataframes:
rich.metadat.geo.roost <- rich.metadat %>%
  filter(Roost=="Yes") %>%
  left_join(., land.cover.roosts)
##################################################################################################

##################################################################################################
#### --- Modeling and model selection --- ####
##################################################################################################

# Create the list of candidate models:
bact.glmms <- list(
  propsav2000 = glmer(asv.richness ~ scale(PropSav2000)+ (1|Site), data=rich.metadat.geo.roost, family=poisson),
  propsug2000 = glmer(asv.richness ~ scale(PropSug2000)+ (1|Site), data=rich.metadat.geo.roost, family=poisson),
  proprur2000 = glmer(asv.richness ~ scale(PropRur2000)+ (1|Site), data=rich.metadat.geo.roost, family=poisson),
  propwat2000 = glmer(asv.richness ~ scale(PropWat2000)+ (1|Site), data=rich.metadat.geo.roost, family=poisson),
  edge2000 = glmer(asv.richness ~ scale(SavEdge2000)+ (1|Site), data=rich.metadat.geo.roost, family=poisson),
  split2000 = glmer(asv.richness ~ scale(SavSplit2000)+ (1|Site), data=rich.metadat.geo.roost, family=poisson),
  propsav120 = glmer(asv.richness ~ scale(PropSav120)+ (1|Site), data=rich.metadat.geo.roost, family=poisson),
  propsug120 = glmer(asv.richness ~ scale(PropSug120)+ (1|Site), data=rich.metadat.geo.roost, family=poisson),
  proprur120 = glmer(asv.richness ~ scale(PropRur120)+ (1|Site), data=rich.metadat.geo.roost, family=poisson),
  edge120 = glmer(asv.richness ~ scale(SavEdge120)+ (1|Site), data=rich.metadat.geo.roost, family=poisson),
  null = glmer(asv.richness ~ 1 + (1|Site), data=rich.metadat.geo.roost, family=poisson)
  
)

# Run models and model selection :
bact.modsel <- model.sel(bact.glmms)
# Null model ranks first, no further analysis
##################################################################################################

##################################################################################################
#### --- End of script, proceed to Script 05 --- ####
##################################################################################################
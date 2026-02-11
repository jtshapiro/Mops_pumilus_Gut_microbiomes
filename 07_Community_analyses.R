##################################################################################################
# Code to perform analysis for                                                               #####
# Shapiro et al.                                                                             #####
# Gut microbiomes of the little free-tailed bat (Mops pumilus) show high                     #####
# prevalence of potential bacterial pathogens but limited responses to land cover            #####
#                                                                                            #####
#                                                                                            #####
# Script 07 : Bacterial community analysis                                                   #####
#                                                                                            #####
#                                                                                            #####
#                                                                                            #####
# Script tested for R version 4.1.2                                                          #####
#                                                                                            #####
#                                                                                            #####
##################################################################################################


##################################################################################################
#### --- Load libraries --- #####
##################################################################################################
library(tidyverse)
library(vegan)
library(geosphere)
library(corrplot)
library(viridis)
library(ggpubr)
##################################################################################################


##################################################################################################
#### --- Data --- ####
##################################################################################################

# Data files created in previous scripts 
tax.filtered2 # Created in Script01
pres.abs.asv.bact # Created in Script 01
rich.metadat.geo.roost # Created in Script04
##################################################################################################


##################################################################################################
#### --- Community analysis set-up --- ####
##################################################################################################

# Create a matrix:
pres.abs.asv.bact.1 <- pres.abs.asv.bact %>%
  column_to_rownames(., var="sample.id")

# Calculate Jaccard distance on the matrix
jac.mat <-vegdist(pres.abs.asv.bact.1, Type = "Jaccard", binary = TRUE) 
# Creates a similarity matrix (Jaccard is the type), based on presence/absence

# Perform NMDS
bact.mds <- monoMDS(jac.mat)

stressplot(bact.mds)

# Extract the first two dimensions from the NMDS
Nmds1.lc<-bact.mds$points[,1] #Extracts dimension 1 as a vector
Nmds2.lc<-bact.mds$points[,2] #Extracts dimension 2 as a vector

Nmds.lc=cbind(Nmds1.lc,Nmds2.lc) #Makes a new sheet with dimension1 and 2
Nmds.lc<-as.data.frame(Nmds.lc) #Creates a new dataframe with dimensions 1 and 2

##################################################################################################
#### --- NMDS plots by host characteristics --- ####
##################################################################################################

# Plots of NMDS (Supplementary Figure 5)
# Visualization based on sex:
sample.MDS.sex <-ggplot(data = Nmds.lc, aes(y=Nmds2.lc, x=Nmds1.lc, color = rich.metadat$Sex)) +
  geom_point(size = 2.2) + 
  scale_color_viridis(name= "Sex",discrete=TRUE) +
  xlab('NMDS Axis 1') +
  ylab("NMDS Axis 2") +
  theme_classic() +
  theme(legend.title=element_text(size=16),
        legend.text=element_text(size=15))

sample.MDS.sex

# Based on age
sample.MDS.age <-ggplot(data = Nmds.lc, aes(y=Nmds2.lc, x=Nmds1.lc, color = rich.metadat$Age)) +
  geom_point(size = 2.2) + 
  scale_color_viridis(name= "Age",discrete=TRUE) +
  xlab('NMDS Axis 1') +
  ylab("NMDS Axis 2") +
  theme_classic() +
  theme(legend.title=element_text(size=16),
        legend.text=element_text(size=15))

sample.MDS.age

# Based on reproductive condition:
sample.MDS.repro <-ggplot(data = Nmds.lc, aes(y=Nmds2.lc, x=Nmds1.lc, color = rich.metadat$ReproStatus)) +
  geom_point(size = 2.2) + 
  scale_color_viridis(name= "Repro. cond.",discrete=TRUE) +
  xlab('NMDS Axis 1') +
  ylab("NMDS Axis 2") +
  theme_classic() +
  theme(legend.title=element_text(size=16),
        legend.text=element_text(size=15))

sample.MDS.repro

# Group all NMDS plots together (Supplementary Figure 5)
nmds.plot <- ggarrange(sample.MDS.sex, sample.MDS.age, sample.MDS.repro, 
                       ncol = 3, nrow = 1, labels = 
                         c("A.","B.", "C."),
                       align='hv',
                       vjust=0, hjust=1) +
  theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"))

# Optional: save figure
ggsave(nmds.plot, filename="nmds.plot.png", width = 17, height =5)
##################################################################################################


##################################################################################################
#### --- Set-up for db-RDA ---####
##################################################################################################

# Create distance matrix (response variable)
RDA.bact.vars <- rich.metadat.geo.roost %>%
  dplyr::select(Bat_ID, starts_with('seq') ) %>%
  column_to_rownames(var="Bat_ID")

pre.RDA.environment.vars=rich.metadat.geo.roost %>%
  dplyr::select(Bat_ID, Site, Sex, Age, ReproStatus, starts_with("Prop"),
                starts_with("Sav"), Latitude, Longitude)

# Correlation analysis and plot between variables to select variable to be used in db-RDA
# Set-up:
rda.corr.dat <- pre.RDA.environment.vars %>%
  select(-BatID, -Site, -Sex, -Age, -ReproStatus, -Latitude, -Longitude) %>%
  distinct()

# Scale
rda.corr.dat.scaled <- rda.corr.dat %>%
  mutate(across(everything(~scale(.x))))

# Correlation:
rda.corr <- cor(rda.corr.dat.scaled)

# Plot (Supplementary Figure 2)
corrplot(rda.corr, method = 'number', type="lower")

# Final environmental variable data frame:  
RDA.environment.vars=rich.metadat.geo.roost %>%
  dplyr::select(Bat_ID, Site, Sex, Age, ReproStatus, starts_with("Prop"),
                starts_with("Sav"), Latitude, Longitude) %>%
  dplyr::select(-SavShape120, -SavShape2000, -SavSplit120, -PropWat120, -PropSav2000, -SavEdge120)

# Make categorical variables factors:
RDA.environment.vars$Sex=as.factor(RDA.environment.vars$Sex) 
RDA.environment.vars$Age=as.factor(RDA.environment.vars$Age)
RDA.environment.vars$ReproStatus=as.factor(RDA.environment.vars$ReproStatus)

# To control for spatial autocorrelation, condition by distance between sites
# To do this:
# Create a distance matrix :
distance_matrix <- distm(rich.metadat.geo.roost[,5060:5059], fun = distVincentyEllipsoid)

# Now will use pcnm to get the distance between all the sites
pcnmr=pcnm(distance_matrix)

# Create object that will be used to condition in the db-RDA
pcnmv=pcnmr$vectors

# Run db-RDA
# Condition based on distances between sites:
cap.all.cond <- dbrda(RDA.bact.vars ~ Sex + Age + ReproStatus + 
                        scale(PropSav120) + scale(PropSug120) + scale(PropRur120)
                      + scale(PropSug2000) + scale(PropRur2000)
                      + scale(PropWat2000) + scale(SavEdge2000) + scale(SavSplit2000)
                      + Condition(pcnmv), 
                      data = RDA.environment.vars,
                      distance = "jaccard") 

# Contribution by term, axis, margin
anova.cca(cap.all.cond)
anova.cca(cap.all.cond, by = "term")
anova.cca(cap.all.cond, by= "axis")
anova.cca(cap.all.cond, by = "margin")

RsquareAdj(cap.all.cond)

summary(cap.all.cond)


# Forward selection of db-RDA
RDA.bact.forward.cond=ordiR2step(dbrda(RDA.bact.vars~1 + Condition(pcnmv), data = RDA.environment.vars,distance="jaccard"),
                                 scope=formula(cap.all.cond), scale=TRUE,Pin=.05, R2scope=F, direction="forward",pstep=1000) 

# Calculate contribution of term, axis, margins
term.aov <- anova.cca(RDA.bact.forward.cond, by="term")
axis.aov <- anova.cca(RDA.bact.forward.cond, by="axis")
margin.aov <- anova.cca(RDA.bact.forward.cond, by="margin")
RsquareAdj(RDA.bact.forward.cond) # R sq = 0.06684044, R sq adjusted = 0.02582627
plot(RDA.bact.forward.cond)

# Plotting the final db-RDA : 
# To make prettier figure, will scale in the data frame and change column names
RDA.environment.vars.pretty <- RDA.environment.vars %>%
  mutate(Bat_ID==as.character(Bat_ID)) %>%
  mutate_if(is.numeric, scale) %>%
  rename("Sug_cov_landscape" = "PropSug2000",
         "Sav_cov_fine" = "PropSav120",
         "Rur_cov_fine"="PropRur120",
         "Sav_split_landscape" = "SavSplit2000")

# Rerun the db-RDA
rda.bact.forward.selected.cond.pretty <- dbrda(RDA.bact.vars ~ Rur_cov_fine + Sug_cov_landscape + 
                                                 + Sav_cov_fine 
                                               +  Sav_split_landscape
                                               + Condition(pcnmv), data = RDA.environment.vars.pretty,
                                               distance = "jaccard") 

# Make and save the plot (Figure 4)
jpeg("pretty.dbrda.vegan20251115.jpeg", width = 8.5, height = 7.5, units = "in", res=900, quality=100)
plot(rda.bact.forward.selected.cond.pretty, type="n")
text(rda.bact.forward.selected.cond.pretty, dis="bp", cex=1.4, col="blue")
points(rda.bact.forward.selected.cond.pretty, display="sites",pch=3, col="red",  cex=1.3)
dev.off()
##################################################################################################


##################################################################################################
#### --- End of script, all analyses complete ---####
##################################################################################################

##################################################################################################
# Code to perform analysis for                                                               #####
# Shapiro et al.                                                                             #####
# Gut microbiomes of the little free-tailed bat (Mops pumilus) show high                     #####
# prevalence of potential bacterial pathogens but limited responses to land cover            #####
#                                                                                            #####
#                                                                                            #####
# Script 02: Calculations of identified ASVs per taxonomic level;                             #####
#           recreates Figure 2 (barplot);                                                    #####
#           calculates rarefaction curve for sampling completeness.                          #####
#                                                                                            #####
#                                                                                            #####
# Script tested for R version 4.1.2                                                          #####
#                                                                                            #####
##################################################################################################


##################################################################################################
#### --- Load libraries --- ####
##################################################################################################

# Install libraries if necessary (first time using script) 
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(vegan)
##################################################################################################


##################################################################################################
#### --- Data --- ####
##################################################################################################
# Data frames created in the previous script:
tax.filtered2 # Created in script 01_DADA2_DataProcessing
rel.abund # Created in script 01_DADA2_DataProcessing
tot.reads.samp.bact # Created in script 01_DADA2_DataProcessing
##################################################################################################


##################################################################################################
#### --- Analysis of ASV taxonomy --- ####
##################################################################################################

# Percent of ASVs ID'ed to each level of classification
# Phylum
phyl <- tax.filtered2 %>%
  filter(Phylum != 'NA')
phyl.id <- nrow(phyl)/nrow(tax.filtered2)*100

# Class
cl <- tax.filtered2 %>%
  filter(Class != 'NA')
cl.id <- nrow(cl)/nrow(tax.filtered2)*100

# Order 
or <- tax.filtered2 %>%
  filter(Order != 'NA')
ord.id <- nrow(or)/nrow(tax.filtered2)*100

# Family
fam <- tax.filtered2 %>%
  filter(Family != 'NA')
fam.id <- nrow(fam)/nrow(tax.filtered2)*100

# Genus
gen <- tax.filtered2 %>%
  filter(Genus != 'NA')
gen.id <- nrow(gen)/nrow(tax.filtered2)*100

# Species
spe <- tax.filtered2 %>%
  filter(Species != 'NA')
sp.id <- nrow(spe)/nrow(tax.filtered2)*100

# Relative abundance for Figure 2 (barplot)
# By Phylum
# Start with : 
rel.abund # Created in script 01_DADA2_DataProcessing

# Pivot the table
rel1 <- t(rel.abund) %>%
  as.data.frame() %>%
  row_to_names(row_number = 1) %>%
  mutate(across(where(is.character), as.numeric))

# Included ASVs
incl.phy <- tax.filtered2 %>%
  dplyr::select(asv, Phylum)

rel2 <- rel1 %>%
  rownames_to_column(., var="asv") %>%
  left_join(., incl.phy) %>%
  dplyr::select(Phylum, everything()) %>%
  dplyr::select(-asv)

# Relabel unknown phyla
rel3 <- rel2 %>%
  mutate(Phylum = ifelse(Phylum=="NA", "Unknown phyla", Phylum)) %>%
  mutate(Phylum=replace_na(Phylum, "Unknown phyla")) %>%
  group_by(Phylum) %>%
  summarise(across(where(is.numeric), ~ sum(.x))) %>%
  pivot_longer(cols = where(is.numeric))

# Count bats per phyla and rename those found in <20 bats as "Minor phyla"
bats.per.phyl <- rel3 %>%
  ungroup() %>%
  filter(value > 0) %>%
  group_by(Phylum) %>%
  mutate(n.bats=n()) %>%
  dplyr::select(-name, -value) %>%
  distinct() %>%
  mutate(phyl.graph = ifelse(n.bats >= 20, Phylum, "Minor phyla"))

# Join the data sets
rel4 <- rel3 %>%
  left_join(., bats.per.phyl, by="Phylum") %>%
  dplyr::select(-Phylum, -n.bats) %>%
  rename(Phylum=phyl.graph)

# Set Phyla as factors and order them for the figure legend
rel4$Phylum<- factor(rel4$Phylum, levels = c("Actinobacteriota","Bacteroidota",
                                             "Firmicutes","Fusobacteriota","Proteobacteria",
                                             "Minor phyla"))

# Make the bar plot for Phyla
phyl.bar <- rel4 %>% 
  ggplot(aes(x=name, y=value, fill=Phylum))+
  geom_bar(position="stack", stat="identity") +
  scale_fill_viridis_d(begin=1, end=0) + 
  ylab("Relative abundance") +
  xlab(" ") +
  theme_classic() +
  coord_cartesian(ylim = c(0,1), expand=F)+
  theme(axis.title=element_text(size=16)) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_text(size=13)) +
  theme(legend.title=element_text(size=16),
        legend.text=element_text(size=14))


# Repeat for Family
rel.abund

# Pivot the table
rel.f1 <- t(rel.abund) %>%
  as.data.frame() %>%
  row_to_names(row_number = 1) %>%
  mutate(across(where(is.character), as.numeric))

incl.fam <- tax.filtered2 %>%
  dplyr::select(asv, Family)

rel.f2 <- rel.f1 %>%
  rownames_to_column(., var="asv") %>%
  left_join(., incl.fam) %>%
  dplyr::select(Family, everything()) %>%
  dplyr::select(-asv)


rel.f3 <- rel.f2 %>%
  mutate(Family = ifelse(Family=="NA", "Unknown families", Family)) %>%
  mutate(Family=replace_na(Family, "Unknown families")) %>%
  group_by(Family) %>%
  summarise(across(where(is.numeric), ~ sum(.x))) %>%
  pivot_longer(cols = where(is.numeric))

bats.per.fam <- rel.f3 %>%
  ungroup() %>%
  filter(value > 0) %>%
  group_by(Family) %>%
  mutate(n.bats=n()) %>%
  dplyr::select(-name, -value) %>%
  distinct() %>%
  mutate(fam.graph = ifelse(n.bats >= 40, Family, "Minor families")) %>%
  mutate(fam.graph = ifelse(Family=="NA", "Unknown family",fam.graph))

rel.f4 <- rel.f3 %>%
  left_join(., bats.per.fam, by="Family") %>%
  dplyr::select(-Family, -n.bats) %>%
  rename(Family=fam.graph)

rel.f4$Family <- factor(rel.f4$Family, levels = c("Aerococcaceae","Clostridiaceae",
                                                  "Corynebacteriaceae","Enterobacteriaceae","Enterococcaceae",
                                                  "Erwiniaceae","Morganellaceae","Mycoplasmataceae","Pasteurellaceae",       
                                                  "Peptostreptococcaceae","Streptococcaceae"     
                                                  ,"Weeksellaceae","Yersiniaceae","Minor families","Unknown families"))
# Create bar plot for family
fam.bar <- rel.f4 %>% 
  ggplot(aes(x=name, y=value, fill=Family))+
  geom_bar(position="stack", stat="identity") +
  scale_fill_viridis_d(begin=1, end=0) + 
  ylab("Relative abundance") +
  xlab("Bat") +
  theme_classic() +
  coord_cartesian(ylim = c(0,1), expand=F)+
  theme(axis.title=element_text(size=16)) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_text(size=13)) +
  theme(legend.title=element_text(size=16),
        legend.text=element_text(size=14))

# Combine Phylum and Family barplots into one figure (recreates Figure 2)
rel.abund.bars <- ggarrange(phyl.bar, fam.bar, 
                            ncol = 1, nrow = 2, labels = 
                              c("A.","B."),
                            align='hv',
                            vjust=0, hjust=1) +
  theme(plot.margin = margin(3,3,3,3, "cm"))

# Optional - save figure
ggsave(rel.abund.bars, filename="yourfilename.png", 
       width= 40, height =30, units = "cm")
##################################################################################################


##################################################################################################
##### --- Species accumulation curves ---####
##################################################################################################

# Start with: 
tot.reads.samp.bact # Created in script 01_DADA2_DataProcessing

# Create an abundance matrix
comm.mat.abund <- tot.reads.samp.bact %>%
  column_to_rownames(., var="sample.id") %>%
  select(-tot.reads)

# Rarefy by individuals (reads per ASV) using abundance matrix above 
# and plot the curve (Supplementary Figure 2)
plot(specaccum(comm.mat.abund, method="rarefaction"), ylim=c(0,2500),
     xlab = "Number bats sampled", ylab = "Number ASVs", cex.axis=0.7, cex.lab=0.8)

# To save figure : 
png(filename="spec.accum.plot.png",width = 12, height=9, units = "cm", res=600)
plot(specaccum(comm.mat.abund, method="rarefaction"), ylim=c(0,2500),
     xlab = "Number bats sampled", ylab = "Number ASVs", cex.axis=0.7, cex.lab=0.8)
dev.off()
##################################################################################################


##################################################################################################
#### --- End of script, proceed to Script 03 --- ####
##################################################################################################
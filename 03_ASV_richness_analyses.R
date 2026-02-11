##################################################################################################
# Code to perform analysis for                                                               #####
# Shapiro et al.                                                                             #####
# Gut microbiomes of the little free-tailed bat (Mops pumilus) show high                     #####
# prevalence of potential bacterial pathogens but limited responses to land cover            #####
#                                                                                            #####
#                                                                                            #####
# Script 03 : T-tests / ANOVAS, correlations analysis for differences in                      #####
#            ASV species richness based on bat characteristics                               #####
#            (age, sex, reproductive condition,body condition)                               #####
#                                                                                            #####
#                                                                                            #####
# Script tested for R version 4.1.2                                                          #####
#                                                                                            #####
##################################################################################################


##################################################################################################
#### --- Load libraries --- ####
##################################################################################################
library(tidyverse)
library(ggpubr)
library(Hmisc)
##################################################################################################


##################################################################################################
#### --- Data --- ####
##################################################################################################

# Read in the data
# NOTE : bat_metadata also read in Script 01
bat_metadata <- read.csv("bat_metadata.csv", header=TRUE)

# Dataframe created in Script 01:
asv.rich
##################################################################################################


##################################################################################################
#### --- Data preparation --- ####
##################################################################################################

# Select relevant columns from asv.rich:
asv.rich2 <- asv.rich %>%
  select(sample.id, asv.richness) %>%
  # Rename sample.id column to Bat_ID 
  rename(Bat_ID = sample.id)

# Join the ASV richness data to the metadata:
rich.metadat <- bat_metadata %>%
  left_join(., asv.rich2)
##################################################################################################


##################################################################################################
##### --- Test for differences in ASV richness based on host characteristics --- ####
##################################################################################################

# Sex (t-test (2 groups))
Sex.test <- t.test(asv.richness ~ Sex, data = rich.metadat)
Sex.test
summary(Sex.test)

# Age (t-test (2 groups))
Age.test <- t.test(asv.richness ~ Age, data = rich.metadat)
Age.test
summary(Age.test)

# Reproductive condition (ANOVA (4 groups))
Repro.anov <- aov(asv.richness ~ ReproStatus, data = rich.metadat)
Repro.anov
summary(Repro.anov)


# Create violin plots for each characteristic (Supplementary Figure 3)
# Sex
Sex.plot<-ggplot(data = rich.metadat, 
                 aes(x = Sex, y = asv.richness)) + 
  labs(x='Sex', y= 'ASV Richness', size=18) + 
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
  stat_summary(fun=mean, geom="point",shape=23, size=2, color="red", fill="red") +
  theme(axis.text = element_text(size=16)) +
  theme_classic()

# Age
Age.plot<-ggplot(data = rich.metadat, 
                 aes(x = Age, y = asv.richness)) + 
  labs(x='Age', y= 'ASV Richness', size=18) + 
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
  stat_summary(fun=mean, geom="point",shape=23, size=2, color="red", fill="red") +
  theme(axis.text = element_text(size=16)) +
  theme_classic()

# Reproductive condition
Repro.plot<-ggplot(data = rich.metadat, 
                   aes(x = ReproStatus, y = asv.richness)) + 
  labs(x='Reproductive Status', y= 'ASV Richness', size=18) + 
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
  stat_summary(fun=mean, geom="point",shape=23, size=2, color="red", fill="red") +
  theme(axis.text = element_text(size=16)) +
  theme_classic()

# Combine all the plots into a single figure (Supplementary Figure 3)
si.anova.plot <- ggarrange(Sex.plot, Age.plot, Repro.plot, 
                           ncol = 1, nrow = 3, labels = 
                             c("A.","B.", "C."),
                           align='hv',
                           vjust=0, hjust=1) +
  theme(plot.margin = margin(1,1,1,1, "cm"))

# Optional: save the plot
ggsave(si.anova.plot, filename="si.anova.plot.png", 
       width= 20, height =25, units = "cm")
##################################################################################################


##################################################################################################
#### --- Correlation between body condition and ASV richness --- ####
##################################################################################################

# First calculate body condition (Weight / Forearm length)
rich.metadat2 <- rich.metadat %>%
  mutate(Weight=as.numeric(Weight),
         Forearm=as.numeric(Forearm)) %>%
  mutate(body.cond = Weight/Forearm)


# Create a scatterplot (Supplementary Figure 4)
asv.bc.plot <- ggplot(rich.metadat2, aes(x=body.cond, y=asv.richness)) + 
  geom_point(size=1.5) +
  labs(x='Body condition', y= 'ASV Richness', size=24) +
  theme(axis.text = element_text(size=20)) +
  theme_classic()

# Optional: save plot
ggsave(asv.bc.plot, filename = "asv.corr.plot.png", 
       width= 15, height =12, units = "cm" )

# Use the rcorr function in Hmisc to calculate correlation between 
# body condition and ASV richness with p values
rcor<-rcorr(rich.metadat3$body.cond, rich.metadat3$asv.richness, type = "spearman")
rcor
##################################################################################################


##################################################################################################
#### --- End of script, proceed to Script 04 --- #### 
##################################################################################################
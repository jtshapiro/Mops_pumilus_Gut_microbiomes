##################################################################################################
# Code to perform analysis for                                                               #####
# Shapiro et al.                                                                             #####
# Gut microbiomes of the little free-tailed bat (Mops pumilus) show high                     #####
# prevalence of potential bacterial pathogens but limited responses to land cover            #####
#                                                                                            #####
#                                                                                            #####
# Script 04 : Calculate phylogenetic diversity                                               #####
#             Perform T-tests / Kruskal-Wallis, correlations analysis                        #####
#             for differences in phylogenetic diversity based on bat characteristics         #####
#            (age, sex, reproductive condition,body condition)                               #####
#                                                                                            #####
#                                                                                            #####
# Script tested for R version 4.1.2                                                          #####
#                                                                                            #####
##################################################################################################


##################################################################################################
#### --- Load libraries --- ####
##################################################################################################
library(devtools)
library(tidyverse)
library(msa)
library(ape)
library(ips)
library(phangorn)
library(DECIPHER)
library("seqinr")
library(phytools)
library(picante)
library(MuMIn)
library(lmerTest)
library(ggpubr)
library(Hmisc)
##################################################################################################


##################################################################################################
#### --- Data --- ####
##################################################################################################
# Start with the fasta file of sequences
big_fasta <- read.fasta(file = "dada_seqs.fa")
length(big_fasta)
names(big_fasta)

# Filter the fasta file to only include sequences that passed all filtering steps
# and are included in the final dadta set
Interesting_OTU_fasta <- big_fasta[names(big_fasta) %in% incl.asv$asv]
Interesting_OTU_fasta <- Interesting_OTU_fasta[order(names(Interesting_OTU_fasta))] #Ensuring they are in the same order

##################################################################################################
#### --- Sequence alignment and phylogenetic tree construction --- ####
##################################################################################################
# Alignn sequences using the ClustalW method
Interesting_OTU_alignment <- msa(Interesting_OTU_fasta, method = "ClustalW", type = "dna", order = "input")

class(Interesting_OTU_alignment)

is(Interesting_OTU_alignment)

# Make a new object in case anything goes wrong
Interesting_OTU_alignment_test <- Interesting_OTU_alignment

# Assign AAMultipleAlignment class to the object
class(Interesting_OTU_alignment_test) <- "AAMultipleAlignment"

# Convert to a seqinr alignment
Interesting_OTU_alignment_seqinr <- msaConvert(Interesting_OTU_alignment_test, 
                                               type = "seqinr::alignment")

# Calculate distances between all sequences
Interesting_OTU_alignment_dist <- seqinr::dist.alignment(Interesting_OTU_alignment_seqinr, 
                                                         matrix = "identity")

# Create a neighbor-joining tree
tree_nj <- nj(Interesting_OTU_alignment_dist)


# Calculate phylogenetic diversity
# Not: Input -- community matrix, the tree, no root needed
phylo_div <- pd(pres.abs.asv.bact.1, tree_nj, include.root=FALSE)

##################################################################################################
#### --- Differences in phylogenetic diversity  --- ####
##################################################################################################
# Join phylo_div to metadata
phylometadat <- phylo_div %>%
  rownames_to_column(var="Bat_ID") %>%
  mutate(Bat_ID = as.integer(Bat_ID)) %>%
  left_join(rich.metadat3)

# Age, Sex, Repro, Body condition
# Sex, t-test
Sex.test.phylo <- t.test(PD ~ Sex, data = phylometadat)
Sex.test.phylo
#summary(Sex.anov)

# Age, t-test
Age.test.phylo <- t.test(PD ~ Age, data = phylometadat)
Age.test.phylo
summary(Age.test.phylo)

# Reproductive condition -- Kruskal-Wallis (4 groups)
Repro.kw.phylo <- kruskal.test(PD ~ ReproStatus, data = phylometadat)
Repro.kw.phylo
summary(Repro.kw.phylo)


# Boxplots for SI
Sex.plot.phylo<-ggplot(data = phylometadat, 
                       aes(x = Sex, y = PD)) + 
  labs(x='Sex', y= 'Phylogenetic diversity', size=18) + 
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
  stat_summary(fun=mean, geom="point",shape=23, size=2, color="red", fill="red") +
  theme(axis.text = element_text(size=16)) +
  theme_classic()

Age.plot.phylo<-ggplot(data = phylometadat, 
                       aes(x = Age, y = PD)) + 
  labs(x='Age', y= 'Phylogenetic diversity', size=18) + 
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
  stat_summary(fun=mean, geom="point",shape=23, size=2, color="red", fill="red") +
  theme(axis.text = element_text(size=16)) +
  theme_classic()

Repro.plot.phylo<-ggplot(data = phylometadat, 
                         aes(x = ReproStatus, y = PD)) + 
  labs(x='Reproductive Status', y= 'Phylogenetic diversity', size=18) + 
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
  stat_summary(fun=mean, geom="point",shape=23, size=2, color="red", fill="red") +
  theme(axis.text = element_text(size=16)) +
  theme_classic()

# Join all plots together as panels in single plot
si.anova.plot.phylo <- ggarrange(Sex.plot.phylo, Age.plot.phylo, Repro.plot.phylo, 
                                 ncol = 1, nrow = 3, labels = 
                                   c("A.","B.", "C."),
                                 align='hv',
                                 vjust=0, hjust=1) +
  theme(plot.margin = margin(1,1,1,1, "cm"))

ggsave(si.anova.plot.phylo, filename="si.anova.plot.phylo.png", 
       width= 20, height =25, units = "cm")

# Correlation between body condition and phylogenetic diversity
# Scatterplot 
phydiv.bc.plot <- ggplot(phylometadat, aes(x=body.cond, y=PD)) + 
  geom_point(size=1.5) +
  labs(x='Body condition', y= 'Phylogenetic diversity', size=24) +
  theme(axis.text = element_text(size=20)) +
  theme_classic()

ggsave(phydiv.bc.plot, filename = "phydiv.corr.plot.png", 
       width= 15, height =12, units = "cm" )

# Calculate correlation with Hmisc package
rcor.phyl<-rcorr(phylometadat$body.cond, phylometadat$PD, type = "pearson")
rcor.phyl


##################################################################################################
#### --- End of script, proceed to Script 05 --- #### 
##################################################################################################
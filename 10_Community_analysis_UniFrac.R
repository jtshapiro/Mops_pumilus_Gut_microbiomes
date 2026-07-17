##################################################################################################
# Code to perform analysis for                                                               #####
# Shapiro et al.                                                                             #####
# Gut microbiomes of the little free-tailed bat (Mops pumilus) show high                     #####
# prevalence of potential bacterial pathogens but limited responses to land cover            #####
#                                                                                            #####
#                                                                                            #####
# Script 10 : Bacterial community analysis, based on UniFrac distances                       #####
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
library(picante)
library(ape)
library(vegan)
library(viridis)
library(ggpubr)
##################################################################################################

##################################################################################################
#### --- Data --- ####
##################################################################################################
# Data files created in previous scripts 
tree_nj # Created in Script05
RDA.bact.vars # Created in Script 08
##################################################################################################


##################################################################################################
#### --- Community analysis set-up - UniFrac distance calculation --- ####
##################################################################################################
# First need to root the phylogenetic tree
# Following this: https://john-quensen.com/r/unifrac-and-tree-roots/
# Custom function
pick_new_outgroup <- function(tree.unrooted){
  require("magrittr")
  require("data.table")
  require("ape") # ape::Ntip
  # tablify parts of tree that we need.
  treeDT <-
    cbind(
      data.table(tree.unrooted$edge),
      data.table(length = tree.unrooted$edge.length)
    )[1:Ntip(tree.unrooted)] %>%
    cbind(data.table(id = tree.unrooted$tip.label))
  # Take the longest terminal branch as outgroup
  new.outgroup <- treeDT[which.max(length)]$id
  return(new.outgroup) }


out.group <- pick_new_outgroup(tree_nj)
out.group ## [1] "seq2456"

rooted.tree <- ape::root(tree_nj, outgroup=out.group, resolve.root=TRUE)
phy_tree(rooted.tree) <- rooted.tree
phy_tree(rooted.tree)


# Create phyloseq OTU table
ps.out.table<- phyloseq::otu_table(RDA.bact.vars, taxa_are_rows = F)

# Join with the tree
ps.object <- merge_phyloseq(ps.out.table, rooted.tree)

# Calculate Unifrac distance
unifrac_dist <- phyloseq::distance(ps.object, method = "unifrac")
##################################################################################################

##################################################################################################
#### --- Run community analyses ---####
##################################################################################################

# Redo db-RDA on unifrac data
cap.all.cond.unifrac <- dbrda(unifrac_dist ~ Sex + Age + ReproStatus + 
                                scale(PropSav120) + scale(PropSug120) + scale(PropRur120)
                              + scale(PropSug2000) + scale(PropRur2000)
                              + scale(PropWat2000) + scale(SavEdge2000) + scale(SavSplit2000)
                              + Condition(pcnmv), 
                              data = RDA.environment.vars,
                              distance = "unifrac") 
anova.cca(cap.all.cond.unifrac)
anova.cca(cap.all.cond.unifrac, by = "term")
anova.cca(cap.all.cond.unifrac, by= "axis")
anova.cca(cap.all.cond.unifrac, by = "margin")
plot(cap.all.cond)
RsquareAdj(cap.all.cond.unifrac)


#Forward selection of dbdrda
RDA.bact.forward.cond.unifrac=ordiR2step(dbrda(unifrac_dist~1 + Condition(pcnmv), data = RDA.environment.vars,distance="unifrac"),
                                         scope=formula(cap.all.cond.unifrac), scale=TRUE,Pin=.05, R2scope=F, direction="forward",pstep=1000) 

term.aov.unifrac <- anova.cca(RDA.bact.forward.cond.unifrac, by="term")
axis.aov.unifrac <- anova.cca(RDA.bact.forward.cond.unifrac, by="axis")
margin.aov.unifrac <- anova.cca(RDA.bact.forward.cond.unifrac, by="margin")
RsquareAdj(RDA.bact.forward.cond.unifrac) # R sq = 0.03270771, R sq adjusted = 0.01322571
plot(RDA.bact.forward.cond)

# To make prettier figure, will scale in the 
RDA.environment.vars.pretty <- RDA.environment.vars %>%
  mutate_if(is.numeric, scale)

rda.bact.forward.selected.cond <- dbrda(RDA.bact.vars ~  Condition(pcnmv) + scale(PropRur120) + scale(PropSug2000) +
                                          scale(PropSav120) + scale(SavSplit2000), data = RDA.environment.vars, distance = "jaccard") 

term.aov <- anova.cca(rda.bact.forward.selected.cond, by="term")
axis.aov <- anova.cca(rda.bact.forward.selected.cond, by="axis")
margin.aov <- anova.cca(rda.bact.forward.selected.cond, by="margin")
RsquareAdj(rda.bact.forward.selected.cond) # R sq = 0.06684044, R sq adjusted = 0.02582627
plot(rda.bact.forward.selected.cond)
##################################################################################################
#### --- End of script, all analyses complete ---####
##################################################################################################
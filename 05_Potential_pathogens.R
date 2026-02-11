##################################################################################################
# Code to perform analysis for                                                               #####
# Shapiro et al.                                                                             #####
# Gut microbiomes of the little free-tailed bat (Mops pumilus) show high                     #####
# prevalence of potential bacterial pathogens but limited responses to land cover            #####
#                                                                                            #####
#                                                                                            #####
# Script 05 : Identification of potential pathogens, prevalence, and analysis of             ##### 
#             associations between host characteristics and their prevalences                #####
#                                                                                            #####
#                                                                                            #####
# Script tested for R version 4.1.2                                                          #####
#                                                                                            #####
##################################################################################################


##################################################################################################
#### --- Load libraries ---####
##################################################################################################
library(tidyverse)
library(ggpubr)
##################################################################################################

##################################################################################################
#### --- Data ---####
##################################################################################################

# Dataframes created in previous scripts:
tax.filtered2 # Created in Script01
rich.metadat # Created in Script 03
rich.metadat.geo.roost # Created in Script04
##################################################################################################

##################################################################################################
#### --- Potential pathogen presence ---####
##################################################################################################

# Filter taxonomy list for the potential pathogens
pathogens <- tax.filtered2 %>%
  filter(Genus=='Bartonella'| Genus== 'Rickettsia'| Genus=='Salmonella'|Genus=='Campylobacter'
         |Genus=='Mycoplasma'|Genus=='Mycobacterium')

# Create dataframe of ASVs for each potential pathogen:
# Bartonella
bart.asvs <- tax.filtered2 %>%
  filter(Genus=='Bartonella')

# Rickettsia
rick.asvs <- tax.filtered2 %>%
  filter(Genus=='Rickettsia')

# Salmonella
salm.asvs <- tax.filtered2 %>%
  filter(Genus=='Salmonella')

# Campylobacter
camp.asvs <- tax.filtered2 %>%
  filter(Genus=='Campylobacter')

# Mycoplasma
mycop.asvs <- tax.filtered2 %>%
  filter(Genus=='Mycoplasma')


# Filter the full data set to keep just the pathogens and create a dataframe for each potential pathogen Genus
# Bartonella
bart.metadat.geo.roost <- rich.metadat.geo.roost %>%
  ungroup() %>%
  dplyr::select(Bat_ID:asv.richness, PropSav120:PropSug2000,
                matches(bart.asvs$asv)) %>%
  mutate(bart.pres = across(starts_with("seq")) %>% rowSums) %>%
  mutate(bart.pres = ifelse(bart.pres > 0, 1, 0))  

# Rickettsia
rick.metadat.geo.roost <- rich.metadat.geo.roost %>%
  ungroup() %>%
  dplyr::select(Bat_ID:asv.richness, PropSav120:PropSug2000,
                matches(rick.asvs$asv)) %>%
  mutate(rick.pres = across(starts_with("seq")) %>% rowSums) %>%
  mutate(rick.pres = ifelse(rick.pres > 0, 1, 0))  

# Salmonella
salm.metadat.geo.roost <- rich.metadat.geo.roost %>%
  ungroup() %>%
  dplyr::select(Bat_ID:asv.richness, PropSav120:PropSug2000,
                matches(salm.asvs$asv)) %>%
  mutate(salm.pres = across(starts_with("seq")) %>% rowSums) %>%
  mutate(salm.pres = ifelse(salm.pres > 0, 1, 0))  

# Campylobacter
camp.metadat.geo.roost <- rich.metadat.geo.roost %>%
  ungroup() %>%
  dplyr::select(Bat_ID:asv.richness, PropSav120:PropSug2000,
                matches(camp.asvs$asv)) %>%
  mutate(camp.pres = across(starts_with("seq")) %>% rowSums) %>%
  mutate(camp.pres = ifelse(camp.pres > 0, 1, 0))  
##################################################################################################


##################################################################################################
#### --- Host characteristics and potential pathogen prevalence --- ####
##################################################################################################
# Reformat the data to run t-tests (ASV richness) and then chi squared test for sex, age, and reproductive condition
# Bartonella
bart.chi.dat <- rich.metadat %>%
  ungroup() %>%
  dplyr::select(Bat_ID:asv.richness,
                matches(bart.asvs$asv)) %>%
  mutate(bart.pres = across(starts_with("seq")) %>% rowSums) %>%
  mutate(bart.pres = ifelse(bart.pres > 0, 1, 0)) %>%
  mutate(bart.factor = ifelse(bart.pres > 0, "Present", "Absent"))

# T-test: ASV richness in bats with / without Bartonella
t.test(bart.chi.dat$asv.richness ~ bart.chi.dat$bart.factor)

# Violin plot
bart.rich.plot <-ggplot(data = bart.chi.dat, 
                        aes(x = bart.factor, y = asv.richness)) + 
  labs(x='Bartonella', y= 'ASV Richness', size=18) + 
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
  stat_summary(fun=mean, geom="point",shape=23, size=4, color="red", fill="red") +
  theme_classic()

# Chi-squared test : 
# Sex
bart.chi.sex <- bart.chi.dat %>%
  group_by(Sex, bart.pres) %>%
  summarise(n=n()) %>%
  pivot_wider(names_from = Sex, values_from = n) %>%
  column_to_rownames(., var="bart.pres")

chisq.test(bart.chi.sex)

# Age
bart.chi.age <- bart.chi.dat %>%
  group_by(Age, bart.pres) %>%
  summarise(n=n()) %>%
  pivot_wider(names_from = Age, values_from = n) %>%
  column_to_rownames(., var="bart.pres") %>%
  mutate(J=ifelse(is.na(J), 0,J))

chisq.test(bart.chi.age)
# Note: No juveniles with Bartonella

# Reproductive condition
bart.chi.repro <- bart.chi.dat %>%
  group_by(ReproStatus, bart.pres) %>%
  summarise(n=n()) %>%
  pivot_wider(names_from = ReproStatus, values_from = n) %>%
  column_to_rownames(., var="bart.pres")

chisq.test(bart.chi.repro)

# Repeat for Rickettsia
rick.chi.dat <- rich.metadat %>%
  ungroup() %>%
  dplyr::select(Bat_ID:asv.richness,
                matches(rick.asvs$asv)) %>%
  mutate(rick.pres = across(starts_with("seq")) %>% rowSums) %>%
  mutate(rick.pres = ifelse(rick.pres > 0, 1, 0)) %>%
  mutate(rick.factor = ifelse(rick.pres > 0, "Present", "Absent"))

# T-test: ASV richness in bats with / without Rickettsia
t.test(rick.chi.dat$asv.richness ~ rick.chi.dat$rick.factor)

# Violin plot
rick.rich.plot <-ggplot(data = rick.chi.dat, 
                        aes(x = rick.factor, y = asv.richness)) + 
  labs(x='Rickettsia', y= 'ASV Richness', size=18) + 
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
  stat_summary(fun=mean, geom="point",shape=23, size=4, color="red", fill="red") +
  theme_classic()

# Chi-squared test:
# Sex
rick.chi.sex <- rick.chi.dat %>%
  group_by(Sex, rick.pres) %>%
  summarise(n=n()) %>%
  pivot_wider(names_from = Sex, values_from = n) %>%
  column_to_rownames(., var="rick.pres")

chisq.test(rick.chi.sex)

# Age
rick.chi.age <- rick.chi.dat %>%
  group_by(Age, rick.pres) %>%
  summarise(n=n()) %>%
  pivot_wider(names_from = Age, values_from = n) %>%
  column_to_rownames(., var="rick.pres")

chisq.test(rick.chi.age)


# Reproductive condition
rick.chi.repro <- rick.chi.dat %>%
  group_by(ReproStatus, rick.pres) %>%
  summarise(n=n()) %>%
  pivot_wider(names_from = ReproStatus, values_from = n) %>%
  column_to_rownames(., var="rick.pres")

chisq.test(rick.chi.repro)

# Repeat for Salmonella
salm.chi.dat <- rich.metadat %>%
  ungroup() %>%
  dplyr::select(Bat_ID:asv.richness,
                matches(salm.asvs$asv)) %>%
  mutate(salm.pres = across(starts_with("seq")) %>% rowSums) %>%
  mutate(salm.pres = ifelse(salm.pres > 0, 1, 0)) %>%
  mutate(salm.factor = ifelse(salm.pres > 0, "Present", "Absent"))

# T-test: ASV richness in bats with / without Salmonella
t.test(salm.chi.dat$asv.richness ~ salm.chi.dat$salm.factor)

# Violin plot
salm.rich.plot <-ggplot(data = salm.chi.dat, 
                        aes(x = salm.factor, y = asv.richness)) + 
  labs(x='Salmonella', y= 'ASV Richness', size=18) + 
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
  stat_summary(fun=mean, geom="point",shape=23, size=4, color="red", fill="red") +
  theme_classic()

# Chi-squared test:
# Sex
salm.chi.sex <- salm.chi.dat %>%
  group_by(Sex, salm.pres) %>%
  summarise(n=n()) %>%
  pivot_wider(names_from = Sex, values_from = n) %>%
  column_to_rownames(., var="salm.pres")

chisq.test(salm.chi.sex)

# Age
salm.chi.age <- salm.chi.dat %>%
  group_by(Age, salm.pres) %>%
  summarise(n=n()) %>%
  pivot_wider(names_from = Age, values_from = n) %>%
  column_to_rownames(., var="salm.pres")

chisq.test(salm.chi.age)

# Reproductive condition
salm.chi.repro <- salm.chi.dat %>%
  group_by(ReproStatus, salm.pres) %>%
  summarise(n=n()) %>%
  pivot_wider(names_from = ReproStatus, values_from = n) %>%
  column_to_rownames(., var="salm.pres")

chisq.test(salm.chi.repro)

# Repeat for Campylobacter
camp.chi.dat <- rich.metadat3 %>%
  ungroup() %>%
  dplyr::select(Bat_ID:asv.richness,
                matches(camp.asvs$asv)) %>%
  mutate(camp.pres = across(starts_with("seq")) %>% rowSums) %>%
  mutate(camp.pres = ifelse(camp.pres > 0, 1, 0))%>%
  mutate(camp.factor = ifelse(camp.pres > 0, "Present", "Absent"))

# T-test: ASV richness in bats with / without Campylobacter
t.test(camp.chi.dat$asv.richness ~ camp.chi.dat$camp.factor)

# Violin plot: 
camp.rich.plot <-ggplot(data = camp.chi.dat, 
                        aes(x = camp.factor, y = asv.richness)) + 
  labs(x='Campylobacter', y= 'ASV Richness', size=18) + 
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
  stat_summary(fun=mean, geom="point",shape=23, size=4, color="red", fill="red") +
  theme_classic()

# Chi-squared test:
# Sex
camp.chi.sex <- camp.chi.dat %>%
  group_by(Sex, camp.pres) %>%
  summarise(n=n()) %>%
  pivot_wider(names_from = Sex, values_from = n) %>%
  column_to_rownames(., var="camp.pres") %>%
  mutate(M=ifelse(is.na(M), 0,M))

chisq.test(camp.chi.sex)

# Age
camp.chi.age <- camp.chi.dat %>%
  group_by(Age, camp.pres) %>%
  summarise(n=n()) %>%
  pivot_wider(names_from = Age, values_from = n) %>%
  column_to_rownames(., var="camp.pres") %>%
  mutate(J=ifelse(is.na(J), 0,J))

chisq.test(camp.chi.age)
# NOTE: No juveniles with Campylobacter

# Reproductive condition:
camp.chi.repro <- camp.chi.dat %>%
  group_by(ReproStatus, camp.pres) %>%
  summarise(n=n()) %>%
  pivot_wider(names_from = ReproStatus, values_from = n) %>%
  column_to_rownames(., var="camp.pres") %>%
  mutate(P=ifelse(is.na(P), 0,P)) %>%
  select(-R)

chisq.test(camp.chi.repro)

# Repeat for Mycoplasma
# Mycoplasma
mycop.chi.dat <- rich.metadat %>%
  ungroup() %>%
  dplyr::select(Bat_ID:asv.richness,
                matches(mycop.asvs$asv)) %>%
  mutate(mycop.pres = across(starts_with("seq")) %>% rowSums) %>%
  mutate(mycop.pres = ifelse(mycop.pres > 0, 1, 0))%>%
  mutate(mycop.factor = ifelse(mycop.pres > 0, "Present", "Absent"))

# T-test: ASV richness in bats with / without Mycoplasma
t.test(mycop.chi.dat$asv.richness ~ mycop.chi.dat$mycop.factor)

# Violin plot:
mycop.rich.plot <-ggplot(data = mycop.chi.dat, 
                         aes(x = mycop.factor, y = asv.richness)) + 
  labs(x='Mycoplasma', y= 'ASV Richness', size=18) + 
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
  stat_summary(fun=mean, geom="point",shape=23, size=4, color="red", fill="red") +
  theme_classic()

# Chi-squared test:
# Sex
mycop.chi.sex <- mycop.chi.dat %>%
  group_by(Sex, mycop.pres) %>%
  summarise(n=n()) %>%
  pivot_wider(names_from = Sex, values_from = n) %>%
  column_to_rownames(., var="mycop.pres") %>%
  mutate(M=ifelse(is.na(M), 0,M))

chisq.test(mycop.chi.sex)

# Age
mycop.chi.age <- mycop.chi.dat %>%
  group_by(Age, mycop.pres) %>%
  summarise(n=n()) %>%
  pivot_wider(names_from = Age, values_from = n) %>%
  column_to_rownames(., var="mycop.pres") %>%
  mutate(J=ifelse(is.na(J), 0,J))

chisq.test(mycop.chi.age)

# Reproductive condition:
mycop.chi.repro <- mycop.chi.dat %>%
  group_by(ReproStatus, mycop.pres) %>%
  summarise(n=n()) %>%
  pivot_wider(names_from = ReproStatus, values_from = n) %>%
  column_to_rownames(., var="mycop.pres") %>%
  mutate(P=ifelse(is.na(P), 0,P)) #%>%
#select(-R)

chisq.test(mycop.chi.repro)


# Combine all the violin plots comparing ASV richness in bats with / without each potential pathogen:
richness.path.plot <- ggarrange(bart.rich.plot, camp.rich.plot, mycop.rich.plot, rick.rich.plot, salm.rich.plot,
                                ncol = 3, nrow = 2, labels = 
                                  c("A.","B.", "C.", "D.", "E."),
                                align='hv',
                                vjust=0, hjust=1) +
  theme(plot.margin = margin(1,1,1,1, "cm"))

# Optional: save
ggsave(richness.path.plot, filename="richness.path.plot.png", width = 10, height =8)
##################################################################################################

##################################################################################################
#### --- Co-occurrence of potential pathogens --- ####
##################################################################################################
# Create dataframes for the bats positive for each potential pathogen and then combine and add up 
# to determine how many co-occur in each individual:

# For Bartonella: 
bart.yes <- bart.chi.dat %>%
  filter(bart.pres==1) %>%
  select(Bat_ID, Site, bart.pres)

# For Rickettsia:
rick.yes <- rick.chi.dat %>%
  filter(rick.pres==1) %>%
  select(Bat_ID, Site, rick.pres)

# For Salmonella:
salm.yes <- salm.chi.dat %>%
  filter(salm.pres==1) %>%
  select(Bat_ID, Site, salm.pres)

# For Campylobacter:
camp.yes <- camp.chi.dat %>%
  filter(camp.pres==1) %>%
  select(Bat_ID, Site, camp.pres)

# For Mycoplasma:
mycop.yes <- mycop.chi.dat %>%
  filter(mycop.pres==1) %>%
  select(Bat_ID, Site, mycop.pres)

# Join all the dataframes and count the number of potential pathogens for each individual bat:
all.path.list <- bart.yes %>%
  full_join(., rick.yes) %>%
  full_join(., salm.yes) %>%
  full_join(., camp.yes) %>%
  full_join(., mycop.yes) %>%
  mutate(across(contains('pres'), replace_na, 0)) %>%
  mutate(n.path = bart.pres + rick.pres + salm.pres + camp.pres + mycop.pres)

co.occur.summ <- all.path.list %>%
  group_by(n.path) %>%
  summarise(n=n())
##################################################################################################


##################################################################################################
#### --- End of script, proceed to Script 06 ---####
##################################################################################################

##################################################################################################
# Code to perform analysis for                                                               #####
# Shapiro et al.                                                                             #####
# Gut microbiomes of the little free-tailed bat (Mops pumilus) show high                     #####
# prevalence of potential bacterial pathogens but limited responses to land cover            #####
#                                                                                            #####
#                                                                                            #####
# Script 01: Initial data processing to set up all downstream analyses, including            #####
#           processing of reads using DADA2, creation of ASVs creation of ASV table,         #####
#           taxonomic identification, data filtering                                         #####
#                                                                                            #####
#                                                                                            #####
# Script tested for R version 4.1.2                                                          #####
#                                                                                            #####
##################################################################################################

##################################################################################################
#### --- Library installation and loading --- ####
##################################################################################################

# Load necessary libraries 
# Install first if necessary (first time only)
#BiocManager::install(version = '3.14')
library(BiocManager)

#BiocManager::install("dada2")
library(dada2)

library("devtools")

#BiocManager::install('phyloseq')
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")

library(janitor)
library(ggpubr)
##################################################################################################


##################################################################################################
#### --- Set directory to wherever the files are located --- ####
main_path <- "~/YourFolder"
##################################################################################################


##################################################################################################
#### --- Processing of reads using DADA2 ---####
##################################################################################################

##################################################################################################
# 1. Matching the sense and antisense reads  ####
##################################################################################################

path <- file.path(main_path, "DADA2_SS")
fns <- list.files(path)
fastqs <- fns[grepl("fastq$", fns)]
fastqs <- sort(fastqs)
fnFs <- fastqs[grepl("_R1.", fastqs)]
fnRs <- fastqs[grepl("_R2.", fastqs)]
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)
match_path <- file.path(path, "matched")
if(!file_test("-d", match_path)) dir.create(match_path)
filtFs <- file.path(match_path, paste0(sample.names, "_F_matched.fastq.gz"))
filtRs <- file.path(match_path, paste0(sample.names, "_R_matched.fastq.gz"))
for(i in seq_along(fnFs)) {
  fastqPairedFilter(c(fnFs[i], fnRs[i]), c(filtFs[i], filtRs[i]),
                    matchIDs=TRUE)
}


path <- file.path(main_path, "DADA2_AS")
fns <- list.files(path)
fastqs <- fns[grepl("fastq$", fns)]
fastqs <- sort(fastqs)
fnFs <- fastqs[grepl("_R2.", fastqs)] # Reverse direction compared to above for
#the "sense reads". In practice the reads are here complement reversed to be in
#the same orientation as the "sense" reads.
fnRs <- fastqs[grepl("_R1.", fastqs)] # See above
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)
match_path <- file.path(path, "matched")
if(!file_test("-d", match_path)) dir.create(match_path)
filtFs <- file.path(match_path, paste0(sample.names, "_F_matched.fastq.gz"))
filtRs <- file.path(match_path, paste0(sample.names, "_R_matched.fastq.gz"))
for(i in seq_along(fnFs)) {
  fastqPairedFilter(c(fnFs[i], fnRs[i]), c(filtFs[i], filtRs[i]),
                    matchIDs=TRUE)
}
##################################################################################################


##################################################################################################
# 2. Filtering and truncation of the reads ####
##################################################################################################

# Filtering of the matched reads. Using DADA2 the paired reads are now filtered 
# (maximum expected errors 2, maximum number of N's (0), truncate length at 215 for R1, and 205 for R2).

# First filtering of the Sense-reads:
path <- file.path(main_path, "DADA2_SS/matched")
fns <- list.files(path)
fastqs <- fns[grepl("matched.fastq", fns)]
fastqs <- sort(fastqs)
fnFs <- fastqs[grepl("_F_", fastqs)]
fnRs <- fastqs[grepl("_R_", fastqs)]

# Can plot a graph of the error rates, learn when to cut the data
# plotQualityProfile(fnRs)
# plotQualityProfile(fnRs)
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)

# Can plot a graph of the error rates, useful to see where to truncate data
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)
plotQualityProfile(fnRs)
plotQualityProfile(fnRs[1])
filt_path <- file.path(path, "filtered")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
for(i in seq_along(fnFs)) {
  fastqPairedFilter(c(fnFs[i], fnRs[i]), c(filtFs[i], filtRs[i]),
                    minLen=10, maxN=0, maxEE=2, truncQ=2, # was 2
                    compress=TRUE, verbose=TRUE)
}

# Next filtering of anti-sense
path <- file.path(main_path, "DADA2_AS/matched") 
fns <- list.files(path)
fastqs <- fns[grepl("matched.fastq", fns)]
fastqs <- sort(fastqs)
fnFs <- fastqs[grepl("_F_", fastqs)]
fnRs <- fastqs[grepl("_R_", fastqs)]
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)

# Can plot a graph of the error rates, useful to see where to truncate data
# plotQualityProfile(fnRs)
# plotQualityProfile(fnRs)
filt_path <- file.path(path, "filtered")
if(!file_test("-d", filt_path)) dir.create(filt_path) 
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
for(i in seq_along(fnFs)) {
  fastqPairedFilter(c(fnFs[i], fnRs[i]), c(filtFs[i], filtRs[i]), 
                    minLen=10, maxN=0, maxEE=2, truncQ=2, 
                    compress=TRUE, verbose=TRUE)
}
##################################################################################################


##################################################################################################
# 3. Actual DADA2 processing of the reads ####
##################################################################################################

# Processing the set of files containing the forward primer in the R1 reads (the sense reads):
filt_path <- file.path(main_path, "DADA2_SS/matched/filtered") 
fns <- list.files(filt_path)
fastqs <- fns[grepl(".fastq.gz", fns)]
fastqs <- sort(fastqs)
fnFs <- fastqs[grepl("_F_", fastqs)]
fnRs <- fastqs[grepl("_R_", fastqs)]
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)
filtFs <- file.path(filt_path, fnFs)
filtRs <- file.path(filt_path, fnRs)

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE,minOverlap = 5)
seqtab_SS <- makeSequenceTable(mergers[names(mergers)])
seqtab.nochim_SS <- removeBimeraDenovo(seqtab_SS, verbose=TRUE)
stSS <- file.path(main_path,"seqtab_SS")
stnsSS <- file.path(main_path,"seqtab.nochim_SS")
saveRDS(seqtab_SS,stSS)
saveRDS(seqtab.nochim_SS,stnsSS)

# Repeat with the anti-sense reads:
filt_path <- file.path(main_path, "DADA2_AS/matched/filtered") 
fns <- list.files(filt_path)
fastqs <- fns[grepl(".fastq.gz", fns)]
fastqs <- sort(fastqs)
fnFs <- fastqs[grepl("_F_", fastqs)]
fnRs <- fastqs[grepl("_R_", fastqs)]
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)
filtFs <- file.path(filt_path, fnFs)
filtRs <- file.path(filt_path, fnRs)

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE,minOverlap = 5)
seqtab_AS <- makeSequenceTable(mergers[names(mergers)])
seqtab.nochim_AS <- removeBimeraDenovo(seqtab_AS, verbose=TRUE)
stAS <- file.path(main_path,"seqtab_AS")
stnsAS <- file.path(main_path,"seqtab.nochim_AS")
saveRDS(seqtab_AS,stAS)
saveRDS(seqtab.nochim_AS,stnsAS)


# Merge the resulting tables of the "sense" and the "antisense" analyses 

# First, define a function for combining two or more tables, collapsing samples with similar names:
sumSequenceTables <- function(table1, table2, ..., orderBy = "abundance") {
  # Combine passed tables into a list
  tables <- list(table1, table2)
  tables <- c(tables, list(...))
  # Validate tables
  if(!(all(sapply(tables, dada2:::is.sequence.table)))) {
    stop("At least two valid sequence tables, and no invalid objects, are expected.")
  }
  sample.names <- rownames(tables[[1]])
  for(i in seq(2, length(tables))) {
    sample.names <- c(sample.names, rownames(tables[[i]]))
  }
  seqs <- unique(c(sapply(tables, colnames), recursive=TRUE))
  sams <- unique(sample.names)
  # Make merged table
  rval <- matrix(0L, nrow=length(sams), ncol=length(seqs))
  rownames(rval) <- sams
  colnames(rval) <- seqs
  for(tab in tables) {
    rval[rownames(tab), colnames(tab)] <- rval[rownames(tab), colnames(tab)] + tab
  }
  # Order columns
  if(!is.null(orderBy)) {
    if(orderBy == "abundance") {
      rval <- rval[,order(colSums(rval), decreasing=TRUE),drop=FALSE]
    } else if(orderBy == "nsamples") {
      rval <- rval[,order(colSums(rval>0), decreasing=TRUE),drop=FALSE]
    }
  }
  rval
}


stAS <- file.path(main_path,"seqtab_AS")
stnsAS <- file.path(main_path,"seqtab.nochim_AS")
stSS <- file.path(main_path,"seqtab_SS")
stnsSS <- file.path(main_path,"seqtab.nochim_SS")
seqtab.nochim_AS <- readRDS(stnsAS)
seqtab.nochim_SS <- readRDS(stnsSS)
seqtab_AS <- readRDS(stAS)
seqtab_SS <- readRDS(stSS)
Plant_sumtable <- sumSequenceTables(seqtab_SS,seqtab_AS)
Plant_nochim_sumtable <- sumSequenceTables(seqtab.nochim_SS,seqtab.nochim_AS)
stBoth <- file.path(main_path,"seqtab_Both")
stnsBoth <- file.path(main_path,"seqtab.nochim_Both")
saveRDS(Plant_sumtable,stBoth)
saveRDS(Plant_nochim_sumtable,stnsBoth)
##################################################################################################


##################################################################################################
# 4. Create a fasta file of sequence variants ####
##################################################################################################
?otu_table
seqs <- colnames(Plant_nochim_sumtable)
otab <- otu_table(Plant_nochim_sumtable, taxa_are_rows=FALSE)
colnames(otab) <- paste0("seq", seq(ncol(otab)))
otab = t(otab)
write.table(seqs, "dada_seqs.txt",quote=FALSE)
write.table(otab, "dada_table.txt",quote=FALSE,sep="\t")

#Can run this in Terminal to create a fasta file with all the sequences and OTU in.
#grep -v '^x' dada_seqs.txt | awk '{print ">seq"$1"\n"$2}' > dada_seqs.fa; rm dada_seqs.txt
##################################################################################################


##################################################################################################
# 5. Assign taxonomy to ASVs  ####
##################################################################################################
genus.species <- assignSpecies(Plant_nochim_sumtable, "silva_species_assignment_v138.fa.gz")
taxa <- assignTaxonomy(Plant_nochim_sumtable, "silva_nr_v138_train_set.fa.gz", multithread=TRUE)
unname(head(taxa))
##################################################################################################


##################################################################################################
# 6. Combine taxonomy with abundance table  ####
##################################################################################################
t.Plant_nochim_sumtable <- t(Plant_nochim_sumtable)#transpose the table
observation <- 1:nrow(t.Plant_nochim_sumtable)
observation <- sub("^", "SeqVar", observation)
#t.seqtab.nochim <- cbind('#OTUID' = rownames(t.seqtab.nochim), t.seqtab.nochim)#Add '#OTUID' to the header (required by biom)
t.Plant_nochim_sumtable <- cbind('#OTUID' = observation, t.Plant_nochim_sumtable)#Add '#OTUID' to the header (required by biom)
Taxonomy <- paste(taxa[,1], taxa[,2], taxa[,3], taxa[,4], taxa[,5], taxa[,6], genus.species[,2], sep=";")
t.Plant_nochim_sumtable.taxa <- data.frame(t.Plant_nochim_sumtable, Taxonomy)
write.table(t.Plant_nochim_sumtable.taxa, "dada2_seq_table.txt", sep='\t', row.names=FALSE, quote=FALSE)
##################################################################################################


##################################################################################################
# 7. Export Taxonomy to check on specific species
##################################################################################################
tax.02032025 <- t.Plant_nochim_sumtable.taxa %>%
  remove_rownames() %>%
  select(X.OTUID, Taxonomy) %>%
  separate(Taxonomy, c('Kingdom','Phylum','Class',
                       'Order','Family','Genus','Species'), sep=";") 
##################################################################################################

##################################################################################################
#### --- End of DADA2 ---####
##################################################################################################

##################################################################################################
#### --- Further processing --- ####
##################################################################################################

##################################################################################################
#### --- Match PCR replicates to bat IDs ---####
##################################################################################################
# Join sample name to PCR ID
tags_to_samp <- read.csv("tags_to_samples.csv", header=T)
otab2 <- t(otab) %>%
  as.data.frame(.)%>%
  rownames_to_column(., var="pcr.id") %>%
  right_join(., tags_to_samp, by="pcr.id") %>%
  select(sample.id, pcr.id, everything())

# Filter out NAs from table
otab3 <- otab2 %>%
  drop_na()
##################################################################################################


##################################################################################################
#### --- ASV filtering --- ####
##################################################################################################

# Identify ASVs found in blanks to filter them out from data set:
blank.otu <- otab3 %>%
  filter(grepl('B|blank', sample.id))

blank.pcr.id <- blank.otu$pcr.id

t.blank.otu <- as.data.frame(t(blank.otu)) %>%
  filter(!row_number() %in% c(1)) %>%
  janitor::row_to_names(1) %>%
  mutate_if(is.character,as.numeric) %>%
  mutate(sumrow = rowSums(.)) %>%
  select(sumrow, everything()) %>%
  filter(sumrow > 0)

blank.seqs <- t.blank.otu %>%
  rownames_to_column(., var = "blank.asv") %>%
  select(blank.asv)

blank.seqs.vec <- as.vector(blank.seqs)

blank.seqs.list <- as.list(blank.seqs)

# Filter blanks out of the ASV table:
otab4 <- otab3 %>%
  filter(!(pcr.id %in% blank.pcr.id)) %>%
  # Remove ASV sequences found in blanks:
  select(-blank.seqs$blank.asv) %>%
  arrange(sample.id)


# Filter to keep only ASV's found in at least 2 PCRs per bat (2/3 triplicates  2/2 duplicates)
# Note : This is going to get a bit comoplicated
otab5 <- otab4 %>%
  select(-pcr.id) 

otab6 <- otab5 %>%
  mutate(across(starts_with('seq'), ~ifelse(.==0,FALSE,TRUE)))

otab6 <- otab5 %>%
  group_by(sample.id) %>%
  arrange(across(starts_with("seq"), .by_group=TRUE)) %>%
  filter(n() == 3) %>%
  arrange(sample.id)

otab7 <- otab6 %>%
  group_by(sample.id) %>%
  
# Separate samples with only 2 replicates:
  otab.rep2 <- otab5 %>% 
  group_by(sample.id) %>% 
  filter(n() == 2)

rep.2.mixed <- otab.rep2 %>%
  group_by(sample.id) %>%
  arrange(., by_group=TRUE) %>%
  mutate(across(starts_with('seq'), ~ifelse(first(.) == 0, 0, lag(.)))) 
# Returns nulls in cases where first row of group is not 0

# Use coalesce to change the NAs (where there are no zeros) back to the original value
rep.2.mixed.1 <- rep.2.mixed %>% 
  group_by(sample.id) %>%
  coalesce(., otab.rep2)

# Mutate each
otab7b <- otab6 %>%
  group_by(sample.id) %>%
  arrange(., by_group=TRUE) %>%
  mutate_each(funs(sum(.==0))) %>%
  arrange(., sample.id)

# Change values depending on how many 0's:
# 3 = all 0's
# 2 = all 0's
# 1 = all NA's
# 0 = all NA's
# Then coalesce with original data set to replace NA's with original read counts

otab7c <- otab7b %>%
  mutate_at(., vars(-sample.id), ~na_if(., 1)) %>%
  mutate_at(., vars(-sample.id), ~na_if(., 0)) 

otab7d <- otab7c %>%
  mutate(across(starts_with('seq'), ~ifelse((.) == 3, 0, .))) %>%
  mutate(across(starts_with('seq'), ~ifelse((.) == 2, 0, .)))

otab7e <- otab7d %>%
  coalesce(., otab6)

# Bind the tables with three and two replicates for complete table
otab8 <- otab7e %>%
  bind_rows(., rep.2.mixed.1) %>%
  arrange(., sample.id)

# Filter to only include samples from Chaerephon / Mops pumilus : 
# First read in metadata file with all Chaerephon / Mops pumilus samples:
bat_metadata <- read.csv("bat_metadata.csv", header=TRUE)

# Then filter
tot.reads.samp <- otab8 %>%
  # Filter to keep relevant
  filter(sample.id %in% bat_metadata$Bat_ID) %>%
  group_by(sample.id) %>%
  # Calculate total reads per sample
  mutate(across(starts_with('seq'), ~sum(.))) %>%
  distinct() %>%
  rowwise() %>% 
  mutate(tot.reads = sum(c_across(starts_with("seq")))) %>%
  dplyr::select(sample.id, tot.reads, everything())

# For next steps, convert read counts to presence / absence
pres.abs.asv <- tot.reads.samp %>%
  dplyr::select(-tot.reads) %>%
  mutate_if(is.numeric, ~1 * (. > 0))

# Make sure no columns represent taxa found in 0 bats
bats.per.asv <- pres.abs.asv %>%
  ungroup() %>%
  dplyr::select(starts_with('seq')) %>%
  t(.) %>%
  as.data.frame(.) %>%
  mutate(num.bats = rowSums(.)) %>%
  dplyr::select(num.bats, everything())

# Find ASVs appearing in no bats
asvs.no.bats <- bats.per.asv %>%
  filter(num.bats==0)

# Edit the taxonomy file and filter out all the seq vars found in 0 bats
tax.03012026.filtered <- tax.02032025 %>% 
  rename("asv"="X.OTUID") %>%
  mutate(asv = str_replace_all(asv,"SeqVar", "seq")) %>%
  filter(asv %in% asvs.bats$asv)

# Filter out the eukaryotes
tax.filtered2 <- tax.03012026.filtered %>%
  filter(Kingdom == 'Bacteria')


incl.asv <- tax.filtered2 %>%
  dplyr::select(asv) 

# Filter out Eukaryotes and Archaea
tot.reads.samp.bact <- tot.reads.samp %>%
  select(sample.id, tot.reads, all_of(incl.asv$asv))
##################################################################################################


##################################################################################################
#### --- Counts / description of ASVs --- ####
##################################################################################################

# Relative abundance of ASVs
rel.abund <- tot.reads.samp.bact %>%
  transmute_at(vars(starts_with("seq")), funs(. / tot.reads))

# Presence / absence of ASVs
pres.abs.asv.bact <- tot.reads.samp.bact %>%
  dplyr::select(-tot.reads) %>%
  mutate_if(is.numeric, ~1 * (. > 0))

# Calculate ASV richness per bat
asv.rich <- pres.abs.asv.bact %>%
  rowwise() %>%
  mutate(asv.richness = sum(across(starts_with("seq")))) %>%
  dplyr::select(sample.id,asv.richness, everything()) %>%
  mutate(sample.id=as.integer(sample.id))

summary(asv.rich$asv.richness)

# Calculate number of bats per ASV
asvs.bats <- bats.per.asv %>%
  filter(num.bats>0) %>%
  rownames_to_column(., var="asv")

# Number of ASVs found in only 1 bat
only1bat <- asvs.bats %>%
  filter(num.bats ==1)
##################################################################################################


##################################################################################################
#### --- End of script, proceed to Script 02 --- ####
##################################################################################################
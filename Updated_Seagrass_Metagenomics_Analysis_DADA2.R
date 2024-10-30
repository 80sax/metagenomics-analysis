#### Updated Seagrass Microbiome Analysis Code ####
############## Author: Ha Lam Hai #################
############## Date: 12 Oct 2024 ##################

#### Loading packages to be used ####

library(dada2)
library(phyloseq)
library(kableExtra)
library(Biostrings)
library(ggplot2)
library(taxonomizr)
library(rBLAST)
library(GEOquery)
library(ggpubr)
library(knitr)
library(dplyr)
library(fantaxtic)
library(microbiome)
library(grid)
library(ggtext)
library(stringr)
library(ggh4x)
library(ggnested)
library(ggsci)
library(ggalluvial)
library(gridExtra)
library(RColorBrewer)
library(ggsignif)
library(tidyverse)
library(broom)
library(DESeq2)
library(vegan)
library(Polychrome)
library(texreg)
library(AICcmodavg)
library(ggstatsplot)

#### Loading training data for taxonomix assignment ####

taxonomy_training_data_wSpecies <- "silva_nr99_v138.1_wSpecies_train_set.fa.gz"

#### Loading data files from previous runs ####

readRDS("seqtab_nochim_combined.RDS") ## combined seqtab without chimeras
readRDS("taxa_combined_all.RDS") ## combined taxas

#### Loading updated metadata ####

metadata_path <- "metadata24.csv" ## path to metadata
meta <- read.csv(metadata_path, header = TRUE) ## read in the metadata
rownames(meta) <- meta$SampleID ## rename rows to sample ID
meta <- meta[order(rownames(meta)),] ## order metadata by sample ID

##----------------- Organizing and Troubleshooting DNA Reads -------------------

#### Pre-processing of data files ####

sample_path <- "N2416594_30-1058248066_Meta_2024-08-12/MS240808-2565" ## input path from working directory to folder containing fastq files
files <- list.files(sample_path) ## regex to search
newname <- sub('DNA-','', files) ## new name
file.rename(file.path(sample_path,files), file.path(sample_path, newname)) ##rename files

for (file in list.files(path = sample_path, pattern = ".fq.gz", full.names = TRUE))
{
  gunzip(file) ## unzip all files
}

for (filename in list.files(path = sample_path, pattern = ".fq", full.names = TRUE))
{
  base <- tools::file_path_sans_ext(filename)
  newname <- paste(sep="",base,".fastq") ## rename .fq files to .fastq (dada2 preferred input file format)
  file.rename(filename, newname)
}

#### Organizing forward and reverse data files ####

fwd_reads_path <- paste0(sample_path, '/fwd') ## name path for forward reads
rev_reads_path <- paste0(sample_path, '/rev') ## name path for reverse reads

dir.create(fwd_reads_path) # create paths
dir.create(rev_reads_path)

for (filename in list.files(path = sample_path, pattern = "1.fastq"))
{
  file.rename(from = paste0(sample_path, "/", filename), to = paste0(fwd_reads_path, "/", filename)) ## transfer files with forward reads to forward path
}

for (filename in list.files(path = sample_path, pattern = "2.fastq"))
{
  file.rename(from = paste0(sample_path, "/", filename), to = paste0(rev_reads_path, "/", filename)) ## transfer files with reverse reads to reverse path
}

#### Troubleshooting incomplete and erroneous data files ####

# In this step, reads from individual files will be checked to determine
# any non-matching read lengths and read quality line lengths. Reads can
# be found as the second line of every four lines in the fastq format.
# The quality line is the fourth line. For dada2, it is important to make
# sure the lengths of these two lines are the same for every read in the
# single file.

faulty_files <- list()

Fs_path <- sort(list.files(fwd_reads_path, full.names = TRUE))
Rs_path <- sort(list.files(rev_reads_path, full.names = TRUE))

paths <- c(Fs_path, Rs_path)

for (path in paths)
{
  for (file in path)
  {
    print(file)
    foo <- readLines(file)
    seqlens <- sapply(foo[seq(2,length(foo),4)], nchar)
    quallens <- sapply(foo[seq(4,length(foo),4)], nchar)
    faulty_files[[file]] <- as.character(which(seqlens != quallens))
  }
}

View(faulty_files)

# For each file listed in the faulty_files dataframe, the number of faulty
# lines in each file is indicated. 0 means there is no faulty lines. The
# faulty lines are also highlighted for each file. Individual inspections 
# and corrections of these files are necessary before proceeding with the
# rest of the codes. These files could be edited using any text editor
# (e.g. I personally used VS Code)

# Proceed ONLY when all data is error-free. This can be checked by running 
# the troubleshoot again until no errors are found.

##-------------------- Cleaning and quantifying reads --------------------------

#### QC of a random sample of reads ####

set.seed(1234)

len <- length(Fs_path) ## total number of files
n <- 12 ## number of sample profiles, n <= len

plotQualityProfile(Fs_path[sample(1:len, n, replace = FALSE)])
plotQualityProfile(Rs_path[sample(1:len, n, replace = FALSE)])

# Quality profiles of each sampled file will be shown as a plot. Y-axis
# is average quality score across reads, while x-axis is the nucleotide
# position. Set the trimming parameters in the following step accordingly.
# However, beware that too much shortening will result in failure to 
# merge the forward and reverse reads in a later step. If merging fails,
# return to this step to adjust the trimming parameters.

fwd_trim <- 0 ## how many nucleotides to trim from the front
rev_trim <- 0 ## and from the back

#### Filtering and trimming ####

Fs_path_filtered <- file.path(fwd_reads_path, "filtered_Fs") ## create paths for filtered reads
Rs_path_filtered <- file.path(rev_reads_path, "filtered_Rs")

sample_names <- str_replace(string = basename(Fs_path), 
                            pattern = "_1\\.fastq",
                            replacement = "") ## obtain all sample names (filenames)

sample_names <- sample_names[sample_names != 'filtered_Fs']

out <- filterAndTrim(fwd = Fs_path, filt = Fs_path_filtered, rev = Rs_path,
                     filt.rev = Rs_path_filtered, matchIDs = TRUE, maxEE = c(2,2), 
                     truncQ = 2, maxN = 0, rm.phix = TRUE, verbose = TRUE, 
                     compress = TRUE, trimLeft = fwd_trim, trimRight = rev_trim) ## filter and trim reads according to pre-set parameters

Fs_filt <- list.files(Fs_path_filtered, full.names = TRUE, pattern = "fastq") ## list filtered file names
Rs_filt <- list.files(Rs_path_filtered, full.names = TRUE, pattern = "fastq")

names(Fs_filt) <- sample_names ## rename with sample names
names(Rs_filt) <- sample_names

#### Predicting error rates ####

errF <- learnErrors(Fs_filt, multithread = TRUE)
errR <- learnErrors(Rs_filt, multithread = TRUE)

#### Creating dada data sets ####

dadaFs <- dada(Fs_filt, err = errF, multithread = TRUE)
dadaRs <- dada(Rs_filt, err = errR, multithread = TRUE)

saveRDS(dadaFs, file = "./dadaFs.RDS") ## temporarily save these files as a pitstop
saveRDS(dadaRs, file = "./dadaRs.RDS")

#### Merging paired reads and making sequence table ####

mergers <- mergePairs(dadaFs, Fs_path_filtered, dadaRs, Rs_path_filtered, verbose = TRUE)
seqtab <- makeSequenceTable(mergers)

#### Removing chimeras ####

seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)

df <- as.dataFrame(rownames(seqtab_nochim)) # to extract sample names out to update metadata
write.csv(df, "rownames.csv", row.names = F)

seqtab_nochim_combined <- mergeSequenceTables(seqtab_nochim_combined, seqtab_nochim6) ## merge the new seqtab_nochim with existing combined seqtab_nochim
saveRDS(seqtab_nochim_combined, file = "./seqtab_nochim_combined.RDS") ## save file for the next batch of data

num_chim_removed <- 1 - (sum(seqtab_nochim/sum(seqtab)))
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab_nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample_names
kableExtra::kable(track)

##------------------------ Taxonomy assignment ---------------------------------

#### Assigning taxa to reads (to species level) ####

# Process described here is broken into every 5000 reads because this 
# step is computationally heavy. If the code is run on a computer with 
# a large enough RAM, one could process every at once, using:

# taxa_combined <- assignTaxonomy(seqtab_nochim, taxonomy_training_data_wSpecies, verbose = TRUE)

taxa_combined <- assignTaxonomy(seqtab_nochim[1,1:5000], taxonomy_training_data_wSpecies, verbose = TRUE) ## initialize the assignment process with the first 5000 reads

for (i in seq(5001, dim(seqtab_nochim)[2], by = 5001)) {
  if (i + 5000 <= dim(seqtab_nochim)[2]) {
    tax <- assignTaxonomy(seqtab_nochim[1,i:(i + 5000)], taxonomy_training_data_wSpecies, verbose = TRUE)
  }
  else {
    tax <- assignTaxonomy(seqtab_nochim[1,i:dim(seqtab_nochim6)[2]], taxonomy_training_data_wSpecies, verbose = TRUE)
  }
  taxa_combined <- rbind(taxa_combined, tax)
}

taxa_combined_all <- rbind(taxa_combined_all, taxa_combined) ## combine the existing taxa_table with the new taxa_table
taxa_combined_all <- taxa_combined_all[!duplicated(rownames(taxa_combined_all)),] ## remove rows containing duplicated reads
saveRDS(taxa_combined_all, file = "./taxa_combined_all.RDS") ## save file for the next batch of data

#### Creating subsets of taxas of interest ####

list_pathogens <- c("Shewanella", "Arcobacter", "Vibrio", "Tenacibaculum", 
                    "Enterococcus", "Streptococcus", "Pseudomonas", "Photobacterium", 
                    "Lactococcus", "Flavobacterium", "Rickettsia", "Pseudoalteromonas", 
                    "Clostridium", "Thalassomonas", "Bacillus", "Chryseobacterium", 
                    "Phormidium", "Corynebacterium") ## list of pathogenic genera of interest

list_biofilm <- c("Gemmatimonadetes", "Calditrichaceae",
                  "Candidatus Entotheonella", "Microtrichales", "Nitrosococcaceae", "Lentisphaera") ## list of biofilm-making bacteria of interest

pathogen_taxa_all <- subset(taxa_combined_all, taxa_combined_all[,"Genus"] %in% list_pathogens)
biofilm_taxa_all <- subset(taxa_combined_all, taxa_combined_all[, "Phylum"] %in% list_biofilm | taxa_combined_all[, "Genus"] %in% list_biofilm | taxa_combined_all[, "Class"] %in% list_biofilm | taxa_combined_all[, "Order"] %in% list_biofilm | taxa_combined_all[, "Family"] %in% list_biofilm)
bdellovibrio_all <- subset(taxa_combined_all, taxa_combined_all[,"Genus"] == "Bdellovibrio")

##---------------------- Constructing phyloseq objects -------------------------

# For all taxas:
ps_all <- phyloseq(otu_table(seqtab_nochim_combined, taxa_are_rows=FALSE), sample_data(meta), tax_table(taxa_combined_all))
saveRDS(ps_all, file = "./ps_all.RDS")

# For taxa subsets:
ps_pathogens <- phyloseq(otu_table(seqtab_nochim_combined, taxa_are_rows=FALSE), sample_data(meta), tax_table(pathogen_taxa_all))
ps_biofilm <- phyloseq(otu_table(seqtab_nochim_combined, taxa_are_rows=FALSE), sample_data(meta), tax_table(biofilm_taxa_all))
ps_bdellovibrio <- phyloseq(otu_table(seqtab_nochim_combined, taxa_are_rows=FALSE), sample_data(meta), tax_table(bdellovibrio_all))
saveRDS(ps_pathogens, file = "./ps_pathogens.RDS")
saveRDS(ps_biofilm, file = "./ps_biofilm.RDS")
saveRDS(ps_bdellovibrio, file = "./ps_bdellovibrio.RDS")

##------------------ Plotting and statistical inference ------------------------

#### Setting up Basic Dataframes ####

# Alpha diversity (Shannon/Simpson)

alpha_diversity.all <- estimate_richness(ps_all, measures = c("shannon", "simpson"))

meta_withH <- merge(meta, alpha_diversity.all, by = 0) ## merging with metadata
rownames(meta_withH) <- meta_withH$SampleID

# Separate different studies

ps_all.8site <- subset_samples(ps_all, Study == "enterodat21")
ps_all.removal <- subset_samples(ps_all, Study == "removal_exp_v3")
ps_all.meso <- subset_samples(ps_all, Study == "mesocosm")
ps_all.rel <- microbiome::transform(ps_all, "compositional") ## calculating relative abundances
ps_all.rel.8site <- subset_samples(ps_all.rel, Study == "enterodat21")
ps_all.rel.removal <- subset_samples(ps_all.rel, Study == "removal_exp_v3")
ps_all.rel.meso <- subset_samples(ps_all.rel, Study == "mesocosm")

ps_pathogens.8site <- subset_samples(ps_pathogens, Study == "enterodat21")
ps_pathogens.removal <- subset_samples(ps_pathogens, Study == "removal_exp_v3")
ps_pathogens.meso <- subset_samples(ps_pathogens, Study == "mesocosm")

ps_biofilm.8site <- subset_samples(ps_biofilm, Study == "enterodat21")
ps_biofilm.removal <- subset_samples(ps_biofilm, Study == "removal_exp_v3")
ps_biofilm.meso <- subset_samples(ps_biofilm, Study == "mesocosm")

ps_bdellovibrio.8site <- subset_samples(ps_bdellovibrio, Study == "enterodat21")
ps_bdellovibrio.removal <- subset_samples(ps_bdellovibrio, Study == "removal_exp_v3")
ps_bdellovibrio.meso <- subset_samples(ps_bdellovibrio, Study == "mesocosm")

#### PERMANOVA ####

calculate_permanova <- function(pseq.rel, core_detec, core_prev, filename) {
  otu <- abundances(pseq.rel)
  meta <- meta(pseq.rel)
  permanova <- adonis2(t(otu) ~ Treatment,
                       data = meta, permutations=99, method = "bray")
  print(permanova)
  capture.output(permanova, file = paste(filename, "all permanova results.txt", sep = " "))
  pseq.core <- core(pseq.rel, detection = core_detec, prevalence = core_prev) ## obtain core microbiome
  otu <- abundances(pseq.core)
  meta <- meta(pseq.core)
  permanova <- adonis2(t(otu) ~ Treatment,
                       data = meta, permutations = 99, method = "bray") ## check for NAs in pseq.core if this returns error.
  print(permanova)
  capture.output(permanova, file = paste(filename, "core permanova results.txt", sep = " "))
}

ps_all.rel.8site.subset <- subset_samples(ps_all.rel.8site, SampleID != "CH11-LEK7958_L1" & SampleID != "CH21-LEK7959_L1" & SampleID !="CH31-LEK7960_L1") ## nothing detected in core

calculate_permanova(ps_all.rel.8site.subset, 0, .3, "8 site") ## View txt files created using Notepad
calculate_permanova(ps_all.rel.removal, 0, .3, "removal")

#### Rarefaction Curve ####

rarefaction_curve <- function(pseq, n_colours, steps) {
  comm_data <- otu_table(pseq) ## obtain community data from phyloseq object
  comm_data <- as(comm_data, "matrix")
  
  plot_colours <- createPalette(n_colours,  c("#ff0000","#ffff00","#00ff00","#0000ff")) ## create a large enough colour palett
  
  S <- specnumber(comm_data) ## observed number of species
  raremax <- min(rowSums(comm_data))
  Srare <- rarefy(comm_data, raremax)
  plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
  abline(0, 1)
  
  rcurve <- rarecurve(comm_data, step = steps, sample = raremax, col = plot_colours, cex = 0.6, ylab = "OTU", label = FALSE, tidy = TRUE)
  colnames(rcurve)[1] <- "Sample_ID"
  colnames(rcurve)[2] <- "Sample_size"
  colnames(rcurve)[3] <- "OTU"
  
  rcurve_ggplot_wLegends <- ggplot(data = rcurve, aes(x = Sample_size, y = OTU, group = Sample_ID)) + geom_line(aes(color = Sample_ID), show.legend = TRUE) + xlab("Sample Size") + ylab("OTU") + theme_classic()
  rcurve_ggplot_woLegends <- ggplot(data = rcurve, aes(x = Sample_size, y = OTU, group = Sample_ID)) + geom_line(aes(color = Sample_ID), show.legend = FALSE) + xlab("Sample Size") + ylab("OTU") + theme_classic()
  grid.draw(rcurve_ggplot_woLegends)
  legend <- cowplot::get_legend(rcurve_ggplot_wLegends)
  grid.newpage()
  grid.draw(legend) ## print out a separate page for legends
}

rarefaction_curve(ps_all, 151, 5000) ## try a larger n_colours if more colours are needed

#### Core Microbiome ####

# Following function converts treatment from 0's and 1's to "+" and "-"
convert_treatment <- function(pseq) {
  for (row in rownames(sample_data(pseq))) {
    if (sample_data(pseq)[row, "Treatment"] == 0) {
      sample_data(pseq)[row, "Treatment"] <- "+"
    }
    else {
      sample_data(pseq)[row, "Treatment"] <- "-"
    }
  }
  return(pseq)
}

# For the 8 sites:
ps_all.rel.wSeagrass.8site <- subset_samples(ps_all.rel, Treatment == 0 & Study == "enterodat21") ## separating samples with and without seagrass
ps_all.rel.woSeagrass.8site <- subset_samples(ps_all.rel, Treatment == 1 & Study == "enterodat21")
pseq.core.wSeagrass <- core(ps_all.rel.wSeagrass.8site, detection = 0, prevalence = .3) ## getting core microbiome
pseq.core.woSeagrass <- core(ps_all.rel.woSeagrass.8site, detection = 0, prevalence = .3)
pseq.core.wSeagrass.merged <- merge_samples(pseq.core.wSeagrass, "Location") ## merging by location
pseq.core.woSeagrass.merged <- merge_samples(pseq.core.woSeagrass, "Location")

pseq.core.wSeagrass.merged <- convert_treatment(pseq.core.wSeagrass.merged)
pseq.core.woSeagrass.merged <- convert_treatment(pseq.core.woSeagrass.merged)

# Following function processes the microbiome data obtained
process_microbiome <- function(pseq) {
  pseq <- subset_taxa(pseq, !is.na(Genus))
  pseq <- tax_glom(pseq, "Genus")
  pseq <- transform_sample_counts(pseq, function(x) x/sum(x))
  pseq.df <- psmelt(pseq)
  pseq.df$Genus <- as.character(pseq.df$Genus)
  pseq.df$Genus[pseq.df$Abundance < 0.05] <- "Others"
  return(pseq.df)
}

pseq.core.wSeagrass.merged.rel.df <- process_microbiome(pseq.core.wSeagrass.merged)
pseq.core.woSeagrass.merged.rel.df <- process_microbiome(pseq.core.woSeagrass.merged)

core.means <- pseq.core.wSeagrass.merged.rel.df %>% group_by(Genus) %>% 
  summarise(Mean.Abundance=mean(Abundance),
            .groups = 'drop')
core.means <- core.means[order(core.means$Mean.Abundance, decreasing = TRUE),] ## ordering the plot by core microbiome relative abundances
order.core <- core.means$Genus
order.core <- order.core[!order.core == "Others"]
order.core.pos <- append("Others", order.core)

plotting_core_microbiome <- function(df, order, x, y, title) {
  ggplot(df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = order))) +
    geom_bar(position = "stack", stat = "identity", width = 1) +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 90, hjust=1,vjust=0.3), axis.text.y = element_text(vjust=0.4)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    ylab(substitute(paste(bold(y)))) + scale_fill_brewer("Genus", palette = "Paired") +
    xlab(substitute(paste(bold(x)))) +
    theme(axis.title = element_text(size = 10)) + ggtitle(title) + theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5))
}

plt.core.pos <- plotting_core_microbiome(pseq.core.wSeagrass.merged.rel.df, order.core.pos, "Relative Abundance", "Location", "Seagrass+")

core.means <- pseq.core.woSeagrass.merged.rel.df %>% group_by(Genus) %>% 
  summarise(Mean.Abundance=mean(Abundance),
            .groups = 'drop')
core.means <- core.means[order(core.means$Mean.Abundance, decreasing = TRUE),]
order.core <- core.means$Genus
order.core <- order.core[!order.core == "Others"]
order.core.neg <- append("Others", order.core) ## to compare this order with the seagrass+ order. 
# Maintaining consistencies with the seagrass+ plot
order.core.neg[2] <- "HIMB11"
order.core.neg[3] <- "Cyanobium PCC-6307"
order.core.neg[4] <- "Synechococcus CC9902"

plt.core.neg <- plotting_core_microbiome(pseq.core.woSeagrass.merged.rel.df, order.core.neg, "Relative Abundance", "Location", "Seagrass-")

require(grid)
plt_combined <- ggarrange(plt.core.neg, plt.core.pos  + rremove("ylab"), ## remove axis labels from plots
                          labels = NULL,
                          ncol = 2, nrow = 1,
                          common.legend = TRUE, legend = "right",
                          align = "v", 
                          font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"), heights = c(1,1.5))

plt_combined ## combined plot

# For the removal plots:

ps_all.rel.removal.noNA <- subset_samples(ps_all.rel.removal, !is.na(sample_data(ps_all.rel.removal)$Plot)) ## removing plots without plot numbers
sample_data(ps_all.rel.removal.noNA)[,7] <- paste0(sample_data(ps_all.rel.removal.noNA)$Plot, sample_data(ps_all.rel.removal.noNA)$Day, sample_data(ps_all.rel.removal.noNA)$Month, sample_data(ps_all.rel.removal.noNA)$Year)
colnames(sample_data(ps_all.rel.removal.noNA))[7] <- "PlotDate"

ps_all.rel.removal.control <- subset_samples(ps_all.rel.removal.noNA, Treatment == 0)
ps_all.rel.removal.treatment <- subset_samples(ps_all.rel.removal.noNA, Treatment == 1)

pseq.core.removal.control <- core(ps_all.rel.removal.control, detection = 0, prevalence = .3)
pseq.core.removal.treatment <- core(ps_all.rel.removal.treatment, detection = 0, prevalence = .3)

pseq.core.removal.control.merged <- merge_samples(pseq.core.removal.control, "PlotDate")
pseq.core.removal.treatment.merged <- merge_samples(pseq.core.removal.treatment, "PlotDate")

pseq.core.removal.control.merged <- convert_treatment(pseq.core.removal.control.merged)
pseq.core.removal.treatment.merged <- convert_treatment(pseq.core.removal.treatment.merged)

pseq.core.removal.control.merged.rel.df <- process_microbiome(pseq.core.removal.control.merged)
pseq.core.removal.treatment.merged.rel.df <- process_microbiome(pseq.core.removal.treatment.merged)

core.means <- pseq.core.removal.control.merged.rel.df %>% group_by(Genus) %>% 
  summarise(Mean.Abundance=mean(Abundance),
            .groups = 'drop')
core.means <- core.means[order(core.means$Mean.Abundance, decreasing = TRUE),]
order.core <- core.means$Genus
order.core <- order.core[!order.core == "Others"]
order.core1 <- append("Others", order.core)

plt.core.control <- plotting_core_microbiome(pseq.core.removal.control.merged.rel.df, order.core1, "Relative Abundance", "Plot Number", "Control")
plt.core.treatment <- plotting_core_microbiome(pseq.core.removal.treatment.merged.rel.df, order.core1, "Relative Abundance", "Plot Number", "Treatment")

require(grid)
plt_combined <- ggarrange(plt.core.control, plt.core.treatment + rremove("ylab"), # remove axis labels from plots
                          labels = NULL,
                          ncol = 2, nrow = 1,
                          common.legend = TRUE, legend = "right",
                          align = "v", 
                          font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"), heights = c(1,1.5))

plt_combined

#### Relative Abundance Plots ####

# For 8 sites:

merge_by_loctreat <- function(pseq.8site) { ## merging data based on Location and Treatment
  sample_data(pseq.8site)[,"LocTreat"] <- paste0(sample_data(pseq.8site)$Location, sample_data(pseq.8site)$Treatment)
  pseq.8site.merged <- merge_samples(pseq.8site, "LocTreat")
  for (row1 in rownames(sample_data(pseq.8site.merged))) {
    for (row2 in rownames(sample_data(pseq.8site))) {
      if (sample_data(pseq.8site)[row2, "LocTreat"] == row1) {
        sample_data(pseq.8site.merged)[row1, "Location"] <- sample_data(pseq.8site)[row2, "Location"]
      }
      next
    }
  }
  return(pseq.8site.merged)
}  

ps_all.8site.merged <- merge_by_loctreat(ps_all.8site)
ps_all.8site.merged <- convert_treatment(ps_all.8site.merged)
ps_all.8site.merged.rel.df <- process_microbiome(ps_all.8site.merged)

means.8site <- ps_all.8site.merged.rel.df %>% group_by(Genus) %>% 
  summarise(Mean.Abundance=mean(Abundance),
            .groups = 'drop')
means.8site <- means.8site[order(means.8site$Mean.Abundance, decreasing = TRUE),]
order.8site <- means.8site$Genus
order.8site <- order.8site[!order.8site == "Others"]
order.8site <- append("Others", order.8site)

plt.8site <- ggplot(ps_all.8site.merged.rel.df, aes(x = Location, y = Abundance, fill = factor(Genus, levels = order.8site))) +
  geom_bar(position = "stack", stat = "identity", width = 1) +
  facet_grid2(~Treatment, scales = "free") + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0))) +
  ylab(substitute(paste(bold("Relative abundance")))) + scale_fill_brewer("Genus", palette = "Paired") +
  xlab(substitute(paste(bold("Location")))) + theme(strip.text.x = element_text(colour = "black", size = 10))

plt.8site

# For removal plots:
merge_by_plotdate <- function(pseq.removal) { ## merging data based on Plot Number and Date
  sample_data(pseq.removal)[,"PlotDate"] <- paste0(sample_data(pseq.removal)$Plot, sample_data(pseq.removal)$Day, sample_data(pseq.removal)$Month, sample_data(pseq.removal)$Year)
  pseq.removal.merged <- merge_samples(pseq.removal, "PlotDate")
}

ps_all.rel.removal.noNA <- subset_samples(ps_all.rel.removal, !is.na(sample_data(ps_all.rel.removal)$Plot) & !is.na(sample_data(ps_all.rel.removal)$Month))
ps_all.rel.removal.noNA.merged <- merge_by_plotdate(ps_all.rel.removal.noNA)

for (row in rownames(sample_data(ps_all.rel.removal.noNA.merged.rel))) {
  if (sample_data(ps_all.rel.removal.noNA.merged.rel)[row, "Treatment"] == 0) {
    sample_data(ps_all.rel.removal.noNA.merged.rel)[row, "Treatment"] <- "Control"
  }
  else {
    sample_data(ps_all.rel.removal.noNA.merged.rel)[row, "Treatment"] <- "Treatment"
  }
}

ps_all.rel.removal.noNA.merged.df <- process_microbiome(ps_all.rel.removal.noNA.merged.rel)

means.removal <- ps_all.rel.removal.noNA.merged.df %>% group_by(Genus) %>% 
  summarise(Mean.Abundance=mean(Abundance),
            .groups = 'drop')
means.removal <- means.removal[order(means.removal$Mean.Abundance, decreasing = TRUE),]
order.removal <- means.removal$Genus
order.removal <- order.removal[!order.removal == "Others"]
order.removal <- append("Others", order.removal)

pair_expanded <- colorRampPalette(brewer.pal(9, "Paired"))(17)
plt.removal <- ggplot(ps_all.rel.removal.noNA.merged.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = order.removal))) +
  geom_bar(position = "stack", stat = "identity", width = 1) +
  facet_nested(~Treatment + Plot, nest_line = FALSE, scales = "free", strip = ggh4x::strip_split(c("top", "bottom"))) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(panel.spacing.x = unit(0,"line")) +
  ylab(substitute(paste(bold("Relative abundance")))) + scale_fill_manual("Genus", values = pair_expanded) +
  xlab(substitute(paste(bold("Plot number")))) + theme(strip.text.x = element_text(colour = "black", size = 10)) + theme(panel.spacing.x = unit(c(rep(0, 3), 1, rep(0, 3)), "lines"))

plt.removal

# For pathogenic species (8 sites - relative abundance):

ps_pathogens.8site.merged <- merge_by_loctreat(ps_pathogens.8site)

ps_pathogens.8site.merged <- convert_treatment(ps_pathogens.8site.merged)

ps_pathogens.8site.merged.noNA <- subset_taxa(ps_pathogens.8site.merged, !is.na(Genus))
ps_pathogens.8site.merged.noNA.glom <- tax_glom(ps_pathogens.8site.merged.noNA, "Genus")
ps_pathogens.8site.merged.noNA.glom.rel <- transform_sample_counts(ps_pathogens.8site.merged.noNA.glom, function(x) x/sum(x))
ps_pathogens.8site.merged.noNA.glom.rel.df <- psmelt(ps_pathogens.8site.merged.noNA.glom.rel)
ps_pathogens.8site.merged.noNA.glom.rel.df$Genus <- as.character(ps_pathogens.8site.merged.noNA.glom.rel.df$Genus)

ps_pathogens.8site.means <- ps_pathogens.8site.merged.noNA.glom.rel.df %>% group_by(Genus) %>% 
  summarise(Mean.Abundance=mean(Abundance),
            .groups = 'drop')
ps_pathogens.8site.means <- ps_pathogens.8site.means[order(ps_pathogens.8site.means$Mean.Abundance, decreasing = TRUE),]
order.ps_pathogens.8site <- ps_pathogens.8site.means$Genus

pair_expanded2 <- colorRampPalette(brewer.pal(12, "Paired"))(14) ## creating a large enough colour palette

ps_pathogens.8site.merged.noNA.glom.rel.df.pos <- ps_pathogens.8site.merged.noNA.glom.rel.df[ps_pathogens.8site.merged.noNA.glom.rel.df$Treatment == "+",]
ps_pathogens.8site.merged.noNA.glom.rel.df.pos$Sample <- substr(ps_pathogens.8site.merged.noNA.glom.rel.df.pos$Sample, 1, nchar(ps_pathogens.8site.merged.noNA.glom.rel.df.pos$Sample)-1)
ps_pathogens.8site.merged.noNA.glom.rel.df.neg <- ps_pathogens.8site.merged.noNA.glom.rel.df[ps_pathogens.8site.merged.noNA.glom.rel.df$Treatment == "-",]
ps_pathogens.8site.merged.noNA.glom.rel.df.neg$Sample <- substr(ps_pathogens.8site.merged.noNA.glom.rel.df.neg$Sample, 1, nchar(ps_pathogens.8site.merged.noNA.glom.rel.df.neg$Sample)-1)

plotting_pathogens_abundance <- function(pseq.df, order, title) {
  ggplot(pseq.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = order))) +
    geom_bar(position = "stack", stat = "identity", width = 0.8) +
    facet_grid(~Location, scales = "free") + 
    theme(axis.ticks.x = element_blank()) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    ylab(substitute(paste(bold("Relative abundance")))) + scale_fill_manual("Genus", values = pair_expanded2) +
    xlab(substitute(paste(bold("Location")))) + theme(strip.text.x = element_blank()) +
    theme(axis.title = element_text(size = 11), axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90)) + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11), legend.title = element_text(face = "bold"))
}

plt.ps_pathogens.8site.pos <- plotting_pathogens_abundance(ps_pathogens.8site.merged.noNA.glom.rel.df.pos, order.ps_pathogens.8site, "Seagrass+")
plt.ps_pathogens.8site.neg <- plotting_pathogens_abundance(ps_pathogens.8site.merged.noNA.glom.rel.df.neg, order.ps_pathogens.8site, "Seagrass-")

plt_combined <- ggarrange(plt.ps_pathogens.8site.pos, plt.ps_pathogens.8site.neg + rremove("ylab"), # remove axis labels from plots
                          labels = NULL, 
                          ncol = 2, nrow = 1, widths = c(1.5, 1),
                          common.legend = TRUE, legend = "right",
                          align = "h", 
                          font.label = list(size = 14, color = "black", face = "bold", family = NULL, position = "top"), heights = c(1,1.5))

plt_combined

# For pathogenic species (8 sites - changes in relative abundance):

ps_pathogens.8site.merged.bySeagrass <- merge_samples(ps_pathogens.8site, "Treatment")
ps_pathogens.8site.merged.bySeagrass <- convert_treatment(ps_pathogens.8site.merged.bySeagrass)

ps_pathogens.8site.merged.bySeagrass <- subset_taxa(ps_pathogens.8site.merged.bySeagrass, !is.na(Genus))
ps_pathogens.8site.merged.bySeagrass <- tax_glom(ps_pathogens.8site.merged.bySeagrass, "Genus")
ps_pathogens.8site.merged.bySeagrass <- transform_sample_counts(ps_pathogens.8site.merged.bySeagrass, function(x) x/sum(x))
ps_pathogens.8site.merged.bySeagrass <- psmelt(ps_pathogens.8site.merged.bySeagrass)
ps_pathogens.8site.merged.bySeagrass$Genus <- as.character(ps_pathogens.8site.merged.bySeagrass$Genus)

ps_pathogens.8site.changes <- data.frame(Genus = unique(ps_pathogens.8site.merged.bySeagrass$Genus))
ps_pathogens.8site.changes$Change <- c(rep(0,14))
rownames(ps_pathogens.8site.changes) <- ps_pathogens.8site.changes$Genus
for (g in rownames(ps_pathogens.8site.changes)) {
  change <- ps_pathogens.8site.merged.bySeagrass[ps_pathogens.8site.merged.bySeagrass$Treatment == "-" & ps_pathogens.8site.merged.bySeagrass$Genus == g,]$Abundance - ps_pathogens.8site.merged.bySeagrass[ps_pathogens.8site.merged.bySeagrass$Treatment == "+" & ps_pathogens.8site.merged.bySeagrass$Genus == g,]$Abundance
  ps_pathogens.8site.changes[g,]$Change <- change
}

ps_pathogens.8site.changes$Seagrass <- c(rep("Present", 14))

for (g in rownames(ps_pathogens.8site.changes)) {
  if (ps_pathogens.8site.changes[g,]$Change < 0) {
    ps_pathogens.8site.changes[g,]$Seagrass <- "Absent"
  }
}

relative_abundance_changes_pathogens_plot <- ggplot(ps_pathogens.8site.changes) + geom_bar(aes(x = reorder(Genus, -Change), y = Change, fill = Seagrass), stat = 'identity') + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3, margin = margin(t = 5))) + xlab("Pathogenic Genus") + ylab("Change in Relative Abundance (%)") + theme(axis.title = element_text(face = "bold"))

relative_abundance_changes_pathogens_plot

#### pCOA ####

# For whole microbiome:

ps_all.8site.prop <- transform_sample_counts(ps_all.8site, function(otu) otu/sum(otu))
sample_data(ps_all.8site.prop)$Location <- as.factor(sample_data(ps_all.8site.prop)$Location)

ps_all.8site.prop <- convert_treatment(ps_all.8site.prop)

pCoa.bray.8site <- ordinate(ps_all.8site.prop, method="PCoA", distance="bray", trymax = 2000)

plt_pCoa.8site <- plot_ordination(ps_all.8site.prop, pCoa.bray.8site, color = "Location", shape = "Seagrass") + geom_point(size = 2.25) + theme_classic() +
  stat_ellipse(geom = "polygon", lwd = 0.5, aes(fill = Location), alpha = 0.05)

plt_pCoa.8site

# For core microbiome:

ps_all.rel.wSeagrass.8site <- subset_samples(ps_all.rel, Treatment == 0 & Study == "enterodat21") ## separating samples with and without seagrass
ps_all.rel.woSeagrass.8site <- subset_samples(ps_all.rel, Treatment == 1 & Study == "enterodat21")
pseq.core.wSeagrass <- core(ps_all.rel.wSeagrass.8site, detection = 0, prevalence = .3) ## getting core microbiome
pseq.core.woSeagrass <- core(ps_all.rel.woSeagrass.8site, detection = 0, prevalence = .3)
pseq.core.combined <- merge_phyloseq(pseq.core.wSeagrass, pseq.core.woSeagrass) ## merging the two core microbiome phyloseqs
pseq.core.combined.prop <- transform_sample_counts(pseq.core.combined, function(otu) otu/sum(otu))
sample_data(pseq.core.combined.prop)$Location <- as.factor(sample_data(pseq.core.combined.prop)$Location)

pseq.core.combined.prop <- convert_treatment(pseq.core.combined.prop)

pseq.core.combined.prop.cleaned <- subset_samples(pseq.core.combined.prop, SampleID != "CH11-LEK7958_L1" & SampleID != "CH31-LEK7960_L1" & SampleID != "CH21-LEK7959_L1") ## removing samples with NA's

pCoa.bray.core <- ordinate(pseq.core.combined.prop.cleaned, method="PCoA", distance = "bray", trymax = 2000)

plt_pCoa.8site.core <- plot_ordination(pseq.core.combined.prop.cleaned, pCoa.bray.core, color = "Location", shape = "Seagrass") + geom_point(size = 2.25) + theme_classic() +
  stat_ellipse(geom = "polygon", lwd = 0.5, aes(fill = Location), alpha = 0.05)

#### Procrustes Analysis ####

ps_all.removal.prop.pos <- subset_samples(ps_all.rel.removal, Treatment == 0)
ps_all.removal.prop.neg <- subset_samples(ps_all.rel.removal, Treatment == 1)

samples1 <- as.vector(sample(rownames(sample_data(ps_all.removal.prop.pos)), 11)) ## need to re-iterate this step if the final plot is not as desired.
samples2 <- c()

# An example of a sample that works:
# samples1 <- c("EB73-LFE8401_L0", "EB223-LGA13932_L1", "EB164-LFK4824_L1", "EB229-LGA13933_L1", "EB92-LFJ3169-RL1_L1", "EB236-LGA13931_L1", "EB51-LFE8384_L0", "EB74-LFE8403-RL1_L1", "EB62-LFE8410_L0", "EB67-LFE8413_L1", "EB05-LFE8371_L0")

ps_all.removal.prop.pos.updated <- subset_samples(ps_all.removal.prop.pos, !SampleID %in% samples1)
ps_all.removal.prop.neg.updated <- subset_samples(ps_all.removal.prop.neg, !SampleID %in% samples2)

pCoa.bray.removal.pos <- ordinate(ps_all.removal.prop.pos.updated, method="RDA", distance="bray", trymax = 2000)
pCoa.bray.removal.neg <- ordinate(ps_all.removal.prop.neg.updated, method="RDA", distance="bray", trymax = 2000)

pro <- procrustes(X = pCoa.bray.removal.pos, Y = pCoa.bray.removal.neg, symmetric = FALSE)

pro

ctest <- data.frame(rda1=pro$Yrot[,1],
                    rda2=pro$Yrot[,2],xrda1=pro$X[,1],
                    xrda2=pro$X[,2])
library(scales)
hex <- hue_pal()(2)

ggplot(ctest) +
  geom_point(aes(x=rda1, y=rda2, color = "#F8766D"), size = 3) +
  geom_point(aes(x=xrda1, y=xrda2, color = "#00BFC4"), size = 3) + 
  labs(color = "Treatment\n") +
  scale_color_manual(labels = c("Negative", "Positive"), values = c("#F8766D", "#00BFC4")) +
  geom_segment(aes(x=rda1,y=rda2,xend=xrda1,yend=xrda2), arrow = arrow(length = unit(0.2, "cm"))) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black")) + xlab("PCoa 1") + ylab("PCoa 2")

protest <- protest(X = pCoa.bray.removal.pos, Y = pCoa.bray.removal.neg, permutations = 100000)

protest

densityplot(permustats(protest))

#### Flow Charts ####

# For biofilm species: 

ps_biofilm.removal.noNA <- subset_samples(ps_biofilm.removal, !is.na(sample_data(ps_biofilm.removal)$Plot) & !is.na(sample_data(ps_biofilm.removal)$Month))
ps_biofilm.removal.noNA.merged <- merge_by_plotdate(ps_biofilm.removal.noNA)
ps_biofilm.removal.noNA.merged <- transform_sample_counts(ps_biofilm.removal.noNA.merged, function(otu) otu/sum(otu))

ps_biofilm.removal.noNA.merged <- convert_treatment(ps_biofilm.removal.noNA.merged)

biofilm_bact <- data.frame(tax_level = c("Class", "Family", "Genus", "Order", "Family", "Genus"), name = list_biofilm) ## creating a dataframe containing the biofilm bacteria names and the associated taxa levels
rownames(biofilm_bact) <- biofilm_bact$name

bact.abun.list <- list() ## creating a list of the biofilm taxas and their relative abundance at each plot number and date
for (bact in rownames(biofilm_bact)) {
  print(bact)
  ps_biofilm.removal.noNA.merged.glom <- tax_glom(ps_biofilm.removal.noNA.merged, taxrank = biofilm_bact[bact, "tax_level"])
  ps_biofilm.removal.noNA.merged.glom <- psmelt(ps_biofilm.removal.noNA.merged.glom)
  bact.abun <- ps_biofilm.removal.noNA.merged.glom[,colnames(ps_biofilm.removal.noNA.merged.glom) %in% c("Treatment", "Plot", "Month", "Abundance", "Class", "Order", "Family", "Genus", "PlotMonth")]
  bact.abun$Plot <- as.factor(bact.abun$Plot)
  bact.abun$Month <- as.factor(bact.abun$Month)
  bact.abun <- bact.abun[bact.abun[,biofilm_bact[bact, "tax_level"]] == bact,]
  bact.abun.list <- append(bact.abun.list, list(list(bact, bact.abun)))
}

plt.list <- list() ## creating a list of plots
for (bact in bact.abun.list) {
  plt <- ggplot(data = bact[[2]], 
                aes(x = Month, y = Abundance, alluvium = Plot)) + 
    geom_alluvium(aes(fill = Treatment, colour = Treatment), alpha = .75, decreasing = FALSE) +
    theme_minimal() +
    ylab("Relative<br>abundance") +
    theme(axis.text=element_text(size=8),
          axis.title=element_text(size=10,face="bold")) + theme(legend.key.size = unit(0.4, "cm")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(axis.line = element_line(colour = "black")) +
    theme(legend.title=element_text(size=9, face = "bold"), legend.text=element_text(face = "italic")) + 
    theme(axis.title.y = element_markdown(), legend.title = element_markdown(), axis.text.x = element_markdown(), legend.title.align = 0.5) +
    ggtitle(bact[[1]]) + theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold")) + theme(axis.text.y = element_blank())
  plt.list <- append(plt.list, list(plt))
}

require(grid)

plt.all.biofilm <- ggarrange(plt.list[[1]], plt.list[[2]] + rremove("ylab"), plt.list[[4]] + rremove("ylab"), plt.list[[5]], plt.list[[6]] + rremove("ylab"),
                             labels = NULL,
                             ncol = 3, nrow = 2,
                             common.legend = TRUE, legend = "right",
                             align = "h", 
                             font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"), heights = c(1,1,1.5))

plt.all.biofilm ## combined plot

# For pathogenic species

ps_pathogens.removal.noNA <- subset_samples(ps_pathogens.removal, !is.na(sample_data(ps_pathogens.removal)$Plot) & !is.na(sample_data(ps_pathogens.removal)$Month))
ps_pathogens.removal.noNA.merged <- merge_by_plotdate(ps_pathogens.removal.noNA)
ps_pathogens.removal.noNA.merged <- transform_sample_counts(ps_pathogens.removal.noNA.merged, function(otu) otu/sum(otu))

ps_pathogens.removal.noNA.merged <- convert_treatment(ps_pathogens.removal.noNA.merged)

ps_pathogens.removal.noNA.merged.noNA <- subset_samples(ps_pathogens.removal.noNA.merged, rownames(otu_table(ps_pathogens.removal.noNA.merged)) != "71442023")
ps_pathogens.removal.noNA.merged.glom <- tax_glom(ps_pathogens.removal.noNA.merged.noNA, taxrank = "Genus")
ps_pathogens.removal.noNA.merged.glom <- psmelt(ps_pathogens.removal.noNA.merged.glom)
bact.abun.path <- ps_pathogens.removal.noNA.merged.glom[,colnames(ps_pathogens.removal.noNA.merged.glom) %in% c("Treatment", "Plot", "Month", "Abundance", "Genus", "PlotDate")]
bact.abun.path$Plot <- as.factor(bact.abun.path$Plot)
bact.abun.path$Month <- as.factor(bact.abun.path$Month)

plt.list.pathogen <- list()

list_pathogens_present <- list("Pseudoalteromonas", "Pseudomonas", "Tenacibaculum", "Vibrio", "Photobacterium", "Shewanella", "Bacillus")

for (bact in list_pathogens_present) {
  abun.data <- bact.abun.path[bact.abun.path$Genus == bact,]
  plt <- ggplot(data = abun.data, 
                aes(x = Month, y = Abundance, alluvium = Plot)) + 
    geom_alluvium(aes(fill = Treatment, colour = Treatment), alpha = .75, decreasing = FALSE) +
    theme_minimal() +
    ylab("Relative<br>abundance") +
    theme(axis.text=element_text(size=8),
          axis.title=element_text(size=10,face="bold")) + theme(legend.key.size = unit(0.4, "cm")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(axis.line = element_line(colour = "black")) +
    theme(legend.title=element_text(size=9, face = "bold"), legend.text=element_text(face = "italic")) + 
    theme(axis.title.y = element_markdown(), legend.title = element_markdown(), axis.text.x = element_markdown(), legend.title.align = 0.5) +
    ggtitle(bact) + theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold")) + theme(axis.text.y = element_blank())
  plt.list.pathogen <- append(plt.list.pathogen, list(plt))
}


require(grid)

plt.all.pathogen <- ggarrange(plt.list.pathogen[[1]], plt.list.pathogen[[2]] + rremove("ylab"), plt.list.pathogen[[3]] + rremove("ylab"), plt.list.pathogen[[4]] + rremove("ylab"),
                              plt.list.pathogen[[5]], plt.list.pathogen[[6]] + rremove("ylab"), plt.list.pathogen[[7]] + rremove("ylab"),# remove axis labels from plots
                              ncol = 4, nrow = 2,
                              common.legend = TRUE, legend = "right",
                              align = "h", 
                              font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"), heights = c(1,1,1.5))

plt.all.pathogen ## combined plot

#### ANOVA ####

calculate_anova <- function(pseq, filename) {
  pseq.glom <- tax_glom(pseq, "Genus")
  pseq.df <- psmelt(transform_sample_counts(pseq.glom, function(otu) otu/sum(otu)))
  
  abundance <- function() {
    genera <- unique(pseq.df$Genus)
    for (gi in genera) {
      print(paste0("Significance results for ", gi))
      abundance.data <- pseq.df[pseq.df$Genus == gi,]
      abundance.data <- abundance.data[,c("Treatment", "Abundance")]
      abundance.aov <- aov(Abundance ~ Treatment, data = abundance.data)
      print(summary(abundance.aov))
    }
  }
  capture.output(abundance(), file = filename)
}

calculate_anova(ps_pathogens.removal, "Relative Abundance Pathogens Anova Results.txt")
calculate_anova(ps_biofilm.removal, "Relative Abundance Biofilm Anova Results.txt")

#### Heat Map ####

processing_heatmap_removal <- function(pseq.removal) {
  pseq.removal.merged <- merge_by_plotdate(pseq.removal)
  pseq.removal.merged <- convert_treatment(pseq.removal.merged)
  pseq.removal.merged <- aggregate_rare(pseq.removal.merged, level = "Genus", prevalence = 0, detection = 10)
  pseq.removal.merged <- microbiome::transform(pseq.removal.merged, "compositional")
  pseq.removal.merged.df <- psmelt(pseq.removal.merged)
}

plotting_heatmap_removal <- function(pseq.df, x, y) {
  order.pseq <- sort(unique(pseq.df$OTU))
  order.pseq <- order.pseq[order.pseq != "Other"]
  order.pseq <- append(order.pseq, "Other")
  plt <- ggplot(pseq.df, aes(Sample, OTU)) + facet_nested(~Treatment + Plot, nest_line = TRUE, scales = "free", switch = "x") +
    geom_tile(aes(fill = Abundance)) +
    scale_fill_gradient(limits=c(0,1), breaks=seq(0, 1,by = +0.25)) +  
    scale_fill_continuous(high = "#00008B", low = "#FFFFFF", limits=c(0,1), breaks=seq(0, 1,by = +0.25)) + 
    ylab(substitute(paste(bold(y)))) + 
    xlab(substitute(paste(bold(x)))) + 
    guides(fill=guide_legend(title="Relative abundance")) + 
    theme(axis.title.x = element_text(vjust=+3), axis.title.y = element_text(vjust= +1)) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
          panel.background = element_blank(), axis.ticks = element_blank()) + 
    scale_y_discrete(limits = rev(order.pseq)) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(panel.spacing.x = unit(0,"line")) + 
    theme(strip.text.x = element_text(colour = "black", size = 6), strip.background = element_blank()) +
    theme(axis.title = element_text(size = 10)) +
    scale_x_discrete(expand = c(0,0))
  
  grid.draw(plt)
}

# For biofilm species:

ps_biofilm.removal.noNA <- subset_samples(ps_biofilm.removal, !is.na(sample_data(ps_biofilm.removal)$Plot) & !is.na(sample_data(ps_biofilm.removal)$Month))
ps_biofilm.removal.noNA.heatmap <- processing_heatmap_removal(ps_biofilm.removal.noNA)
plotting_heatmap_removal(ps_biofilm.removal.noNA.heatmap, "Plot", "Known biofilm-making bacteria")

# For pathogenic species:

ps_pathogens.removal.noNA <- subset_samples(ps_pathogens.removal, !is.na(sample_data(ps_pathogens.removal)$Plot) & !is.na(sample_data(ps_pathogens.removal)$Month))
ps_pathogens.removal.noNA.heatmap <- processing_heatmap_removal(ps_pathogens.removal.noNA)
plotting_heatmap_removal(ps_pathogens.removal.noNA.heatmap, "Plot", "Pathogenic bacteria")

#### Boxplots ####

# Pathogens at 8 sites:

ps_all.rel.8site.noNA <- subset_taxa(ps_all.rel.8site, !is.na(Genus))
ps_all.rel.8site.noNA.pathogens <- subset_taxa(ps_all.rel.8site.noNA, Genus %in% list_pathogens)
ps_all.rel.8site.noNA.pathogens.glom <- tax_glom(ps_all.rel.8site.noNA.pathogens, "Kingdom")

ps_all.rel.8site.noNA.pathogens.glom <- convert_treatment(ps_all.rel.8site.noNA.pathogens.glom)

ps_all.rel.8site.noNA.pathogens.glom.log <- transform_sample_counts(ps_all.rel.8site.noNA.pathogens.glom, log)
ps_all.rel.8site.noNA.pathogens.glom.log.df <- psmelt(ps_all.rel.8site.noNA.pathogens.glom.log)
ps_all.rel.8site.noNA.pathogens.glom.log.df$Treatment <- factor(ps_all.rel.8site.noNA.pathogens.glom.log.df$Treatment, levels = c("+", "-"), ordered = TRUE)

pathogen_boxplot <- ggplot(ps_all.rel.8site.noNA.pathogens.glom.log.df) + geom_boxplot(aes(x = Treatment, y = Abundance, fill = Treatment)) + theme_classic() + xlab("Seagrass") + ylab("Log-transformed Relative Abundance") + theme(axis.title = element_text(face = "bold")) + guides(fill=guide_legend(title="Seagrass")) + scale_fill_manual(values = c("#00BFC4", "#F8766D"))

pathogen_boxplot

wilcox.test(Abundance ~ Treatment,
            data = ps_all.rel.8site.noNA.pathogens.glom.log.df)
wilcox.test(Abundance ~ Treatment,
            data = ps_all.rel.8site.noNA.pathogens.glom.df)

#### SIMPER ####

# For 8 sites:

ps_pathogens.8site.merged <- merge_by_loctreat(ps_pathogens.8site)
ps_pathogens.8site.merged <- convert_treatment(ps_pathogens.8site.merged)

ps_pathogens.8site.merged.noNA <- subset_taxa(ps_pathogens.8site.merged, !is.na(Genus))
ps_pathogens.8site.merged.noNA.glom <- tax_glom(ps_pathogens.8site.merged.noNA, "Genus")
ps_pathogens.8site.merged.noNA.glom.rel <- transform_sample_counts(ps_pathogens.8site.merged.noNA.glom, function(x) x/sum(x))

ps_pathogens.8site.abundance.matrix <- as(otu_table(ps_pathogens.8site.merged.noNA.glom.rel), "matrix")
colnames(ps_pathogens.8site.abundance.matrix) <- tax_table(ps_pathogens.8site.merged.noNA.glom.rel)[,"Genus"]

sim <- with(sample_data(ps_pathogens.8site.merged.noNA.glom.rel), simper(ps_pathogens.8site.abundance.matrix, Treatment))
summary(sim)


# For core microbiomes at 8 sites:

pseq.core.8site.glom <- tax_glom(pseq.core.8site, "Genus")

pseq.core.8site.glom <- convert_treatment(pseq.core.8site.glom)

pseq.core.8site.glom.rel <- transform_sample_counts(pseq.core.8site.glom, function(x) x/sum(x))
pseq.core.8site.abundance.matrix <- as(otu_table(pseq.core.8site.glom.rel), "matrix")
colnames(pseq.core.8site.abundance.matrix) <- tax_table(pseq.core.8site.glom.rel)[,"Genus"]

sim <- with(sample_data(pseq.core.8site.glom.rel), simper(pseq.core.8site.abundance.matrix, Treatment))
summary(sim)

#### Alpha Diversity ####

# Calculating the alpha diversities

alpha_diversity.all <- estimate_richness(ps_all, measures = c("shannon", "simpson")) ## calculating diversity indices 
meta_withH <- merge(meta, alpha_diversity.all, by = 0) ## merging with the metadata
rownames(meta_withH) <- meta_withH$SampleID ## rename each row with sampleID
for (sample in rownames(meta_withH)) {
  if (is.na(meta_withH[sample,]$Treatment)) {
    next
  }
  if (meta_withH[sample,]$Treatment == 0) {
    meta_withH[sample,]$Treatment <- "+"
  }
  else {
    meta_withH[sample,]$Treatment <- "-"
  }
}

# Dot plot for removal plots

meta_withH.removal <- meta_withH[meta_withH$Study == "removal_exp_v3",]
meta_withH.removal.noNA <- meta_withH.removal[!is.na(meta_withH.removal$Plot) && !is.na(meta_withH.removal$Month),]
meta_withH.removal.noNA$PlotDate <- paste0(meta_withH.removal.noNA$Plot, meta_withH.removal.noNA$Day, meta_withH.removal.noNA$Month, meta_withH.removal.noNA$Year)

data_summary <- function(data, varname, groupnames){ ## this function gives a summary of the data table
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE)),
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

shannon_data.removal <- meta_withH.removal.noNA[,c("Shannon", "Plot", "Day", "Month", "Year", "Treatment", "PlotDate")]
simpson_data.removal <- meta_withH.removal.noNA[,c("Simpson", "Plot", "Day", "Month", "Year", "Treatment", "PlotDate")]

shannon_data.summary <- data_summary(shannon_data.removal, varname="Shannon", 
                                     groupnames=c("Plot", "Month", "Treatment"))
shannon_data.summary$Month <- as.factor(shannon_data.summary$Month)

simpson_data.summary <- data_summary(simpson_data.removal, varname="Simpson", 
                                     groupnames=c("Plot", "Month", "Treatment"))
simpson_data.summary$Month <- as.factor(simpson_data.summary$Month)

plt_shannon.removal <- ggplot(shannon_data.summary, aes(x=Month, y=Shannon, color=Treatment), group = Plot) +
  geom_point(size = 2, position=position_dodge(0.05)) + geom_errorbar(aes(ymin=Shannon-sd, ymax=Shannon+sd), width=.05, position=position_dodge(0.05)) +
  theme_classic()

plt_simpson.removal <- ggplot(simpson_data.summary, aes(x=Month, y=Simpson, color=Treatment), group = Plot) +
  geom_point(size = 2, position=position_dodge(0.05)) + geom_errorbar(aes(ymin=Simpson-sd, ymax=Simpson+sd), width=.05, position=position_dodge(0.05)) +
  theme_classic()

plt_shannon.removal
plt_simpson.removal

# Dot plot for 8 sites

meta_withH.8site <- meta_withH[meta_withH$Study == "enterodat21",]
meta_withH.8site.noNA <- meta_withH.8site[!is.na(meta_withH.8site$Location) && !is.na(meta_withH.8site$Treatment),]
meta_withH.8site.noNA$LocTreat <- paste0(meta_withH.8site.noNA$Location, meta_withH.8site.noNA$Treatment)

shannon_data.8site <- meta_withH.8site.noNA[,c("Shannon", "Location", "Treatment", "LocTreat")]
simpson_data.8site <- meta_withH.8site.noNA[,c("Simpson", "Location", "Treatment", "LocTreat")]

shannon_data.summary.8site <- data_summary(shannon_data, varname="Shannon", 
                                       groupnames=c("Location", "Treatment"))

plt_shannon.8site <- ggplot(shannon_data.summary.8site, aes(x=Location, y=Shannon, color=Treatment)) +
  geom_point(size = 2, position=position_dodge(0.3)) + geom_errorbar(aes(ymin=Shannon-sd, ymax=Shannon+sd), width=0, position=position_dodge(0.3)) +
  theme_classic()

simpson_data.summary.8site <- data_summary(simpson_data, varname="Simpson", 
                                       groupnames=c("Location", "Treatment"))

plt_simpson.8site <- ggplot(simpson_data.summary.8site, aes(x=Location, y=Simpson, color=Treatment)) +
  geom_point(size = 2, position=position_dodge(0.3)) + geom_errorbar(aes(ymin=Simpson-sd, ymax=Simpson+sd), width=0, position=position_dodge(0.3)) +
  theme_classic()

plt_shannon.8site
plt_simpson.8site

# Violin plot for 8 sites

violin.plt_shannon.8site <-
  ggbetweenstats(
    plot.type = "violin", 
    data = shannon_data,
    x = LocTreat,
    y = Shannon, bf.message = FALSE,
    ggsignif.args = list(textsize = 3.5, step_increase = 0.02, tip_length = 0.01)) + 
  scale_fill_brewer(palette = "Paired") + theme(axis.ticks = element_blank(),
                                                axis.line = element_line(colour = "grey50"),
                                                panel.grid = element_line(color = "#b4aea9"),
                                                panel.grid.minor = element_blank(),
                                                panel.grid.major.x = element_blank(),
                                                panel.grid.major.y = element_line(linetype = "dashed")) + 
  ylab("Shannon index") + ylim(c(3.5, 8))

violin.plt_simpson.8site <- ggbetweenstats(
  plot.type = "violin", 
  data = simpson_data,
  x = Location,
  y = Simpson, bf.message = FALSE, 
  ggsignif.args = list(textsize = 3.5, step_increase = 0.08, tip_length = 0.01)) + 
  scale_fill_npg() + scale_color_npg() + theme(axis.ticks = element_blank(),
                                               axis.line = element_line(colour = "grey50"),
                                               panel.grid = element_line(color = "#b4aea9"),
                                               panel.grid.minor = element_blank(),
                                               panel.grid.major.x = element_blank(),
                                               panel.grid.major.y = element_line(linetype = "dashed")) + 
  ylab("Simpson index") + ylim(c(0, 1.5))

violin.plt_shannon.8site
violin.plt_simpson.8site

# Box plot for Shannon Index over time - 8 sites

meta_withH.removal.noNA$Month <- as.character(meta_withH.removal.noNA$Month)
meta_withH.removal.noNA$Plot <- as.character(meta_withH.removal.noNA$Plot)

for (sample in rownames(meta_withH.removal.noNA)) {
  if (meta_withH.removal.noNA[sample,]$Treatment == "+") {
    meta_withH.removal.noNA[sample,]$Treatment <- "Control"
  }
  else {
    meta_withH.removal.noNA[sample,]$Treatment <- "Removal"
  }
}

ggplot(meta_withH.removal.noNA, aes(x=Month, y=Shannon)) + facet_nested_wrap(~ Treatment + Plot, nrow = 2, ncol = 4, nest_line = element_line(linetype = 1)) +
  geom_boxplot(width = 0.5) +
  geom_dotplot(aes(fill = Plot), dotsize = 1, binaxis='y', stackdir='center', position = "jitter") +
  theme_classic() + theme(strip.background = element_blank(),
                          ggh4x.facet.nestline = element_line(colour = "black"))

# ANOVA

sh.aov <- aov(Shannon ~ Treatment, data = meta_withH.removal)
capture.output(summary(sh.aov), file = "./Shannon Removal Anova Results.txt")

meta_withH.8site.noNA <- meta_withH.8site[!is.na(meta_withH.8site$Shoots),]
meta_withH.8site.noNA <- meta_withH.8site.noNA[meta_withH.8site.noNA$Treatment != 1,]
sh.aov <- aov(MPN ~ Location, data = meta_withH.8site.noNA)
capture.output(summary(sh.aov), TukeyHSD(sh.aov), file = "MPN Anova Results.txt")

# Linear regressions

shoots_shannon_corr <- lm(Shannon~Shoots, data = meta_withH.8site) ## shannon vs shoots
summary(shoots_shannon_corr)
plot(meta_withH.8site$Shoots,  meta_withH.8site$Shannon, pch = 16, cex = 1.3, col = "blue", main = "Shannon Index Against Shoot Density", xlab = "Shoot Density", ylab = "Shannon Index")
abline(shoots_shannon_corr)

shoots_simpson_corr <- lm(Simpson~Shoots, data = meta_withH.8site) ## simpson vs shoots
summary(shoots_simpson_corr)
plot(meta_withH.8site$Shoots,  meta_withH.8site$Simpson, pch = 16, cex = 1.3, col = "blue", main = "Simpson Index Against Shoot Density", xlab = "Shoot Density", ylab = "Simpson Index")
abline(shoots_simpson_corr)

mpn_shannon_corr <- lm(MPN~Shannon, data = meta_withH.8site) ## MPN vs shannon
summary(mpn_shannon_corr)
plot(meta_withH.8site$Shannon,  meta_withH.8site$MPN, pch = 16, cex = 1.3, col = "blue", main = "MPN of Enterococcus Against Shannon Index", xlab = "Shannon Index", ylab = "MPN of Enterococcus")
abline(mpn_shannon_corr)

mpn_simpson_corr <- lm(MPN~Simpson, data = meta_withH.8site) ## MPN vs simpson
summary(mpn_simpson_corr)
plot(meta_withH.8site$Simpson,  meta_withH.8site$MPN, pch = 16, cex = 1.3, col = "blue", main = "MPN of Enterococcus Against Simpson Index", xlab = "Simpson Index", ylab = "MPN of Enterococcus")
abline(mpn_simpson_corr)

shoots_MPN_corr <- lm(MPN~Shoots, data = meta_withH.8site) ## MPN vs shoots
summary(shoots_MPN_corr)
plot(meta_withH.8site$Shoots,  meta_withH.8site$MPN, pch = 16, cex = 1.3, col = "blue", main = "MPN Against Shoot Density", xlab = "Shoot Density", ylab = "MPN")
abline(shoots_MPN_corr)

# PERMANOVA

meta_withH.8site.noNA <- meta_withH.8site[!is.na(meta_withH.8site$MPN),] ## for MPN across locations
meta_withH.8site.noNA$MPN <- as.numeric(meta_withH.8site.noNA$MPN)
meta_withH.8site.noNA <- meta_withH.8site.noNA[meta_withH.8site.noNA$MPN != 0,]
MPN_data <- as.data.frame(meta_withH.8site.noNA$MPN)
rownames(MPN_data) <- rownames(meta_withH.8site.noNA)
MPN_data.noNA <- MPN_data[!is.na(MPN_data$MPN),]
permanova <- adonis2(MPN_data ~ Location,
                     data = meta_withH.8site.noNA, permutations=99, method = "bray")
MPN_loc <- as.data.frame(meta_withH.8site.noNA %>% group_by(Location) %>% summarise(across(c(MPN), mean)))
permanova <- adonis2(t(MPN_loc$MPN) ~ MPN_loc$Location,
                     data = MPN_loc, permutations=99, method = "bray")
print(as.data.frame(permanova$aov.tab)["Location", "Pr(>F)"])

meta_withH.8site.noNA <- meta_withH.8site[!is.na(meta_withH.8site$Shoots),] ## for shoots across locations
meta_withH.8site.noNA$Shoots <- as.numeric(meta_withH.8site.noNA$Shoots)
meta_withH.8site.noNA <- meta_withH.8site.noNA[meta_withH.8site.noNA$Shoots != 0,]
Shoots_data <- as.data.frame(meta_withH.8site.noNA$Shoots)
rownames(Shoots_data) <- rownames(meta_withH.8site.noNA)
permanova <- adonis2(Shoots_data ~ Location,
                     data = meta_withH.8site.noNA, permutations=99, method = "bray")

permanova

# GLM

pseq.core.8site <- core(ps_all.rel.8site, detection = 0, prevalence = .3)
pseq.core.8site <- subset_samples(pseq.core.8site, SampleID != "CH11-LEK7958_L1" & SampleID != "CH21-LEK7959_L1" & SampleID !="CH31-LEK7960_L1") #Nothing detected in the core
alpha_diversity.core <- estimate_richness(pseq.core.8site, measures = c("shannon", "simpson"))

meta_withH_core <- merge(meta, alpha_diversity.core, by = 0)
rownames(meta_withH_core) <- meta_withH_core$SampleID
meta_withH_core$MPN <- as.numeric(meta_withH_core$MPN)

m1 <- glm(MPN ~ Shoots, data = meta_withH_core)
m2 <- glm(MPN ~ Shannon, data = meta_withH_core)
m3 <- glm(MPN ~ Shoots + Shannon, data = meta_withH_core)
m4 <- glm(MPN ~ Shoots*Shannon, data = meta_withH_core)

ms <- list(MPN_Shoots = m1, MPN_Shannon = m2, MPN_Shoots_plus_Shannon = m3, MPN_Shoots_times_Shannon = m4)
tbl <- aictab(ms)

bind_rows(lapply(ms, tidy), .id="key")

screenreg(ms)

# Shoots Mean

meta_withH.8site.noNA <- meta_withH.8site[!is.na(meta_withH.8site$Shoots),]
meta_withH.8site.noNA <- meta_withH.8site.noNA[meta_withH.8site.noNA$Treatment != 1,]
summarize(group_by(meta_withH.8site.noNA, Location), mean_shoots= mean(Shoots))

#### Low Abundance Species ####

ps_all.rel.8site.merged <- merge_by_loctreat(ps_all.rel.8site)
ps_all.rel.8site.merged.df <- psmelt(ps_all.rel.8site.merged)
ps_all.rel.8site.merged.df$Genus <- as.character(ps_all.rel.8site.merged.df$Genus)
ps_all.rel.8site.merged.df.lowabun <- ps_all.rel.8site.merged.df[ps_all.rel.8site.merged.df$Abundance < 0.05,]
ps_all.rel.8site.merged.df.lowabun <- ps_all.rel.8site.merged.df.lowabun[ps_all.rel.8site.merged.df.lowabun$Abundance > 0.035,]

lowabun.8site.means <- ps_all.rel.8site.merged.df.lowabun %>% group_by(Genus) %>% 
  summarise(Mean.Abundance=mean(Abundance),
            .groups = 'drop')
lowabun.8site.means <- lowabun.8site.means[order(lowabun.8site.means$Mean.Abundance, decreasing = TRUE),]
order.lowabun.8site <- lowabun.8site.means$Genus
order.lowabun.8site <- order.lowabun.8site[!order.lowabun.8site == "Others"]
order.lowabun.8site <- append("Others", order.lowabun.8site)

plt.8site.lowabun <- ggplot(ps_all.rel.8site.merged.df.lowabun, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = order.lowabun.8site))) +
  geom_bar(position = "stack", stat = "identity", width = 1) +
  facet_nested(~Location + Treatment, nest_line = TRUE, scales = "free", switch = "x") + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(panel.spacing.x = unit(0,"line")) +
  ylab(substitute(paste(bold("Relative abundance")))) + scale_fill_brewer("Genus", palette = "Paired") +
  xlab(substitute(paste(bold("Location")))) + theme(strip.text.x = element_text(colour = "black", size = 6), strip.background = element_blank()) +
  theme(axis.title = element_text(size = 10)) + scale_x_discrete(expand = c(0,0))

plt.8site.lowabun

#### Combined Plot for Pathogenic Taxas (Archaic) ####

# This code is outdated, due to changes in the data structures and naming
# conventions. However, this could be salvaged, knowing that:
# ps2.entero = ps_pathogens.8site
# toplot = most of the data can be found in meta_withH

library(stringr) 
toplot <- data.frame(LocTreat = shannon_toplot$LocTreat, Location = NA, Treatment = NA,
                     H = shannon_toplot$mean, Hmin = shannon_toplot$min, 
                     Hmax = shannon_toplot$max,
                     meanShoot = c(55.39285714,204.35,0,18.6,35.18181818,0,54.25,0,7.769230769,0,0.783783784,0,33.60344828,0),
                     sdShoot = c(51.47911427,102.7061598,0,11.44582442,30.66905458,0,49.38264028,0,9.435340006,0,2.439662262,0,30.2358795,0))
toplot$Location <- str_sub(toplot$LocTreat, end = -2)
toplot$Treatment <- str_sub(toplot$LocTreat, start = -1)
toplot <- toplot[toplot$LocTreat != "Cyrene-",] #Cyrene- no data for pathogen abundance

sample_data(ps2.entero)$LocTreat <- paste0(sample_data(ps2.entero)$Location, sample_data(ps2.entero)$Treatment)

ps2.merged <- merge_samples(ps2.entero, c("LocTreat"))
ps2.merged <- transform_sample_counts(ps2.merged, function(otu) otu/sum(otu))
ps2.merged <- subset_samples(ps2.merged, rownames(sample_data(ps2.merged)) != "Cyrene1") #Cyrene- no data for pathogen abundance
ps2.merged <- tax_glom(ps2.merged, "Genus")
ps2.merged <- psmelt(ps2.merged)

for (row in rownames(ps2.merged)) {
  if (ps2.merged[row, "Treatment"] == 0) {
    ps2.merged[row, "Treatment"] <- "+"
  }
  else {
    ps2.merged[row, "Treatment"] <- "-"
  }
}
ps2.merged$Location <- str_sub(ps2.merged$Sample, end = -2)
ps2.merged$LocTreat <- paste0(ps2.merged$Location, ps2.merged$Treatment)
pathogen.abundance <- ps2.merged[,colnames(ps2.merged) %in% c("LocTreat", "Abundance", "Genus")]
toplot2 <- merge(toplot, pathogen.abundance, by = "LocTreat")
toplot2 <- toplot2 %>% arrange(desc(meanShoot))
toplot2.pos <- toplot2[toplot2$Treatment == "+",]
toplot2 <- toplot2 %>% arrange(meanShoot)
toplot2.neg <- toplot2[toplot2$Treatment == "-",]

desired_order <- c('Changi', 'Bendera Bay', 'Eagle Bay', 'Cyrene', 'Tanah Merah', 'Chek Jawa', 'Semakau', 'Sentosa')
desired_order_neg <- c('Changi', 'Eagle Bay', 'Tanah Merah', 'Semakau', 'Sentosa')
pos.sample <- c(1, 15, 29, 43, 57, 71, 85, 99)
neg.sample <- c(1, 15, 29, 43, 57)
p1 <- ggplot(toplot2.neg[neg.sample,]) +
  geom_bar(aes(x=Location, y=meanShoot), stat="identity",fill="skyblue",alpha=0.7) +
  geom_errorbar(aes(x=Location,y=meanShoot,ymin=meanShoot-sdShoot,ymax=meanShoot+sdShoot), width=0.1, size=0.1) + 
  theme_minimal() + theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + 
  ylab("Mean<br>shoot density") + theme(axis.title.y = element_markdown()) + 
  scale_x_discrete(limits = desired_order_neg) +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10,face="bold")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black")) + ylim(-5,350) +
  geom_text(aes(x=Location, label=round(meanShoot,2), y=meanShoot+sdShoot), vjust = -1.1, size=2.5)

p2 <- ggplot(data=toplot2.neg, aes(x=Location, y=Abundance, group=Genus, fill=Genus)) + 
  geom_line(aes(col = Genus), linewidth=0.1) + 
  geom_point(aes(col = Genus), size = 1.5, alpha = 0.5) + scale_color_hue() + theme_minimal() +
  ylab("Pathogen<br>abundance") + xlab("Location") +
  scale_x_discrete(limits = desired_order_neg) +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10,face="bold")) + theme(legend.key.size = unit(0.4, "cm")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  theme(legend.title=element_text(size=9, face = "bold"), legend.text=element_text(face = "italic")) + 
  theme(axis.title.y = element_markdown(), legend.title = element_markdown(), axis.text.x = element_markdown(), legend.title.align = 0.5)

p3 <- ggplot(toplot2.neg[neg.sample,]) +
  geom_pointrange(aes(x=Location,y=H,ymin=Hmin, ymax=Hmax), size = 0.3, color = "darkred") + 
  theme_minimal() + theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + 
  ylab("Bacterial<br>diversity") + 
  scale_x_discrete(limits = desired_order_neg) +
  theme(axis.title.y = element_markdown()) +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10,face="bold")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black"))

require(grid)
plt_combined <- ggarrange(p3 + rremove("xlab"), p2, # remove axis labels from plots
                          labels = NULL,
                          ncol = 1, nrow = 2,
                          common.legend = TRUE, legend = "right",
                          align = "v", 
                          font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"), heights = c(1,1.5))

plt_combined
dev.off() 
options("devEval/args/path"=file.path("."))
devEval("tiff", name="plot1", width=120, height=800, {
  plt_combined;
})


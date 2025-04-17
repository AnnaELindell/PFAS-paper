library(dada2)

#######################################################
### Non-quantitative analysis of all samples  #########
#######################################################


path <- './data_nq/'
list.files(path)

# Forward and reverse fastq filenames 
fnFs <- sort(list.files(path, pattern="_1_sequence.txt.gz",
                        full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2_sequence.txt.gz",
                        full.names = TRUE))

# Extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 5)

## Inspect read quality profiles
# Forward reads
plotQualityProfile(fnFs[1]) # over 30 is good quality, green line is mean quality score
# Reverse reads
plotQualityProfile(fnRs[3])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

## Change filter parameters for different primer set, default 16S from EMBL
## would be 515F - 806R, therefore the 'truncLen = c(240, 160)' is recommended
## 
## Here the test dataset has a peak around ~410 bp (unknown primer set), took 
## ~20s to run with 3 samples on Mac with Apple M1 Pro chip, 16Gb
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(230,210),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
# On Windows set multithread=FALSE
head(out)

write.table(out, "output_nq/M4_16S_out.txt", quote = F)


# Learn the error rate, took ~90s to run for each 
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)


## Sample inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]


## merge paired reads 
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers[[1]])
write.table(mergers, "output_nq/M4_16S_mergers.txt", sep = "\t", dec = ".", na = "NA",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

## construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
saveRDS(seqtab, "output_nq/M4_16S_seqtab.rds")
# should be 3 2025 using 240, 210 as trunc parameters

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))



# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus",
                                    multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
saveRDS(seqtab.nochim, "output_nq/M4_16S_seqtab.nochim.rds")

sum(seqtab.nochim)/sum(seqtab)



# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN),
               sapply(dadaRs, getN),
               sapply(mergers, getN),
               rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls:
# e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR",
                     "merged", "nonchim")
rownames(track) <- sample.names
head(track) # numbers in output are numbers of reads

write.table(track, "output_nq/M4_16S_track.txt", quote = F)

track <- read.table("output_nq/M4_16S_track.txt", sep = " ", header = TRUE, check.names = FALSE)
sum(track$input)
mean(track$input)


## Hand off to phyloseq
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")

theme_set(theme_bw())







# Assign taxonomy with custom database, ~207s to run, ~3.5 min
seqtab.nochim <- readRDS("output_nq/M4_16S_seqtab.nochim.rds")

taxa <- assignTaxonomy(seqtab.nochim,
                       "./community_data/Low_High_sonja.txt",
                       multithread=TRUE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print) 
length(taxa.print)

write.table(taxa, "output_nq/M4_16S_taxa.txt", sep = "\t")




# read in metadata:
metadata <- read.table("m4_16s_metadata.txt", sep = "\t", header = TRUE, check.names = FALSE)
rownames(metadata) <- metadata$"16S_ID"
metadata$time.point.d <- as.character(metadata$time.point.d)


# create phyloseq object linking tables
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(metadata), 
               tax_table(taxa))



#ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
saveRDS(ps, "output_nq/M4_16S_ps.rds")

library(tidyverse)
library(dplyr)

## Plot alpha-diversity
plot_richness(ps, x="mouse",
              measures=c("Shannon", "Observed"), color="group", shape = "tissue")
ggsave("figures_nq/alpha_diversity.png",
       width = 35, height = 10, unit = "cm")




# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

plot_ordination(ps.prop, ord.nmds.bray, color="group", shape = "tissue", title="Bray NMDS")
ggsave("figures_nq/bray_curtis_distance.png",
       width = 15, height = 10, unit = "cm")




## Taxonomic profile
top30 <- names(sort(taxa_sums(ps), decreasing=TRUE)) #[1:30]
ps.top30 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top30 <- prune_taxa(top30, ps.top30)

my_plot_bar = function (physeq, x = "Sample", y = "Abundance", fill = NULL, title = NULL, 
                        facet_grid = NULL) {
  mdf = psmelt(physeq)
  p = ggplot(mdf, aes_string(x = x, y = y, fill = fill))
  p = p + geom_bar(stat = "identity", position = "stack")
  p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  if (!is.null(facet_grid)) {
    p <- p + facet_grid(facet_grid)
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}


my_plot_bar(ps.top30, x="time_tissue", fill="Species")  + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.title.x = element_blank()) +
  facet_wrap(~group, scales="free")
ggsave("figures_nq/abundance_by_treatment_species.png",
       width = 15, height = 15, unit = "cm")

my_plot_bar(ps.top30, x="time_tissue", fill="Family")  + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.title.x = element_blank()) +
  facet_wrap(~group, scales="free")
ggsave("figures_nq/abundance_by_treatment_family.png",
       width = 15, height = 15, unit = "cm")

my_plot_bar(ps.top30, x="time_tissue", fill="Phylum")  + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.title.x = element_blank()) +
  facet_wrap(~group, scales="free")
ggsave("figures_nq/abundance_by_treatment_phylum.png",
       width = 15, height = 15, unit = "cm")





### Analysis of PFNA samples

data <- readRDS("output_nq/M4_16S_ps.rds")

mdf = psmelt(data)
mdf[is.na(mdf)] <- "not assigned"


total_count <- mdf %>%
  group_by(mouse, time_tissue) %>%
  summarize(total_count = sum(Abundance))

mdf <- merge(mdf, total_count, by = c("mouse", "time_tissue"))

mdf <- mdf %>%
  mutate(Count = Abundance,
         Abundance = Count / total_count)
mdf.sum <- mdf %>%
  group_by(mouse, group, time_tissue, Kingdom, Phylum, Class, Order, Family, Genus, Species, total_count) %>%
  summarise(abs.Abundance = sum(Count),
            rel.Abundance = sum(Abundance))

mdf.sum$Species <- as.factor(mdf.sum$Species)
levels(mdf.sum$Species)
mdf.sum$Species <- factor(mdf.sum$Species, 
                          levels = c("Akkermansia muciniphila", "Collinsella aerofaciens", "Enterocloster bolteae", "Escherichia coli", "Ruminococcus gnavus",
                                     "Bacteroides fragilis", "Bacteroides thetaiotaomicron", "Bacteroides uniformis", "Lacrimispora saccharolytica", "Phocaeicola vulgatus", "not assigned"))


# Check germfree mice

gf <- mdf.sum %>% filter(group == "GF untreated")
ggplot(gf, aes_string(x = "time_tissue", y = "abs.Abundance", color = "mouse"))  + 
  geom_jitter(width = 0.1, height = 0, shape = 1, size = 3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        text = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 10, angle = 90, , hjust = 1, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 10)) +
  facet_wrap(~Species, scales="free_y", ncol = 5)
ggsave("figures_nq/absolute_abundance_gf.png",
       width = 25, height = 15, unit = "cm")

ggplot(gf, aes_string(x = "time_tissue", y = "rel.Abundance", color="mouse"))  + 
  geom_jitter(width = 0.1, height = 0, shape = 1, size = 3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        text = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 10, angle = 90, , hjust = 1, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 10)) +
  facet_wrap(~Species, scales="free_y", ncol = 5)
ggsave("figures_nq/relative_abundance_gf.png",
       width = 25, height = 15, unit = "cm")

# All GF samples contain A. muciniphila
# A118 SI contaminated with B. theta and P. vulgatus
# A119 and A118 contaminated with R. gnavus



# LC and HC small intestine

mdf.sum1 <- mdf.sum %>% 
  filter(group %in% c("HC + PFNA", "LC + PFNA"), 
         time_tissue == "3_Small intestine")

ggplot(mdf.sum1, aes_string(x = "group", y = "abs.Abundance"))  + 
  geom_jitter(width = 0.1, height = 0, shape = 1, size = 3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        text = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 10, angle = 90, , hjust = 1, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 10)) +
  facet_wrap(~Species, scales="free_y", ncol = 5)
ggsave("figures_nq/absolute_abundance_si.png",
       width = 25, height = 15, unit = "cm")

ggplot(mdf.sum1, aes_string(x = "group", y = "rel.Abundance"))  + 
  geom_jitter(width = 0.1, height = 0, shape = 1, size = 3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        text = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 10, angle = 90, , hjust = 1, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 10)) +
  facet_wrap(~Species, scales="free_y", ncol = 5)
ggsave("figures_nq/relative_abundance_si.png",
       width = 25, height = 15, unit = "cm")


# All HC samples clean apart from A. muciniphila, which is present in all samples.
# One LC sample contaminated with B. uniformis to low levels



# Feces & colon

mdf.sum2 <- mdf.sum %>% 
  filter(group %in% c("HC + PFNA", "LC + PFNA"), 
         time_tissue != "3_Small intestine")

ggplot(mdf.sum2, aes_string(x = "time_tissue", y = "abs.Abundance", color = "group"))  + 
  geom_jitter(width = 0.1, height = 0, shape = 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        text = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 10, angle = 90, , hjust = 1, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 10)) +
  facet_wrap(~Species, scales="free_y", ncol = 5)
ggsave("figures_nq/absolute_abundance_feces.png",
       width = 21, height = 15, unit = "cm")

ggplot(mdf.sum2, aes_string(x = "time_tissue", y = "rel.Abundance", color = "group"))  + 
  geom_jitter(width = 0.1, height = 0, shape = 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        text = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 10, angle = 90, , hjust = 1, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 10)) +
  facet_wrap(~Species, scales="free_y", ncol = 5)
ggsave("figures_nq/relative_abundance_feces.png",
       width = 21, height = 15, unit = "cm")





# PFNA over time feces & colon

mdf.sum3 <- mdf.sum %>% 
  filter(group %in% c("HC + PFNA", "LC + PFNA"), 
         time_tissue != "3_Small intestine") %>%
  separate(time_tissue, c("time.point.d", "tissue"), sep = "_", remove = FALSE)
mdf.sum3$time.point.d <- as.numeric(mdf.sum3$time.point.d)

ggplot(mdf.sum3, aes_string(x = "time.point.d", y = "abs.Abundance"))  + 
  geom_line(aes(color = mouse)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        text = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 10, angle = 90, , hjust = 1, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 10)) +
  facet_wrap(~Species, scales="free_y", ncol = 5)
ggsave("figures_nq/absolute_abundance_feces_over_time.png",
       width = 25, height = 15, unit = "cm")

ggplot(mdf.sum3, aes_string(x = "time.point.d", y = "rel.Abundance"))  + 
  geom_line(aes(color = mouse)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        text = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 10, angle = 90, , hjust = 1, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 10)) +
  facet_wrap(~Species, scales="free_y", ncol = 5)
ggsave("figures_nq/relative_abundance_feces_over_time.png",
       width = 25, height = 15, unit = "cm")



# Check which strains colonised how many mice

mdf.sum3.count <- mdf.sum3 %>%
  filter(abs.Abundance != 0, Species != "not assigned") %>%
  group_by(time_tissue, time.point.d, tissue) %>%
  count(Species)

# C. carofaciens only found in 7-9 mice (low abundance)




# plot relative abundance without spike-in
lchc <- mdf.sum %>% 
  filter(group %in% c("HC + PFNA", "LC + PFNA")) %>%
  group_by(group, time_tissue, Species) %>%
  summarize(rel.Abundance = sum(rel.Abundance) / 9)

ggplot(lchc, aes_string(x = "time_tissue", y = "rel.Abundance", fill="Species"))  + 
  geom_bar(stat = "identity", position = "stack") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        text = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 10, angle = 90, , hjust = 1, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 10)) +
  facet_wrap(~group, scales="free")
ggsave("figures_nq/relative_abundance.pdf",
       width = 20, height = 15, unit = "cm")

write.table(lchc, file = "figures_nq/relative_abundance.txt", 
            sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE)







#######################################################
### Quantitative analysis of all samples  #########
#######################################################


path <- './data_q/'
list.files(path)

# Forward and reverse fastq filenames 
fnFs <- sort(list.files(path, pattern="_1_sequence.txt.gz",
                        full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2_sequence.txt.gz",
                        full.names = TRUE))

# Extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 5)

## Inspect read quality profiles
# Forward reads
plotQualityProfile(fnFs[1]) # over 30 is good quality, green line is mean quality score
# Reverse reads
plotQualityProfile(fnRs[3])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

## Change filter parameters for different primer set, default 16S from EMBL
## would be 515F - 806R, therefore the 'truncLen = c(240, 160)' is recommended
## 
## Here the test dataset has a peak around ~410 bp (unknown primer set), took 
## ~20s to run with 3 samples on Mac with Apple M1 Pro chip, 16Gb
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(230,210),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
# On Windows set multithread=FALSE
head(out)

write.table(out, "output_q/M4_16S_out.txt", quote = F)


# Learn the error rate, took ~90s to run for each 
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)


## Sample inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]


## merge paired reads 
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers[[1]])
write.table(mergers, "output_q/M4_16S_mergers.txt", sep = "\t", dec = ".", na = "NA",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

## construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
saveRDS(seqtab, "output_q/M4_16S_seqtab.rds")
# should be 3 2025 using 240, 210 as trunc parameters

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))



# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus",
                                    multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
saveRDS(seqtab.nochim, "output_q/M4_16S_seqtab.nochim.rds")

sum(seqtab.nochim)/sum(seqtab)



# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN),
               sapply(dadaRs, getN),
               sapply(mergers, getN),
               rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls:
# e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR",
                     "merged", "nonchim")
rownames(track) <- sample.names
head(track) # numbers in output are numbers of reads

write.table(track, "output_q/M4_16S_track.txt", quote = F)



## Hand off to phyloseq
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")

theme_set(theme_bw())







# Assign taxonomy with custom database, ~207s to run, ~3.5 min
seqtab.nochim <- readRDS("output_q/M4_16S_seqtab.nochim.rds")

taxa <- assignTaxonomy(seqtab.nochim,
                       "./community_data/Low_High_afabarum_spikein_sonja.txt",
                       multithread=TRUE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print) 
length(taxa.print)

write.table(taxa, "output_q/M4_16S_taxa.txt", sep = "\t")




# read in metadata:
metadata <- read.table("m4_16s_metadata.txt", sep = "\t", header = TRUE, check.names = FALSE)
rownames(metadata) <- metadata$"16S_ID"
metadata$time.point.d <- as.character(metadata$time.point.d)


# create phyloseq object linking tables
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(metadata), 
               tax_table(taxa))



#ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
saveRDS(ps, "output_q/M4_16S_ps.rds")

library(tidyverse)
library(dplyr)

## Plot alpha-diversity
plot_richness(ps, x="mouse",
              measures=c("Shannon", "Observed"), color="group", shape = "tissue")
ggsave("figures_q/alpha_diversity.png",
       width = 35, height = 10, unit = "cm")




# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

plot_ordination(ps.prop, ord.nmds.bray, color="group", shape = "tissue", title="Bray NMDS")
ggsave("figures_q/bray_curtis_distance.png",
       width = 15, height = 10, unit = "cm")




## Taxonomic profile
top30 <- names(sort(taxa_sums(ps), decreasing=TRUE)) #[1:30]
ps.top30 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top30 <- prune_taxa(top30, ps.top30)

my_plot_bar = function (physeq, x = "Sample", y = "Abundance", fill = NULL, title = NULL, 
                        facet_grid = NULL) {
  mdf = psmelt(physeq)
  p = ggplot(mdf, aes_string(x = x, y = y, fill = fill))
  p = p + geom_bar(stat = "identity", position = "stack")
  p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  if (!is.null(facet_grid)) {
    p <- p + facet_grid(facet_grid)
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}


my_plot_bar(ps.top30, x="time_tissue", fill="Species")  + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.title.x = element_blank()) +
  facet_wrap(~group, scales="free")
ggsave("figures_q/composition_with_spikein/abundance_by_treatment_species.png",
       width = 15, height = 15, unit = "cm")

my_plot_bar(ps.top30, x="time_tissue", fill="Family")  + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.title.x = element_blank()) +
  facet_wrap(~group, scales="free")
ggsave("figures_q/composition_with_spikein/abundance_by_treatment_family.png",
       width = 15, height = 15, unit = "cm")

my_plot_bar(ps.top30, x="time_tissue", fill="Phylum")  + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.title.x = element_blank()) +
  facet_wrap(~group, scales="free")
ggsave("figures_q/composition_with_spikein/abundance_by_treatment_phylum.png",
       width = 15, height = 15, unit = "cm")





### Analysis of PFNA samples

data <- readRDS("output_q/M4_16S_ps.rds")

mdf = psmelt(data)
mdf[is.na(mdf)] <- "not assigned"


# Check A. fabarum spike in

total_count_af <- mdf %>%
  group_by(mouse, time_tissue) %>%
  summarize(total_count = sum(Abundance))
mdf.af <- merge(mdf, total_count_af, by = c("mouse", "time_tissue"))
mdf.af <- mdf.af %>%
  mutate(Count = Abundance,
         Abundance = Count / total_count)
mdf.af <- mdf.af %>%
  group_by(mouse, group, time_tissue, Kingdom, Phylum, Class, Order, Family, Genus, Species, total_count) %>%
  summarise(abs.Abundance = sum(Count),
            rel.Abundance = sum(Abundance))

mdf.af$Species <- as.factor(mdf.af$Species)
levels(mdf.af$Species)
mdf.af$Species <- factor(mdf.af$Species, 
                         levels = c("Akkermansia muciniphila", "Collinsella aerofaciens", "Enterocloster bolteae", "Escherichia coli", "Ruminococcus gnavus",
                                    "Bacteroides fragilis", "Bacteroides thetaiotaomicron", "Bacteroides uniformis", "Lacrimispora saccharolytica", "Phocaeicola vulgatus", "not assigned", "Acetobacter fabarum"))


ggplot(subset(mdf.af, Species == "Acetobacter fabarum"), 
       aes_string(x = "group", y = "abs.Abundance", color="group"))  + 
  geom_jitter(width = 0.1, shape = 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        text = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 10, angle = 90, , hjust = 1, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 10)) +
  facet_wrap(~time_tissue, scales="free_x", nrow = 1)
ggsave("figures_q/composition_with_spikein/afabarum_absolute_abundance.png",
       width = 20, height = 7, unit = "cm")

ggplot(subset(mdf.af, Species == "Acetobacter fabarum"), 
       aes_string(x = "group", y = "rel.Abundance", color="group"))  + 
  geom_jitter(width = 0.1, shape = 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        text = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 10, angle = 90, , hjust = 1, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 10)) +
  facet_wrap(~time_tissue, scales="free_x", nrow = 1)
ggsave("figures_q/composition_with_spikein/afabarum_replative_abundance.png",
       width = 20, height = 7, unit = "cm")

af <- mdf.af %>%
  filter(Species == "Acetobacter fabarum",
         group %in% c("HC + PFNA", "LC + PFNA"))
ggplot(af, aes_string(x = "group", y = "rel.Abundance"))  + 
  geom_jitter(width = 0.1, height = 0, shape = 1, size = 3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        text = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 10, angle = 90, , hjust = 1, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 10))
ggsave("figures_q/composition_with_spikein/afabarum_replative_abundance_1.png",
       width = 7, height = 10, unit = "cm")
write.table(af, file = "figures_q/composition_with_spikein/afabarum_replative_abundance_1.txt", 
            sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE)



lc <- mdf.af %>% filter(group == "LC + PFNA", Species == "Acetobacter fabarum")
hc <- mdf.af %>% filter(group == "HC + PFNA", Species == "Acetobacter fabarum")
t.test(lc$abs.Abundance, hc$abs.Abundance, alternative = "two.sided") # p = 1.749e-08
t.test(lc$rel.Abundance, hc$rel.Abundance, alternative = "two.sided") # p = 3.929e-15
median(hc$rel.Abundance) / median(lc$rel.Abundance) # 2-fold change

# HC has higher absolute and relative amount of A. fabarum spike in, suggesting lower bacterial biomass in HC colonised mice.
# Unfortunately not all samples have detectable levels of A. fabarum, making quantitative analysis not possible.



# Calculate relative abundances without the A. fabarum spike in

mdf <- mdf %>% filter(Species != "Acetobacter fabarum")

total_count <- mdf %>%
  group_by(mouse, time_tissue) %>%
  summarize(total_count = sum(Abundance))

mdf <- merge(mdf, total_count, by = c("mouse", "time_tissue"))

mdf <- mdf %>%
  mutate(Count = Abundance,
         Abundance = Count / total_count)
mdf.sum <- mdf %>%
  group_by(mouse, group, time_tissue, Kingdom, Phylum, Class, Order, Family, Genus, Species, total_count) %>%
  summarise(abs.Abundance = sum(Count),
            rel.Abundance = sum(Abundance))

mdf.sum$Species <- as.factor(mdf.sum$Species)
levels(mdf.sum$Species)
mdf.sum$Species <- factor(mdf.sum$Species, 
                          levels = c("Akkermansia muciniphila", "Collinsella aerofaciens", "Enterocloster bolteae", "Escherichia coli", "Ruminococcus gnavus",
                                     "Bacteroides fragilis", "Bacteroides thetaiotaomicron", "Bacteroides uniformis", "Lacrimispora saccharolytica", "Phocaeicola vulgatus", "not assigned"))


# Check germfree mice and blank sample

blank <- mdf.sum %>% filter(group == "Blank")
# Blank contains 166 reads of A. muciniphila

gf <- mdf.sum %>% filter(group == "GF untreated")
ggplot(gf, aes_string(x = "time_tissue", y = "abs.Abundance", color = "mouse"))  + 
  geom_jitter(width = 0.1, height = 0, shape = 1, size = 3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        text = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 10, angle = 90, , hjust = 1, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 10)) +
  facet_wrap(~Species, scales="free_y", ncol = 5)
ggsave("figures_q/composition_without_spikein/absolute_abundance_gf.png",
       width = 25, height = 15, unit = "cm")

ggplot(gf, aes_string(x = "time_tissue", y = "rel.Abundance", color="mouse"))  + 
  geom_jitter(width = 0.1, height = 0, shape = 1, size = 3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        text = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 10, angle = 90, , hjust = 1, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 10)) +
  facet_wrap(~Species, scales="free_y", ncol = 5)
ggsave("figures_q/composition_without_spikein/relative_abundance_gf.png",
       width = 25, height = 15, unit = "cm")

# All GF samples contain A. muciniphila




# LC and HC small intestine

mdf.sum1 <- mdf.sum %>% 
  filter(group %in% c("HC + PFNA", "LC + PFNA"), 
         time_tissue == "3_Small intestine")

ggplot(mdf.sum1, aes_string(x = "group", y = "abs.Abundance"))  + 
  geom_jitter(width = 0.1, height = 0, shape = 1, size = 3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        text = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 10, angle = 90, , hjust = 1, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 10)) +
  facet_wrap(~Species, scales="free_y", ncol = 5)
ggsave("figures_q/composition_without_spikein/absolute_abundance_si.png",
       width = 25, height = 15, unit = "cm")

ggplot(mdf.sum1, aes_string(x = "group", y = "rel.Abundance"))  + 
  geom_jitter(width = 0.1, height = 0, shape = 1, size = 3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        text = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 10, angle = 90, , hjust = 1, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 10)) +
  facet_wrap(~Species, scales="free_y", ncol = 5)
ggsave("figures_q/composition_without_spikein/relative_abundance_si.png",
       width = 25, height = 15, unit = "cm")




# Feces & colon

mdf.sum2 <- mdf.sum %>% 
  filter(group %in% c("HC + PFNA", "LC + PFNA"), 
         time_tissue != "3_Small intestine")

ggplot(mdf.sum2, aes_string(x = "time_tissue", y = "abs.Abundance", color = "group"))  + 
  geom_jitter(width = 0.1, height = 0, shape = 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        text = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 10, angle = 90, , hjust = 1, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 10)) +
  facet_wrap(~Species, scales="free_y", ncol = 5)
ggsave("figures_q/composition_without_spikein/absolute_abundance_feces.png",
       width = 21, height = 15, unit = "cm")

ggplot(mdf.sum2, aes_string(x = "time_tissue", y = "rel.Abundance", color = "group"))  + 
  geom_jitter(width = 0.1, height = 0, shape = 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        text = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 10, angle = 90, , hjust = 1, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 10)) +
  facet_wrap(~Species, scales="free_y", ncol = 5)
ggsave("figures_q/composition_without_spikein/relative_abundance_feces.png",
       width = 21, height = 15, unit = "cm")





# PFNA over time feces & colon

mdf.sum3 <- mdf.sum %>% 
  filter(group %in% c("HC + PFNA", "LC + PFNA"), 
         time_tissue != "3_Small intestine") %>%
  separate(time_tissue, c("time.point.d", "tissue"), sep = "_", remove = FALSE)
mdf.sum3$time.point.d <- as.numeric(mdf.sum3$time.point.d)

ggplot(mdf.sum3, aes_string(x = "time.point.d", y = "abs.Abundance"))  + 
  geom_line(aes(color = mouse)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        text = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 10, angle = 90, , hjust = 1, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 10)) +
  facet_wrap(~Species, scales="free_y", ncol = 5)
ggsave("figures_q/composition_without_spikein/absolute_abundance_feces_over_time.png",
       width = 25, height = 15, unit = "cm")

ggplot(mdf.sum3, aes_string(x = "time.point.d", y = "rel.Abundance"))  + 
  geom_line(aes(color = mouse)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        text = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 10, angle = 90, , hjust = 1, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 10)) +
  facet_wrap(~Species, scales="free_y", ncol = 5)
ggsave("figures_q/composition_without_spikein/relative_abundance_feces_over_time.png",
       width = 25, height = 15, unit = "cm")





# plot relative abundance without spike-in
lchc <- mdf.sum %>% 
  filter(group %in% c("HC + PFNA", "LC + PFNA")) %>%
  group_by(group, time_tissue, Species) %>%
  summarize(rel.Abundance = sum(rel.Abundance) / 9)

ggplot(lchc, aes_string(x = "time_tissue", y = "rel.Abundance", fill="Species"))  + 
  geom_bar(stat = "identity", position = "stack") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        text = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 10, angle = 90, , hjust = 1, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 10)) +
  facet_wrap(~group, scales="free")
ggsave("figures_q/composition_without_spikein/relative_abundance.png",
       width = 20, height = 15, unit = "cm")


lchc <- mdf.af %>% 
  filter(group %in% c("HC + PFNA", "LC + PFNA")) %>%
  group_by(group, time_tissue, Species) %>%
  summarize(rel.Abundance = sum(rel.Abundance) / 9)

ggplot(lchc, aes_string(x = "time_tissue", y = "rel.Abundance", fill="Species"))  + 
  geom_bar(stat = "identity", position = "stack") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        text = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 10, angle = 90, , hjust = 1, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 10)) +
  facet_wrap(~group, scales="free")
ggsave("figures_q/composition_with_spikein/relative_abundance.png",
       width = 20, height = 15, unit = "cm")






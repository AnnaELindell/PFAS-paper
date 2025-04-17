library(dada2)
path <- './data/'
list.files(path)

# Forward and reverse fastq filenames 
fnFs <- sort(list.files(path, pattern="_1_sequence.txt.gz",
                        full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2_sequence.txt.gz",
                        full.names = TRUE))

# Extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 6)

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

write.table(out, "output/20230912_16S_out.txt", quote = F)


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
write.table(mergers, "output/20230912_16S_mergers.txt", sep = "\t", dec = ".", na = "NA",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

## construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
saveRDS(seqtab, "output/20230912_16S_seqtab.rds")
# should be 3 2025 using 240, 210 as trunc parameters

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus",
                                    multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
saveRDS(seqtab.nochim, "output/20230912_16S_seqtab.nochim.rds")

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

write.table(track, "output/20230912_16S_track.txt", quote = F)

track <- read.table("output/20230912_16S_track.txt", sep = " ", header = TRUE, check.names = FALSE)
track$'16S_ID' <- rownames(track)
metadata <- read.table("samples-metadata.txt", sep = "\t", header = TRUE, check.names = FALSE)
track <- merge(track, metadata, by = "16S_ID")
track <- track %>% filter(treatment %in% c("PFNA", "PFNA GF"))
sum(track$input)
mean(track$input)


## Hand off to phyloseq
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")

theme_set(theme_bw())



# Assign taxonomy, ~207s to run, ~3.5 min
taxa <- assignTaxonomy(seqtab.nochim,
                       "./tax/silva_nr99_v138.1_wSpecies_train_set.fa.gz",
                       multithread=TRUE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print) 
length(taxa.print)



# Assign taxonomy with custom database, ~207s to run, ~3.5 min
taxa <- assignTaxonomy(seqtab.nochim,
                       "./tax/Com20_Ecoli-spike-in.txt",
                       multithread=TRUE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print) 
length(taxa.print)

write.table(taxa, "output/20230912_16S_taxa.txt", sep = "\t")




# read in metadata:
metadata <- read.table("samples-metadata.txt", sep = "\t", header = TRUE, check.names = FALSE)
rownames(metadata) <- metadata$"16S_ID"
metadata$time_point <- as.character(metadata$time_point)


# create phyloseq object linking tables for all compounds
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(metadata), 
               tax_table(taxa))

# create phyloseq object linking tables for PFNA only
metadata.pfna <- metadata %>%
  filter(treatment %in% c("DMSO", "PFNA", "PFNA GF"))
samples.pfna <- metadata.pfna$'16S_ID'
seqtab.nochim.pfna <- seqtab.nochim[samples.pfna,]
ps.pfna <- phyloseq(otu_table(seqtab.nochim.pfna, taxa_are_rows=FALSE), 
               sample_data(metadata.pfna), 
               tax_table(taxa))



#ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
saveRDS(ps, "output/20230912_16S_ps.rds")

dna.pfna <- Biostrings::DNAStringSet(taxa_names(ps.pfna))
names(dna.pfna) <- taxa_names(ps.pfna)
ps.pfna <- merge_phyloseq(ps.pfna, dna.pfna)
taxa_names(ps.pfna) <- paste0("ASV", seq(ntaxa(ps.pfna)))
ps.pfna
saveRDS(ps.pfna, "output/20230912_16S_ps_pfna.rds")

library(tidyverse)

## Plot alpha-diversity for all compounds
plot_richness(ps, x="mouse_number",
              measures=c("Shannon", "Observed"), color="treatment", shape = "specimen")
ggsave("figures/alpha_diversity_all.png",
       width = 35, height = 10, unit = "cm")

## Plot alpha-diversity only for DMSO and PFNA treated groups
plot_richness(ps.pfna, x="mouse_number",
              measures=c("Shannon", "Observed"), color="treatment", shape = "specimen")
ggsave("figures/alpha_diversity_pfna.png",
       width = 25, height = 10, unit = "cm")


# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

plot_ordination(ps.prop, ord.nmds.bray, color="treatment", shape = "specimen", title="Bray NMDS")
ggsave("figures/bray_curtis_distance_all.png",
       width = 15, height = 10, unit = "cm")

ps.prop.pfna <- transform_sample_counts(ps.pfna, function(otu) otu/sum(otu))
ord.nmds.bray.pfna <- ordinate(ps.prop.pfna, method="NMDS", distance="bray")

plot_ordination(ps.prop.pfna, ord.nmds.bray.pfna, color="treatment", shape = "specimen", title="Bray NMDS")
ggsave("figures/bray_curtis_distance_pfna.png",
       width = 15, height = 10, unit = "cm")





## Taxonomic profile for all compounds
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


my_plot_bar(ps.top30, x="time_specimen", fill="Species")  + 
  facet_wrap(~treatment, scales="free")
ggsave("figures/composition_with_spikein_all/species_abundance_by_treatment.png",
       width = 25, height = 15, unit = "cm")

my_plot_bar(ps.top30, x="time_specimen", fill="Family")  + 
  facet_wrap(~treatment, scales="free")
ggsave("figures/composition_with_spikein_all/family_abundance_by_treatment.png",
       width = 25, height = 15, unit = "cm")

my_plot_bar(ps.top30, x="time_specimen", fill="Phylum")  + 
  facet_wrap(~treatment, scales="free")
ggsave("figures/composition_with_spikein_all/phylum_abundance_by_treatment.png",
       width = 25, height = 15, unit = "cm")



## Taxonomic profile for pfna
top30.pfna <- names(sort(taxa_sums(ps.pfna), decreasing=TRUE)) #[1:30]
ps.pfna.top30 <- transform_sample_counts(ps.pfna, function(OTU) OTU/sum(OTU))
ps.pfna.top30 <- prune_taxa(top30.pfna, ps.pfna.top30)

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


my_plot_bar(ps.pfna.top30, x="time_specimen", fill="Species")  + 
  facet_wrap(~treatment, scales="free")
ggsave("figures/composition_with_spikein_pfna/species_abundance_by_treatment.png",
       width = 25, height = 15, unit = "cm")

my_plot_bar(ps.pfna.top30, x="time_specimen", fill="Family")  + 
  facet_wrap(~treatment, scales="free")
ggsave("figures/composition_with_spikein_pfna/family_abundance_by_treatment.png",
       width = 25, height = 15, unit = "cm")

my_plot_bar(ps.pfna.top30, x="time_specimen", fill="Phylum")  + 
  facet_wrap(~treatment, scales="free")
ggsave("figures/composition_with_spikein_pfna/phylum_abundance_by_treatment.png",
       width = 25, height = 15, unit = "cm")





library(tidyverse)



### Analysis of PFNA samples

data <- readRDS("output/20230912_16S_ps_pfna.rds")

mdf = psmelt(data)
mdf[is.na(mdf)] <- "not assigned"


# Check E. coli spike in

total_count_ecoli <- mdf %>%
  group_by(mouse_number, time_specimen) %>%
  summarize(total_count = sum(Abundance))
mdf.ecoli <- merge(mdf, total_count_ecoli, by = c("mouse_number", "time_specimen"))
mdf.ecoli <- mdf.ecoli %>%
  mutate(Count = Abundance,
         Abundance = Count / total_count)
mdf.ecoli.sum <- mdf.ecoli %>%
  group_by(mouse_number, treatment, time_specimen, Kingdom, Phylum, Class, Order, Family, Genus, Species, total_count) %>%
  summarise(abs.Abundance = sum(Count),
            rel.Abundance = sum(Abundance))

ggplot(subset(mdf.ecoli.sum, Species == "Escherichia coli"), 
       aes_string(x = "time_specimen", y = "abs.Abundance", color="treatment"))  + 
  geom_jitter(width = 0.1, shape = 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        text = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 10, angle = 90, , hjust = 1, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 10)) +
  facet_wrap(~Species, scales="free_y")
ggsave("figures/composition_with_spikein_pfna/ecoli_absolute_abundance1_pfna_dmso_gf.png",
       width = 10, height = 10, unit = "cm")
ggplot(subset(mdf.ecoli.sum, Species == "Escherichia coli"), 
       aes_string(x = "time_specimen", y = "rel.Abundance", color="treatment"))  + 
  geom_jitter(width = 0.1, shape = 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        text = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 10, angle = 90, , hjust = 1, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 10)) +
  facet_wrap(~Species, scales="free")
ggsave("figures/composition_with_spikein_pfna/ecoli_relative_abundance1_pfna_dmso_gf.png",
       width = 10, height = 10, unit = "cm")

ggplot(subset(mdf.ecoli.sum, Species == "Escherichia coli" & treatment != "PFNA GF"), 
       aes_string(x = "time_specimen", y = "abs.Abundance", color="treatment"))  + 
  geom_jitter(width = 0.1, shape = 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        text = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 10, angle = 90, , hjust = 1, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 10)) +
  facet_wrap(~Species, scales="free_y")
ggsave("figures/composition_with_spikein_pfna/ecoli_absolute_abundance1_pfna_dmso.png",
       width = 10, height = 10, unit = "cm")
ggplot(subset(mdf.ecoli.sum, Species == "Escherichia coli" & treatment != "PFNA GF"), 
       aes_string(x = "time_specimen", y = "rel.Abundance", color="treatment"))  + 
  geom_jitter(width = 0.1, shape = 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        text = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 10, angle = 90, , hjust = 1, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 10)) +
  facet_wrap(~Species, scales="free")
ggsave("figures/composition_with_spikein_pfna/ecoli_relative_abundance1_pfna_dmso.png",
       width = 10, height = 10, unit = "cm")

ggplot(subset(mdf.ecoli.sum, Species == "Escherichia coli" & treatment != "PFNA GF"), 
       aes_string(x = "time_specimen", y = "abs.Abundance", color="treatment"))  + 
  geom_jitter(width = 0.1, shape = 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        text = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 10, angle = 90, , hjust = 1, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 10)) +
  facet_wrap(~Species) +
  ylim(0, 300)
ggsave("figures/composition_with_spikein_pfna/ecoli_absolute_abundance1_pfna_dmso_crop.png",
       width = 10, height = 10, unit = "cm")
ggplot(subset(mdf.ecoli.sum, Species == "Escherichia coli" & treatment != "PFNA GF"), 
       aes_string(x = "time_specimen", y = "rel.Abundance", color="treatment"))  + 
  geom_jitter(width = 0.1, shape = 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        text = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 10, angle = 90, , hjust = 1, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 10)) +
  facet_wrap(~Species) +
  ylim(0, 0.025)
ggsave("figures/composition_with_spikein_pfna/ecoli_relative_abundance1_pfna_dmso_crop.png",
       width = 10, height = 10, unit = "cm")

# no difference in E. coli spike in between DMSO and PFNA treated mice suggests same amount of bacterial biomass.




# Calculate relative abundances without the E. coli spike in

mdf <- mdf %>% filter(Species != "Escherichia coli")

total_count <- mdf %>%
  group_by(mouse_number, time_specimen) %>%
  summarize(total_count = sum(Abundance))

mdf <- merge(mdf, total_count, by = c("mouse_number", "time_specimen"))

mdf <- mdf %>%
  mutate(Count = Abundance,
         Abundance = Count / total_count)
mdf.sum <- mdf %>%
  group_by(mouse_number, treatment, time_specimen, Kingdom, Phylum, Class, Order, Family, Genus, Species, total_count) %>%
  summarise(abs.Abundance = sum(Count),
            rel.Abundance = sum(Abundance))




# Check germfree mice

ggplot(subset(mdf.sum, treatment == "PFNA GF"), aes_string(x = "time_specimen", y = "abs.Abundance", color="treatment"))  + 
  geom_jitter(width = 0.1, shape = 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        text = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 10, angle = 90, , hjust = 1, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 10)) +
  facet_wrap(~Species, scales="free_y")
ggsave("figures/composition_without_spikein_pfna/absolute_abundance_gf.png",
       width = 28, height = 15, unit = "cm")

ggplot(subset(mdf.sum, treatment == "PFNA GF"), aes_string(x = "time_specimen", y = "rel.Abundance", color="treatment"))  + 
  geom_jitter(width = 0.1, shape = 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        text = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 10, angle = 90, , hjust = 1, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 10)) +
  facet_wrap(~Species, scales="free_y")
ggsave("figures/composition_without_spikein_pfna/relative_abundance_gf.png",
       width = 28, height = 15, unit = "cm")

# one small intestine sample contaminated (probably happened after sample collection, as colon and feces samples from the same mouse look clean)



# PFNA vs. DMSO treatment small intestine

mdf.sum1 <- mdf.sum %>% filter(treatment != "PFNA GF", time_specimen == "3_Small intestine")

ggplot(mdf.sum1, aes_string(x = "treatment", y = "abs.Abundance"))  + 
  geom_jitter(width = 0.1, shape = 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        text = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 10, angle = 90, , hjust = 1, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 10)) +
  facet_wrap(~Species, scales="free_y")
ggsave("figures/composition_without_spikein_pfna/absolute_abundance_dmso_pfna_si.png",
       width = 25, height = 15, unit = "cm")

ggplot(mdf.sum1, aes_string(x = "treatment", y = "rel.Abundance"))  + 
  geom_jitter(width = 0.1, shape = 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        text = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 10, angle = 90, , hjust = 1, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 10)) +
  facet_wrap(~Species, scales="free_y")
ggsave("figures/composition_without_spikein_pfna/relative_abundance_pfna_dmso_si.png",
       width = 25, height = 15, unit = "cm")

ttest.si <- mdf.sum1 %>%
  group_by(time_specimen, Species) %>%
  summarize(p.value = t.test(rel.Abundance ~ treatment, alternative = "two.sided", paired = FALSE, var.equal =TRUE)$p.value) 
write.table(ttest.si,
            file = "results/m1_16s_si_relabundance_pvalues.txt",
            sep = "\t", dec = ".", na = "NA",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
# no significant differences between DMSO and PFNA treatemnt for SI samples



# PFNA vs. DMSO treatment feces & colon

mdf.sum2 <- mdf.sum %>% filter(treatment != "PFNA GF", time_specimen != "3_Small intestine")

ggplot(mdf.sum2, aes_string(x = "time_specimen", y = "abs.Abundance", color = "treatment"))  + 
  geom_jitter(width = 0.1, shape = 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        text = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 10, angle = 90, , hjust = 1, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 10)) +
  facet_wrap(~Species, scales="free_y")
ggsave("figures/composition_without_spikein_pfna/absolute_abundance_dmso_pfna_feces.png",
       width = 28, height = 15, unit = "cm")

ggplot(mdf.sum2, aes_string(x = "time_specimen", y = "rel.Abundance", color = "treatment"))  + 
  geom_jitter(width = 0.1, shape = 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        text = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 10, angle = 90, , hjust = 1, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 10)) +
  facet_wrap(~Species, scales="free_y")
ggsave("figures/composition_without_spikein_pfna/relative_abundance_pfna_dmso_feces.png",
       width = 28, height = 15, unit = "cm")

ttest.feces <- mdf.sum2 %>%
  group_by(time_specimen, Species) %>%
  summarize(p.value = t.test(rel.Abundance ~ treatment, alternative = "two.sided", paired = FALSE, var.equal =TRUE)$p.value) 
write.table(ttest.feces,
            file = "results/m1_16s_feces_relabundance_pvalues.txt",
            sep = "\t", dec = ".", na = "NA",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
subset(ttest.feces, p.value < 0.05)
# Only B. fragilis significantly different in Colon sample (DMSO higher than PFNA treatment)




# PFNA over time feces & colon

mdf.sum3 <- mdf.sum %>% 
  filter(treatment == "PFNA", time_specimen != "3_Small intestine") %>%
  separate(time_specimen, c("time.point.d", "tissue"), sep = "_", remove = FALSE)
mdf.sum3$time.point.d <- as.numeric(mdf.sum3$time.point.d)

ggplot(mdf.sum3, aes_string(x = "time.point.d", y = "abs.Abundance"))  + 
  geom_line(aes(color = mouse_number)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        text = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 10, angle = 90, , hjust = 1, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 10)) +
  facet_wrap(~Species, scales="free_y")
ggsave("figures/composition_without_spikein_pfna/absolute_abundance_pfna_feces.png",
       width = 28, height = 15, unit = "cm")

ggplot(mdf.sum3, aes_string(x = "time.point.d", y = "rel.Abundance"))  + 
  geom_line(aes(color = mouse_number)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        text = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 10, angle = 90, , hjust = 1, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 10)) +
  facet_wrap(~Species, scales="free_y")
ggsave("figures/composition_without_spikein_pfna/relative_abundance_pfna_feces.png",
       width = 28, height = 15, unit = "cm")



# Check which strains colonised how many mice

mdf.sum3.count <- mdf.sum3 %>%
  filter(abs.Abundance != 0, Species != "not assigned") %>%
  group_by(time_specimen, time.point.d, tissue) %>%
  count(Species)

# B. fragilis, C. perfringens, C. saccharolyticum only found in a few samples





### get abundance of accumulators and non-accumulators in PFNA treated mice

ggplot(mdf, aes_string(x = "time_specimen", y = "Abundance", fill="Species"))  + 
  geom_bar(stat = "identity", position = "stack") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        text = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 10, angle = 90, , hjust = 1, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 10)) +
  facet_wrap(~treatment, scales="free")
ggsave("figures/composition_without_spikein_pfna/relative_abundance.png",
       width = 20, height = 15, unit = "cm")

plot1.data <- mdf %>%
  group_by(treatment, time_specimen, Species) %>%
  summarize(Abundance = sum(Abundance) / 9)
ggplot(plot1.data, aes_string(x = "time_specimen", y = "Abundance", fill="Species"))  + 
  geom_bar(stat = "identity", position = "stack") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        text = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 10, angle = 90, , hjust = 1, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 10)) +
  labs(y = "Relative Abundance") +
  facet_wrap(~treatment, scales="free")
ggsave("figures/composition_without_spikein_pfna/relative_abundance_normalised.pdf",
       width = 20, height = 15, units = "cm")
write.table(plot1.data, "figures/composition_without_spikein_pfna/relative_abundance_normalised.txt", 
            sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE)



pfna_class <- read.table("tax/pfna-accumulation-classification.txt", sep = "\t", header = T)
mdf <- merge(mdf, pfna_class, by = "Species")

ggplot(mdf, aes_string(x = "time_specimen", y = "Abundance", fill="pfna_accumulation"))  + 
  geom_bar(stat = "identity", position = "stack") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        text = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 10, angle = 90, , hjust = 1, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 10)) +
  facet_wrap(~treatment, scales="free")
ggsave("figures/composition_without_spikein_pfna/relative_abundance_high_low.png",
       width = 20, height = 15, unit = "cm")
plot2.data <- mdf %>%
  group_by(treatment, time_specimen, pfna_accumulation) %>%
  summarize(Abundance = sum(Abundance) / 9)
ggplot(plot2.data, aes_string(x = "time_specimen", y = "Abundance", fill="pfna_accumulation"))  + 
  geom_bar(stat = "identity", position = "stack") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        text = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 10, angle = 90, , hjust = 1, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 10)) +
  labs(y = "Relative Abundance") +
  facet_wrap(~treatment, scales="free")
ggsave("figures/composition_without_spikein_pfna/relative_abundance_high_low_norm.pdf",
       width = 20, height = 15, unit = "cm")
write.table(plot2.data, "figures/composition_without_spikein_pfna/relative_abundance_high_low_norm.txt", 
            sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE)




# read in feces PFNA concentration
pfna_conc <- read.table("m1_feces_data_pfna.txt", sep = "\t", header = T)
pfna_conc <- pfna_conc %>%
  filter(treatment == "Com20 + PFNA") %>%
  unite(time_specimen, time.point.d, tissue) %>%
  select(mouse_number, time_specimen, conc2)

# subset mdf data
mdf_pfna <- mdf %>%
  filter(treatment == "PFNA") %>%
  select(mouse_number, time_specimen, Abundance, Species, pfna_accumulation)

# merge data
mdf_pfna <- merge(mdf_pfna, pfna_conc, by = c("mouse_number", "time_specimen"))
mdf_pfna_sum <- mdf_pfna %>%
  group_by(mouse_number, time_specimen, pfna_accumulation) %>%
  summarise(Abundance = sum(Abundance),
            conc2 = mean(conc2))
  
ggplot(subset(mdf_pfna_sum, pfna_accumulation == "High-accumulating"), aes(x = Abundance, y = conc2, colour = mouse_number)) +
  geom_point(size = 3, alpha = 0.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        text = element_text(size=10),
        axis.text.x = element_text(size=10, color = "black"),
        axis.text.y = element_text(size=10, color = "black"),
        strip.background = element_blank()) +
  facet_wrap(~time_specimen, scales = "free") +
  labs(y = "PFNA concentration [ug/g feces]", x = "Abundance [%]") +
  ylim(0,NA)

mdf_pfna_sum %>%
  filter(pfna_accumulation == "High-accumulating",
         time_specimen %in% c("1_Feces", "2_Feces", "3_Colon")) %>%
  ggplot(aes(x = Abundance, y = conc2, colour = mouse_number)) +
  geom_point(size = 3, alpha = 0.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        text = element_text(size=10),
        axis.text.x = element_text(size=10, color = "black"),
        axis.text.y = element_text(size=10, color = "black"),
        strip.background = element_blank())




################## Plot the phyla profile of top 200 ASVs

################## Repeat the practice with GTDB taxonomic assignment

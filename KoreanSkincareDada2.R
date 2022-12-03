# Install bioconductor controller, BiocManager
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#Install dada2 from bioconductor
#BiocManager::install(version = "3.15")
#BiocManager::install("dada2", version = "3.15")

save.image(file = "120322.Rdata")
#load("120322.Rdata")


# Load dada2 and check the version
library("dada2"); packageVersion("dada2")

# UPDATED FROM TUTORIAL
# Make a variable that lists all the files in the KoreanSkincareTrimmedReads_g folder and then list them.
path <- "/cloud/project/KoreanSkincareTrimmedReads_g" # using option g right now
list.files(path)

# UPDATED FROM TUTORIAL
# Forward and reverse fastq filenames end in _1.fastq.gz and _2.fastq.gz
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))
# UPDATED FROM TUTORIAL
# Just changed 1 to 2 to get the info after the first underscore
# Extract sample names, assuming filenames start with Trimmed_SAMPLENAME
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 2)

# Plot quality profile of fastq files
plotQualityProfile(fnFs[1])
plotQualityProfile(fnRs[1])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

# Assign sample names to the column names of the filtered data
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Trim reads
# Updated truncLen based on how our reads look post trimming
# Trying to keep as many reads as possible while cutting off the very bad ones
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(275,230),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)

# These are test files to check error file and read file to get error in each read
errF <- learnErrors(filtFs, multithread=FALSE)
errR <- learnErrors(filtRs, multithread=FALSE)

# Plot the errors
plotErrors(errF, nominalQ=TRUE)

# Actually run dada
dadaFs <- dada(filtFs, err=errF, multithread=FALSE)
dadaRs <- dada(filtRs, err=errR, multithread=FALSE)
dadaFs[[1]]

# Perform quality control and trimming of the reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Change the dataframe to a sequence table and print dimensions
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# Remove bimers and print dimensions
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

# Return number of sequences in file
getN <- function(x) sum(getUniques(x))

# Make table of number of reads before and after processing
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

# Name columns
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)


# Now begin Part 2 (Assigning taxonomy)

# READS MIGHT NOT BE APPROPRIATELY ASSIGNED: WILL DISCUSS DURING MEETING
# Need to upload the silva_nr_v132_train_set.fa.gz file into the KoreanSkincareTrimmedReads_g folder
# Assign taxonomy using SILVA database
taxa <- assignTaxonomy(seqtab.nochim, "/cloud/project/KoreanSkincareTrimmedReads_g/silva_nr_v132_train_set.fa.gz", multithread=FALSE)

# Optional addition
# Need to upload silva_species_assignment_v132.fa.gz file into the KoreanSkincareTrimmedReads_g folder
taxa <- addSpecies(taxa, "/cloud/project/KoreanSkincareTrimmedReads_g/silva_species_assignment_v132.fa.gz")

# Create new datafrom from old one
taxa.print <- taxa 
# Remove row names
rownames(taxa.print) <- NULL
# Display new dataframe to print
head(taxa.print)

# Just to see what we're working with
rownames(seqtab.nochim)

# I BELIEVE THIS CODE IS ONLY FOR THE TUTORIAL
# WE DO NOT HAVE A MOCK COMMUNITY TO EVALUATE DADA2'S ACCURACY AGAINST
# Subset to mock samples
# unqs.mock <- seqtab.nochim["Mock",]
# unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock, sort decreasing
# Print number of exact matches
# cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
# Load mock reference
# mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
# Count number of unique mock species in ref
# match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
# Print number of exact matches
# cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")

# THE CODE BELOW THIS IS PRETTY SPECIFIC TO THE TUTORIAL
# OUR FILES ARE NAMED VERY DIFFERENTLY, SO WE'D NEED TO GENERATE PLOTS IN A DIFFERENT WAY
# Load packages
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")

# Setting black/white color theme
theme_set(theme_bw())

# Code that extracts subject, gender, day and time (early/late) from sample names
metadata <- read.table(file = 'Metadata.txt', sep = ',', header = TRUE)
View(metadata)


id <- metadata$Run
week <- substr(metadata$Submitter_Id, 4,5)
samdf <- data.frame(ID = id, Week = week)


# Setting phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                sample_data(samdf), 
                tax_table(taxa))
# ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample

# Make a DNAStringSet object from taxa names in phyloseq
# dna <- Biostrings::DNAStringSet(taxa_names(ps))
# Change the names of the DNA sequences to match the names of the species
# names(dna) <- taxa_names(ps)
#merge phyloseq object with dna sequence data
# ps <- merge_phyloseq(ps, dna)
# assign names to the taxa of a phyloseq object
# taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
# ps

# Graph features
# plot_richness(ps, x="Day", measures=c("Shannon", "Simpson"), color="When")

# Transform data to proportions as appropriate for Bray-Curtis distances
# ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
# ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

# Set more graph features
# plot_ordination(ps.prop, ord.nmds.bray, color="When", title="Bray NMDS")

# This code sorts the OTUs by the decreasing value of the taxa_sums function. This function will return the sum of the OTUs for each taxa. 
# The names of the OTUs are then displayed using the first 20 numbers of the list of OTUs.
# top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
# The next step is to transform the sample counts of the OTUs. The transform_sample_counts function will take the OTUs and divide it by the sum of all of the OTUs.
# The next step is to prune the taxa. The prune_taxa function will only display the top 20 OTUs. 
# ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
# ps.top20 <- prune_taxa(top20, ps.top20)
# The final code will plot the bar chart. In this chart, the x axis is the day, the fill is the family, and the facet will wrap the when part of the data. The scales are set to be "free" on the x axis.
# plot_bar(ps.top20, x="Day", fill="Family") + facet_wrap(~When, scales="free_x")

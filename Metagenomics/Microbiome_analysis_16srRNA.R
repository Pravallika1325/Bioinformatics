library('dada2')
packageVersion('dada2')
library('phyloseq')
path <- "F:\R dada2\new r dada2/MiSeq_SOP"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

fnFs
fnRs

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

sample.names


#plot qulatiy
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])



# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered3", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered3", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
filtFs
filtRs


# trimming the reads:

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)
tail(out)


# learning the error rates:

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ = TRUE)


# sample inference:

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]


# merge sequence pairs:
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)


# Inspect the merger data.frame from the first sample
head(mergers[[1]])


# construct sequence table:
seqtab <- makeSequenceTable(mergers)
dim(seqtab)


# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))



# removing chimeric sequences:
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)



# track reads through pipeline:

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))


# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)





# Save the file
write.csv(track,  "readinfor_table.csv")


#Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "E:/R dada2/silva_nr_v138_train_set.fa.gz", multithread=FALSE,tryRC=TRUE)

taxa
#Assign species
taxa <- addSpecies(taxa, "C:/Users/Amrita/Downloads/16S_Metagenomics-20230213T143924Z-001/16S_Metagenomics/silva_species_assignment_v138.fa.gz",tryRC=TRUE)


#Inspect the taxonomic assignments
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
               tax_table(taxa))

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps


# preparing the ASV table (merging the taxonomic table and abundance table)
df1 <- t(otu_table(ps))
df2 <- tax_table(ps)
df3 <- cbind(df2, df1)


# Save the file
write.csv(df3, "file_ASV.csv")
View(df2)
View(dna)
View(track)










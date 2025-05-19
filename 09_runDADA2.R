# DADA2 -------------------------------------------------------------------

rm(list = ls())

# Setup -------------------------------------------------------------------

getwd()
setwd("~/sme_rhizobiome/")
getwd()

path <- paste0(getwd(), "/analysis/dada2/")
path

library("tidyverse")    
library("dada2")        
library("DECIPHER")
library("phangorn")
library("phyloseq")
library("ggtree")    

# Functions ---------------------------------------------------------------

tax.dfying <- function(tax.matrix) {
  
  if (colnames(tax.matrix) %>% str_count(., "Species") %>% sum() > 1) {
    
    tax.matrix %>% 
      data.frame() %>% 
      mutate(across(everything(), ~ str_remove_all(., "^[kdpcofgs]__"))) %>% 
      mutate(across(everything(), ~ str_remove_all(., "_[:digit:]+$"))) %>%
      mutate(across(everything(), ~ str_remove_all(., "_[:LETTER:]+$"))) %>%
      mutate(seq = rownames(.)) %>%
      select(seq, Kingdom:Genus, Species = Species.1) %>% 
      mutate(Species = paste0(Genus, "_", Species)) %>%
      mutate(Species = str_remove_all(Species, ".*NA$")) %>%
      mutate(Species = str_remove_all(Species, "\\(.*$")) %>%
      mutate(across(Kingdom:Species, ~ gsub(" ", "_", .))) %>%
      mutate(across(Kingdom:Species, ~ gsub("^$", NA, .))) %>%
      mutate(across(Kingdom:Species, ~ gsub("^_[:letter:]", NA, .))) -> tax.df
    
  } else {
    
    tax.matrix %>% 
      data.frame() %>% 
      mutate(across(everything(), ~ str_remove_all(., "^[kdpcofgs]__"))) %>% 
      mutate(across(everything(), ~ str_remove_all(., "_[:digit:]+$"))) %>%
      mutate(across(everything(), ~ str_remove_all(., "_[:LETTER:]+$"))) %>%
      mutate(seq = rownames(.)) %>%
      select(seq, Kingdom:Species) %>% 
      mutate(Species = paste0(Genus, "_", Species)) %>%
      mutate(Species = str_remove_all(Species, ".*NA$")) %>%
      mutate(Species = str_remove_all(Species, "\\(.*$")) %>%
      mutate(across(Kingdom:Species, ~ gsub(" ", "_", .))) %>%
      mutate(across(Kingdom:Species, ~ gsub("^$", NA, .))) %>%
      mutate(across(Kingdom:Species, ~ gsub("^_[:letter:]", NA, .))) -> tax.df
    
  }
  
  return(tax.df)
  
}

summary.tax <- function(ordered_names, ...){
  
  tmp_names <- ordered_names
  tmp_list <- list(...)
  
  tmp_summary1 <- data.frame(check.names = F, 
                             Database = NA, 
                             Kingdom = NA, 
                             Phylum = NA, 
                             Class = NA, 
                             Order = NA,
                             Family = NA, 
                             Genus = NA, 
                             Species = NA)
  
  for (i in seq(length(tmp_list))) {
    
    tmp_db <- tmp_list[[i]] %>% select(Kingdom:Species)
    tmp_name <- tmp_names[[i]]
    
    tmp_ksum <- sum(table(tmp_db["Kingdom"]))
    tmp_psum <- sum(table(tmp_db["Phylum"]))
    tmp_csum <- sum(table(tmp_db["Class"]))
    tmp_osum <- sum(table(tmp_db["Order"]))
    tmp_fsum <- sum(table(tmp_db["Family"]))
    tmp_gsum <- sum(table(tmp_db["Genus"]))
    tmp_ssum <- sum(table(tmp_db["Species"]))
    
    tmp_summary1 <- tmp_summary1 %>% add_row(Database = tmp_name, 
                                             Kingdom = tmp_ksum, 
                                             Phylum = tmp_psum, 
                                             Class = tmp_csum, 
                                             Order = tmp_osum,
                                             Family = tmp_fsum, 
                                             Genus = tmp_gsum, 
                                             Species = tmp_ssum)
    
  }
  
  tmp_summary1 <- tmp_summary1[-1,]
  
  tmp_summary2 <- data.frame(check.names = F, 
                             Database = NA,
                             UnqKingdom = NA, 
                             UnqPhylum = NA, 
                             UnqClass = NA, 
                             UnqOrder = NA,
                             UnqFamily = NA, 
                             UnqGenus = NA, 
                             UnqSpecies = NA)
  
  for (i in seq(length(tmp_list))) {
    
    tmp_db <- tmp_list[[i]] %>% select(Kingdom:Species)
    tmp_name <- tmp_names[[i]]
    
    tmp_ksum_u <- sum(table(unique(tmp_db["Kingdom"])))
    tmp_psum_u <- sum(table(unique(tmp_db["Phylum"])))
    tmp_csum_u <- sum(table(unique(tmp_db["Class"])))
    tmp_osum_u <- sum(table(unique(tmp_db["Order"])))
    tmp_fsum_u <- sum(table(unique(tmp_db["Family"])))
    tmp_gsum_u <- sum(table(unique(tmp_db["Genus"])))
    tmp_ssum_u <- sum(table(unique(tmp_db["Species"])))
    
    tmp_summary2 <- tmp_summary2 %>% add_row(Database = tmp_name, 
                                             UnqKingdom = tmp_ksum_u, 
                                             UnqPhylum = tmp_psum_u, 
                                             UnqClass = tmp_csum_u, 
                                             UnqOrder = tmp_osum_u, 
                                             UnqFamily = tmp_fsum_u, 
                                             UnqGenus = tmp_gsum_u, 
                                             UnqSpecies = tmp_ssum_u)
    
  }
  
  tmp_summary2 <- tmp_summary2[-1,]
  
  return(list(tmp_summary1, tmp_summary2))
  
}

summary.tax.2 <- function(count_table, ordered_names, ...){
  
  tmp_count <- count_table
  tmp_names <- ordered_names
  tmp_list <- list(...)
  
  tmp_summary1 <- data.frame(check.names = F, 
                             Database = NA, 
                             Kingdom = NA, 
                             Phylum = NA, 
                             Class = NA, 
                             Order = NA,
                             Family = NA, 
                             Genus = NA, 
                             Species = NA)
  
  for (i in seq(length(tmp_list))) {
    
    tmp_db <- merge(tmp_list[[i]], tmp_count, by = "seq") %>% mutate(RowSum = rowSums(select(., where(is.numeric)))) %>% select(Kingdom:Species, RowSum)
    tmp_name <- tmp_names[[i]]
    
    tmp_ksum <- tmp_db %>% filter(!is.na(Kingdom)) %>% pull(RowSum) %>% sum()
    tmp_psum <- tmp_db %>% filter(!is.na(Phylum)) %>% pull(RowSum) %>% sum()
    tmp_csum <- tmp_db %>% filter(!is.na(Class)) %>% pull(RowSum) %>% sum()
    tmp_osum <- tmp_db %>% filter(!is.na(Order)) %>% pull(RowSum) %>% sum()
    tmp_fsum <- tmp_db %>% filter(!is.na(Family)) %>% pull(RowSum) %>% sum()
    tmp_gsum <- tmp_db %>% filter(!is.na(Genus)) %>% pull(RowSum) %>% sum()
    tmp_ssum <- tmp_db %>% filter(!is.na(Species)) %>% pull(RowSum) %>% sum()
    
    tmp_summary1 <- tmp_summary1 %>% add_row(Database = tmp_name, 
                                             Kingdom = tmp_ksum, 
                                             Phylum = tmp_psum, 
                                             Class = tmp_csum, 
                                             Order = tmp_osum,
                                             Family = tmp_fsum, 
                                             Genus = tmp_gsum, 
                                             Species = tmp_ssum)
    
  }
  
  tmp_summary1 <- tmp_summary1[-1,]
  
  tmp_summary2 <- data.frame(check.names = F, 
                             Database = NA, 
                             Kingdom = NA, 
                             Phylum = NA, 
                             Class = NA, 
                             Order = NA,
                             Family = NA, 
                             Genus = NA, 
                             Species = NA)
  
  for (i in seq(length(tmp_list))) {
    
    tmp_db <- merge(tmp_list[[i]], tmp_count, by = "seq") %>% mutate(RowSum = rowSums(select(., where(is.numeric)))) %>% select(Kingdom:Species, RowSum)
    tmp_name <- tmp_names[[i]]
    
    tmp_ksum <- round((tmp_db %>% filter(!is.na(Kingdom)) %>% pull(RowSum) %>% sum() / tmp_db %>% pull(RowSum) %>% sum()) * 100, 2)
    tmp_psum <- round((tmp_db %>% filter(!is.na(Phylum)) %>% pull(RowSum) %>% sum() / tmp_db %>% pull(RowSum) %>% sum()) * 100, 2)
    tmp_csum <- round((tmp_db %>% filter(!is.na(Class)) %>% pull(RowSum) %>% sum() / tmp_db %>% pull(RowSum) %>% sum()) * 100, 2)
    tmp_osum <- round((tmp_db %>% filter(!is.na(Order)) %>% pull(RowSum) %>% sum() / tmp_db %>% pull(RowSum) %>% sum()) * 100, 2)
    tmp_fsum <- round((tmp_db %>% filter(!is.na(Family)) %>% pull(RowSum) %>% sum() / tmp_db %>% pull(RowSum) %>% sum()) * 100, 2)
    tmp_gsum <- round((tmp_db %>% filter(!is.na(Genus)) %>% pull(RowSum) %>% sum() / tmp_db %>% pull(RowSum) %>% sum()) * 100, 2)
    tmp_ssum <- round((tmp_db %>% filter(!is.na(Species)) %>% pull(RowSum) %>% sum() / tmp_db %>% pull(RowSum) %>% sum()) * 100, 2)
    
    tmp_summary2 <- tmp_summary2 %>% add_row(Database = tmp_name, 
                                             Kingdom = tmp_ksum, 
                                             Phylum = tmp_psum, 
                                             Class = tmp_csum, 
                                             Order = tmp_osum,
                                             Family = tmp_fsum, 
                                             Genus = tmp_gsum, 
                                             Species = tmp_ssum)
    
  }
  
  tmp_summary2 <- tmp_summary2[-1,]
  
  return(list(tmp_summary1, tmp_summary2))
  
}

fix.tax <- function(tax) {
  
  tmp <- as.data.frame(check.names = F, tax) %>% dplyr::mutate(across(Kingdom:Species, ~ str_replace_na(., "Unclassified")))
  
  tmp[2] <- str_remove_all(tmp[[2]],"k__")
  tmp[3] <- str_remove_all(tmp[[3]],"p__")
  tmp[4] <- str_remove_all(tmp[[4]],"c__")
  tmp[5] <- str_remove_all(tmp[[5]],"o__")
  tmp[6] <- str_remove_all(tmp[[6]],"f__")
  tmp[7] <- str_remove_all(tmp[[7]],"g__")
  tmp[8] <- str_remove_all(tmp[[8]],"s__")
  
  tmp[2] <- paste0("k__",tmp[[2]])
  tmp[3] <- paste0("p__",tmp[[3]])
  tmp[4] <- paste0("c__",tmp[[4]])
  tmp[5] <- paste0("o__",tmp[[5]])
  tmp[6] <- paste0("f__",tmp[[6]])
  tmp[7] <- paste0("g__",tmp[[7]])
  tmp[8] <- paste0("s__",tmp[[8]]) 
  
  names(tmp) <- c("ID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  
  tmp[grep("_Unclassified$", tmp$Species), 8] <- "s__Unclassified"
  
  tmp$Species <- sub("s__\\(.*", NA, tmp[,8])
  tmp$Species <- gsub("\\(.*", "", tmp[,8])
  
  return(tmp)
  
}

# ASV infering ------------------------------------------------------------

reads <- list.files("analysis/final/", pattern = ".fastq", full.names=TRUE)
reads

sample.names <- paste0(sapply(strsplit(basename(reads), "_final.fastq"), `[`, 1))
sample.names

reads.filt <- file.path(path, "filt", paste0(sample.names, ".fastq"))
names(reads.filt) <- sample.names
reads.filt

qc <- filterAndTrim(fwd = reads, 
                    filt = reads.filt, 
                    truncLen = 0, 
                    maxEE = 4, 
                    maxN = 0, 
                    rm.phix = T, 
                    compress = F, 
                    multithread = T)
qc

errors.model <- learnErrors(reads.filt, multithread = T) 

plotErrors(errors.model, nominalQ=TRUE)

reads.derep <- derepFastq(reads.filt) 

reads.dada <- dada(reads.derep, err = errors.model, pool = T, multithread = T) 

asvs <- makeSequenceTable(reads.dada)
dim(asvs)

asvs.nochim <- removeBimeraDenovo(asvs, method = "pooled", multithread = T)
dim(asvs.nochim)

getN <- function(x) sum(getUniques(x))

count.dada <- data.frame(check.names = F, cbind(sample.names, qc[,2], sapply(reads.dada, getN), rowSums(asvs.nochim)))
count.dada

colnames(count.dada) <- c("ID", "filterAndTrim", "Denoised", "ASVs")
rownames(count.dada) <- NULL

count.preproc <- read.table("./misc/counts.tsv", header = T, sep = "\t") 
count.preproc

count.proc <- merge(count.preproc, count.dada, by = "ID") %>% mutate_at(-1, as.numeric)
count.proc

# Taxonomic Assignment ----------------------------------------------------

gc()

# # RefRDP (NCBI RefSeq 16S rRNA database supplemented by RDP):
tax.refseq.rdp <- assignTaxonomy(asvs.nochim, "~/work_lamoroso/db/16S/RefSeq_16S_6-11-20_RDPv16_fullTaxo.fa.gz", multithread = T)
tax.refseq.rdp <- addSpecies(tax.refseq.rdp, "~/work_lamoroso/db/16S/RefSeq_16S_6-11-20_RDPv16_Species.fa.gz", allowMultiple = F)
gc()

# # RDP (Ribosomal Database Project):
tax.rdp.18 <- assignTaxonomy(asvs.nochim, "~/work_lamoroso/db/16S/rdp_train_set_18.fa.gz", multithread = T)
tax.rdp.18 <- addSpecies(tax.rdp.18, "~/work_lamoroso/db/16S/rdp_species_assignment_18.fa.gz", allowMultiple = F)
gc()

tax.rdp.19 <- assignTaxonomy(asvs.nochim, "~/work_lamoroso/db/16S/rdp_19_toSpecies_trainset.fa.gz", multithread = T)
gc()

# SILVA: 
tax.silva.138.1 <- assignTaxonomy(asvs.nochim, "~/work_lamoroso/db/16S/silva_nr99_v138.1_train_set.fa.gz", multithread = T)
tax.silva.138.1 <- addSpecies(tax.silva.138.1, "~/work_lamoroso/db/16S/silva_species_assignment_v138.1.fa.gz", allowMultiple = F)
gc()

tax.silva.138.2 <- assignTaxonomy(asvs.nochim, "~/work_lamoroso/db/16S/silva_nr99_v138.2_toGenus_trainset.fa.gz", multithread = T)
tax.silva.138.2 <- addSpecies(tax.silva.138.2, "~/work_lamoroso/db/16S/silva_v138.2_assignSpecies.fa.gz", allowMultiple = F)
gc()


# GTDB
tax.gtdb.202 <- assignTaxonomy(asvs.nochim, "~/work_lamoroso/db/16S/GTDB_bac120_arc122_ssu_r202_fullTaxo.fa.gz", multithread = T)
tax.gtdb.202 <- addSpecies(tax.gtdb.202, "~/work_lamoroso/db/16S/GTDB_bac120_arc122_ssu_r202_Species.fa.gz", allowMultiple = F)
gc()

tax.gtdb.220 <- assignTaxonomy(asvs.nochim, "~/work_lamoroso/db/16S/GTDB_bac120_arc53_ssu_r220_fullTaxo.fa.gz", multithread = T)
tax.gtdb.220 <- addSpecies(tax.gtdb.220, "~/work_lamoroso/db/16S/GTDB_bac120_arc53_ssu_r220_species.fa.gz", allowMultiple = F)
gc()

# GG
tax.gg2 <- assignTaxonomy(asvs.nochim, "~/work_lamoroso/db/16S/gg2_2024_09_toSpecies_trainset.fa.gz", multithread = T)
gc()

# Tax DFying

tax.refseq.rdp.df <- tax.dfying(tax.refseq.rdp)
tax.rdp.18.df <- tax.dfying(tax.rdp.18)
tax.rdp.19.df <- tax.dfying(tax.rdp.19)
tax.silva.138.1.df <- tax.dfying(tax.silva.138.1)
tax.silva.138.2.df <- tax.dfying(tax.silva.138.2)
tax.gtdb.202.df <- tax.dfying(tax.gtdb.202)
tax.gtdb.220.df <- tax.dfying(tax.gtdb.220)
tax.gg2.df <- tax.dfying(tax.gg2)

counts.df <- data.frame(check.names = F, t(asvs.nochim)) %>% 
  mutate(seq = rownames(.)) %>% 
  select(seq, everything())

# Filtering

list.tax.df <- list("RefSeq-RDP" = tax.refseq.rdp.df, 
                    "RDPv18" = tax.rdp.18.df,
                    "RDPv19" = tax.rdp.19.df,
                    "Silvav138-1" = tax.silva.138.1.df,
                    "Silvav138-2" = tax.silva.138.2.df,
                    "GTDBv202" = tax.gtdb.202.df,
                    "GTDBv220" = tax.gtdb.220.df,
                    "GG2" = tax.gg2.df)


for (i in seq_along(names(list.tax.df))) {
  
  tmp.tax.df <- 
    list.tax.df[[i]] %>% 
    mutate(!!names(list.tax.df)[[i]] := paste(Kingdom, Phylum, Class, Order, Family, Genus, Species, sep = "; ") %>% str_remove_all(., "NA")) %>% 
    select(seq, any_of(names(list.tax.df)[[i]]))
  
  if (i == 1) {
    
    tax.all.df <- tmp.tax.df
    
  } else {
    
    tax.all.df <- merge(tax.all.df, tmp.tax.df, by = "seq")
    
  }
  
}

tax.all.df <- 
  merge(tax.all.df, counts.df, by = "seq") %>% 
  mutate(RowSum = rowSums(select(., where(is.numeric)))) %>%
  arrange(-RowSum) %>%
  select(-RowSum)

to.remove <- c("Chloroplast;", "Mitochondria;", "^Unclassified;", "^Eukaryota;", "^;")

seqs.to.keep <- 
  tax.all.df %>% 
  select(seq, any_of(c("RefSeq-RDP", "RDPv19", "Silvav138-2", "GG2"))) %>% 
  filter_at(vars(-seq), all_vars(str_detect(., paste(to.remove, collapse = "|"), negate = T))) %>% 
  pull(seq)

tax.refseq.rdp.df.filt <- tax.refseq.rdp.df %>% filter(seq %in% seqs.to.keep)
tax.rdp.18.df.filt <- tax.rdp.18.df %>% filter(seq %in% seqs.to.keep)
tax.rdp.19.df.filt <- tax.rdp.19.df %>% filter(seq %in% seqs.to.keep)
tax.silva.138.1.df.filt <- tax.silva.138.1.df %>% filter(seq %in% seqs.to.keep)
tax.silva.138.2.df.filt <- tax.silva.138.2.df %>% filter(seq %in% seqs.to.keep)
tax.gtdb.202.df.filt <- tax.gtdb.202.df %>% filter(seq %in% seqs.to.keep)
tax.gtdb.220.df.filt <- tax.gtdb.220.df %>% filter(seq %in% seqs.to.keep)
tax.gg2.df.filt <- tax.gg2.df %>% filter(seq %in% seqs.to.keep)
counts.df.filt <- counts.df %>% filter(seq %in% seqs.to.keep)

list.tax.df <- list("RefSeq-RDP" = tax.refseq.rdp.df.filt, 
                    "RDPv18" = tax.rdp.18.df.filt,
                    "RDPv19" = tax.rdp.19.df.filt,
                    "Silvav138-1" = tax.silva.138.1.df.filt,
                    "Silvav138-2" = tax.silva.138.2.df.filt,
                    "GTDBv202" = tax.gtdb.202.df.filt,
                    "GTDBv220" = tax.gtdb.220.df.filt,
                    "GG2" = tax.gg2.df.filt)

for (i in seq_along(names(list.tax.df))) {
  
  tmp.tax.df <- 
    list.tax.df[[i]] %>% 
    mutate(!!names(list.tax.df)[[i]] := paste(Kingdom, Phylum, Class, Order, Family, Genus, Species, sep = "; ") %>% str_remove_all(., "NA")) %>% 
    select(seq, any_of(names(list.tax.df)[[i]]))
  
  if (i == 1) {
    
    tax.all.df.filt <- tmp.tax.df
    
  } else {
    
    tax.all.df.filt <- merge(tax.all.df.filt, tmp.tax.df, by = "seq")
    
  }
  
}

tax.all.df.filt <- 
  merge(tax.all.df.filt, counts.df.filt, by = "seq") %>% 
  mutate(RowSum = rowSums(select(., where(is.numeric)))) %>%
  arrange(-RowSum) %>%
  select(-RowSum)

# Summary

tax.class <- summary.tax(ordered_names = names(list.tax.df), 
                         tax.refseq.rdp.df.filt, 
                         tax.rdp.18.df.filt, 
                         tax.rdp.19.df.filt,
                         tax.silva.138.1.df.filt,
                         tax.silva.138.2.df.filt,
                         tax.gtdb.202.df.filt,
                         tax.gtdb.220.df.filt,
                         tax.gg2.df.filt)

tax.class.asvs <- tax.class[[1]]
tax.class.asvs

tax.class.uniques <- tax.class[[2]]
tax.class.uniques

tax.class.2 <- summary.tax.2(count_table = counts.df.filt,
                             ordered_names = names(list.tax.df), 
                             tax.refseq.rdp.df.filt, 
                             tax.rdp.18.df.filt, 
                             tax.rdp.19.df.filt,
                             tax.silva.138.1.df.filt,
                             tax.silva.138.2.df.filt,
                             tax.gtdb.202.df.filt,
                             tax.gtdb.220.df.filt,
                             tax.gg2.df.filt)

tax.class2.abs <- tax.class.2[[1]]
tax.class2.abs

tax.class2.rel <- tax.class.2[[2]]
tax.class2.rel

count.proc$`ASVs (Clean)` <- counts.df.filt %>% select(2:ncol(.)) %>% colSums()

master.df <- merge(tax.gg2.df.filt, counts.df.filt) %>%
  mutate("Count (Sum)" = rowSums(.[9:ncol(.)])) %>%
  mutate("Count (Prevalence)" = rowSums(.[9:ncol(.)] != 0) - 1) %>%
  arrange(-`Count (Sum)`) %>% 
  filter(`Count (Prevalence)` > 2) %>% 
  mutate(ID = sprintf("%s%04d", "ASV", 1:n())) %>%
  select(ID, seq, everything())

tax.class.perc <- data.frame(check.names = F, 
                             "Tax" = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                             "Perc" = c(
                               round(sum(master.df %>% filter(!is.na(Kingdom)) %>% pull(`Count (Sum)`)) / sum(master.df$`Count (Sum)`) * 100, 2),
                               round(sum(master.df %>% filter(!is.na(Phylum)) %>% pull(`Count (Sum)`)) / sum(master.df$`Count (Sum)`) * 100, 2),
                               round(sum(master.df %>% filter(!is.na(Class)) %>% pull(`Count (Sum)`)) / sum(master.df$`Count (Sum)`) * 100, 2),
                               round(sum(master.df %>% filter(!is.na(Order)) %>% pull(`Count (Sum)`)) / sum(master.df$`Count (Sum)`) * 100, 2),
                               round(sum(master.df %>% filter(!is.na(Family)) %>% pull(`Count (Sum)`)) / sum(master.df$`Count (Sum)`) * 100, 2),
                               round(sum(master.df %>% filter(!is.na(Genus)) %>% pull(`Count (Sum)`)) / sum(master.df$`Count (Sum)`)  * 100, 2),
                               round(sum(master.df %>% filter(!is.na(Species)) %>% pull(`Count (Sum)`)) / sum(master.df$`Count (Sum)`)  * 100, 2)))
tax.class.perc

count.proc <- 
  merge(count.proc, master.df %>% select(10:ncol(.)) %>% colSums() %>% data.frame(check.names = F), by.x = "ID", by.y = 0, all.x = T) %>% 
  dplyr::rename("ASVs (Filt)" = ".") %>%
  mutate(`Usable (%)` = round((`ASVs (Filt)` / Raw) * 100, 2))

count.proc

# Phylogenetic tree -------------------------------------------------------

pt.seqs <- master.df$seq
names(pt.seqs) <- master.df$ID
pt.alignment <- AlignSeqs(DNAStringSet(pt.seqs), anchor = NA)
pt.phang.align <- phyDat(as(pt.alignment, "matrix"), type = "DNA")
pt.dm <- dist.ml(pt.phang.align)
pt.nj <- NJ(pt.dm) # Note, tip order != sequence order
pt.fit <- pml(pt.nj, data = pt.phang.align)

pt.fit.gtr <- update(pt.fit, k = 4, inv = 0.2)
# pt.fit.gtr <- optim.pml(pt.fit.gtr, model = "GTR", optInv = TRUE, optGamma = TRUE, rearrangement = "stochastic", control = pml.control(trace = 0))
# 
# pt.bs <- bootstrap.pml(pt.fit.gtr, bs = 100, optNni = TRUE, multicore = TRUE, mc.cores = 12, control = pml.control(trace = 0), verbose = T)
# pt.tree.bs <- plotBS(pt.fit.gtr$tree, "none")

# Microbiome analysis tables ----------------------------------------------

tb.metadata <- 
  read.delim(file = "misc/metadata.tsv") %>% 
  mutate(across(where(is.character), ~str_trim(.))) %>% 
  select(`#NAME` = ID, everything())

tb.count <- 
  master.df %>% 
  select(ID, tb.metadata$`#NAME`)

tb.taxonomy <- 
  master.df %>% 
  select(ID, Kingdom:Species) %>% 
  fix.tax() %>% 
  dplyr::rename("#TAXONOMY" = "ID")

write.table(tb.count, file = "results/tables/TB01_counts.txt", quote = F, sep = "\t", col.names = T, row.names = F)
write.table(tb.metadata, file = "results/tables/TB02_metadata.txt", quote = F, sep = "\t", col.names = T, row.names = F)
write.table(tb.taxonomy, file = "results/tables/TB03_taxonomy.txt", quote = F, sep = "\t", col.names = T, row.names = F)
write.tree(pt.fit.gtr$tree, file = "results/tables/TB04_tree.nwk")

master.df <- master.df %>% select(ID:Species, tb.metadata$`#NAME`, `Count (Sum)`, `Count (Prevalence)`)

# Reads Tracking by Group -------------------------------------------------

count.proc.group <- 
  merge(tb.metadata %>% select(`#NAME`, Group), count.proc, by.x = "#NAME", by.y = "ID") %>%
  group_by(Group) %>%
  summarise_if(is.numeric, list(`(mean)` = mean, `(sd)` = sd), na.rm = TRUE) %>%
  data.frame(check.names = F) %>%
  mutate_if(is.numeric, round, 2) %>%
  select(Group, 
         starts_with("Raw"),
         starts_with("filterAndTrim"),
         starts_with("Denoised"),
         starts_with("ASVs"),
         starts_with("ASVs (Clean)"),
         starts_with("ASVs (Filt)"))

proc.tables <- list("Master" = master.df, 
                    "Tracking" = count.proc,
                    "Tracking (Grouped)" = count.proc.group,
                    "TaxClass (ASVs)" = tax.class.asvs,
                    "TaxClass (Uniques)" = tax.class.uniques,
                    "TaxClass (%)" = tax.class2.rel,
                    "All TaxDBs (Contaminated)" = tax.all.df,
                    "All TaxDBs (Clean)" = tax.all.df.filt)

openxlsx::write.xlsx(x = proc.tables, file ="results/tables/Processing.xlsx", rowNames = F)

save.image("analysis/dada2/dada2.RData")

# BLAST -------------------------------------------------------------------

rm(list = ls())

load("analysis/dada2/dada2.RData")

library("tidyverse")
library("dada2")        
library("rBLAST")
library("R.utils")
library("rjson")
library("XML")
library("rentrez")
library("openxlsx")

Sys.setenv(PATH = paste0(Sys.getenv("PATH"), ":/usr/local/bioinfo/"))

blast.work.db <- blast("~/work_lamoroso/db/NCBI/DBs/amplicons/proka.16SrRNA.fna")
blast.headers <- read.delim("~/work_lamoroso/db/NCBI/DBs/amplicons/headers-amplicons.tsv") %>% dplyr::rename(sseqid = 1, description = 2)
blast.my.seqs <- master.df %>% select(ID, seq) %>% column_to_rownames("ID")
blast.results.table <- data.frame()

for (n in seq(1, nrow(blast.my.seqs))) {
  
  tmp.seq <- blast.my.seqs[n,]
  names(tmp.seq) <- row.names(blast.my.seqs)[n]
  tmp.seq <- DNAStringSet(tmp.seq, use.names = T)
  
  tmp.results.table <- 
    tryCatch(predict(blast.work.db, tmp.seq, BLAST_args = paste0('-max_target_seqs ', 5)), 
             error = function(e) {}, 
             warning = function(w) {})
  
  if (nrow(tmp.results.table)>0) {
    
    tmp.results.table.tax <- 
      merge(tmp.results.table, blast.headers, sort = F) %>%
      select(qseqid, sseqid, description, everything())
    
    print(tmp.results.table.tax)
    
    blast.results.table <- bind_rows(blast.results.table, tmp.results.table.tax)
    
    rm(tmp.results.table)
    
  }
  
}

blast.results.table.bkp <- blast.results.table

blast.results.table <- blast.results.table %>% 
  select(qseqid, sseqid, description, bitscore, everything())

rm(list = ls(pattern = "^tmp."))

View(blast.results.table)

blast.results.table.full <- full_join(x = blast.my.seqs %>% rownames_to_column("qseqid"), y = blast.results.table)

blast.results.table.best <-
  blast.results.table.full %>% 
  group_by(qseqid) %>% 
  top_n(., 1, bitscore) %>% 
  top_n(., 1, pident) %>% 
  ungroup() %>% 
  data.frame()

blast.tables <- list("BLAST (Best hit)" = blast.results.table.best, 
                     "BLAST (Top 5 hits)" = blast.results.table.full)

write.xlsx(x = blast.tables, file ="results/tables/BLAST.xlsx", rowNames = F)

# # Feature data and Seq table ----------------------------------------------

featureTable <- master.df %>% select(ID, tb.metadata$`#NAME`) %>% dplyr::rename("#OTU ID" = "ID")
seqTable <- master.df %>% select(ID, seq) %>% mutate(ID = str_replace(ID, "^", ">"))
seqTable <- as.vector(t(seqTable))

write.table(featureTable, "picrust/feature-table.tsv", quote = F, row.names = F, col.names = T, sep = "\t")
write.table(seqTable, "picrust/dna-sequences.fasta", quote = F, row.names = F, col.names = F, sep = "\t")

save.image("analysis/dada2/dada2.RData")

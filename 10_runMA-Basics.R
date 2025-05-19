# Microbiome analysis -----------------------------------------------------

rm(list = ls())

# Setup -------------------------------------------------------------------

library(rstatix)
library(tidyverse)
library(phyloseq)
library(metagMisc)
library(microbiome)
library(ampvis2)
library(artyfarty)
library(scales)
library(svglite)
library(ggpubr)
library(ggdendro)
library(reshape2)
library(vegan)
library(ape)
library(jcolors)
library(xlsx)
library(ggplot2) 
library(grid)
library(gridExtra)
library(ggbeeswarm)
library(corncob)
library(DESeq2)
library(patchwork)
library(cowplot)
library(agricolae)
library(pairwiseAdonis)
library(picante)
library(ggforestplot)
library(ggvenn)
library(ggtext)
library(extrafont) 
library(UpSetR)

#font_import()
loadfonts(device = "all")

set.seed("123")

# Data Wrangling ----------------------------------------------------------

getwd()
setwd("~/sme_rhizobiome/")
getwd()

path <- paste0(getwd(), "/results/")
path

list.files(path)

tb.metad <- 
  read.delim(file = "misc/metadata.tsv") %>% 
  mutate(across(where(is.character), ~str_trim(.))) %>% 
  mutate(across(where(is.character), ~str_replace_all(., "Spirulina maxima", "*Spirulina maxima*"))) %>% 
  mutate(across(where(is.character), ~str_replace_all(., "-1\\)", "<sup>-1</sup>)"))) %>% 
  mutate(across(where(is.character), ~str_replace_all(., "-", "−"))) %>% 
  mutate(FullGroup = paste0(Group, ": ", Treat2)) %>% 
  select(-oldID) %>% 
  select(`#NAME` = ID, everything()) %>% 
  dplyr::rename("ID" = "#NAME")

tb.count <- 
  read.delim(file = paste0(path, "tables/TB01_counts.txt"), header = T, sep = "\t", stringsAsFactors = F, check.names = F)

names(tb.count) <- c("ID", as.character(tb.metad$ID))

tb.metad <- 
  tb.metad %>% 
  mutate(across(.cols = everything(), .fns = ~factor(., unique(.))))

tb.count <- tb.count %>% select(ID, tb.metad$ID)

tb.taxon <- read.delim(file = paste0(path, "tables/TB03_taxonomy.txt"), header = T, sep = "\t", stringsAsFactors = F, check.names = F) %>% 
  dplyr::rename("ID" = "#TAXONOMY")

tb.tree <- read.tree(file = paste0(path, "tables/TB04_tree.nwk"))

row.names(tb.count) <- tb.count$ID
row.names(tb.taxon) <- tb.taxon$ID
row.names(tb.metad) <- tb.metad$ID

tb.count <- tb.count %>% select(-ID)
tb.taxon <- tb.taxon %>% select(-ID)

tb.count <- as.matrix(tb.count)
tb.taxon <- as.matrix(tb.taxon)

ps.count <- otu_table(tb.count, taxa_are_rows = TRUE)
ps.taxon <- tax_table(tb.taxon)
ps.metad <- sample_data(tb.metad)
ps.tree <- phy_tree(tb.tree)

ps.table <- phyloseq(ps.count, ps.taxon, ps.metad, ps.tree)

# Setup graphs ------------------------------------------------------------

tr.colors <- c("#525252FF", "#8cc084FF", "#38a3a5FF", "#0466c8FF")
names(tr.colors) <- unique(tb.metad$FullGroup)

addSmallLegend <- function(myPlot, pointSize = 3, textSize = 8, spaceLegend = 0.8) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_markdown(size = textSize + 2, face = "bold"),
          legend.text  = element_markdown(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}

my_theme <-   theme(legend.position = "right", 
                    legend.text = element_markdown(size = 8, colour = 'black'),
                    legend.title = element_blank(),
                    legend.justification = c(0.5, 0.5),
                    
                    strip.text = element_markdown(size = 10, colour = 'black', face = "bold"),
                    strip.background = element_rect(fill = NA, color = NA), 
                    
                    panel.background = element_rect(fill = "white"),
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.border = element_rect(color = "black", linewidth = 1, fill = NA),
                    
                    axis.ticks = element_line(color = "black", linewidth = .5),
                    axis.ticks.length = unit(2, "pt"),
                    
                    axis.title.x = element_markdown(size = 9, face = 'bold', colour = 'black'),
                    axis.title.y = element_markdown(size = 9, face = 'bold', colour = 'black'),
                    axis.title.y.right = element_markdown(size = 9, face = 'bold', colour = 'black'),
                    axis.text.x = element_markdown(size = 7, colour = 'black'),
                    axis.text.y = element_markdown(size = 7, colour = 'black'),
                    
                    plot.title = element_markdown(size = 10, face = 'bold', colour = 'black'),
                    plot.subtitle = element_markdown(size = 8, colour = 'black'),
                    plot.margin = unit(c(0.5, 0.75, 0.5, 0.75), 'cm'))

# Rarefaction curve -------------------------------------------------------

amp.count <- data.frame(check.names = F, 'OTU' = row.names(otu_table(ps.table)), otu_table(ps.table))
amp.taxon <- data.frame(check.names = F, 'OTU' = row.names(tax_table(ps.table)), tax_table(ps.table))
amp.metad <- data.frame(check.names = F, sample_data(ps.table))

amp.count.taxon <- merge(amp.count, amp.taxon, by = 'OTU')
amp.table <- amp_load(otutable = amp.count.taxon, metadata = amp.metad)

gg.rarecurve <- 
  amp_rarecurve(amp.table, stepsize = 1000, color_by = "FullGroup") + 
  geom_point(x = NA) + 
  theme_bw(base_family = "Liberation Sans") + 
  geom_vline(xintercept = min(sample_sums(ps.table)), color='black', lty='dashed', lwd=0.5) +
  facet_wrap( ~ Group, scales = "free", nrow = 2, ncol = 2) +
  scale_colour_manual(name="", values = tr.colors) +
  scale_x_continuous(expand = expansion(mult = .01, add = 0), 
                     limits = c(0, max(sample_sums(ps.table))),
                     breaks = seq(0, max(sample_sums(ps.table)), 15000),
                     labels = scales::number_format(accuracy = 1, scale = 1/1000, suffix = "k", big.mark = ",")) +
  scale_y_continuous(expand = expansion(mult = c(0, .05), add = c(0, 0))) +
  ylab('Number of observed ASVs<br>') +
  xlab('<br>Sequencing depth') +
  guides(colour = 
           guide_legend(order = 1, 
                        override.aes = list(shape = 22, 
                                            fill = tr.colors, 
                                            color = "black", 
                                            size = 3,
                                            linetype = 0),
                        ncol = 1)) +
  my_theme + 
  theme(legend.position = "right", 
        legend.box.background = element_rect(fill = "white", colour = "black", linewidth = .75))

gg.rarecurve

ggsave(filename = paste0(path, "figures/RareCurve.svg"), 
       plot = gg.rarecurve, 
       dpi = 1200,
       width = 5000 * 1.75,
       height = 5000 * 1.15,
       units = "px")

ggsave(filename = paste0(path, "figures/RareCurve.png"), 
       plot = gg.rarecurve, 
       dpi = 1200,
       width = 5000 * 1.75,
       height = 5000 * 1.15,
       units = "px")

rm(amp.count, amp.taxon, amp.metad, amp.count.taxon, amp.table, gg.rarecurve)

# Rarefaction -------------------------------------------------------------

ps.rarefied <- rarefy_even_depth(ps.table, rngseed = 123, sample.size = min(sample_sums(ps.table)), replace = F, trimOTUs = F)

ps.count.rare <- 
  data.frame(otu_table(ps.rarefied), check.names = F) %>% 
  rownames_to_column("ID") %>%
  arrange(ID) %>%
  column_to_rownames("ID") %>% 
  as.matrix() %>% 
  otu_table(taxa_are_rows = TRUE)

ps.rarefied <- phyloseq(ps.count.rare, ps.taxon, ps.metad, ps.tree)
ps.rarefied <- prune_taxa(taxa_sums(ps.rarefied) > 0, ps.rarefied)
ps.rarefied.tss <- phyloseq_standardize_otu_abundance(physeq = ps.rarefied, method = 'total')

tb.master <- 
  openxlsx::read.xlsx(xlsxFile = paste0(path, "tables/Processing.xlsx"), sheet = 1, check.names = F)

tb.master.raw <- tb.master %>% 
  mutate(across(where(is.character), ~str_replace_all(., " ", "_"))) %>% 
  mutate(across(where(is.character), ~str_replace_all(., "_sp$", "_sp."))) %>% 
  select(ID:Species, tb.metad$ID) %>% 
  mutate("Count (Sum)" = rowSums(.[10:ncol(.)])) %>%
  mutate("Count (Prevalence)" = rowSums(.[10:ncol(.)] != 0) - 1)

tb.master.raref <- tb.master %>% 
  mutate(across(where(is.character), ~str_replace_all(., " ", "_"))) %>% 
  mutate(across(where(is.character), ~str_replace_all(., "_sp$", "_sp."))) %>% 
  select(ID:Species) %>% 
  filter(ID %in% rownames(otu_table(ps.rarefied)))

tb.master.raref <- 
  base::merge(tb.master.raref, data.frame(otu_table(ps.rarefied), check.names = F), by.x = "ID", by.y = 0, all.x = T) %>% 
  mutate("Count (Sum)" = rowSums(.[10:ncol(.)])) %>%
  mutate("Count (Prevalence)" = rowSums(.[10:ncol(.)] != 0) - 1) %>% 
  arrange(-`Count (Sum)`)

count.proc <-
  openxlsx::read.xlsx(xlsxFile = paste0(path, "tables/Processing.xlsx"), sheet = 2, check.names = F)

count.proc <-
  merge(count.proc, tb.master.raref %>% select(10:ncol(.)) %>% colSums() %>% data.frame(check.names = F), by.x = "ID", by.y = 0, all.x = T) %>% 
  dplyr::rename("ASVs (Raref)" = ".") 

names(count.proc) <- names(count.proc) %>% str_replace_all(., "\\.", " ")

count.proc.group <- 
  merge(tb.metad %>% select(ID, Group), count.proc, by = "ID") %>%
  group_by(Group) %>%
  summarise_if(is.numeric, list(`(mean)` = mean, `(sd)` = sd), na.rm = TRUE) %>%
  data.frame(check.names = F) %>%
  mutate_if(is.numeric, round, 2) %>% 
  select(-contains("Usable"))

names(count.proc.group) <- names(count.proc.group) %>% str_replace_all(., "_", " ")

count.proc.group <- 
  count.proc.group %>% 
  select(Group = Group, 
         starts_with("Raw"),
         starts_with("wAdapt"),
         starts_with("QC"),
         starts_with("Merged"),
         starts_with("Trunc"),
         starts_with("filterAndTrim"),
         starts_with("Denoised"),
         starts_with("ASVs (mean)"),
         starts_with("ASVs (sd)"),
         starts_with("ASVs (Clean)"),
         starts_with("ASVs (Filt)"),
         starts_with("ASVs (Raref)")) %>% 
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.), NA, .)))

perc.class <- function(master) {
  
  perc <- data.frame(check.names = F, 
                     "Rank" = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                     "Perc" = c(
                       round(sum(master %>% filter(!str_detect(Kingdom, "Unclassified")) %>% pull(`Count (Sum)`)) / sum(master$`Count (Sum)`) * 100, 2),
                       round(sum(master %>% filter(!str_detect(Phylum, "Unclassified")) %>% pull(`Count (Sum)`)) / sum(master$`Count (Sum)`) * 100, 2),
                       round(sum(master %>% filter(!str_detect(Class, "Unclassified")) %>% pull(`Count (Sum)`)) / sum(master$`Count (Sum)`) * 100, 2),
                       round(sum(master %>% filter(!str_detect(Order, "Unclassified")) %>% pull(`Count (Sum)`)) / sum(master$`Count (Sum)`) * 100, 2),
                       round(sum(master %>% filter(!str_detect(Family, "Unclassified")) %>% pull(`Count (Sum)`)) / sum(master$`Count (Sum)`) * 100, 2),
                       round(sum(master %>% filter(!str_detect(Genus, "Unclassified")) %>% pull(`Count (Sum)`)) / sum(master$`Count (Sum)`)  * 100, 2),
                       round(sum(master %>% filter(!str_detect(Species, "Unclassified")) %>% pull(`Count (Sum)`)) / sum(master$`Count (Sum)`)  * 100, 2)))
  return(perc)
}

perc.raw <- perc.class(tb.master.raw)
perc.raref <- perc.class(tb.master.raref)

ls.ann <- list("Master (Raw)" = tb.master.raw, 
               "%Classified (Raw)" = perc.raw, 
               "Master (Rarefied)" = tb.master.raref, 
               "%Classified (Rarefied)" = perc.raref, 
               "Tracking" = count.proc, 
               "Tracking (Grouped)" = count.proc.group)

openxlsx::write.xlsx(x = ls.ann, file = paste0(path, "tables/MasterTable.xlsx"), rowNames = F)

# New MA tables

ps.rarefied

otu_table(ps.rarefied) %>% as.data.frame() %>% rownames_to_column("#NAME") %>% arrange(`#NAME`)
tax_table(ps.rarefied) %>% as.data.frame() %>% rownames_to_column("#TAXONOMY") %>% arrange(`#TAXONOMY`)
phy_tree(ps.rarefied)

write.table(otu_table(ps.rarefied) %>% as.data.frame() %>% rownames_to_column("#NAME") %>% arrange(`#NAME`), file = "results/tables/rTB01_counts.txt", quote = F, sep = "\t", col.names = T, row.names = F)
write.table(sample_data(ps.rarefied), file = "results/tables/rTB02_metadata.txt", quote = F, sep = "\t", col.names = T, row.names = F)
write.table(tax_table(ps.rarefied) %>% as.data.frame() %>% rownames_to_column("#TAXONOMY") %>% arrange(`#TAXONOMY`), file = "results/tables/rTB03_taxonomy.txt", quote = F, sep = "\t", col.names = T, row.names = F)
write.tree(phy_tree(ps.rarefied), file = "results/tables/rTB04_tree.nwk")

# New PICRUSt2 tables

featureTable <- tb.master.raref %>% select(ID, tb.metad$ID) %>% dplyr::rename("#OTU ID" = "ID")
seqTable <- tb.master.raref %>% select(ID, seq) %>% mutate(ID = str_replace(ID, "^", ">"))
seqTable <- as.vector(t(seqTable))

write.table(featureTable, "picrust/feature-table.tsv", quote = F, row.names = F, col.names = T, sep = "\t")
write.table(seqTable, "picrust/dna-sequences.fasta", quote = F, row.names = F, col.names = F, sep = "\t")

rm(tb.count, tb.taxon, tb.tree, ps.count, ps.taxon, ps.metad, ps.tree, ps.table)
rm(ps.count.rare, tb.master, tb.master.raw, tb.master.raref, perc.class, perc.raw, perc.raref, ls.ann, count.proc, count.proc.group)
rm(featureTable, seqTable)

# Taxonomic agglomerate tables --------------------------------------------

fun.tb.glom.ind <- function(my.ps.table, my.ps.table.tss) {
  
  tb.glom.k.abs <- tax_glom(my.ps.table, taxrank="Kingdom"); tb.glom.k.abs <- merge(tax_table(tb.glom.k.abs), otu_table(tb.glom.k.abs), by = "row.names", all = T)[-1] %>% mutate(RowSum = rowSums(select(., where(is.numeric)))) %>% arrange(-RowSum) %>% select(-RowSum)
  tb.glom.p.abs <- tax_glom(my.ps.table, taxrank="Phylum"); tb.glom.p.abs <- merge(tax_table(tb.glom.p.abs), otu_table(tb.glom.p.abs), by = "row.names", all = T)[-1] %>% mutate(RowSum = rowSums(select(., where(is.numeric)))) %>% arrange(-RowSum) %>% select(-RowSum)
  tb.glom.c.abs <- tax_glom(my.ps.table, taxrank="Class"); tb.glom.c.abs <- merge(tax_table(tb.glom.c.abs), otu_table(tb.glom.c.abs), by = "row.names", all = T)[-1] %>% mutate(RowSum = rowSums(select(., where(is.numeric)))) %>% arrange(-RowSum) %>% select(-RowSum)
  tb.glom.o.abs <- tax_glom(my.ps.table, taxrank="Order"); tb.glom.o.abs <- merge(tax_table(tb.glom.o.abs), otu_table(tb.glom.o.abs), by = "row.names", all = T)[-1] %>% mutate(RowSum = rowSums(select(., where(is.numeric)))) %>% arrange(-RowSum) %>% select(-RowSum)
  tb.glom.f.abs <- tax_glom(my.ps.table, taxrank="Family"); tb.glom.f.abs <- merge(tax_table(tb.glom.f.abs), otu_table(tb.glom.f.abs), by = "row.names", all = T)[-1] %>% mutate(RowSum = rowSums(select(., where(is.numeric)))) %>% arrange(-RowSum) %>% select(-RowSum)
  tb.glom.g.abs <- tax_glom(my.ps.table, taxrank="Genus"); tb.glom.g.abs <- merge(tax_table(tb.glom.g.abs), otu_table(tb.glom.g.abs), by = "row.names", all = T)[-1] %>% mutate(RowSum = rowSums(select(., where(is.numeric)))) %>% arrange(-RowSum) %>% select(-RowSum)
  tb.glom.s.abs <- tax_glom(my.ps.table, taxrank="Species"); tb.glom.s.abs <- merge(tax_table(tb.glom.s.abs), otu_table(tb.glom.s.abs), by = "row.names", all = T)[-1] %>% mutate(RowSum = rowSums(select(., where(is.numeric)))) %>% arrange(-RowSum) %>% select(-RowSum)
  
  tb.glom.k <- tax_glom(my.ps.table.tss, taxrank="Kingdom"); tb.glom.k <- merge(tax_table(tb.glom.k), otu_table(tb.glom.k), by = "row.names", all = T)[-1] %>% mutate(RowSum = rowSums(select(., where(is.numeric)))) %>% arrange(-RowSum) %>% select(-RowSum)
  tb.glom.p <- tax_glom(my.ps.table.tss, taxrank="Phylum"); tb.glom.p <- merge(tax_table(tb.glom.p), otu_table(tb.glom.p), by = "row.names", all = T)[-1] %>% mutate(RowSum = rowSums(select(., where(is.numeric)))) %>% arrange(-RowSum) %>% select(-RowSum)
  tb.glom.c <- tax_glom(my.ps.table.tss, taxrank="Class"); tb.glom.c <- merge(tax_table(tb.glom.c), otu_table(tb.glom.c), by = "row.names", all = T)[-1] %>% mutate(RowSum = rowSums(select(., where(is.numeric)))) %>% arrange(-RowSum) %>% select(-RowSum)
  tb.glom.o <- tax_glom(my.ps.table.tss, taxrank="Order"); tb.glom.o <- merge(tax_table(tb.glom.o), otu_table(tb.glom.o), by = "row.names", all = T)[-1] %>% mutate(RowSum = rowSums(select(., where(is.numeric)))) %>% arrange(-RowSum) %>% select(-RowSum)
  tb.glom.f <- tax_glom(my.ps.table.tss, taxrank="Family"); tb.glom.f <- merge(tax_table(tb.glom.f), otu_table(tb.glom.f), by = "row.names", all = T)[-1] %>% mutate(RowSum = rowSums(select(., where(is.numeric)))) %>% arrange(-RowSum) %>% select(-RowSum)
  tb.glom.g <- tax_glom(my.ps.table.tss, taxrank="Genus"); tb.glom.g <- merge(tax_table(tb.glom.g), otu_table(tb.glom.g), by = "row.names", all = T)[-1] %>% mutate(RowSum = rowSums(select(., where(is.numeric)))) %>% arrange(-RowSum) %>% select(-RowSum)
  tb.glom.s <- tax_glom(my.ps.table.tss, taxrank="Species"); tb.glom.s <- merge(tax_table(tb.glom.s), otu_table(tb.glom.s), by = "row.names", all = T)[-1] %>% mutate(RowSum = rowSums(select(., where(is.numeric)))) %>% arrange(-RowSum) %>% select(-RowSum)
  
  ls.ann <- 
    list(
      "Kingdom" = tb.glom.k.abs[, c(1:1,8:dim(tb.glom.k.abs)[2]) ],
      "Phylum" = tb.glom.p.abs[, c(1:2,8:dim(tb.glom.p.abs)[2]) ],
      "Class" = tb.glom.c.abs[, c(1:3,8:dim(tb.glom.c.abs)[2]) ],
      "Order" = tb.glom.o.abs[, c(1:4,8:dim(tb.glom.o.abs)[2]) ],
      "Family" = tb.glom.f.abs[, c(1:5,8:dim(tb.glom.f.abs)[2]) ],
      "Genus" = tb.glom.g.abs[, c(1:6,8:dim(tb.glom.g.abs)[2]) ],
      "Species" = tb.glom.s.abs[, c(1:7,8:dim(tb.glom.s.abs)[2]) ],
      "%Kingdom" = tb.glom.k[, c(1:1,8:dim(tb.glom.k)[2]) ],
      "%Phylum" = tb.glom.p[, c(1:2,8:dim(tb.glom.p)[2]) ],
      "%Class" = tb.glom.c[, c(1:3,8:dim(tb.glom.c)[2]) ],
      "%Order" = tb.glom.o[, c(1:4,8:dim(tb.glom.o)[2]) ],
      "%Family" = tb.glom.f[, c(1:5,8:dim(tb.glom.f)[2]) ],
      "%Genus" = tb.glom.g[, c(1:6,8:dim(tb.glom.g)[2]) ],
      "%Species" = tb.glom.s[, c(1:7,8:dim(tb.glom.s)[2]) ]
    )
  
  openxlsx::write.xlsx(x = ls.ann, file = paste0(path, "tables/TaxCountSample", ".xlsx"), rowNames = F)
  
}

fun.tb.glom.gp <- function(my.ps.table, group) {
  
  tb.glom.k.abs.gp <- merge(data.frame(check.names = F, tax_table(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Kingdom")))[1:1],
                            data.frame(check.names = F, t(otu_table(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Kingdom")))), by = 0)[-1] %>% mutate(RowSum = rowSums(select(., where(is.numeric)))) %>% arrange(-RowSum) %>% select(-RowSum)
  tb.glom.p.abs.gp <- merge(data.frame(check.names = F, tax_table(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Phylum")))[1:2],
                            data.frame(check.names = F, t(otu_table(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Phylum")))), by = 0)[-1] %>% mutate(RowSum = rowSums(select(., where(is.numeric)))) %>% arrange(-RowSum) %>% select(-RowSum)
  tb.glom.c.abs.gp <- merge(data.frame(check.names = F, tax_table(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Class")))[1:3],
                            data.frame(check.names = F, t(otu_table(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Class")))), by = 0)[-1] %>% mutate(RowSum = rowSums(select(., where(is.numeric)))) %>% arrange(-RowSum) %>% select(-RowSum)
  tb.glom.o.abs.gp <- merge(data.frame(check.names = F, tax_table(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Order")))[1:4],
                            data.frame(check.names = F, t(otu_table(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Order")))), by = 0)[-1] %>% mutate(RowSum = rowSums(select(., where(is.numeric)))) %>% arrange(-RowSum) %>% select(-RowSum)
  tb.glom.f.abs.gp <- merge(data.frame(check.names = F, tax_table(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Family")))[1:5],
                            data.frame(check.names = F, t(otu_table(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Family")))), by = 0)[-1] %>% mutate(RowSum = rowSums(select(., where(is.numeric)))) %>% arrange(-RowSum) %>% select(-RowSum)
  tb.glom.g.abs.gp <- merge(data.frame(check.names = F, tax_table(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Genus")))[1:6],
                            data.frame(check.names = F, t(otu_table(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Genus")))), by = 0)[-1] %>% mutate(RowSum = rowSums(select(., where(is.numeric)))) %>% arrange(-RowSum) %>% select(-RowSum)
  tb.glom.s.abs.gp <- merge(data.frame(check.names = F, tax_table(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Species")))[1:7],
                            data.frame(check.names = F, t(otu_table(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Species")))), by = 0)[-1] %>% mutate(RowSum = rowSums(select(., where(is.numeric)))) %>% arrange(-RowSum) %>% select(-RowSum)
  
  tb.glom.k.gp <- merge(data.frame(check.names = F, tax_table(phyloseq_standardize_otu_abundance(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Kingdom"))))[1:1],
                        data.frame(check.names = F, t(otu_table(phyloseq_standardize_otu_abundance(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Kingdom"))))), by = 0)[-1] %>% mutate(RowSum = rowSums(select(., where(is.numeric)))) %>% arrange(-RowSum) %>% select(-RowSum)
  tb.glom.p.gp <- merge(data.frame(check.names = F, tax_table(phyloseq_standardize_otu_abundance(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Phylum"))))[1:2],
                        data.frame(check.names = F, t(otu_table(phyloseq_standardize_otu_abundance(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Phylum"))))), by = 0)[-1] %>% mutate(RowSum = rowSums(select(., where(is.numeric)))) %>% arrange(-RowSum) %>% select(-RowSum)
  tb.glom.c.gp <- merge(data.frame(check.names = F, tax_table(phyloseq_standardize_otu_abundance(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Class"))))[1:3],
                        data.frame(check.names = F, t(otu_table(phyloseq_standardize_otu_abundance(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Class"))))), by = 0)[-1] %>% mutate(RowSum = rowSums(select(., where(is.numeric)))) %>% arrange(-RowSum) %>% select(-RowSum)
  tb.glom.o.gp <- merge(data.frame(check.names = F, tax_table(phyloseq_standardize_otu_abundance(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Order"))))[1:4],
                        data.frame(check.names = F, t(otu_table(phyloseq_standardize_otu_abundance(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Order"))))), by = 0)[-1] %>% mutate(RowSum = rowSums(select(., where(is.numeric)))) %>% arrange(-RowSum) %>% select(-RowSum)
  tb.glom.f.gp <- merge(data.frame(check.names = F, tax_table(phyloseq_standardize_otu_abundance(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Family"))))[1:5],
                        data.frame(check.names = F, t(otu_table(phyloseq_standardize_otu_abundance(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Family"))))), by = 0)[-1] %>% mutate(RowSum = rowSums(select(., where(is.numeric)))) %>% arrange(-RowSum) %>% select(-RowSum)
  tb.glom.g.gp <- merge(data.frame(check.names = F, tax_table(phyloseq_standardize_otu_abundance(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Genus"))))[1:6],
                        data.frame(check.names = F, t(otu_table(phyloseq_standardize_otu_abundance(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Genus"))))), by = 0)[-1] %>% mutate(RowSum = rowSums(select(., where(is.numeric)))) %>% arrange(-RowSum) %>% select(-RowSum)
  tb.glom.s.gp <- merge(data.frame(check.names = F, tax_table(phyloseq_standardize_otu_abundance(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Species"))))[1:7],
                        data.frame(check.names = F, t(otu_table(phyloseq_standardize_otu_abundance(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Species"))))), by = 0)[-1] %>% mutate(RowSum = rowSums(select(., where(is.numeric)))) %>% arrange(-RowSum) %>% select(-RowSum)
  
  ls.ann <- 
    list(
      "Kingdom" = tb.glom.k.abs.gp,
      "Phylum" = tb.glom.p.abs.gp,
      "Class" = tb.glom.c.abs.gp,
      "Order" = tb.glom.o.abs.gp,
      "Family" = tb.glom.f.abs.gp,
      "Genus" = tb.glom.g.abs.gp,
      "Species" = tb.glom.s.abs.gp,
      "%Kingdom" = tb.glom.k.gp,
      "%Phylum" = tb.glom.p.gp,
      "%Class" = tb.glom.c.gp,
      "%Order" = tb.glom.o.gp,
      "%Family" = tb.glom.f.gp,
      "%Genus" = tb.glom.g.gp,
      "%Species" = tb.glom.s.gp
    )
  
  openxlsx::write.xlsx(x = ls.ann, file = paste0(path, "tables/TaxCountGrouped_", str_to_title(group), ".xlsx"), rowNames = F)
  
}

fun.tb.taxclass.gp <- function(my.ps.table, group) {
  
  tb.taxclass.sum <- cbind(
    
    merge(data.frame(check.names = F, tax_table(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Kingdom")))["Kingdom"],
          data.frame(check.names = F, t(otu_table(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Kingdom")))),
          by = 0)[-1] %>% filter(!str_detect(Kingdom, "Unclassified$")) %>% select(-1) %>% colSums() %>% as.data.frame(check.names = F),
    
    merge(data.frame(check.names = F, tax_table(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Phylum")))["Phylum"],
          data.frame(check.names = F, t(otu_table(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Phylum")))),
          by = 0)[-1] %>% filter(!str_detect(Phylum, "Unclassified$")) %>% select(-1) %>% colSums() %>% as.data.frame(check.names = F),
    
    merge(data.frame(check.names = F, tax_table(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Class")))["Class"],
          data.frame(check.names = F, t(otu_table(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Class")))),
          by = 0)[-1] %>% filter(!str_detect(Class, "Unclassified$")) %>% select(-1) %>% colSums() %>% as.data.frame(check.names = F),
    
    merge(data.frame(check.names = F, tax_table(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Order")))["Order"],
          data.frame(check.names = F, t(otu_table(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Order")))),
          by = 0)[-1] %>% filter(!str_detect(Order, "Unclassified$")) %>% select(-1) %>% colSums() %>% as.data.frame(check.names = F),
    
    merge(data.frame(check.names = F, tax_table(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Family")))["Family"],
          data.frame(check.names = F, t(otu_table(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Family")))),
          by = 0)[-1] %>% filter(!str_detect(Family, "Unclassified$")) %>% select(-1) %>% colSums() %>% as.data.frame(check.names = F),
    
    merge(data.frame(check.names = F, tax_table(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Genus")))["Genus"],
          data.frame(check.names = F, t(otu_table(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Genus")))),
          by = 0)[-1] %>% filter(!str_detect(Genus, "Unclassified$")) %>% select(-1) %>% colSums() %>% as.data.frame(check.names = F),
    
    merge(data.frame(check.names = F, tax_table(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Species")))["Species"],
          data.frame(check.names = F, t(otu_table(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Species")))),
          by = 0)[-1] %>% filter(!str_detect(Species, "Unclassified$")) %>% select(-1) %>% colSums() %>% as.data.frame(check.names = F)
    
  )
  
  tb.taxclass.unq <- cbind(
    
    merge(data.frame(check.names = F, tax_table(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Kingdom")))["Kingdom"],
          data.frame(check.names = F, t(otu_table(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Kingdom")))),
          by = 0)[-1] %>% filter(!str_detect(Kingdom, "Unclassified$")) %>% select(-1) %>% mutate_all(~ ifelse(. > 0, 1, .)) %>% colSums() %>% as.data.frame(check.names = F),
    
    merge(data.frame(check.names = F, tax_table(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Phylum")))["Phylum"],
          data.frame(check.names = F, t(otu_table(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Phylum")))),
          by = 0)[-1] %>% filter(!str_detect(Phylum, "Unclassified$")) %>% select(-1) %>% mutate_all(~ ifelse(. > 0, 1, .)) %>% colSums() %>% as.data.frame(check.names = F),
    
    merge(data.frame(check.names = F, tax_table(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Class")))["Class"],
          data.frame(check.names = F, t(otu_table(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Class")))),
          by = 0)[-1] %>% filter(!str_detect(Class, "Unclassified$")) %>% select(-1) %>% mutate_all(~ ifelse(. > 0, 1, .)) %>% colSums() %>% as.data.frame(check.names = F),
    
    merge(data.frame(check.names = F, tax_table(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Order")))["Order"],
          data.frame(check.names = F, t(otu_table(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Order")))),
          by = 0)[-1] %>% filter(!str_detect(Order, "Unclassified$")) %>% select(-1) %>% mutate_all(~ ifelse(. > 0, 1, .)) %>% colSums() %>% as.data.frame(check.names = F),
    
    merge(data.frame(check.names = F, tax_table(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Family")))["Family"],
          data.frame(check.names = F, t(otu_table(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Family")))),
          by = 0)[-1] %>% filter(!str_detect(Family, "Unclassified$")) %>% select(-1) %>% mutate_all(~ ifelse(. > 0, 1, .)) %>% colSums() %>% as.data.frame(check.names = F),
    
    merge(data.frame(check.names = F, tax_table(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Genus")))["Genus"],
          data.frame(check.names = F, t(otu_table(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Genus")))),
          by = 0)[-1] %>% filter(!str_detect(Genus, "Unclassified$")) %>% select(-1) %>% mutate_all(~ ifelse(. > 0, 1, .)) %>% colSums() %>% as.data.frame(check.names = F),
    
    merge(data.frame(check.names = F, tax_table(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Species")))["Species"],
          data.frame(check.names = F, t(otu_table(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Species")))),
          by = 0)[-1] %>% filter(!str_detect(Species, "Unclassified$")) %>% select(-1) %>% mutate_all(~ ifelse(. > 0, 1, .)) %>% colSums() %>% as.data.frame(check.names = F)
    
  )
  
  names(tb.taxclass.sum) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  names(tb.taxclass.unq) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  
  ls.ann <- 
    list(
      "Total" = tb.taxclass.sum,
      "Uniques" = tb.taxclass.unq
    )
  
  openxlsx::write.xlsx(x = ls.ann, file = paste0(path, "tables/TaxCountSummary_", str_to_title(group), ".xlsx"), rowNames = T)
  
}

fun.tb.glom.ind(ps.rarefied, ps.rarefied.tss)

fun.tb.glom.gp(ps.rarefied, "Group")

fun.tb.taxclass.gp(ps.rarefied, "Group")

rm(fun.tb.glom.ind, fun.tb.glom.gp, fun.tb.taxclass.gp)

# Taxabar (Samples) -------------------------------------------------------

tb.glom.p <- tax_glom(ps.rarefied.tss, taxrank="Phylum"); tb.glom.p <- merge(tax_table(tb.glom.p), otu_table(tb.glom.p), by = "row.names", all = T)[-1]

ranks <- c("Phylum")

fun_bestax <- function(df) {
  
  for (n_taxa in 1:nrow(df)) {
    for (n_rank in 1:7) {
      if (str_detect(df[[n_taxa, n_rank]], "Unclassified")) {
        
        df[[n_taxa, n_rank]] <- 
          str_replace(df[[n_taxa, n_rank]], "Unclassified", df[[n_taxa, (n_rank-1)]] %>% str_remove(., ".*__") %>% paste0(., " (Un.)")) %>% 
          str_replace(., " \\(Un\\.\\).*", " (Un.)")
        
      }
    }
  }
  return(df)
}

tb.taxabar <- list(
  fun_bestax(tb.glom.p[, c(1:2, 8:ncol(tb.glom.p))])
)

gg.bp.plots <- list()

for (i in seq_along(ranks)) {
  
  tmp.rank <- ranks[[i]]
  tmp.tb.glom <- tb.taxabar[[i]]
  
  tmp.tb.glom <- tmp.tb.glom %>% 
    unite("taxa", 2:(i+1), sep = "; ") %>% 
    select(-1) %>% 
    mutate("countSum" = rowSums(.[2:ncol(.)]),
           "taxa" = make.unique(taxa),
           "taxa" = str_replace_all(taxa, pattern = "__", replacement = ": "),
           "taxa" = str_replace_all(taxa, pattern = "_", replacement = " ")) %>% 
    arrange(., desc(countSum))
  
  top_in <- 10 ; top_out = top_in + 1
  
  tmp.tb.glom.top <- rbind(tmp.tb.glom[1:top_in, c(-ncol(tmp.tb.glom))], 
                           c("Others", colSums(tmp.tb.glom[top_out:nrow(tmp.tb.glom), 2:c(ncol(tmp.tb.glom)-1)]))) %>% 
    filter(complete.cases(.)) %>% 
    mutate("taxa" = factor(taxa, levels = taxa)) %>% 
    mutate_at(2:ncol(.), as.numeric) %>% 
    merge(x = melt(.), y = tb.metad, by.x = "variable", by.y = "row.names")
  
  tmp.colors <- colorRampPalette(c(
    
    paletteer::paletteer_d("futurevisions::pegasi")
    
  ))(length(unique(tmp.tb.glom.top$taxa))-1)
  
  names(tmp.colors) <- unique(tmp.tb.glom.top$taxa)[1:top_in]
  
  tmp.colors["Others"] <- "#E6E6E6FF"
  
  gg.bp <- 
    ggplot(tmp.tb.glom.top, aes(fill = taxa, x = ID, y = value)) + 
    geom_bar(position = position_stack(reverse = F), stat = 'summary', fun = 'mean', width = .95) +
    theme_bw(base_family = "Liberation Sans") + 
    facet_grid(. ~ Group, scales = "free_x", space = "free_x") +
    scale_y_continuous(expand = expansion(mult = 0.0025, add = 0), 
                       labels = scales::percent_format(accuracy = 1)) +
    scale_x_discrete(expand = expansion(mult = 0, add = 0)) +
    scale_fill_manual(values = tmp.colors) +  
    labs(fill = tmp.rank) + xlab("") + ylab("Relative abundance (%)<br>") +
    my_theme + 
    theme(legend.position = "none",
          axis.ticks.length = unit(2.5, "pt"),
          axis.text.x = element_markdown(size = 7, colour = 'black', angle = 90, hjust = 0, vjust = .5))
  
  gg.bp <-
    gg.bp +
    plot_grid(
      ggpubr::get_legend(addSmallLegend(gg.bp + theme(legend.position = "right")) + 
                           theme(legend.justification = c(0, 0.5), 
                                 legend.box.margin = margin(0, 0, 0, 0),
                                 legend.box.background = element_rect(fill = "white", colour = "black", linewidth = .75))),
      ncol = 2, align = "hv") + plot_layout(widths = c(.5, 1.5))
  
  gg.bp.plots[[i]] <- gg.bp
  
}

gg.bp.p <- gg.bp.plots[[1]]; gg.bp.p

ggsave(filename = paste0(path, "figures/TaxProfileSample_1-Phylum.png"), plot = gg.bp.p, width = 5000*3.5, height = 5000 * 1.05, units = 'px', dpi = 1200)
ggsave(filename = paste0(path, "figures/TaxProfileSample_1-Phylum.svg"), plot = gg.bp.p, width = 5000*3.5, height = 5000 * 1.05, units = 'px', dpi = 1200)

# Taxabar (Group) ---------------------------------------------------------

gg.bp.plots <- list()

for (i in seq_along(ranks)) {
  
  tmp.rank <- ranks[[i]]
  tmp.tb.glom <- tb.taxabar[[i]]
  
  tmp.tb.glom <- tmp.tb.glom %>% 
    unite("taxa", 2:(i+1), sep = "; ") %>% 
    select(-1) %>% 
    mutate("countSum" = rowSums(.[2:ncol(.)]),
           "taxa" = make.unique(taxa),
           "taxa" = str_replace_all(taxa, pattern = "__", replacement = ": "),
           "taxa" = str_replace_all(taxa, pattern = "_", replacement = " ")) %>% 
    arrange(., desc(countSum))
  
  top_in <- 10 ; top_out = top_in + 1
  
  tmp.tb.glom.top <- rbind(tmp.tb.glom[1:top_in, c(-ncol(tmp.tb.glom))], 
                           c("Others", colSums(tmp.tb.glom[top_out:nrow(tmp.tb.glom), 2:c(ncol(tmp.tb.glom)-1)]))) %>% 
    filter(complete.cases(.)) %>% 
    mutate("taxa" = factor(taxa, levels = taxa)) %>% 
    mutate_at(2:ncol(.), as.numeric) %>% 
    merge(x = melt(.), y = tb.metad, by.x = "variable", by.y = "row.names")
  
  tmp.colors <- colorRampPalette(c(
    
    paletteer::paletteer_d("futurevisions::pegasi")
    
  ))(length(unique(tmp.tb.glom.top$taxa))-1)
  
  names(tmp.colors) <- unique(tmp.tb.glom.top$taxa)[1:top_in]
  
  tmp.colors["Others"] <- "#E6E6E6FF"
  
  gg.bp <- 
    ggplot(tmp.tb.glom.top, aes(fill = taxa, x = Group, y = value)) + 
    geom_bar(position = position_stack(reverse = F), stat = 'summary', fun = 'mean', width = .95) +
    theme_bw(base_family = "Liberation Sans") + 
    facet_grid(. ~ Group, scales = "free_x", space = "free_x") +
    scale_y_continuous(expand = expansion(mult = 0.0025, add = 0), 
                       labels = scales::percent_format(accuracy = 1)) +
    scale_x_discrete(expand = expansion(mult = 0, add = 0)) +
    scale_fill_manual(values = tmp.colors) +  
    labs(fill = tmp.rank) + xlab("") + ylab("Relative abundance (%)<br>") +
    my_theme + 
    theme(legend.position = "none",
          axis.ticks.length.x = unit(0, "pt"),
          axis.ticks.length.y = unit(2.5, "pt"),
          axis.text.x = element_blank())
  
  gg.bp <-
    gg.bp +
    plot_grid(
      ggpubr::get_legend(addSmallLegend(gg.bp + theme(legend.position = "right")) + 
                           theme(legend.justification = c(0, 0.5), 
                                 legend.box.margin = margin(0, 0, 0, 0),
                                 legend.box.background = element_rect(fill = "white", colour = "black", linewidth = .75))),
      ncol = 2, align = "hv") + plot_layout(widths = c(.25, 1.5))
  
  gg.bp.plots[[i]] <- gg.bp
  
}

gg.bp.p <- gg.bp.plots[[1]]; gg.bp.p

ggsave(filename = paste0(path, "figures/TaxProfileGroup_1-Phylum.png"), plot = gg.bp.p, width = 5000*3.5, height = 5000 * 1, units = 'px', dpi = 1200)
ggsave(filename = paste0(path, "figures/TaxProfileGroup_1-Phylum.svg"), plot = gg.bp.p, width = 5000*3.5, height = 5000 * 1, units = 'px', dpi = 1200)

# Taxabar (Both) ----------------------------------------------------------

gg.bp.plots <- list()

for (i in seq_along(ranks)) {
  
  tmp.rank <- ranks[[i]]
  tmp.tb.glom <- tb.taxabar[[i]]
  
  tmp.tb.glom <- tmp.tb.glom %>% 
    unite("taxa", 2:(i+1), sep = "; ") %>% 
    select(-1) %>% 
    mutate("countSum" = rowSums(.[2:ncol(.)]),
           "taxa" = make.unique(taxa),
           "taxa" = str_replace_all(taxa, pattern = "__", replacement = ": "),
           "taxa" = str_replace_all(taxa, pattern = "_", replacement = " ")) %>% 
    arrange(., desc(countSum))
  
  top_in <- 10 ; top_out = top_in + 1
  
  tmp.tb.glom.top <- rbind(tmp.tb.glom[1:top_in, c(-ncol(tmp.tb.glom))], 
                           c("Others", colSums(tmp.tb.glom[top_out:nrow(tmp.tb.glom), 2:c(ncol(tmp.tb.glom)-1)]))) %>% 
    filter(complete.cases(.)) %>% 
    mutate("taxa" = factor(taxa, levels = taxa)) %>% 
    mutate_at(2:ncol(.), as.numeric) %>% 
    merge(x = melt(.), y = tb.metad, by.x = "variable", by.y = "row.names")
  
  tmp.colors <- colorRampPalette(c(
    
    paletteer::paletteer_d("futurevisions::pegasi")
    
  ))(length(unique(tmp.tb.glom.top$taxa))-1)
  
  names(tmp.colors) <- unique(tmp.tb.glom.top$taxa)[1:top_in]
  
  tmp.colors["Others"] <- "#E6E6E6FF"
  
  gg.bp.1 <- 
    ggplot(tmp.tb.glom.top, aes(fill = taxa, x = ID, y = value)) + 
    geom_bar(position = position_stack(reverse = F), stat = 'summary', fun = 'mean', width = .95) +
    theme_bw(base_family = "Liberation Sans") + 
    facet_grid(. ~ Group, scales = "free_x", space = "free_x") +
    scale_y_continuous(expand = expansion(mult = 0.0025, add = 0), 
                       labels = scales::percent_format(accuracy = 1)) +
    scale_x_discrete(expand = expansion(mult = 0, add = 0)) +
    scale_fill_manual(values = tmp.colors) +  
    labs(fill = tmp.rank) + xlab("") + ylab("Relative abundance (%)<br>") +
    my_theme + 
    theme(legend.position = "none",
          axis.ticks.length.y = unit(2.5, "pt"),
          axis.text.x = element_markdown(size = 7, colour = 'black', angle = 90, hjust = 0, vjust = .5))
  
  
  gg.bp.2 <- 
    ggplot(tmp.tb.glom.top, aes(fill = taxa, x = Group, y = value)) + 
    geom_bar(position = position_stack(reverse = F), stat = 'summary', fun = 'mean', width = .95) +
    theme_bw(base_family = "Liberation Sans") + 
    facet_grid(. ~ Group, scales = "free_x", space = "free_x") +
    scale_y_continuous(expand = expansion(mult = 0.0025, add = 0), 
                       labels = scales::percent_format(accuracy = 1)) +
    scale_x_discrete(expand = expansion(mult = 0, add = 0)) +
    scale_fill_manual(values = tmp.colors) +  
    labs(fill = tmp.rank) + xlab("") + ylab("Relative abundance (%)<br>") +
    my_theme + 
    theme(legend.position = "none",
          axis.ticks.length.x = unit(0, "pt"),
          axis.ticks.length.y = unit(2.5, "pt"),
          axis.text.x = element_blank())
  
  gg.bp <-
    gg.bp.1 + gg.bp.2 +
    plot_grid(
      ggpubr::get_legend(addSmallLegend(gg.bp.2 + theme(legend.position = "right")) + 
                           theme(legend.justification = c(0, 0.5), 
                                 legend.box.margin = margin(0, 0, 0, 0),
                                 legend.box.background = element_rect(fill = "white", colour = "black", linewidth = .75))),
      ncol = 2, align = "hv") + plot_layout(widths = c(.5, .25, 1))
  
  gg.bp.plots[[i]] <- gg.bp
  
}

gg.bp.p <- gg.bp.plots[[1]]; gg.bp.p

ggsave(filename = paste0(path, "figures/TaxProfileBoth_1-Phylum.png"), plot = gg.bp.p, width = 5000*3.5, height = 5000 * 1, units = 'px', dpi = 1200)
ggsave(filename = paste0(path, "figures/TaxProfileBoth_1-Phylum.svg"), plot = gg.bp.p, width = 5000*3.5, height = 5000 * 1, units = 'px', dpi = 1200)

rm(list = ls(pattern = "gg.bp"))
rm(list = ls(pattern = "tb.glom"))
rm(tb.taxabar, i, ranks, tmp.colors, tmp.rank, top_in, top_out, fun_bestax)

# Heatmap -----------------------------------------------------------------

library("ComplexHeatmap")

# Top curves

get.topcurve <- function(my.ps.table, rank, my.top) {
  
  return(
    
    merge(data.frame(check.names = F, tax_table(tax_glom(my.ps.table, taxrank= rank))) %>% select(1:any_of(rank)),
          data.frame(check.names = F, otu_table(tax_glom(my.ps.table, taxrank= rank))), by = 0)[-1] %>%
      select(Phylum, any_of(rank), where(is.numeric)) %>% 
      mutate(RowSum = rowSums(select(., where(is.numeric)))) %>% 
      arrange(-RowSum) %>%
      mutate(RowSum = RowSum/sum(RowSum)) %>% 
      select(RowSum) %>% 
      cumsum() %>% 
      mutate(top = 1:nrow(.)) %>% 
      add_row(top = 0, .before = 1) %>% 
      mutate(across(everything(), ~ replace_na(., 0))) %>% 
      
      ggplot(., aes(x = top, y = RowSum)) +
      geom_line() +
      geom_vline(xintercept = my.top, color = 'red', lty = 'dashed', lwd = 0.5) +
      theme_bw(base_family = "Liberation Sans") + 
      scale_x_continuous(expand = c(0.01, 0.05, 0.01, 0.05)) +
      scale_y_continuous(expand = c(0.01, 0.001, 0.01, 0.05), labels = scales::percent_format(accuracy = 1)) +
      ylab('Cumulative abundance<br>') +
      xlab('<br>Top taxa (n)') +
      labs(title = paste0("Rank: ", rank, " (top ", my.top, ")")) +
      my_theme
    
  )
  
}

hp.topcurve.c <- get.topcurve(ps.rarefied, "Class", 20)
hp.topcurve.o <- get.topcurve(ps.rarefied, "Order", 30)
hp.topcurve.f <- get.topcurve(ps.rarefied, "Family", 50)

png(filename = paste0(path, "figures/TopCurve.png"), width = 5000 * 2.5, height = 5000 * 0.75, units = "px", res = 1200, type = "cairo")

hp.topcurve.c + hp.topcurve.o + hp.topcurve.f

dev.off()

# Top tables

fun_bestax <- function(df) {
  
  for (n_taxa in 1:nrow(df)) {
    for (n_rank in 1:7) {
      if (str_detect(df[[n_taxa, n_rank]], "Unclassified")) {
        
        df[[n_taxa, n_rank]] <- 
          str_replace(df[[n_taxa, n_rank]], "Unclassified", df[[n_taxa, (n_rank-1)]] %>% str_remove(., ".*__") %>% paste0(., " (Un.)")) %>% 
          str_replace(., " \\(Un\\.\\).*", " (Un.)")
        
      }
    }
  }
  return(df)
}

get.top.fun.names <- function(my.ps.table, rank, top) {
  
  return(
    
    merge(data.frame(check.names = F, tax_table(tax_glom(my.ps.table, taxrank= rank))) %>% select(1:any_of(rank)),
          data.frame(check.names = F, otu_table(tax_glom(my.ps.table, taxrank= rank))), by = 0)[-1] %>%
      fun_bestax(.) %>% 
      mutate(RowSum = rowSums(select(., where(is.numeric)))) %>% 
      arrange(-RowSum) %>% 
      mutate(!!rank := str_replace_all(get(rank), pattern = "__", replacement = ": "),
             !!rank := str_replace_all(get(rank), pattern = "_", replacement = " ")) %>% 
      head(top) %>% 
      select(-RowSum) %>% 
      pull(rank)
    
  )
  
}

get.top.fun.table <- function(my.ps.table, rank, top) {
  
  return(
    
    merge(data.frame(check.names = F, tax_table(tax_glom(my.ps.table, taxrank= rank))) %>% select(1:any_of(rank)),
          data.frame(check.names = F, otu_table(tax_glom(my.ps.table, taxrank= rank))), by = 0)[-1] %>%
      fun_bestax(.) %>% 
      select(Phylum, any_of(rank), where(is.numeric)) %>% 
      mutate(RowSum = rowSums(select(., where(is.numeric)))) %>% 
      arrange(-RowSum) %>% 
      mutate(Phylum = str_replace_all(Phylum, pattern = "__", replacement = ": "),
             Phylum = str_replace_all(Phylum, pattern = "_", replacement = " "),
             Taxa = str_replace_all(get(rank), pattern = "__", replacement = ": "),
             Taxa = str_replace_all(Taxa, pattern = "_", replacement = " ")) %>% 
      head(top) %>% 
      select(-RowSum) %>% 
      select(Phylum, Taxa, where(is.numeric)) %>% 
      mutate(mean = rowMeans(select(., where(is.numeric)))) 
    
  )
  
}

hp.top.names.p <- get.top.fun.names(ps.rarefied.tss, "Phylum", 10)

hp.top.tb.c <- get.top.fun.table(ps.rarefied.tss, "Class", 20)
hp.top.tb.o <- get.top.fun.table(ps.rarefied.tss, "Order", 30)
hp.top.tb.f <- get.top.fun.table(ps.rarefied.tss, "Family", 50)

hp.top.tb.c <- hp.top.tb.c %>% mutate(Phylum = ifelse(Phylum %in% hp.top.names.p, Phylum, "Others"))
hp.top.tb.o <- hp.top.tb.o %>% mutate(Phylum = ifelse(Phylum %in% hp.top.names.p, Phylum, "Others"))
hp.top.tb.f <- hp.top.tb.f %>% mutate(Phylum = ifelse(Phylum %in% hp.top.names.p, Phylum, "Others"))

# HPs

# Class

hp.top.tb <- hp.top.tb.c

# Anno. Col
hp.anno.col <- tb.metad %>% select("Treatment        " = Group)
hp.anno.col.col <- c("#525252FF", "#8cc084FF", "#38a3a5FF", "#0466c8FF")
names(hp.anno.col.col) <- unique(hp.anno.col$`Treatment        `)

# Anno. Row
hp.anno.row <- 
  hp.top.tb %>% 
  select(Taxa, Phylum, `Rel. abundance (%)        ` = mean) %>% 
  mutate(`Rel. abundance (%)        ` = `Rel. abundance (%)        ` * 100) %>% 
  mutate(Phylum = factor(Phylum, levels = c(hp.top.names.p, "Others"))) %>% 
  mutate(Taxa = paste0(Taxa, "        ")) %>% 
  column_to_rownames("Taxa")

hp.anno.row.col.p <- colorRampPalette(c(
  
  paletteer::paletteer_d("futurevisions::pegasi")
  
))(length(unique(hp.top.names.p)))

names(hp.anno.row.col.p) <- unique(hp.top.names.p)

hp.anno.row.col.p["Others"] <- "#E6E6E6FF"

hp.anno.col.col.ls <- 
  list(
    `Rel. abundance (%)        ` = c("white", "black"),
    Phylum = hp.anno.row.col.p,
    `Treatment        ` = hp.anno.col.col
  )

# Mat. prep.
hp.values <- 
  hp.top.tb %>% 
  select(Taxa, tb.metad$ID) %>% 
  mutate(Taxa = paste0(Taxa, "        ")) %>% 
  column_to_rownames("Taxa") %>%
  compositions::clr()

hp.breaks <- max(abs(min(hp.values, na.rm = T)), max(hp.values, na.rm = T))

# Plot
ComplexHeatmap::pheatmap(hp.values, name = "CLR norm. abundance",
                         scale = "none",
                         breaks = c(seq(-ceiling(hp.breaks), ceiling(hp.breaks), .01)),
                         cluster_row = T, 
                         cluster_cols = T, 
                         fontsize_col = 10, 
                         fontsize_row = 10,
                         border_color = NA,
                         annotation_names_col = F,
                         annotation_names_row = F,
                         annotation_col = hp.anno.col, 
                         annotation_row = hp.anno.row, 
                         annotation_colors = hp.anno.col.col.ls,
                         color = colorRampPalette(c("#364B9AFF", "#4A7BB7FF", "#6EA6CDFF", "#98CAE1FF", "#FFFFFFFF", "#FDB366FF", "#F67E4BFF", "#DD3D2DFF", "#A50026FF"))(length(seq(-ceiling(hp.breaks), ceiling(hp.breaks), .01))-1),
                         clustering_method = "ward.D2") -> hp.plot

hp.plot

png(filename = paste0(path, "figures/Heatmap_2-Class.png"), width = 5000 * 2.5, height = 5000 * 1.35, units = "px", res = 1200, type = "cairo")

hp.plot

dev.off()

# Order

hp.top.tb <- hp.top.tb.o

# Anno. Col
hp.anno.col <- tb.metad %>% select("Treatment        " = Group)
hp.anno.col.col <- c("#525252FF", "#8cc084FF", "#38a3a5FF", "#0466c8FF")
names(hp.anno.col.col) <- unique(hp.anno.col$`Treatment        `)

# Anno. Row
hp.anno.row <- 
  hp.top.tb %>% 
  select(Taxa, Phylum, `Rel. abundance (%)        ` = mean) %>% 
  mutate(`Rel. abundance (%)        ` = `Rel. abundance (%)        ` * 100) %>% 
  mutate(Phylum = factor(Phylum, levels = c(hp.top.names.p, "Others"))) %>% 
  mutate(Taxa = paste0(Taxa, "        ")) %>% 
  column_to_rownames("Taxa")

hp.anno.row.col.p <- colorRampPalette(c(
  
  paletteer::paletteer_d("futurevisions::pegasi")
  
))(length(unique(hp.top.names.p)))

names(hp.anno.row.col.p) <- unique(hp.top.names.p)

hp.anno.row.col.p["Others"] <- "#E6E6E6FF"

hp.anno.col.col.ls <- 
  list(
    `Rel. abundance (%)        ` = c("white", "black"),
    Phylum = hp.anno.row.col.p,
    `Treatment        ` = hp.anno.col.col
  )

# Mat. prep.
hp.values <- 
  hp.top.tb %>% 
  select(Taxa, tb.metad$ID) %>% 
  mutate(Taxa = paste0(Taxa, "        ")) %>% 
  column_to_rownames("Taxa") %>%
  compositions::clr()

hp.breaks <- max(abs(min(hp.values, na.rm = T)), max(hp.values, na.rm = T))

# Plot
ComplexHeatmap::pheatmap(hp.values, name = "CLR norm. abundance",
                         scale = "none",
                         breaks = c(seq(-ceiling(hp.breaks), ceiling(hp.breaks), .01)),
                         cluster_row = T, 
                         cluster_cols = T, 
                         fontsize_col = 10, 
                         fontsize_row = 10,
                         border_color = NA,
                         annotation_names_col = F,
                         annotation_names_row = F,
                         annotation_col = hp.anno.col, 
                         annotation_row = hp.anno.row, 
                         annotation_colors = hp.anno.col.col.ls,
                         color = colorRampPalette(c("#364B9AFF", "#4A7BB7FF", "#6EA6CDFF", "#98CAE1FF", "#FFFFFFFF", "#FDB366FF", "#F67E4BFF", "#DD3D2DFF", "#A50026FF"))(length(seq(-ceiling(hp.breaks), ceiling(hp.breaks), .01))-1),
                         clustering_method = "ward.D2") -> hp.plot

hp.plot

png(filename = paste0(path, "figures/Heatmap_3-Order.png"), width = 5000 * 2.5, height = 5000 * 1.75, units = "px", res = 1200, type = "cairo")

hp.plot

dev.off()

# Family

hp.top.tb <- hp.top.tb.f

# Anno. Col
hp.anno.col <- tb.metad %>% select("Treatment        " = Group)
hp.anno.col.col <- c("#525252FF", "#8cc084FF", "#38a3a5FF", "#0466c8FF")
names(hp.anno.col.col) <- unique(hp.anno.col$`Treatment        `)

# Anno. Row
hp.anno.row <- 
  hp.top.tb %>% 
  select(Taxa, Phylum, `Rel. abundance (%)        ` = mean) %>% 
  mutate(`Rel. abundance (%)        ` = `Rel. abundance (%)        ` * 100) %>% 
  mutate(Phylum = factor(Phylum, levels = c(hp.top.names.p, "Others"))) %>% 
  mutate(Taxa = paste0(Taxa, "        ")) %>% 
  column_to_rownames("Taxa")

hp.anno.row.col.p <- colorRampPalette(c(
  
  paletteer::paletteer_d("futurevisions::pegasi")
  
))(length(unique(hp.top.names.p)))

names(hp.anno.row.col.p) <- unique(hp.top.names.p)

hp.anno.row.col.p["Others"] <- "#E6E6E6FF"

hp.anno.col.col.ls <- 
  list(
    `Rel. abundance (%)        ` = c("white", "black"),
    Phylum = hp.anno.row.col.p,
    `Treatment        ` = hp.anno.col.col
  )

# Mat. prep.
hp.values <- 
  hp.top.tb %>% 
  select(Taxa, tb.metad$ID) %>% 
  mutate(Taxa = paste0(Taxa, "        ")) %>% 
  column_to_rownames("Taxa") %>%
  compositions::clr()

hp.breaks <- max(abs(min(hp.values, na.rm = T)), max(hp.values, na.rm = T))

# Plot
ComplexHeatmap::pheatmap(hp.values, name = "CLR norm. abundance",
                         scale = "none",
                         breaks = c(seq(-ceiling(hp.breaks), ceiling(hp.breaks), .01)),
                         cluster_row = T, 
                         cluster_cols = T, 
                         fontsize_col = 10, 
                         fontsize_row = 10,
                         border_color = NA,
                         annotation_names_col = F,
                         annotation_names_row = F,
                         annotation_col = hp.anno.col, 
                         annotation_row = hp.anno.row, 
                         annotation_colors = hp.anno.col.col.ls,
                         color = colorRampPalette(c("#364B9AFF", "#4A7BB7FF", "#6EA6CDFF", "#98CAE1FF", "#FFFFFFFF", "#FDB366FF", "#F67E4BFF", "#DD3D2DFF", "#A50026FF"))(length(seq(-ceiling(hp.breaks), ceiling(hp.breaks), .01))-1),
                         clustering_method = "ward.D2") -> hp.plot

hp.plot

png(filename = paste0(path, "figures/Heatmap_4-Family.png"), width = 5000 * 2.5, height = 5000 * 2.65, units = "px", res = 1200, type = "cairo")

hp.plot

dev.off()

rm(list = ls(pattern = "hp."))
rm(get.top.fun.names, get.top.fun.table, fun_bestax, get.topcurve)

# Venn Diagram ------------------------------------------------------------

fun.venn <- function(my.ps.table, my.title) {
  
  ps.venn <- merge_samples(my.ps.table, "Group", fun=sum)
  tb.venn <- data.frame(check.names = F, t(ps.venn@otu_table)) %>% rownames_to_column(var = "ID")
  names(tb.venn)
  
  ls.venn <- list(
    `CTRL` = tb.venn %>% filter(`CTRL` != 0) %>% pull(ID),
    `SM05` = tb.venn %>% filter(`SM05` != 0) %>% pull(ID),
    `SM10` = tb.venn %>% filter(`SM10` != 0) %>% pull(ID),
    `SM20` = tb.venn %>% filter(`SM20` != 0) %>% pull(ID)
  )
  
  gg.venn <- 
    ggvenn(ls.venn, 
           fill_color = unname(tr.colors), 
           stroke_color = NA, 
           stroke_size = 0.5, 
           set_name_size = 0, 
           fill_alpha = .5, 
           text_size = 2.75,  
           show_percentage = F) + 
    labs(title = my.title) + 
    theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5)); gg.venn
  
  return(gg.venn)
  
}

gg.venn.a <- fun.venn(ps.rarefied, "ASVs")
gg.venn.p <- fun.venn(tax_glom(ps.rarefied, taxrank="Phylum"), "Phylum")
gg.venn.c <- fun.venn(tax_glom(ps.rarefied, taxrank="Class"), "Class")
gg.venn.o <- fun.venn(tax_glom(ps.rarefied, taxrank="Order"), "Order")
gg.venn.f <- fun.venn(tax_glom(ps.rarefied, taxrank="Family"), "Family")
gg.venn.g <- fun.venn(tax_glom(ps.rarefied, taxrank="Genus"), "Genus")

gg.venn <- gg.venn.a + gg.venn.p + gg.venn.c + gg.venn.o + gg.venn.f + gg.venn.g + plot_layout(nrow = 1); gg.venn

ggsave(filename = paste0(path, "figures/VennDiagram.png"), plot = gg.venn, bg = "white", dpi = 1200, width = 5000*4.25, height = 5000*1.25, units = "px")
ggsave(filename = paste0(path, "figures/VennDiagram.svg"), plot = gg.venn, bg = "white", dpi = 1200, width = 5000*4.25, height = 5000*1.25, units = "px")

rm(list = ls(pattern = "venn"))

# Alpha diversity ---------------------------------------------------------

adiv.table <- 
  microbiome::alpha(ps.rarefied) %>% 
  select(observed, diversity_shannon, diversity_gini_simpson)
names(adiv.table) <- c("Richness", "Shannon's diversity", "Gini-Simpson")

adiv.table <- merge(sample_data(ps.rarefied), adiv.table, by = "row.names", sort = F)[-1]

adiv.table.longer <- 
  adiv.table %>% 
  mutate(across(where(is.character), ~factor(., unique(.)))) %>% 
  pivot_longer(cols = Richness:`Gini-Simpson`, names_to = "index", values_to = "values") %>% 
  data.frame(check.names = F) %>% 
  mutate(across(where(is.character), ~factor(., unique(.))))

adiv.table.summary <-
  adiv.table.longer %>% 
  group_by(Group, FullGroup, index) %>% 
  summarise(temp = list(c(psych::describe(values), quantile(values, na.rm = T), IQR = iqr(values, na.rm = T)))) %>% 
  unnest_wider(temp) %>% 
  select(-c(vars, trimmed, mad)) %>% 
  mutate(across(where(is.numeric), ~case_when(is.nan(.) ~ NA_integer_, is.infinite(.) ~ NA_integer_,.default = .))) %>% 
  data.frame(check.names = F)

get_kruskal <- function(tb) {
  
  tb %>%
    group_by(index) %>%
    nest() %>%
    mutate(kruskal = map(.x = data, ~kruskal.test(values ~ Group, data = .x) %>% broom::tidy())) %>%
    unnest(kruskal) %>%
    mutate(across(where(is.numeric), ~round(., 4))) -> p
  
  return(p)
  
}

adiv.kruskal <- get_kruskal(adiv.table.longer)

adiv.kruskal.2 <-
  adiv.kruskal %>% 
  mutate(p.full = paste0("Kruskal-Wallis, *p*: ", ifelse(p.value < 0.001, "<0.001", sprintf("%1.3f", round(p.value, 3)))))

adiv.get.posthoc.lsd <- function(tb, tb.summary, tb.stat, ls.invert = NULL) {
  
  tb.compm <- data.frame(index = "", Group = "", group = "")[-1,]
  
  for (i in unique(tb$index)) {
    
    tmp.tb <- tb %>% filter(index == i)
    
    if (i %in% ls.invert) {
      
      tmp.tb$values <- -1 * tmp.tb$values
      
    }
    
    tmp.lsd <- with(tmp.tb, agricolae::kruskal(tmp.tb$values, tmp.tb$Group, alpha = 0.05))
    tmp.lsd <-
      cbind(Group = row.names(tmp.lsd[["groups"]]), tmp.lsd[["groups"]]) %>% 
      remove_rownames() %>% 
      mutate(index = i, group = str_trim(groups)) %>% 
      select(index, Group, group)
    
    if (tb.stat %>% filter(index == i) %>% pull(p.value) > 0.05) {
      
      tmp.lsd <- mutate(tmp.lsd, group = NA_character_)
      
    }
    
    tb.compm <- bind_rows(tb.compm, tmp.lsd)
    
  }
  
  tb.compm <- merge(tb.summary, tb.compm, sort = F, all.x = T)
  
  return(tb.compm)
  
} 

adiv.compm.lsd <- adiv.get.posthoc.lsd(adiv.table.longer, adiv.table.summary, adiv.kruskal) %>% arrange(index)

adiv.gg <- 
  ggplot() +
  facet_wrap( ~ index, scales = "free", nrow = 1, ncol = 5) +
  geom_boxplot(data = adiv.table.longer, aes(x = FullGroup, y = values, fill = FullGroup, color = FullGroup), 
               width = .5, outliers = F, show.legend = F, alpha = .5) + 
  geom_point(data = adiv.table.longer, aes(x = FullGroup, y = values, color = FullGroup),
             position = position_jitterdodge(dodge.width=0.9), size = 1.75, shape = 16, alpha = .75) +
  geom_richtext(mapping = aes(label = p.full, x = 2.5, y = Inf, vjust = 2.5),
                data = adiv.kruskal.2, position = position_dodge(width = 0),
                size = 2.5, hjust = 0.5, color = "black", inherit.aes = T, 
                label.padding = unit(0, "pt"), fill = NA, label.color = NA) +
  geom_richtext(mapping = aes(label = group, x = as.numeric(FullGroup), y = max), 
                data = adiv.compm.lsd,
                size = 2.5, vjust = -1, hjust = 0.5, color = "black", inherit.aes = T,
                label.padding = unit(0, "pt"), fill = NA, label.color = NA) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.5)), breaks = scales::extended_breaks(n = 6)) +
  scale_x_discrete(expand = expansion(mult = c(0.25, 0.25))) +
  scale_fill_manual(values = tr.colors, name = NULL) +
  scale_color_manual(values = tr.colors, name = NULL) +
  ylab('') + xlab('') +
  guides(color = 
           guide_legend(order = 1, 
                        override.aes = list(shape = 22, 
                                            fill = tr.colors, 
                                            color = "black", 
                                            size = 3,
                                            linetype = 0),
                        ncol = 2)) +
  my_theme + 
  theme(legend.position = "bottom", 
        legend.box.background = element_rect(fill = "white", colour = "black", linewidth = .75),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

adiv.gg

ggsave(filename = paste0(path, "figures/AlphaDiv.png"), plot = adiv.gg, dpi = 1200, width = 5000*1.75, height = 5000*.85, units = "px")
ggsave(filename = paste0(path, "figures/AlphaDiv.svg"), plot = adiv.gg, dpi = 1200, width = 5000*1.75, height = 5000*.85, units = "px")

ls.adiv <- 
  
  list(
    
    "Alpha (Measures)" = 
      adiv.table %>% 
      select(ID, Group, Richness:`Gini-Simpson`) %>%
      mutate(across(where(is.factor), ~as.character(.))),
    
    "Alpha (Kruskal-Wallis)" = 
      adiv.kruskal.2 %>%
      select(Index = index, statistic, p = p.value, method) %>% 
      mutate(across(where(is.factor), ~as.character(.))),
    
    "Alpha (Summary)" =
      adiv.compm.lsd %>% 
      select(Group, FullGroup, Index = index, everything()) %>% 
      mutate(across(where(is.factor), ~as.character(.))) %>% 
      mutate(across(where(is.character), ~str_remove_all(., "<.?sup>"))) %>% 
      mutate(across(where(is.character), ~str_remove_all(., "\\*")))
    
  )

openxlsx::write.xlsx(x = ls.adiv, file = paste0(path, "tables/AlphaDiv.xlsx"), rowNames = F)

rm(list = ls(pattern = "adiv"))

# Beta diversity ----------------------------------------------------------

bdiv.indexes <- c("bray", "jaccard", "wunifrac", "unifrac")
names(bdiv.indexes) <- c("Bray-Curtis", "Jaccard", "UniFrac (Weighted)", "UniFrac (Unweighted)")

bdiv.plots <- list()

for (i in seq_along(bdiv.indexes)) {
  
  print(names(bdiv.indexes)[i])
  
  # Calculo das distancias
  bdiv.dist <- phyloseq::distance(ps.rarefied.tss, bdiv.indexes[[i]])
  
  bdiv.hc <- hclust(bdiv.dist, "complete")
  bdiv.dd <- as.dendrogram(bdiv.hc)
  bdiv.dd <- dendro_data(bdiv.dd, type = "rectangle")
  
  bdiv.dd.gg <- 
    ggplot() + 
    geom_segment(data = bdiv.dd$segments, aes(x = x, y = y, xend = xend, yend = yend)) + 
    geom_point(data = merge(bdiv.dd$labels, tb.metad, by.x = "label", by.y = "ID"),
               aes(y = y, x = x, fill = FullGroup), color = "black", shape = 21, size = 3) + 
    geom_text(data = bdiv.dd$labels,
              aes(y = y, x = x, label = label), color = "black", size = 2.5, angle = 90, hjust = 1.25) + 
    scale_y_continuous(expand = expansion(mult = c(0.35, 0.1)), breaks = seq(0, max(bdiv.dd$segments$yend), .1)) +
    scale_fill_manual(values = tr.colors) + 
    ylab(paste(names(bdiv.indexes)[i], "distance")) + xlab("") +
    my_theme + 
    theme(legend.position = "none", 
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
  bdiv.dd.gg
  
  # Análise PERMANOVA do indice de Bray-Curtis 
  
  bdiv.pval <- adonis2(bdiv.dist ~ Group, data = sample_data(ps.rarefied.tss) %>% data.frame(), method = bdiv.indexes[[i]], permutations = 9999)
  
  # Há diferença entre os grupos?
  bdiv.perma.tb <- 
    bdiv.pval %>% 
    as.data.frame() %>% 
    rownames_to_column() %>% 
    dplyr::rename("Source" = 1, "DF" = 2, "SS" = 3, "F" = 5, "P" = 6) %>% 
    mutate(Source = str_replace(Source, "Group", "Treatment"),
           Source = str_to_title(Source),
           sig = case_when(P <= 0.001 ~ "***", P <= 0.01 ~ "**", P <= 0.05 ~ "*", P <= 0.1 ~ ".", .default = " "))
  
  bdiv.posthoc.c1 <- pairwiseAdonis::pairwise.adonis(bdiv.dist, sample_data(ps.rarefied.tss)$Group, p.adjust.m = "holm", perm = 9999) %>%
    dplyr::rename("Pairs" = 1, "DF" = 2, "SS" = 3, "F" = 4, "P" = 6, "P(adj)" = 7) %>%  select(1, 2, 3, 5, 4, 6, 7, 8) %>%
    mutate(sig = case_when(`P(adj)` <= 0.001 ~ "***", `P(adj)` <= 0.01 ~ "**", `P(adj)` <= 0.05 ~ "*", `P(adj)` <= 0.1 ~ ".", .default = " "))
  
  # Ordenação 
  bdiv.ord <- ordinate(ps.rarefied.tss, "PCoA", distance = bdiv.dist)
  
  # Gráfico
  
  bdiv.pcoa <- 
    plot_ordination(ps.rarefied.tss, bdiv.ord, axes = c(1, 2)) +
    geom_vline(xintercept = 0, size = 0.25, linetype = 2) +
    geom_hline(yintercept = 0, size = 0.25, linetype = 2) +
    stat_ellipse(geom = "polygon", 
                 aes(
                   fill = sample_data(ps.rarefied.tss)$FullGroup,
                   color = sample_data(ps.rarefied.tss)$FullGroup),
                 linetype = 1, alpha = 0, show.legend = F) +
    geom_point(size = 2, shape = 16, aes(color = sample_data(ps.rarefied.tss)$FullGroup), inherit.aes = T) +
    scale_fill_manual(values = tr.colors) +
    scale_color_manual(values = tr.colors) +
    xlab(paste0("<br>PCoA 1 (",round(bdiv.ord[["values"]][["Relative_eig"]][1]*100, 2),"%)")) +
    ylab(paste0("PCoA 2 (",round(bdiv.ord[["values"]][["Relative_eig"]][2]*100, 2),"%)<br>")) +
    labs(title = names(bdiv.indexes)[i],
         subtitle = paste0("PERMANOVA, *p*: ", ifelse(bdiv.perma.tb[1, "P"] < 0.001, "<0.001", sprintf("%1.3f", round(bdiv.perma.tb[1, "P"], 3))))) +
    my_theme + 
    theme(legend.position = "none",
          plot.title = element_markdown(size = 10, face = 'bold', colour = 'black'),
          plot.subtitle = element_markdown(size = 7, colour = 'black'))
  
  bdiv.pcoa
  
  bdiv.disp <- 
    betadisper(d = bdiv.dist, group = tb.metad$FullGroup, type = "centroid")$distances %>% 
    data.frame() %>% rownames_to_column("ID") %>% dplyr::rename("values" = ".") %>% 
    merge(sample_data(ps.rarefied.tss) %>% data.frame(), .)
  
  bdiv.disp.summary <-
    bdiv.disp %>% 
    group_by(FullGroup) %>% 
    summarise(temp = list(c(psych::describe(values), quantile(values, na.rm = T), IQR = iqr(values, na.rm = T)))) %>% 
    unnest_wider(temp) %>% 
    select(-c(vars, trimmed, mad)) %>% 
    mutate(across(where(is.numeric), ~case_when(is.nan(.) ~ NA_integer_, is.infinite(.) ~ NA_integer_,.default = .))) %>% 
    data.frame(check.names = F)
  
  get_kruskal <- function(tb) {
    
    tb %>%
      nest() %>%
      mutate(kruskal = map(.x = data, ~kruskal.test(values ~ Group, data = .x) %>% broom::tidy())) %>%
      unnest(kruskal) %>%
      mutate(across(where(is.numeric), ~round(., 4))) %>% 
      select(-data) %>% 
      data.frame(check.names = F) -> p
    
    return(p)
    
  }
  
  bdiv.disp.kruskal <- get_kruskal(bdiv.disp)
  
  bdiv.disp.kruskal.2 <-
    bdiv.disp.kruskal %>% 
    mutate(p.full = paste0("Kruskal-Wallis, *p*: ", ifelse(p.value < 0.001, "<0.001", sprintf("%1.3f", round(p.value, 3))))) %>% 
    mutate(Index = "Betadisp") %>% 
    select(Index, everything())
  
  bdiv.disp.get.posthoc.lsd <- function(tb, tb.summary, tb.stat) {
    
    tb.compm <- data.frame(FullGroup = "", group = "")[-1,]
    
    tmp.tb <- tb 
    
    tmp.lsd <- with(tmp.tb, agricolae::kruskal(tmp.tb$values, tmp.tb$FullGroup, alpha = 0.05))
    tmp.lsd <-
      cbind(FullGroup = row.names(tmp.lsd[["groups"]]), tmp.lsd[["groups"]]) %>% 
      remove_rownames() %>% 
      mutate(group = str_trim(groups)) %>% 
      select(FullGroup, group)
    
    if (tb.stat %>% pull(p.value) > 0.05) {
      
      tmp.lsd <- mutate(tmp.lsd, group = NA_character_)
      
    }
    
    tb.compm <- bind_rows(tb.compm, tmp.lsd)
    
    tb.compm <- merge(tb.summary, tb.compm, sort = F, all.x = T)
    
    return(tb.compm)
    
  } 
  
  bdiv.disp.compm.lsd <- bdiv.disp.get.posthoc.lsd(bdiv.disp, bdiv.disp.summary, bdiv.disp.kruskal)
  
  bdiv.disp.gg <- 
    ggplot() +
    geom_boxplot(data = bdiv.disp, aes(x = FullGroup, y = values, fill = FullGroup, color = FullGroup), 
                 width = .5, outliers = F, show.legend = F, alpha = .5) + 
    geom_point(data = bdiv.disp, aes(x = FullGroup, y = values, color = FullGroup),
               position = position_jitterdodge(dodge.width=0.9), size = 1.75, shape = 16, alpha = .75) +
    geom_richtext(mapping = aes(label = p.full, x = 2.5, y = Inf, vjust = 2.5),
                  data = bdiv.disp.kruskal.2, position = position_dodge(width = 0),
                  size = 2.5, hjust = 0.5, color = "black", inherit.aes = T, 
                  label.padding = unit(0, "pt"), fill = NA, label.color = NA) +
    geom_richtext(mapping = aes(label = group, x = as.numeric(FullGroup), y = max), 
                  data = bdiv.disp.compm.lsd,
                  size = 2.5, vjust = -1, hjust = 0.5, color = "black", inherit.aes = T,
                  label.padding = unit(0, "pt"), fill = NA, label.color = NA) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.5)), breaks = scales::extended_breaks(n = 6)) +
    scale_x_discrete(expand = expansion(mult = c(0.25, 0.25))) +
    scale_fill_manual(values = tr.colors, name = NULL) +
    scale_color_manual(values = tr.colors, name = NULL) +
    ylab('') + xlab('') + labs(subtitle = "Beta dispersion") +
    guides(color = 
             guide_legend(order = 1, 
                          override.aes = list(shape = 22, 
                                              fill = tr.colors, 
                                              color = "black", 
                                              size = 3,
                                              linetype = 0),
                          ncol = 2)) +
    my_theme + 
    theme(legend.position = "none", 
          legend.box.background = element_rect(fill = "white", colour = "black", linewidth = .75),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(), 
          plot.subtitle = element_markdown(face = "bold", color = "black", hjust = 0.5, size = 10),
          plot.margin = unit(c(0.5, 0.75, 0.5, 0), 'cm'))
  
  bdiv.plot <- bdiv.pcoa + bdiv.disp.gg + plot_layout(widths = c(1.15, 1))
  bdiv.plot
  
  ggsave(filename = paste0(path, "figures/BetaDiv_", str_to_title(bdiv.indexes[[i]]), ".png"), plot = bdiv.plot, dpi = 1200, width = 5000 * 1.75, height = 5000 * 1, units = "px")
  ggsave(filename = paste0(path, "figures/BetaDiv_", str_to_title(bdiv.indexes[[i]]), ".svg"), plot = bdiv.plot, dpi = 1200, width = 5000 * 1.75, height = 5000 * 1, units = "px")
  
  ls.bdiv <- 
    
    list(
      "Beta PERMANOVA" = bdiv.perma.tb %>% 
        mutate(across(where(is.character), ~str_remove_all(., "\\*"))),
      
      "Beta post-hoc" = bdiv.posthoc.c1,
      
      "Beta Dispersion (Measures)" = 
        bdiv.disp %>% 
        select(ID, Group = Group, Betadisp = values) %>%
        mutate(across(where(is.factor), ~as.character(.))) %>% 
        mutate(across(where(is.character), ~str_remove_all(., "<.?sup>"))),
      
      "Beta Dispersion (Kruskal)" = 
        bdiv.disp.kruskal.2 %>%
        select(Index, statistic, p = p.value, method) %>% 
        mutate(across(where(is.factor), ~as.character(.))) %>% 
        mutate(across(where(is.character), ~str_remove_all(., "\\*"))),
      
      "Beta Dispersion (Summary)" =
        bdiv.disp.compm.lsd %>% 
        dplyr::rename(Group = FullGroup) %>% 
        mutate(across(where(is.factor), ~as.character(.))) %>% 
        mutate(across(where(is.character), ~str_remove_all(., "<.?sup>"))) %>% 
        mutate(across(where(is.character), ~str_remove_all(., ":.*")))
      
    )
  
  openxlsx::write.xlsx(x = ls.bdiv, file = paste0(path, "tables/BetaDiv_", str_to_title(bdiv.indexes[[i]]), ".xlsx"), rowNames = F)
  
}

rm(list = ls(pattern = "bdiv"))
rm(i)

# DAA ---------------------------------------------------------------------

# TSS table

daa.glom.gp.fun <- function(my.ps.table, group) {
  
  return(
    
    rbind(
      
      merge(data.frame(check.names = F, tax_table(phyloseq_standardize_otu_abundance(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Kingdom"))))["Kingdom"],
            data.frame(check.names = F, t(otu_table(phyloseq_standardize_otu_abundance(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Kingdom"))))), by = 0)[-1] %>% 
        dplyr::mutate(taxa = Kingdom, rank = "Kingdom") %>% select(taxa, rank, unique(tb.metad[[group]])),
      
      merge(data.frame(check.names = F, tax_table(phyloseq_standardize_otu_abundance(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Phylum"))))["Phylum"],
            data.frame(check.names = F, t(otu_table(phyloseq_standardize_otu_abundance(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Phylum"))))), by = 0)[-1] %>% 
        dplyr::mutate(taxa = Phylum, rank = "Phylum") %>% select(taxa, rank, unique(tb.metad[[group]])),
      
      merge(data.frame(check.names = F, tax_table(phyloseq_standardize_otu_abundance(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Class"))))["Class"],
            data.frame(check.names = F, t(otu_table(phyloseq_standardize_otu_abundance(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Class"))))), by = 0)[-1] %>% 
        dplyr::mutate(taxa = Class, rank = "Class") %>% select(taxa, rank, unique(tb.metad[[group]])),
      
      merge(data.frame(check.names = F, tax_table(phyloseq_standardize_otu_abundance(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Order"))))["Order"],
            data.frame(check.names = F, t(otu_table(phyloseq_standardize_otu_abundance(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Order"))))), by = 0)[-1] %>% 
        dplyr::mutate(taxa = Order, rank = "Order") %>% select(taxa, rank, unique(tb.metad[[group]])),
      
      merge(data.frame(check.names = F, tax_table(phyloseq_standardize_otu_abundance(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Family"))))["Family"],
            data.frame(check.names = F, t(otu_table(phyloseq_standardize_otu_abundance(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Family"))))), by = 0)[-1] %>% 
        dplyr::mutate(taxa = Family, rank = "Family") %>% select(taxa, rank, unique(tb.metad[[group]])),
      
      merge(data.frame(check.names = F, tax_table(phyloseq_standardize_otu_abundance(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Genus"))))["Genus"],
            data.frame(check.names = F, t(otu_table(phyloseq_standardize_otu_abundance(tax_glom(merge_samples(my.ps.table, group, fun=sum), taxrank="Genus"))))), by = 0)[-1] %>% 
        dplyr::mutate(taxa = Genus, rank = "Genus") %>% select(taxa, rank, unique(tb.metad[[group]]))
      
    ) %>% 
      
      filter(!str_detect(taxa, "Unclassified"))
    
  )
  
}


daa.glom0 <- daa.glom.gp.fun(ps.rarefied, "ID")
daa.glom1 <- daa.glom.gp.fun(ps.rarefied, "Group")

daa.glom.tss <- 
  daa.glom0 %>% pivot_longer(cols = tb.metad$ID, names_to = "ID", values_to = "rel_abund") %>% 
  merge(tb.metad, ., by = "ID") %>% 
  mutate(taxa = factor(taxa, levels = daa.glom0$taxa %>% unique())) %>% 
  mutate(rank = factor(rank, levels = daa.glom0$rank %>% unique())) %>% 
  arrange(taxa, Group, ID)

daa.glom.tss.summary <-
  daa.glom.tss %>% 
  select(Group, taxa, rank, rel_abund) %>% 
  group_by(Group, taxa, rank) %>% 
  mutate(mean = mean(rel_abund), sd = sd(rel_abund)) %>% 
  ungroup() %>% 
  select(-rel_abund) %>% 
  distinct() %>% 
  mutate()

# -------------------------------------------------------------------------

daa.glom.tss %>%
  nest(data = -c(taxa, rank)) %>%
  mutate(kruskal = map(.x=data, ~kruskal.test(rel_abund ~ Group, data = .x) %>% broom::tidy())) %>%
  unnest(c(kruskal)) %>% 
  select(-data) %>% 
  data.frame() -> tb.kruskal.tmp

bind_rows(tb.kruskal.tmp %>% select(taxa, rank) %>% mutate(group1 = "CTRL", group2 = "SM05"),
          tb.kruskal.tmp %>% select(taxa, rank) %>% mutate(group1 = "CTRL", group2 = "SM10"),
          tb.kruskal.tmp %>% select(taxa, rank) %>% mutate(group1 = "CTRL", group2 = "SM20")) %>% 
  arrange(taxa) %>% 
  data.frame() -> tb.wilcox.frame

daa.glom.tss %>%
  nest(data = -c(taxa, rank)) %>%
  mutate(pwwilcox = map(.x=data, ~pairwise.wilcox.test(.x$rel_abund, .x$Group) %>% broom::tidy())) %>%
  unnest(c(pwwilcox)) %>% 
  select(taxa, rank, group0 = group2, group2 = group1, p = p.value) %>% dplyr::rename(group1 = group0) %>% 
  mutate(p.signif = case_when(p < 0.0001 ~ "****", p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", p >= 0.05 ~ NA_character_)) %>% 
  mutate(method = "Wilcoxon") %>% 
  data.frame() -> tb.wilcox.tmp

base::merge(tb.wilcox.frame, tb.wilcox.tmp, by = c("taxa", "rank", "group1", "group2"), all.x = T, sort = T) %>% 
  mutate(p.signif = str_replace_na(replacement = "ns", p.signif)) %>% 
  mutate(method = str_replace_na(replacement = "Wilcoxon", method)) -> tb.wilcox.tmp

tb.kruskal.tmp
tb.wilcox.tmp

for (i.taxa in tb.kruskal.tmp$taxa) {
  
  if ((tb.kruskal.tmp %>% filter(taxa == i.taxa) %>% pull(p.value))>= 0.05 | 
      is.na(tb.kruskal.tmp %>% filter(taxa == i.taxa) %>% pull(p.value))) {
    
    tb.wilcox.tmp[which(tb.wilcox.tmp$taxa == i.taxa), "p"] <- NA_real_
    tb.wilcox.tmp[which(tb.wilcox.tmp$taxa == i.taxa), "p.signif"] <- NA_character_
    
  }
  
}

tb.wilcox.tmp %>% 
  merge(., table(tb.metad$Group) %>% data.frame(), by.x = "group1", by.y = "Var1") %>% 
  merge(., table(tb.metad$Group) %>% data.frame(), by.x = "group2", by.y = "Var1") %>% 
  select(taxa, rank, group1, group2, n1 = Freq.x, n2 = Freq.y, everything()) %>% 
  mutate(group1 = factor(group1, levels = levels(tb.metad$Group))) %>% 
  mutate(group2 = factor(group2, levels = levels(tb.metad$Group))) %>% 
  arrange(taxa, group1) -> tb.wilcox.tmp

for (i.n in 1:nrow(tb.wilcox.tmp)) {
  
  tb.wilcox.tmp[i.n, "mean1"] <- daa.glom.tss.summary %>% filter(taxa == tb.wilcox.tmp[i.n, "taxa"] & Group == tb.wilcox.tmp[i.n, "group1"]) %>% pull(mean)
  tb.wilcox.tmp[i.n, "sd1"] <- daa.glom.tss.summary %>% filter(taxa == tb.wilcox.tmp[i.n, "taxa"] & Group == tb.wilcox.tmp[i.n, "group1"]) %>% pull(sd)
  tb.wilcox.tmp[i.n, "mean2"] <- daa.glom.tss.summary %>% filter(taxa == tb.wilcox.tmp[i.n, "taxa"] & Group == tb.wilcox.tmp[i.n, "group2"]) %>% pull(mean)
  tb.wilcox.tmp[i.n, "sd2"] <- daa.glom.tss.summary %>% filter(taxa == tb.wilcox.tmp[i.n, "taxa"] & Group == tb.wilcox.tmp[i.n, "group2"]) %>% pull(sd)
  
  tb.wilcox.tmp[i.n, "difference"] <- tb.wilcox.tmp[i.n, "mean1"] - tb.wilcox.tmp[i.n, "mean2"]
  tb.wilcox.tmp[i.n, "SE"] <- sqrt((tb.wilcox.tmp[i.n, "sd1"]^2 / tb.wilcox.tmp[i.n, "n1"]) + (tb.wilcox.tmp[i.n, "sd2"]^2 / tb.wilcox.tmp[i.n, "n2"]))
  
  tb.wilcox.tmp[i.n, "CI_lower"] <- tb.wilcox.tmp[i.n, "difference"] - (1.96 * tb.wilcox.tmp[i.n, "SE"])
  tb.wilcox.tmp[i.n, "CI_upper"] <- tb.wilcox.tmp[i.n, "difference"] + (1.96 * tb.wilcox.tmp[i.n, "SE"])
  
}

tb.kruskal <- tb.kruskal.tmp
tb.wilcox <- tb.wilcox.tmp

# STAMP -------------------------------------------------------------------

tr.colors <- c("#525252FF", "#8cc084FF", "#38a3a5FF", "#0466c8FF")
names(tr.colors) <- unique(tb.metad$Group)

gg.stamp.fun <- function(my.tb.wilcox, g1, g2) {
  
  daa.cont <- 
    my.tb.wilcox %>% 
    filter(group1 == g1 & group2 == g2 & rank != "Kingdom" & !is.na(p) & !is.na(p.signif)) %>% 
    filter(p < 0.05)
  
  daa.cont.bp <- 
    daa.cont %>% 
    mutate(taxa = str_replace(taxa, "__", ": "),
           taxa = str_replace_all(taxa, "_", " "),
           taxa = fct_reorder(taxa, difference)) %>% 
    arrange(rank, taxa, group1) %>% 
    select(taxa:group2, p.signif, starts_with("mean"), starts_with("sd"), difference, starts_with("CI"), p) %>% 
    pivot_longer(., cols = !c(taxa, rank, difference, starts_with("CI"), p), names_pattern = "(.*)(1|2)$", names_to = c(".value", "names")) %>% 
    filter(!is.na(group))
  
  daa.ranks <- unique(daa.cont.bp$rank)
  
  daa.bp.ls <- list()
  daa.dp.ls <- list()
  daa.tp.ls <- list()
  daa.hg.ls <- list()
  
  for (i in seq_along(daa.ranks)) {
    
    daa.cont.bp %>% 
      filter(rank == daa.ranks[i]) %>% 
      select(taxa) %>% 
      distinct() %>% 
      nrow() -> gg.da.hg
    
    daa.cont.bp %>% 
      filter(rank == daa.ranks[i]) %>% 
      ggplot(., aes(x = mean, y = taxa, fill = group)) + 
      geom_stripes(odd = "#88888833", even = "#FFFFFF00") +
      geom_errorbar(width = .25, position=position_dodge(.75), stat='identity', aes(xmin = 0, xmax = mean + sd)) +
      geom_bar(width = .75, position='dodge', stat='identity', linewidth = .25) +
      theme_forest() +
      scale_fill_manual(values = tr.colors) +
      scale_x_continuous(labels = scales::percent_format(), expand = expansion(mult = c(0, 0.1)), trans='sqrt') +
      labs(y = paste0(daa.ranks[i], "\n"), x = "Mean proportion\n(%, square root transformed)") +
      theme(legend.position="none",
            axis.title.x = element_text(size = 8, face = "bold"),
            axis.text.x = element_text(size = 6),
            axis.title.y = element_text(size = 10, face = "bold"),
            axis.text.y = element_text(size = 8),
            axis.line.x = element_line(colour = "black", linewidth = 0.5, linetype = "solid"),
            axis.line.y = element_line(colour = "black", linewidth = 0.5, linetype = "solid")) -> gg.da.bp
    
    daa.cont.bp %>% 
      filter(rank == daa.ranks[i]) %>%
      mutate(my.fill = ifelse(difference > 0, g1, g2)) %>% 
      select(taxa:CI_upper, my.fill) %>% 
      distinct() %>% 
      ggplot(., aes(x = difference, y = taxa, fill = my.fill)) + 
      geom_stripes(odd = "#88888833", even = "#FFFFFF00") +
      geom_errorbar(width=.15, position=position_dodge(.9), stat='identity', aes(xmin = CI_lower, xmax = CI_upper)) +
      geom_point(size = 2.5, shape = 21, color = "black") +
      geom_vline(xintercept = 0, linewidth = 0.25, linetype = 2) +
      theme_forest() +
      scale_x_continuous(limits = ~ c(-1, 1) * max(abs(.x))) +
      scale_fill_manual(values = tr.colors) +
      labs(y = daa.ranks[i], x = "Difference in mean proportions\n(%, 95% confidence intervals)") +
      theme(legend.position="none",
            axis.title.x = element_text(size = 8, face = "bold"),
            axis.text.x = element_text(size = 6),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.line.x = element_line(colour = "black", linewidth = 0.5, linetype = "solid")) -> gg.da.dp
    
    daa.cont.bp %>% 
      filter(rank == daa.ranks[i]) %>% 
      ggplot(., aes(y = taxa)) +
      geom_richtext(aes(x = 0, label = paste0("*p*: ", scales::scientific(p))), 
                    size = 2.5, color = "black", 
                    label.padding = unit(0, "pt"), fill = NA, label.color = NA) +
      theme_void() -> gg.da.tp
    
    daa.bp.ls[[i]] <- gg.da.bp
    daa.dp.ls[[i]] <- gg.da.dp
    daa.tp.ls[[i]] <- gg.da.tp
    daa.hg.ls[[i]] <- gg.da.hg
    
  }
  
  return(list(daa.bp.ls, daa.dp.ls, daa.tp.ls, daa.hg.ls))
  
}

gg.stamp.cont.1 <- gg.stamp.fun(tb.wilcox, g1 = "CTRL", g2 = "SM05")
gg.stamp.cont.2 <- gg.stamp.fun(tb.wilcox, g1 = "CTRL", g2 = "SM10")
gg.stamp.cont.3 <- gg.stamp.fun(tb.wilcox, g1 = "CTRL", g2 = "SM20")

(((gg.stamp.cont.1[[1]][[1]] + gg.stamp.cont.1[[2]][[1]] + gg.stamp.cont.1[[3]][[1]])) + plot_layout(widths = c(1, 1, .25))) /
  (((gg.stamp.cont.3[[1]][[1]] + gg.stamp.cont.3[[2]][[1]] + gg.stamp.cont.3[[3]][[1]])) + plot_layout(widths = c(1, 1, .25))) +
  plot_layout(heights = c(gg.stamp.cont.1[[4]][[1]], gg.stamp.cont.3[[4]][[1]])) -> gg.daa.p; gg.daa.p

(((gg.stamp.cont.1[[1]][[2]] + gg.stamp.cont.1[[2]][[2]] + gg.stamp.cont.1[[3]][[2]])) + plot_layout(widths = c(1, 1, .25))) /
  (((gg.stamp.cont.2[[1]][[1]] + gg.stamp.cont.2[[2]][[1]] + gg.stamp.cont.2[[3]][[1]])) + plot_layout(widths = c(1, 1, .25))) /
  (((gg.stamp.cont.3[[1]][[2]] + gg.stamp.cont.3[[2]][[2]] + gg.stamp.cont.3[[3]][[2]])) + plot_layout(widths = c(1, 1, .25))) + 
  plot_layout(heights = c(gg.stamp.cont.1[[4]][[2]], gg.stamp.cont.2[[4]][[1]], gg.stamp.cont.3[[4]][[2]])) -> gg.daa.c; gg.daa.c

(((gg.stamp.cont.1[[1]][[3]] + gg.stamp.cont.1[[2]][[3]] + gg.stamp.cont.1[[3]][[3]])) + plot_layout(widths = c(1, 1, .25))) /
  (((gg.stamp.cont.2[[1]][[2]] + gg.stamp.cont.2[[2]][[2]] + gg.stamp.cont.2[[3]][[2]])) + plot_layout(widths = c(1, 1, .25))) /
  (((gg.stamp.cont.3[[1]][[3]] + gg.stamp.cont.3[[2]][[3]] + gg.stamp.cont.3[[3]][[3]])) + plot_layout(widths = c(1, 1, .25))) + 
  plot_layout(heights = c(gg.stamp.cont.1[[4]][[3]], gg.stamp.cont.2[[4]][[2]], gg.stamp.cont.3[[4]][[3]])) -> gg.daa.o; gg.daa.o

(((gg.stamp.cont.1[[1]][[4]] + gg.stamp.cont.1[[2]][[4]] + gg.stamp.cont.1[[3]][[4]])) + plot_layout(widths = c(1, 1, .25))) /
  (((gg.stamp.cont.2[[1]][[3]] + gg.stamp.cont.2[[2]][[3]] + gg.stamp.cont.2[[3]][[3]])) + plot_layout(widths = c(1, 1, .25))) /
  (((gg.stamp.cont.3[[1]][[4]] + gg.stamp.cont.3[[2]][[4]] + gg.stamp.cont.3[[3]][[4]])) + plot_layout(widths = c(1, 1, .25))) + 
  plot_layout(heights = c(gg.stamp.cont.1[[4]][[4]], gg.stamp.cont.2[[4]][[3]], gg.stamp.cont.3[[4]][[4]])) -> gg.daa.f; gg.daa.f

ggsave(filename = paste0(path, "figures/DAA_1-Phylum.png"), plot = gg.daa.p, dpi = 1200, width = 5000 * 1.75, height = 5000 * .65, units = "px")
ggsave(filename = paste0(path, "figures/DAA_1-Phylum.svg"), plot = gg.daa.p, dpi = 1200, width = 5000 * 1.75, height = 5000 * .65, units = "px")

# ggsave(filename = paste0(path, "figures/DAA_2-Class.png"), plot = gg.daa.c, dpi = 1200, width = 5000 * 1.8, height = 5000 * 1.25, units = "px")
# ggsave(filename = paste0(path, "figures/DAA_2-Class.svg"), plot = gg.daa.c, dpi = 1200, width = 5000 * 1.8, height = 5000 * 1.25, units = "px")
# 
# ggsave(filename = paste0(path, "figures/DAA_3-Order.png"), plot = gg.daa.o, dpi = 1200, width = 5000 * 1.75, height = 5000 * 1.5, units = "px")
# ggsave(filename = paste0(path, "figures/DAA_3-Order.svg"), plot = gg.daa.o, dpi = 1200, width = 5000 * 1.75, height = 5000 * 1.5, units = "px")

ggsave(filename = paste0(path, "figures/DAA_4-Family.png"), plot = gg.daa.f, dpi = 1200, width = 5000 * 2, height = 5000 * 2.75, units = "px")
ggsave(filename = paste0(path, "figures/DAA_4-Family.svg"), plot = gg.daa.f, dpi = 1200, width = 5000 * 2, height = 5000 * 2.75, units = "px")

ls.daa <- list("Kruskal" = tb.kruskal, "Wilcoxon (post-hoc)" = tb.wilcox)

openxlsx::write.xlsx(x = ls.daa, file = paste0(path, "tables/DAAtaxa.xlsx"), rowNames = F)

# BoxPlot -----------------------------------------------------------------

tb.kruskal.bp <- 
  tb.kruskal %>% 
  filter(p.value < 0.05 & rank %in% c("Phylum", "Family"))%>% 
  select(taxa, rank, p.value) %>% 
  mutate(taxa = str_replace(taxa, "__", ": "), taxa = str_replace_all(taxa, "_", " ")) %>% 
  mutate(p.fancy = paste0("Kruskal-Wallis, *p*: ", ifelse(p.value < 0.001, "<0.001", sprintf('%1.3f', round(p.value, 3)))))

tb.wilcox.bp <-
  tb.wilcox %>% 
  mutate(taxa = str_replace(taxa, "__", ": "), taxa = str_replace_all(taxa, "_", " ")) %>% 
  filter(taxa %in% tb.kruskal.bp$taxa) %>% 
  group_by(taxa) %>%
  filter(sum(p.signif == "ns") < 3) %>%
  ungroup() %>% 
  select(taxa, rank, group1, group2, p.signif) %>% 
  mutate(x_start = as.numeric(group1), 
         x_end = as.numeric(group2))

tb.daa.bp <- 
  daa.glom.tss %>% 
  mutate(taxa = str_replace(taxa, "__", ": "), taxa = str_replace_all(taxa, "_", " ")) %>% 
  filter(taxa %in% tb.kruskal.bp$taxa)

tb.wilcox.bp <-
  tb.daa.bp %>% 
  group_by(taxa, rank) %>% 
  mutate(my_max = max(rel_abund)) %>% 
  select(taxa, rank, my_max) %>% 
  distinct() %>% 
  merge(tb.wilcox.bp, ., by = c("taxa", "rank")) %>% 
  arrange(taxa) %>% 
  mutate(y_start =
           case_when(group1 == "CTRL" & group2 == "SM05" ~ my_max * 1.15,
                     group1 == "CTRL" & group2 == "SM10" ~ my_max * 1.30,
                     group1 == "CTRL" & group2 == "SM20" ~ my_max * 1.45)) %>% 
  mutate(p.signif = str_replace_all(p.signif, "\\*", "\\\\*"))

for (i.taxa in unique(tb.wilcox.bp$taxa)) {
  
  ggsave(
    tb.daa.bp %>% 
      filter(taxa == i.taxa) %>% 
      ggplot(., aes(x = Group, y = rel_abund, fill = Group, color = Group)) +
      facet_grid(. ~ taxa) + 
      
      geom_boxplot(width = 0.75, outliers = F, alpha = .25) + 
      geom_point(position = position_jitterdodge(dodge.width=0.9), size = 2, shape = 16, alpha = .75) +
      geom_richtext(aes(x = (x_start + x_end)/2, y = y_start * 1.05, label = p.signif), size = 2.5,
                    colour = "black", data = tb.wilcox.bp %>% filter(taxa == i.taxa), inherit.aes = F,
                    label.padding = unit(0, "pt"), fill = NA, label.color = NA) +
      geom_richtext(aes(x = 2.5, y = Inf, label = p.fancy), position = position_dodge(width = 0), size = 3, 
                    hjust = 0.5, vjust = 3, colour = "black", data = tb.kruskal.bp %>% filter(taxa == i.taxa), inherit.aes = F,
                    label.padding = unit(0, "pt"), fill = NA, label.color = NA) +
      geom_segment(aes(x = x_start, xend = x_end, y = y_start), colour = "black", data = tb.wilcox.bp %>% filter(taxa == i.taxa), inherit.aes = F) +
      geom_segment(aes(x = x_start, xend = x_start, y = y_start*.99, yend = y_start), colour = "black", data = tb.wilcox.bp %>% filter(taxa == i.taxa), inherit.aes = F) +
      geom_segment(aes(x = x_end, xend = x_end, y = y_start, yend = y_start*.99), colour = "black", data = tb.wilcox.bp %>% filter(taxa == i.taxa), inherit.aes = F) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) + 
      ylab('Rel. abundance (%)<br>') + xlab("") + 
      scale_fill_manual(values = tr.colors, name = NULL) +
      scale_color_manual(values = tr.colors, name = NULL) +
      my_theme + 
      theme(legend.position = "none", 
            strip.text = element_markdown(face = "bold", color = "black", size = 10)),
    filename = paste0(path, "figures/DAAsingle/DAA-", i.taxa %>% 
                        str_replace_all(., "/", "-") %>% 
                        str_replace_all(., " ", "_") %>% 
                        str_replace_all(., ":", "_"),".png"), dpi = 1200, width = 4500 * 1, height = 4500 * 1.25, units = "px"
  )
  
}

rm(list = ls())

# End ---------------------------------------------------------------------

# Microbiome Analysis - Network -------------------------------------------

# Setup -------------------------------------------------------------------

rm(list = ls())

# install.packages("psych", dependencies = TRUE)
# install.packages("igraph", dependencies = TRUE)
# install.packages("qgraph", dependencies = TRUE)
# install.packages("vegan", dependencies = TRUE)
# install.packages("MCL", dependencies = TRUE)
# install.packages("devtools", dependencies = TRUE); library("devtools")
# devtools::install_github("zdk123/SpiecEasi", dependencies = TRUE)

# flattenCorrMatrix
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    from = rownames(cormat)[row(cormat)[ut]],
    to = rownames(cormat)[col(cormat)[ut]],
    r  = cormat[ut],
    p = pmat[ut]
  )
}

## Helper functions
# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

reorder_cor_and_p <- function(cormat, pmat){
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
  pmat <- pmat[hc$order, hc$order]
  list(r = cormat, p = pmat)
}

# Load libraries
library("ggrepel")
library("patchwork")
library("ape")
library("psych")
library("tidyverse")
library("reshape2")
library("igraph")
library("qgraph")
library("vegan")
library("MCL")
library("SpiecEasi")
library("phyloseq")
library("metagMisc")
library("xlsx")
library("extrafont") 

#font_import()
loadfonts(device = "all")

set.seed("123")

getwd()
setwd("~/sme_rhizobiome/")
getwd()

path <- paste0(getwd(), "/results/")
path

list.files(path)

load("analysis/network.RData")

# Data Wrangling ----------------------------------------------------------

tb.metad <- 
  read.delim(file = paste0(path, "tables/rTB02_metadata.txt"))

tb.count <- 
  read.delim(file = paste0(path, "tables/rTB01_counts.txt"), header = T, sep = "\t", stringsAsFactors = F, check.names = F)

names(tb.count) <- c("ID", as.character(tb.metad$ID))

tb.metad <- 
  tb.metad %>% 
  mutate(across(.cols = everything(), .fns = ~factor(., unique(.))))

tb.count <- tb.count %>% select(ID, tb.metad$ID)

tb.taxon <- read.delim(file = paste0(path, "tables/rTB03_taxonomy.txt"), header = T, sep = "\t", stringsAsFactors = F, check.names = F) %>% 
  dplyr::rename("ID" = "#TAXONOMY")

tb.tree <- read.tree(file = paste0(path, "tables/rTB04_tree.nwk"))

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

ps.rarefied <- phyloseq(ps.count, ps.taxon, ps.metad, ps.tree)
ps.rarefied.tss <- phyloseq_standardize_otu_abundance(physeq = ps.rarefied, method = 'total')

# -------------------------------------------------------------------------

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

get.taxa.glom.table <- function(my.ps.table, rank, trans = F) {
  
  return(
    
    if (trans) {
      
      merge(data.frame(check.names = F, tax_table(tax_glom(my.ps.table, taxrank= rank))) %>% select(1:any_of(rank)),
            data.frame(check.names = F, t(otu_table(tax_glom(my.ps.table, taxrank= rank)))), by = 0)[-1] %>%
        fun_bestax(.) %>% 
        select(Phylum, any_of(rank), where(is.numeric)) %>% 
        mutate(Phylum = str_replace_all(Phylum, pattern = "__", replacement = ": "),
               Phylum = str_replace_all(Phylum, pattern = "_", replacement = " "),
               Taxa = str_replace_all(get(rank), pattern = "__", replacement = ": "),
               Taxa = str_replace_all(Taxa, pattern = "_", replacement = " ")) %>% 
        select(Phylum, Taxa, where(is.numeric))
      
    } else {
      
      merge(data.frame(check.names = F, tax_table(tax_glom(my.ps.table, taxrank= rank))) %>% select(1:any_of(rank)),
            data.frame(check.names = F, otu_table(tax_glom(my.ps.table, taxrank= rank))), by = 0)[-1] %>%
        fun_bestax(.) %>% 
        select(Phylum, any_of(rank), where(is.numeric)) %>% 
        mutate(Phylum = str_replace_all(Phylum, pattern = "__", replacement = ": "),
               Phylum = str_replace_all(Phylum, pattern = "_", replacement = " "),
               Taxa = str_replace_all(get(rank), pattern = "__", replacement = ": "),
               Taxa = str_replace_all(Taxa, pattern = "_", replacement = " ")) %>% 
        select(Phylum, Taxa, where(is.numeric))
      
    }
    
  )
  
}

ls.net.names.top.p <- get.top.fun.names(ps.rarefied.tss, "Phylum", 10)

tb.net.f.abs <- 
  get.taxa.glom.table(ps.rarefied, "Family") %>% 
  mutate(Phylum = ifelse(Phylum %in% ls.net.names.top.p, Phylum, "Others"))

tb.net.f.rel <- 
  get.taxa.glom.table(phyloseq_standardize_otu_abundance(physeq = merge_samples(ps.rarefied, "Group", fun=sum), method = 'total'), "Family", T) %>% 
  mutate(Phylum = ifelse(Phylum %in% ls.net.names.top.p, Phylum, "Others"))

names(tb.net.f.rel) <- paste0(names(tb.net.f.rel), ".rel")

tb.net.f <- 
  cbind(tb.net.f.abs, tb.net.f.rel) %>% 
  select(Phylum, Taxa, 
         starts_with("CTRL"), starts_with("SM05"),
         starts_with("SM10"), starts_with("SM20")) %>% 
  mutate(RowSum = 
           (rowSums(select(., any_of(tb.metad$ID))))/
           sum(rowSums(select(., any_of(tb.metad$ID)))),
         RowName =
           make.unique(Taxa)
  ) %>% 
  arrange(-RowSum) %>% 
  filter(RowSum > 0.0001) %>% 
  column_to_rownames(., "RowName")

# Subset by treat
tb.count.1 <- tb.net.f %>% select(tb.metad %>% filter(Group == "CTRL") %>% pull(ID)) %>% t() %>% data.frame(check.names = F) %>% select(where(~ sum(.) != 0))
tb.count.2 <- tb.net.f %>% select(tb.metad %>% filter(Group == "SM05") %>% pull(ID)) %>% t() %>% data.frame(check.names = F) %>% select(where(~ sum(.) != 0))
tb.count.3 <- tb.net.f %>% select(tb.metad %>% filter(Group == "SM10") %>% pull(ID)) %>% t() %>% data.frame(check.names = F) %>% select(where(~ sum(.) != 0))
tb.count.4 <- tb.net.f %>% select(tb.metad %>% filter(Group == "SM20") %>% pull(ID)) %>% t() %>% data.frame(check.names = F) %>% select(where(~ sum(.) != 0))

# Correlation tables
otu.corr.pearson.1 <- corr.test(tb.count.1, method = "pearson", use = "pairwise", adjust = "none", alpha = 0.05)
otu.corr.pearson.2 <- corr.test(tb.count.2, method = "pearson", use = "pairwise", adjust = "none", alpha = 0.05)
otu.corr.pearson.3 <- corr.test(tb.count.3, method = "pearson", use = "pairwise", adjust = "none", alpha = 0.05)
otu.corr.pearson.4 <- corr.test(tb.count.4, method = "pearson", use = "pairwise", adjust = "none", alpha = 0.05)

# SPARCC

get.sparcc.fun <- function(tb) {
  
  my.tb.sparcc <- tb
  
  sparcc.out <- sparcc(my.tb.sparcc)
  sparcc.out
  
  rownames(sparcc.out$Cor) <- colnames(my.tb.sparcc)
  colnames(sparcc.out$Cor) <- colnames(my.tb.sparcc)
  rownames(sparcc.out$Cov) <- colnames(my.tb.sparcc)
  colnames(sparcc.out$Cov) <- colnames(my.tb.sparcc)
  
  sparcc.out.portable <- sparcc.out$Cor %>% reorder_cormat %>% get_upper_tri() %>% reshape2::melt() %>% na.omit()
  
  # sparcc.out.plot.1 <- sparcc.out.portable %>% ggplot(aes(x = Var2, y = Var1, fill = value)) + geom_tile() + scale_fill_gradient2() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  # sparcc.out.plot.1
  
  sparcc.out.boot <- sparccboot(my.tb.sparcc, R = 50, ncpus = 10)
  sparcc.out.boot.p <- pval.sparccboot(sparcc.out.boot)
  
  sparcc.out.boot.cors <- sparcc.out.boot.p$cors
  sparcc.out.boot.pvals <- sparcc.out.boot.p$pvals
  
  sparcc.out.pcors <- diag(0.5, nrow = dim(sparcc.out$Cor)[1], ncol = dim(sparcc.out$Cor)[1])
  sparcc.out.pcors[upper.tri(sparcc.out.pcors, diag=FALSE)] <- sparcc.out.boot.cors
  sparcc.out.pcors <- sparcc.out.pcors + t(sparcc.out.pcors)
  
  sparcc.out.pval <- diag(0.5, nrow = dim(sparcc.out$Cor)[1], ncol = dim(sparcc.out$Cor)[1])
  sparcc.out.pval[upper.tri(sparcc.out.pval, diag=FALSE)] <- sparcc.out.boot.pvals
  sparcc.out.pval <- sparcc.out.pval + t(sparcc.out.pval)
  
  rownames(sparcc.out.pcors) <- colnames(my.tb.sparcc)
  colnames(sparcc.out.pcors) <- colnames(my.tb.sparcc)
  rownames(sparcc.out.pval) <- colnames(my.tb.sparcc)
  colnames(sparcc.out.pval) <- colnames(my.tb.sparcc)
  
  sparcc.reordered.all <- reorder_cor_and_p(sparcc.out.pcors, sparcc.out.pval)
  sparcc.reordered.cor <- sparcc.reordered.all$r
  sparcc.reordered.p <- sparcc.reordered.all$p
  
  sparcc.reordered.cor.processed <- sparcc.reordered.cor  %>% get_upper_tri() %>% reshape2::melt() %>% na.omit() %>% rename(cor = value)
  sparcc.reordered.p.processed <- sparcc.reordered.p  %>% get_upper_tri() %>% reshape2::melt() %>% na.omit() %>% rename(p = value)
  
  sparcc.out.fdr <- 
    left_join(sparcc.reordered.cor.processed, sparcc.reordered.p.processed, by = c("Var1", "Var2")) %>%
    filter(Var1 != Var2) %>% 
    mutate(fdr = p.adjust(p, method = "BH"))
  
  sparcc.out.fdr.sig <- 
    sparcc.out.fdr %>% 
    filter(p < 0.05 & abs(cor) > 0.5) 
  
  sparcc.out.plot.2 <- 
    sparcc.out.fdr %>% 
    ggplot(aes(x = Var2, y = Var1, fill = cor)) + 
    geom_tile() + 
    geom_point(data = sparcc.out.fdr.sig, shape = 1) +
    scale_fill_gradient2() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          aspect.ratio = 1)
  
  print(sparcc.out.plot.2)
  
  return(sparcc.out.fdr)
  
}

otu.corr.sparcc.1 <- get.sparcc.fun(tb.count.1)
otu.corr.sparcc.2 <- get.sparcc.fun(tb.count.2)
otu.corr.sparcc.3 <- get.sparcc.fun(tb.count.3)
otu.corr.sparcc.4 <- get.sparcc.fun(tb.count.4)

otu.corr.pearson.flatten.1 <- flattenCorrMatrix(otu.corr.pearson.1$r, otu.corr.pearson.1$p)
otu.corr.pearson.flatten.2 <- flattenCorrMatrix(otu.corr.pearson.2$r, otu.corr.pearson.2$p)
otu.corr.pearson.flatten.3 <- flattenCorrMatrix(otu.corr.pearson.3$r, otu.corr.pearson.3$p)
otu.corr.pearson.flatten.4 <- flattenCorrMatrix(otu.corr.pearson.4$r, otu.corr.pearson.4$p)

# Nets ----------------------------------------------------------------------------------------

# Net 1
otu.corr.edges.1 <- 
  otu.corr.sparcc.1 %>% 
  filter(abs(cor) >= 0.75 & p <= 0.05) %>%
  dplyr::rename(from = Var1, to = Var2, correlation = cor, p.value = p) %>% 
  mutate(type = if_else(correlation > 0, 2, 1))

otu.corr.nodes.1 <- 
  merge(x = sort(unique(c(otu.corr.edges.1$from, otu.corr.edges.1$to))), y = tb.net.f, by.x = 1, by.y = 0, all.x = T) %>% 
  select(Phylum, Taxa, starts_with("CTRL")) 

# Net 2
otu.corr.edges.2 <- 
  otu.corr.sparcc.2 %>% 
  filter(abs(cor) >= 0.75 & p <= 0.05) %>%
  dplyr::rename(from = Var1, to = Var2, correlation = cor, p.value = p) %>% 
  mutate(type = if_else(correlation > 0, 2, 1))

otu.corr.nodes.2 <- 
  merge(x = sort(unique(c(otu.corr.edges.2$from, otu.corr.edges.2$to))), y = tb.net.f, by.x = 1, by.y = 0, all.x = T) %>% 
  select(Phylum, Taxa, starts_with("SM05")) 

# Net 3
otu.corr.edges.3 <- 
  otu.corr.sparcc.3 %>% 
  filter(abs(cor) >= 0.75 & p <= 0.05) %>%
  dplyr::rename(from = Var1, to = Var2, correlation = cor, p.value = p) %>% 
  mutate(type = if_else(correlation > 0, 2, 1))

otu.corr.nodes.3 <- 
  merge(x = sort(unique(c(otu.corr.edges.3$from, otu.corr.edges.3$to))), y = tb.net.f, by.x = 1, by.y = 0, all.x = T) %>% 
  select(Phylum, Taxa, starts_with("SM10")) 

# Net 4
otu.corr.edges.4 <- 
  otu.corr.sparcc.4 %>% 
  filter(abs(cor) >= 0.75 & p <= 0.05) %>%
  dplyr::rename(from = Var1, to = Var2, correlation = cor, p.value = p) %>% 
  mutate(type = if_else(correlation > 0, 2, 1))

otu.corr.nodes.4 <- 
  merge(x = sort(unique(c(otu.corr.edges.4$from, otu.corr.edges.4$to))), y = tb.net.f, by.x = 1, by.y = 0, all.x = T) %>% 
  select(Phylum, Taxa, starts_with("SM20")) 

# Getting IDs

top.p <- c(ls.net.names.top.p, "Others")

all.ids <- 
  reduce(list(otu.corr.nodes.1, otu.corr.nodes.2, otu.corr.nodes.3, otu.corr.nodes.4), full_join, by = c("Phylum", "Taxa")) %>% 
  mutate(Sum = rowSums(select(., matches(".rel")), na.rm = TRUE)) %>%
  arrange(desc(Sum)) %>% 
  mutate(IDT = seq_len(nrow(.))) %>%
  mutate(Phylum = factor(Phylum, levels = top.p)) %>%
  mutate(IDP = as.numeric(Phylum)) %>% 
  mutate(ID = Taxa) %>% 
  select(ID, IDP, Phylum, IDT, Taxa)

otu.corr.nodes.1 <- merge(all.ids, otu.corr.nodes.1, sort = F) %>% select(ID, IDP, Phylum, IDT, Taxa, everything())
otu.corr.nodes.2 <- merge(all.ids, otu.corr.nodes.2, sort = F) %>% select(ID, IDP, Phylum, IDT, Taxa, everything())
otu.corr.nodes.3 <- merge(all.ids, otu.corr.nodes.3, sort = F) %>% select(ID, IDP, Phylum, IDT, Taxa, everything())
otu.corr.nodes.4 <- merge(all.ids, otu.corr.nodes.4, sort = F) %>% select(ID, IDP, Phylum, IDT, Taxa, everything())

net.1 <- graph_from_data_frame(d=otu.corr.edges.1, vertices=otu.corr.nodes.1, directed = F)
net.2 <- graph_from_data_frame(d=otu.corr.edges.2, vertices=otu.corr.nodes.2, directed = F)
net.3 <- graph_from_data_frame(d=otu.corr.edges.3, vertices=otu.corr.nodes.3, directed = F)
net.4 <- graph_from_data_frame(d=otu.corr.edges.4, vertices=otu.corr.nodes.4, directed = F)

otu.corr.nodes.1 <- 
  merge(otu.corr.nodes.1, 
        reduce(list(data.frame(degree = degree(net.1, mode="all")) %>% rownames_to_column("Taxa"), 
                    data.frame(betweenness = betweenness(net.1, directed=T, weights=NA, normalized = T)) %>% rownames_to_column("Taxa"), 
                    data.frame(closeness = closeness(net.1, mode="all", normalized = T)) %>% rownames_to_column("Taxa")), 
               full_join, by = "Taxa")) %>% 
  mutate(sig = ifelse(degree >= 10 & betweenness >= 0.05, "*", "")) %>% 
  select(ID, everything()) %>% 
  dplyr::rename("rel" = "CTRL.rel")

otu.corr.nodes.2 <- 
  merge(otu.corr.nodes.2, 
        reduce(list(data.frame(degree = degree(net.2, mode="all")) %>% rownames_to_column("Taxa"), 
                    data.frame(betweenness = betweenness(net.2, directed=T, weights=NA, normalized = T)) %>% rownames_to_column("Taxa"), 
                    data.frame(closeness = closeness(net.2, mode="all", normalized = T)) %>% rownames_to_column("Taxa")), 
               full_join, by = "Taxa")) %>% 
  mutate(sig = ifelse(degree >= 10 & betweenness >= 0.05, "*", "")) %>% 
  select(ID, everything()) %>% 
  dplyr::rename("rel" = "SM05.rel")

otu.corr.nodes.3 <- 
  merge(otu.corr.nodes.3, 
        reduce(list(data.frame(degree = degree(net.3, mode="all")) %>% rownames_to_column("Taxa"), 
                    data.frame(betweenness = betweenness(net.3, directed=T, weights=NA, normalized = T)) %>% rownames_to_column("Taxa"), 
                    data.frame(closeness = closeness(net.3, mode="all", normalized = T)) %>% rownames_to_column("Taxa")), 
               full_join, by = "Taxa")) %>% 
  mutate(sig = ifelse(degree >= 10 & betweenness >= 0.05, "*", "")) %>% 
  select(ID, everything()) %>% 
  dplyr::rename("rel" = "SM10.rel")

otu.corr.nodes.4 <- 
  merge(otu.corr.nodes.4, 
        reduce(list(data.frame(degree = degree(net.4, mode="all")) %>% rownames_to_column("Taxa"), 
                    data.frame(betweenness = betweenness(net.4, directed=T, weights=NA, normalized = T)) %>% rownames_to_column("Taxa"), 
                    data.frame(closeness = closeness(net.4, mode="all", normalized = T)) %>% rownames_to_column("Taxa")), 
               full_join, by = "Taxa")) %>% 
  mutate(sig = ifelse(degree >= 10 & betweenness >= 0.05, "*", "")) %>% 
  select(ID, everything()) %>% 
  dplyr::rename("rel" = "SM20.rel")

# Net 1
net.1 <- graph_from_data_frame(d=otu.corr.edges.1, vertices=otu.corr.nodes.1, directed = F)

V(net.1)$name <- seq(1, length(V(net.1)$name)) # Customized layout to avoid nodes overlapping. Convert node label from names to numerical IDs
e.1 <- get.edgelist(net.1) ; class(e.1) <- "numeric"
l.1 <- qgraph.layout.fruchtermanreingold(e.1, vcount=vcount(net.1), area=8*(vcount(net.1)^2.2), repulse.rad=(vcount(net.1)^2.5))
#marks.1 <- lapply(1:length(articulation.points(net.1)), function(x) which(V(net.1)$name == articulation.points(net.1)[[x]]))

# Net 2
net.2 <- graph_from_data_frame(d=otu.corr.edges.2, vertices=otu.corr.nodes.2, directed = F)

V(net.2)$name <- seq(1, length(V(net.2)$name)) # Customized layout to avoid nodes overlapping. Convert node label from names to numerical IDs
e.2 <- get.edgelist(net.2) ; class(e.2) <- "numeric"
l.2 <- qgraph.layout.fruchtermanreingold(e.2, vcount=vcount(net.2), area=8*(vcount(net.2)^2.2), repulse.rad=(vcount(net.2)^2.5))
#marks.2 <- lapply(1:length(articulation.points(net.2)), function(x) which(V(net.2)$name == articulation.points(net.2)[[x]]))

# Net 3
net.3 <- graph_from_data_frame(d=otu.corr.edges.3, vertices=otu.corr.nodes.3, directed = F)

V(net.3)$name <- seq(1, length(V(net.3)$name)) # Customized layout to avoid nodes overlapping. Convert node label from names to numerical IDs
e.3 <- get.edgelist(net.3) ; class(e.3) <- "numeric"
l.3 <- qgraph.layout.fruchtermanreingold(e.3, vcount=vcount(net.3), area=8*(vcount(net.3)^2.2), repulse.rad=(vcount(net.3)^2.5))
#marks.3 <- lapply(1:length(articulation.points(net.3)), function(x) which(V(net.3)$name == articulation.points(net.3)[[x]]))

# Net 4
net.4 <- graph_from_data_frame(d=otu.corr.edges.4, vertices=otu.corr.nodes.4, directed = F)

V(net.4)$name <- seq(1, length(V(net.4)$name)) # Customized layout to avoid nodes overlapping. Convert node label from names to numerical IDs
e.4 <- get.edgelist(net.4) ; class(e.4) <- "numeric"
l.4 <- qgraph.layout.fruchtermanreingold(e.4, vcount=vcount(net.4), area=8*(vcount(net.4)^2.2), repulse.rad=(vcount(net.4)^2.5))
#marks.4 <- lapply(1:length(articulation.points(net.4)), function(x) which(V(net.4)$name == articulation.points(net.4)[[x]]))

# cols

cols.p <- colorRampPalette(c(
  
  paletteer::paletteer_d("futurevisions::pegasi")
  
))(length(unique(top.p))-1)

names(cols.p) <- unique(top.p)[1:(length(unique(top.p))-1)]

cols.p["Others"] <- "#E6E6E6FF"

cols.cor <- c("#9BC3C9FF", "#E88170FF")
names(cols.cor) <- c("Negative", "Positive")

# Main plot function

# SVG graphics device
svg(paste0(path, "figures/Network/Network1.svg"),
    width=25,
    height=25,
    pointsize=50)

png(filename = paste0(path, "figures/Network/Network.png"), width = 5000 * 1.5, height = 5000 * 1.75, units = "px", res = 1200, type = "cairo")

par(mfrow = c(2, 2), mai = c(0, 0, 0.5, 0), bg = 'white')

# Net 1
plot.igraph(net.1, 
            vertex.size = 4 + log2(V(net.1)$degree)*.75,
            vertex.label = V(net.1)$sig,
            vertex.label.cex = 1 + log2(V(net.1)$degree)/100,
            vertex.label.color = "black",
            vertex.frame.color = "black",
            vertex.color = cols.p[V(net.1)$IDP],
            edge.color = cols.cor[E(net.1)$type],
            edge.width = .5,
            edge.arrow.mode = "-", 
            edge.arrow.size=0,
            edge.curved=0.5,
            layout = l.1, 
            rescale = T,
            main = "CTRL", 
            margin = c(0, 0, -0.05, 0))

# Net 2
plot.igraph(net.2, 
            vertex.size = 4 + log2(V(net.2)$degree)*.75,
            vertex.label = V(net.2)$sig,
            vertex.label.cex = 1 + log2(V(net.2)$degree)/100,
            vertex.label.color = "black",
            vertex.frame.color = "black",
            vertex.color = cols.p[V(net.2)$IDP],
            edge.color = cols.cor[E(net.2)$type],
            edge.width = .5,
            edge.arrow.mode = "-", 
            edge.arrow.size=0,
            edge.curved=0.5,
            layout = l.2, 
            rescale = T,
            main = "SM05",
            margin = c(0, 0, -0.05, 0))

# Net 3
plot.igraph(net.3, 
            vertex.size = 4 + log2(V(net.3)$degree)*.75,
            vertex.label = V(net.3)$sig,
            vertex.label.cex = 1 + log2(V(net.3)$degree)/100,
            vertex.label.color = "black",
            vertex.frame.color = "black",
            vertex.color = cols.p[V(net.3)$IDP],
            edge.color = cols.cor[E(net.3)$type],
            edge.width = .5,
            edge.arrow.mode = "-", 
            edge.arrow.size=0,
            edge.curved=0.5,
            layout = l.3, 
            rescale = T,
            main = "SM10",
            margin = c(0, 0, -0.05, 0))

# Net 4
plot.igraph(net.4, 
            vertex.size = 4 + log2(V(net.4)$degree)*.75,
            vertex.label = V(net.4)$sig,
            vertex.label.cex = 1 + log2(V(net.4)$degree)/100,
            vertex.label.color = "black",
            vertex.frame.color = "black",
            vertex.color = cols.p[V(net.4)$IDP],
            edge.color = cols.cor[E(net.4)$type],
            edge.width = .5,
            edge.arrow.mode = "-", 
            edge.arrow.size=0,
            edge.curved=0.5,
            layout = l.4, 
            rescale = T,
            main = "SM20",
            margin = c(0, 0, -0.05, 0))

# Close the graphics device
dev.off() 

# Legendas ------------------------------------------------------------------------------------

library(ggtext)

my_theme <-   theme(legend.position = "right", 
                    legend.text = element_text(size = 8, colour = 'black'),
                    legend.title = element_blank(),
                    legend.justification = c(0.5, 0.5),
                    
                    strip.text = element_text(size = 10, colour = 'black', face = "bold"),
                    strip.background = element_rect(fill = NA, color = NA), 
                    
                    panel.background = element_rect(fill = "white"),
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.border = element_rect(color = "black", linewidth = 1, fill = NA),
                    
                    axis.ticks = element_line(color = "black", linewidth = .5),
                    axis.ticks.length = unit(2, "pt"),
                    
                    axis.title.x = element_text(size = 9, face = 'bold', colour = 'black'),
                    axis.title.y = element_text(size = 9, face = 'bold', colour = 'black'),
                    axis.title.y.right = element_text(size = 9, face = 'bold', colour = 'black'),
                    axis.text.x = element_text(size = 7, colour = 'black'),
                    axis.text.y = element_text(size = 7, colour = 'black'),
                    
                    plot.title = element_text(size = 10, face = 'bold', colour = 'black'),
                    plot.subtitle = element_text(size = 8, colour = 'black'),
                    plot.margin = unit(c(0.5, 0.75, 0.5, 0.75), 'cm'))


gg.leg.p.full <- 
  ggplot(top.p %>% data.frame(), aes(x = factor(., levels = unique(.)), y = 1:length(top.p))) +
  geom_point(aes(fill = factor(., levels = unique(.))), color = "black", shape = 22, size = 3) +
  scale_fill_manual(values = cols.p) +  
  labs(fill = "Phylum") + 
  theme_bw() +
  my_theme + 
  theme(legend.position = "right",
        axis.ticks.length.x = unit(0, "pt"),
        axis.ticks.length.y = unit(2.5, "pt"),
        axis.text.x = element_blank(),
        legend.justification = c(0, .5), 
        legend.box.margin = margin(0, 0, 0, 0),
        legend.title = element_text(size = 7 + 2, face = "bold"),
        legend.text  = element_text(size = 7),
        legend.key.size = unit(0.65, "lines"))
  
gg.leg.p <- ggpubr::get_legend(gg.leg.p.full)

gg.leg.line.full <- 
  ggplot(data.frame(y = factor(rep(c("Positive", "Negative"), 2), levels = c("Positive", "Negative")), x = c(1, 2, 2, 1)), aes(x = x, y = y, group = y, color = y)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = cols.cor, name = "Correlation (SparCC)") +  
  theme_bw() +
  my_theme + 
  theme(legend.position = "right",
        axis.ticks.length.x = unit(0, "pt"),
        axis.ticks.length.y = unit(2.5, "pt"),
        axis.text.x = element_blank(),
        legend.justification = c(0, .5), 
        legend.box.margin = margin(0, 0, 0, 0),
        legend.title = element_text(size = 7 + 2, face = "bold"),
        legend.text  = element_text(size = 7),
        legend.key.size = unit(0.65, "lines"))

gg.leg.line <- ggpubr::get_legend(gg.leg.line.full)

gg.leg.degree.full <- 
  ggplot(data.frame(y = 1:3), aes(x = 1, y = y, size = y)) +
  geom_point() + 
  theme_bw() +
  my_theme + 
  scale_size_area(max_size = 5, breaks=c(1, 2, 3), labels=c("", "", ""), name = "Degree") +
  theme(legend.position = "right",
        axis.ticks.length.x = unit(0, "pt"),
        axis.ticks.length.y = unit(2.5, "pt"),
        axis.text.x = element_blank(),
        legend.justification = c(0, .5), 
        legend.box.margin = margin(0, 0, 0, 0),
        legend.title = element_text(size = 7 + 2, face = "bold"),
        legend.text  = element_text(size = 7),
        legend.key.size = unit(0.65, "lines"))

gg.leg.degree <- ggpubr::get_legend(gg.leg.degree.full)

png(filename = paste0(path, "figures/Network/Network_Legend.png"), width = 5000 * .5, height = 5000 * 1, units = "px", res = 1200, type = "cairo")

ggpubr::as_ggplot(gg.leg.p) + ggpubr::as_ggplot(gg.leg.line) + ggpubr::as_ggplot(gg.leg.degree) + plot_layout(ncol = 1, heights = c(5, 2, 3))

dev.off() 

# Keystone ----------------------------------------------------------------

tb.keystone <-
  bind_rows(
    otu.corr.nodes.1 %>% select(Phylum, Taxa, degree, betweenness, closeness, sig) %>% mutate(Group = "CTRL"),
    otu.corr.nodes.2 %>% select(Phylum, Taxa, degree, betweenness, closeness, sig) %>% mutate(Group = "SM05"),
    otu.corr.nodes.3 %>% select(Phylum, Taxa, degree, betweenness, closeness, sig) %>% mutate(Group = "SM10"),
    otu.corr.nodes.4 %>% select(Phylum, Taxa, degree, betweenness, closeness, sig) %>% mutate(Group = "SM20")) %>% 
  mutate(Group = factor(Group, levels = c("CTRL", "SM05", "SM10", "SM20"))) %>% 
  mutate(my.label = ifelse(sig == "*", Taxa, NA_character_)) 

ggplot(tb.keystone, aes(y = degree, x = betweenness)) +
  facet_wrap( ~ Group, scales = "fixed") +
  geom_point(shape = 21, size = 2.5, alpha = .75, color = "black", aes(fill = Phylum), show.legend = T) +
  geom_hline(yintercept = 10, size = 0.5, linetype = 2, color = "black") +
  geom_vline(xintercept = 0.05, size = 0.5, linetype = 2) +
  geom_label_repel(aes(label = my.label), size = 2, min.segment.length = 0.2, max.overlaps = 999, fill = "white", segment.color = 'black', color = "black", parse = T) +
  scale_y_continuous(limits = c(0, 45)) +
  scale_x_continuous(limits = c(0, 0.2)) +
  ylab('Degree\n') +
  xlab('\nBetweenness centrality') +
  scale_fill_manual(name = "", values = cols.p) +
  guides(fill = 
           guide_legend(order = 1, 
                        override.aes = list(shape = 21, 
                                            fill = cols.p, 
                                            color = "black", 
                                            size = 3,
                                            linetype = 0,
                                            alpha = 1),
                        ncol = 1)) +
  theme_bw() + 
  my_theme +
  theme(legend.position = "right", 
        legend.box.background = element_rect(fill = "white", colour = "black", linewidth = .75),
        legend.justification = c(0, .5), 
        legend.box.margin = margin(0, 0, 0, 0),
        legend.text  = element_text(size = 7),
        legend.key.size = unit(0.65, "lines")) -> a

a

ggsave(filename = paste0(path, "figures/Network/Network_Keystone.png"), 
       plot = a, 
       dpi = 1200,
       width = 5000 * 2.25,
       height = 5000 * 2,
       units = "px")

ggsave(filename = paste0(path, "figures/Network/Network_Keystone.svg"), 
       plot = a, 
       dpi = 1200,
       width = 5000 * 2.25,
       height = 5000 * 2,
       units = "px")

# Comp. Measures ----------------------------------------------------------

tb.keystone

tb.metad <- 
  read.delim(file = paste0(path, "tables/rTB02_metadata.txt")) %>% 
  mutate(across(.cols = everything(), .fns = ~factor(., unique(.))))

tr.colors <- c("#525252FF", "#8cc084FF", "#38a3a5FF", "#0466c8FF")
names(tr.colors) <- unique(tb.metad$FullGroup)

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

net.measure.table <- merge(tb.keystone, tb.metad %>% select(Group, FullGroup) %>% distinct(), by = "Group", sort = F)

net.measure.table.longer <- 
  net.measure.table %>% 
  select(Group, FullGroup, Phylum, Taxa, degree:closeness) %>% 
  pivot_longer(cols = degree:closeness, names_to = "measure", values_to = "values") %>% 
  data.frame(check.names = F) %>% 
  mutate(measure = str_to_title(measure)) %>% 
  mutate(measure = factor(measure, levels = unique(measure)))

net.measure.table.summary <-
  net.measure.table.longer %>% 
  group_by(Group, FullGroup, measure) %>% 
  summarise(temp = list(c(psych::describe(values), quantile(values, na.rm = T), IQR = matrixStats::iqr(values, na.rm = T)))) %>% 
  unnest_wider(temp) %>% 
  select(-c(vars, trimmed, mad)) %>% 
  mutate(across(where(is.numeric), ~case_when(is.nan(.) ~ NA_integer_, is.infinite(.) ~ NA_integer_,.default = .))) %>% 
  data.frame(check.names = F)

get_kruskal <- function(tb) {
  
  tb %>%
    group_by(measure) %>%
    nest() %>%
    mutate(kruskal = map(.x = data, ~kruskal.test(values ~ Group, data = .x) %>% broom::tidy())) %>%
    unnest(kruskal) %>%
    mutate(across(where(is.numeric), ~round(., 4))) -> p
  
  return(p)
  
}

net.measure.kruskal <- get_kruskal(net.measure.table.longer)

net.measure.kruskal.2 <-
  net.measure.kruskal %>% 
  mutate(p.full = paste0("Kruskal-Wallis, *p*: ", ifelse(p.value < 0.001, "<0.001", sprintf("%1.3f", round(p.value, 3)))))

get_posthoc_lsd <- function(tb, tb.summary, tb.stat, ls.invert = NULL) {
  
  tb.compm <- data.frame(measure = "", Group = "", group = "")[-1,]
  
  for (i in unique(tb$measure)) {
    
    tmp.tb <- tb %>% filter(measure == i)
    
    if (i %in% ls.invert) {
      
      tmp.tb$values <- -1 * tmp.tb$values
      
    }
    
    tmp.lsd <- with(tmp.tb, agricolae::kruskal(tmp.tb$values, tmp.tb$Group, alpha = 0.05))
    tmp.lsd <-
      cbind(Group = row.names(tmp.lsd[["groups"]]), tmp.lsd[["groups"]]) %>% 
      remove_rownames() %>% 
      mutate(measure = i, group = str_trim(groups)) %>% 
      select(measure, Group, group)
    
    if (tb.stat %>% filter(measure == i) %>% pull(p.value) > 0.05) {
      
      tmp.lsd <- mutate(tmp.lsd, group = NA_character_)
      
    }
    
    tb.compm <- bind_rows(tb.compm, tmp.lsd)
    
  }
  
  tb.compm <- merge(tb.summary, tb.compm, sort = F, all.x = T)
  
  return(tb.compm)
  
} 

net.measure.compm.lsd <- get_posthoc_lsd(net.measure.table.longer, net.measure.table.summary, net.measure.kruskal) %>% arrange(measure)

net.measure.gg <- 
  ggplot(net.measure.table.longer, aes(x = Group, y = values)) + 
  facet_wrap( ~ measure, scales = "free", nrow = 1, ncol = 5) +
  geom_jitter(aes(fill = FullGroup), width = 0.25, alpha = .25, size = 1, shape = 21) +
  geom_richtext(mapping = aes(label = p.full, x = 2.5, y = Inf, vjust = 2.5),
                data = net.measure.kruskal.2, position = position_dodge(width = 0),
                size = 2.5, hjust = 0.5, color = "black", inherit.aes = T, 
                label.padding = unit(0, "pt"), fill = NA, label.color = NA) +
  geom_richtext(mapping = aes(label = group, x = as.numeric(FullGroup), y = max), 
                data = net.measure.compm.lsd,
                size = 2.5, vjust = -1, hjust = 0.5, color = "black", inherit.aes = T,
                label.padding = unit(0, "pt"), fill = NA, label.color = NA) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.5)), breaks = scales::extended_breaks(n = 6)) +
  scale_x_discrete(expand = expansion(mult = c(0.25, 0.25))) +
  scale_fill_manual(values = tr.colors, name = NULL) +
  ylab('') + xlab('') +
  guides(fill = 
           guide_legend(order = 1, 
                        override.aes = list(shape = 22, 
                                            fill = tr.colors, 
                                            color = "black", 
                                            size = 3,
                                            linetype = 0, 
                                            alpha = 1),
                        ncol = 2)) +
  my_theme + 
  theme(legend.position = "bottom", 
        legend.box.background = element_rect(fill = "white", colour = "black", linewidth = .75),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

net.measure.gg

ggsave(filename = paste0(path, "figures/Network/Network_Measures.png"), plot = net.measure.gg, dpi = 1200, width = 5000*1.5, height = 5000*1, units = "px")
ggsave(filename = paste0(path, "figures/Network/Network_Measures.svg"), plot = net.measure.gg, dpi = 1200, width = 5000*1.5, height = 5000*1, units = "px")

# Network parameters --------------------------------------------------------------------------

nets <- c("net.1", "net.2", "net.3", "net.4")

net.name <- c("CTRL", "SM05", "SM10", "SM20")

net.attr.table <- data.frame("Attributes" = c("n. of nodes", 
                                              "n. of edges", 
                                              "positive edges",
                                              "negative edges",
                                              "n. of clusters",
                                              "size of clusters",
                                              "edge density",
                                              "avg. degree",
                                              "max. degree",
                                              "avg. closeness",
                                              "avg. betweenness",
                                              "max. betweenness",
                                              "n. keystone taxa",
                                              "modularity coef.",
                                              "articulation points",
                                              "clustering coefficient",
                                              "diameter"))


for (i in seq_along(nets)) {
  
  n <- net.name[i]
  
  net <- get(nets[i])
  
  # Getting IDG back to name
  V(net)$name <- paste0(V(net)$IDT, ") ", V(net)$Taxa)
  
  # n. of vertices
  attr01 <- vcount(net)
  
  # n. of edges
  attr02 <- ecount(net)
  
  # n. of positive edges
  attr03 <- paste0(unname(table(E(net)$type)[2]), " (", round((unname(table(E(net)$type)[2])/ecount(net))*100, 2), "%)")
  
  # n. of negative edges
  attr04 <- paste0(unname(table(E(net)$type)[1]), " (", round((unname(table(E(net)$type)[1])/ecount(net))*100, 2), "%)")
  
  # n. of clusters
  # clusters(net)$membership                                  # Membership of vertices
  attr05 <- components(net)$no                                   # Number of clusters
  attr06 <- paste(components(net)$csize, collapse = ", ")        # Sizes of clusters
  
  # Hub detection
  attr07 <- round(edge_density(net, loops=F), 3) # The proportion of present edges from all possible edges in the network
  
  attr08 <- paste0(round(mean(degree(net, mode = "all")), 3),
                   " (±", round(sd(degree(net, mode = "all")), 3), ")") # Average degrees of nodes
  
  attr09 <- max(degree(net, mode="all")) # Node/s with the most degrees
  
  attr10 <- paste0(round(mean(closeness(net, mode = "all", weights = NA, normalized = T)), 3),
                   " (±", round(sd(closeness(net, mode = "all", weights = NA, normalized = T)), 3), ")") # Average closeness of nodes
  
  attr11 <- paste0(round(mean(betweenness(net, directed=T, weights=NA, normalized = T)), 3),
                   " (±", round(sd(betweenness(net, directed=T, weights=NA, normalized = T)), 3), ")") # Average betweenness of nodes
  
  attr12 <- round(max(betweenness(net, directed=T, weights=NA, normalized = T)), 3) # Node/s with the most degrees
  
  # Sort the species based on hubbiness score
  attr13 <- unname(table(V(net)$sig)["*"])
    
  # Modularity
  attr14 <- round(modularity(net, membership(components(net))), 4)
  
  # Articulation.points
  attr15 <- length(articulation_points(net))
  
  # Clustering Coeficient
  attr16 <- round(transitivity(net, type="global"), 3)
  
  # Diameter
  attr17 <- diameter(net)
  
  attr.col <- t(rbind(lapply(ls(pattern = "^attr[0-9]"), get)))
  
  net.attr.table <- cbind(net.attr.table, attr.col)
  names(net.attr.table)[i+1] <- n 
  
  
}

# -------------------------------------------------------------------------

ls.tables <- 
  
  list(
    
    "Net (Summary)" = 
      net.attr.table %>% 
      mutate(Attributes = str_to_title(Attributes) %>% str_replace(., "Of", "of")) %>% 
      mutate(across(everything(), ~as.character(.) %>% str_replace("NA", NA_character_))) %>% 
      data.frame(),
    
    "Net (Measures)" = 
      net.measure.table %>% 
      select(Net = Group, Phylum, Taxa, degree:closeness, keystone = sig) %>%
      mutate(across(where(is.factor), ~as.character(.))),
    
    "Net Measures (Kruskal-Wallis)" = 
      net.measure.kruskal.2 %>%
      select(Measure = measure, statistic, p = p.value, method) %>% 
      mutate(across(where(is.factor), ~as.character(.))) %>% 
      data.frame(),
    
    "Net Measures (Summary)" =
      net.measure.compm.lsd %>% 
      select(Group, FullGroup, Measure = measure, everything()) %>% 
      mutate(across(where(is.factor), ~as.character(.))) %>% 
      mutate(across(where(is.character), ~str_remove_all(., "<.?sup>"))) %>% 
      mutate(across(where(is.character), ~str_remove_all(., "\\*")))
    
  )

openxlsx::write.xlsx(x = ls.tables, file = paste0(path, "tables/Network.xlsx"), rowNames = F)

# End ---------------------------------------------------------------------

save.image("analysis/network.RData")

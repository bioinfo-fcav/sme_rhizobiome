
# Functional --------------------------------------------------------------

rm(list = ls())

# Setup -------------------------------------------------------------------

library(tidyverse)
library(ggh4x)
library(patchwork)
library(cowplot)
library(ggprism)
library(xlsx)
library(DESeq2)
library(pheatmap) 
library(FactoMineR)
library(factoextra)
library(ggtext)
library(ggrepel)
library(ggalluvial)
library(extrafont) 

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

# Data Wrangling ----------------------------------------------------------

tb.metad <- 
  read.delim(file = paste0(path, "tables/rTB02_metadata.txt")) %>% 
  dplyr::rename(SAMPLE_ID = 1)

tb.count.ko <-
  read.delim(file = "picrust/results/KO_pred.tsv") %>% 
  dplyr::rename(ID = 1) %>% 
  column_to_rownames("ID") %>% 
  select(tb.metad$SAMPLE_ID) %>% 
  filter(rowSums(select(., tb.metad$SAMPLE_ID) * 1, na.rm = TRUE) != 0)

tb.ann.ko <- 
  read.delim(file = "~/work_lamoroso/db/OTHERS/KO_info_2024-12-17.tsv") %>% 
  select(KO_ID, D1:D4)

tb.ann.ko.filt <- 
  read.delim(file = "~/work_lamoroso/db/OTHERS/KO_info_2024-12-17.tsv") %>% 
  select(KO_ID, D1:D4) %>% 
  filter(D1 %in% c("Metabolism", "Environmental Information Processing", "Genetic Information Processing", "Cellular Processes"))

tb.ann.pb <- 
  read.delim(file = "~/work_lamoroso/db/OTHERS/pgpm_traits.tsv") %>% 
  select(PGPT_ID, GENE, KO_ID = KO, PGPT = TRAIT) %>% 
  filter(!is.na(KO_ID)) %>% 
  separate_longer_delim(., PGPT, "#")

tb.ann.count.ko <- 
  base::merge(tb.ann.ko, tb.count.ko, by.x = "KO_ID", by.y = 0, all.y = T) %>% 
  rename_with(~ str_replace(., "^D(\\d+)$", "Level \\1"), starts_with("D"))

tb.ann.count.ko.filt <- 
  base::merge(tb.ann.ko.filt, tb.count.ko, by.x = "KO_ID", by.y = 0, all = F) %>% 
  rename_with(~ str_replace(., "^D(\\d+)$", "Level \\1"), starts_with("D"))

unique(tb.ann.count.ko.filt$`Level 1`)
unique(tb.ann.count.ko.filt$`Level 2`)
unique(tb.ann.count.ko.filt$`Level 3`)
unique(tb.ann.count.ko.filt$`Level 4`)

tb.ann.count.pb <- 
  base::merge(tb.ann.pb, tb.count.ko, by.x = "KO_ID", by.y = 0, all.y = T) %>% 
  separate_wider_delim(PGPT, delim = "; ", names = c("Level 1", "Level 2", "Level 3", "Level 4", "Level 5", "Level 6")) %>% 
  mutate(across(starts_with("Level"), ~ifelse(is.na(.x), "NON-PGPT", .x))) %>% 
  mutate(across(starts_with("Level"), ~str_replace_all(., "_", " ")))

unique(tb.ann.count.pb$`Level 1`)
unique(tb.ann.count.pb$`Level 2`)
unique(tb.ann.count.pb$`Level 3`)
unique(tb.ann.count.pb$`Level 4`)

# Graph Setup -------------------------------------------------------------

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

# PCA ---------------------------------------------------------------------

library("FactoMineR")
library("factoextra")

tb.count.ko.longer <-
  tb.count.ko %>% 
  t() %>% data.frame(check.names = F) %>% 
  rownames_to_column("SAMPLE_ID") %>% 
  merge(tb.metad, .) %>% 
  mutate(across(where(is.character), ~factor(., unique(.)))) %>% 
  pivot_longer(cols = -(names(tb.metad)), names_to = "KO_ID", values_to = "Value")

tb.count.ko
tb.count.ko.longer

mv.factor <-
  tb.count.ko.longer %>%
  select(SAMPLE_ID, Group, FullGroup) %>% 
  distinct() %>% 
  data.frame(check.names = F)

mv.matrix <-
  tb.count.ko.longer %>% 
  select(SAMPLE_ID, KO_ID, Value) %>% 
  distinct() %>% 
  pivot_wider(names_from = KO_ID, values_from = Value) %>% 
  data.frame(check.names = F) %>% 
  column_to_rownames("SAMPLE_ID")

tr.colors.1 <- c("#525252FF", "#8cc084FF", "#38a3a5FF", "#0466c8FF")
names(tr.colors.1) <- unique(mv.factor$Group)

tr.colors.2 <- c("#525252FF", "#8cc084FF", "#38a3a5FF", "#0466c8FF")
names(tr.colors.2) <- unique(mv.factor$FullGroup)

# PCA using FactoMineR
res.pca <- PCA(mv.matrix, graph = F, scale.unit = T)

# Screeplot of explained variances
fviz_screeplot(res.pca, addlabels = T, linecolor ="red", barfill="#495057", barcolor ="#000000") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), breaks = scales::extended_breaks(n = 6)) +
  scale_x_discrete(expand = expansion(mult = c(0.1, 0.1))) +
  my_theme -> pca.scree

pca.scree

fviz_contrib(res.pca, choice="var", axes = c(1,2),  top = 10, fill="#495057", color ="#000000") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), breaks = scales::extended_breaks(n = 6)) +
  scale_x_discrete(expand = expansion(mult = c(0.1, 0.1))) +
  my_theme -> pca.contr

pca.contr

# PCA - IND

fviz_pca_ind(res.pca, geom="point") +
  geom_point(size = 2.25, shape = 21, aes(fill = mv.factor$FullGroup), color = "#00000000", inherit.aes = T) +
  xlab(paste0("<br>PC1 (",round(res.pca$eig[1,3], 2), "%)")) +
  ylab(paste0("PC2 (",round(res.pca$eig[2,3]-res.pca$eig[1,3], 2),"%)<br>")) +
  scale_x_continuous(limits = ~ c(-1, 1) * max(abs(.x)), labels = scales::label_number(scale = .01, accuracy = 0.01)) +
  scale_y_continuous(limits = ~ c(-1, 1) * max(abs(.x)), labels = scales::label_number(scale = .01, accuracy = 0.01)) +
  scale_fill_manual(values = tr.colors.2) +
  labs(title = "KEGG orthologs") +
  guides(fill = 
           guide_legend(order = 1, 
                        override.aes = list(shape = 22, 
                                            fill = tr.colors, 
                                            color = "black", 
                                            size = 3,
                                            linetype = 0, 
                                            alpha = 1),
                        ncol = 1)) + 
  theme_bw() + 
  my_theme + 
  theme(legend.position = "right", 
        legend.box.background = element_rect(fill = "white", colour = "black", linewidth = .75),
        legend.justification = c(0.1, 0.5)) -> gg.pca.ko

gg.pca.ko

ggsave(filename = paste0(path, "figures/Functional/PCA_KOs.png"), plot = gg.pca.ko, device = "png", width = 5000 * 1.75, height = 5000 * 1.0, units = "px", dpi = 1200, limitsize = F)
ggsave(filename = paste0(path, "figures/Functional/PCA_KOs.svg"), plot = gg.pca.ko, device = "svg", width = 5000 * 1.75, height = 5000 * 1.0, units = "px", dpi = 1200, limitsize = F)

# Alluvial Plot -----------------------------------------------------------

tb.pb.alluv <- 
  tb.ann.count.pb %>%
  select(KO_ID, `Level 2`:`Level 3`, tb.metad$SAMPLE_ID) %>%
  distinct() %>% 
  select(-KO_ID) %>% 
  mutate(across(.cols = where(is.numeric), ~ ./sum(.x, na.rm = TRUE)))  %>% 
  pivot_longer(cols = tb.metad$SAMPLE_ID, names_to = "SAMPLE_ID", values_to = "Freq") %>%
  mutate(Group = str_remove_all(SAMPLE_ID, "\\.[:digit:]")) %>% 
  select(-SAMPLE_ID) %>% 
  group_by(`Level 2`, `Level 3`, Group) %>% 
  summarise(across(.cols = everything(), ~ sum(.x, na.rm = TRUE))) %>% 
  ungroup() %>% 
  mutate(`Level 2` = factor(`Level 2`, levels = unique(.$`Level 2`))) %>% 
  mutate(`Level 2` = fct_relevel(`Level 2`, "NON-PGPT", after = Inf)) %>% 
  mutate(cat = ifelse(
    `Level 2` %in% c("DIRECT EFFECTS", "INDIRECT EFFECTS"),
    paste0(.$`Level 2`, ": ", .$`Level 3`),
    paste0(.$`Level 2`))) %>% 
  arrange(Group, `Level 2`, -Freq) %>% 
  mutate(across(where(is.character), ~factor(., unique(.)))) %>% 
  mutate(`Level 3` = as.character(`Level 3`)) %>% 
  arrange(-Freq) %>% 
  mutate(`Level 3` = factor(`Level 3`, unique(.$`Level 3`))) %>% 
  arrange(Group, `Level 2`, -Freq)



colors.pb <- 
  c(rev(colorRampPalette(paletteer::paletteer_d("rcartocolor::Emrld")[-c(1)])(length(unique(filter(tb.pb.alluv, `Level 2` == "DIRECT EFFECTS") %>% pull(`Level 3`))))),
    rev(colorRampPalette(paletteer::paletteer_d("rcartocolor::SunsetDark")[-c(1)])(length(unique(filter(tb.pb.alluv, `Level 2` == "INDIRECT EFFECTS") %>% pull(`Level 3`))))),
    colorRampPalette(paletteer::paletteer_d("miscpalettes::grayscale"))(length(unique(filter(tb.pb.alluv, !`Level 2` %in% c("DIRECT EFFECTS", "INDIRECT EFFECTS")) %>% pull(`Level 3`)))))

names(colors.pb) <- 
  c(
    unique(filter(tb.pb.alluv, `Level 2` == "DIRECT EFFECTS") %>% pull(`cat`)),
    unique(filter(tb.pb.alluv, `Level 2` == "INDIRECT EFFECTS") %>% pull(`cat`)),
    unique(filter(tb.pb.alluv, !`Level 2` %in% c("DIRECT EFFECTS", "INDIRECT EFFECTS")) %>% pull(`cat`)))

gg.alluv <- 
  ggplot(data = tb.pb.alluv,
         aes(axis1 = `Level 2`, axis2 = `Level 3`, y = Freq/6, fill = cat)) +
  facet_wrap(.~Group, scales = "fixed", nrow = 1) + 
  geom_alluvium(alpha = 1) +
  geom_stratum(fill = "white", color = "#000000", alpha = 0) +
  #geom_richtext(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("Cat.", "Sub<br>cat."), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), position = "left", name = "Relative abundance (%)<br>",  
                     labels = scales::percent_format(accuracy = 1),
                     sec.axis = sec_axis(trans=~., 
                                         labels = scales::percent_format(accuracy = 1),
                                         name = "")) + 
  scale_fill_manual(values = colors.pb, name = NULL) +
  guides(fill = 
           guide_legend(order = 1, 
                        override.aes = list(shape = 16, 
                                            fill = colors.pb, 
                                            color = "black", 
                                            size = 1,
                                            linetype = 0),
                        ncol = 1)) +
  theme_bw() +
  my_theme + 
  theme(legend.position = "none")

gg.alluv <- 
  gg.alluv +
  plot_grid(
    ggpubr::get_legend(addSmallLegend(gg.alluv + theme(legend.position = "right")) + 
                         theme(legend.justification = c(0, 0.5), 
                               legend.box.margin = margin(0, 0, 0, 0),
                               legend.box.background = element_rect(fill = "white", colour = "black", linewidth = .75))),
    ncol = 2, align = "hv") + plot_layout(widths = c(1, 1))

gg.alluv

ggsave(filename = paste0(path, "figures/Functional/ProfilePLABASE1.png"), plot = gg.alluv, width = 5000*3.5, height = 5000 * 1.05, units = 'px', dpi = 1200)
ggsave(filename = paste0(path, "figures/Functional/ProfilePLABASE1.svg"), plot = gg.alluv, width = 5000*3.5, height = 5000 * 1.05, units = 'px', dpi = 1200)

#

tb.ko.alluv <- 
  tb.ann.count.ko.filt %>%
  select(KO_ID, `Level 1`:`Level 2`, tb.metad$SAMPLE_ID) %>%
  distinct() %>% 
  select(-KO_ID) %>% 
  mutate(across(.cols = where(is.numeric), ~ ./sum(.x, na.rm = TRUE)))  %>% 
  mutate(across(.cols = c(`Level 1`:`Level 2`), ~ toupper(.)))  %>% 
  pivot_longer(cols = tb.metad$SAMPLE_ID, names_to = "SAMPLE_ID", values_to = "Freq") %>%
  mutate(Group = str_remove_all(SAMPLE_ID, "\\.[:digit:]")) %>% 
  select(-SAMPLE_ID) %>% 
  group_by(`Level 1`, `Level 2`, Group) %>% 
  summarise(across(.cols = everything(), ~ sum(.x, na.rm = TRUE))) %>% 
  ungroup() %>% 
  mutate(`Level 1` = factor(`Level 1`, levels = c("METABOLISM", "ENVIRONMENTAL INFORMATION PROCESSING", "GENETIC INFORMATION PROCESSING", "CELLULAR PROCESSES"))) %>% 
  mutate(cat = paste0(.$`Level 1`, ": ", .$`Level 2`)) %>% 
  arrange(Group, `Level 1`, -Freq) %>% 
  mutate(across(where(is.character), ~factor(., unique(.)))) %>% 
  mutate(`Level 2` = as.character(`Level 2`)) %>% 
  arrange(-Freq) %>% 
  mutate(`Level 2` = factor(`Level 2`, unique(.$`Level 2`))) %>% 
  arrange(Group, `Level 1`, -Freq)



colors.ko <- 
  c(
    rev(colorRampPalette(c(paletteer::paletteer_d("RColorBrewer::YlOrRd"), "#541F3FFF"))(length(unique(filter(tb.ko.alluv, `Level 1` == "METABOLISM") %>% pull(`Level 2`))))),
    rev(colorRampPalette(paletteer::paletteer_d("RColorBrewer::Blues")[-c(1:2)])(length(unique(filter(tb.ko.alluv, `Level 1` == "ENVIRONMENTAL INFORMATION PROCESSING") %>% pull(`Level 2`))))),
    rev(colorRampPalette(paletteer::paletteer_d("RColorBrewer::Greens")[-c(1:2)])(length(unique(filter(tb.ko.alluv, `Level 1` == "GENETIC INFORMATION PROCESSING") %>% pull(`Level 2`))))),
    rev(colorRampPalette(paletteer::paletteer_d("RColorBrewer::Purples")[-c(1:2)])(length(unique(filter(tb.ko.alluv, `Level 1` == "CELLULAR PROCESSES") %>% pull(`Level 2`))))))

names(colors.ko) <- 
  c(
    unique(filter(tb.ko.alluv, `Level 1` == "METABOLISM") %>% pull(`cat`)),
    unique(filter(tb.ko.alluv, `Level 1` == "ENVIRONMENTAL INFORMATION PROCESSING") %>% pull(`cat`)),
    unique(filter(tb.ko.alluv, `Level 1` == "GENETIC INFORMATION PROCESSING") %>% pull(`cat`)),
    unique(filter(tb.ko.alluv, `Level 1` == "CELLULAR PROCESSES") %>% pull(`cat`)))

tb.ko.alluv

gg.alluv <- 
  ggplot(data = tb.ko.alluv,
         aes(axis1 = `Level 1`, axis2 = `Level 2`, y = Freq/6, fill = cat)) +
  facet_wrap(.~Group, scales = "fixed", nrow = 1) + 
  geom_alluvium(alpha = 1) +
  geom_stratum(fill = "white", color = "#000000", alpha = 0) +
  #geom_richtext(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("Cat.", "Sub<br>cat."), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), position = "left", name = "Relative abundance (%)<br>",  
                     labels = scales::percent_format(accuracy = 1),
                     sec.axis = sec_axis(trans=~., 
                                         labels = scales::percent_format(accuracy = 1),
                                         name = "")) + 
  scale_fill_manual(values = colors.ko, name = NULL) +
  guides(fill = 
           guide_legend(order = 1, 
                        override.aes = list(shape = 16, 
                                            fill = colors.ko, 
                                            color = "black", 
                                            size = 1,
                                            linetype = 0),
                        ncol = 1)) +
  theme_bw() +
  my_theme + 
  theme(legend.position = "none")

gg.alluv <- 
  gg.alluv +
  plot_grid(
    ggpubr::get_legend(addSmallLegend(gg.alluv + theme(legend.position = "right")) + 
                         theme(legend.justification = c(0, 0.5), 
                               legend.box.margin = margin(0, 0, 0, 0),
                               legend.box.background = element_rect(fill = "white", colour = "black", linewidth = .75))),
    ncol = 2, align = "hv") + plot_layout(widths = c(1, 1))

gg.alluv

ggsave(filename = paste0(path, "figures/Functional/ProfileKO1.png"), plot = gg.alluv, width = 5000*3.5, height = 5000 * 1.05, units = 'px', dpi = 1200)
ggsave(filename = paste0(path, "figures/Functional/ProfileKO1.svg"), plot = gg.alluv, width = 5000*3.5, height = 5000 * 1.05, units = 'px', dpi = 1200)

# Heatmap -----------------------------------------------------------------

library("ComplexHeatmap")

hp.top.tb <- 
  tb.ann.count.pb %>% 
  select(KO_ID, `Level 2`:`Level 4`, tb.metad$SAMPLE_ID) %>%
  distinct() %>% 
  select(-KO_ID) %>% 
  group_by(across(`Level 2`:`Level 4`)) %>% 
  summarise(across(.cols = everything(), ~ sum(.x, na.rm = TRUE))) %>% 
  ungroup() %>% 
  mutate(across(.cols = where(is.numeric), ~ ./sum(.x, na.rm = TRUE))) %>% 
  mutate(cat = ifelse(
    `Level 2` %in% c("DIRECT EFFECTS", "INDIRECT EFFECTS"),
    paste0(.$`Level 2`, ": ", .$`Level 3`),
    paste0(.$`Level 2`))) %>% 
  select(cat, `Level 2`:`Level 4`, tb.metad$SAMPLE_ID) %>% 
  mutate(mean = rowMeans(select(., where(is.numeric)))) 

hp.anno.col <- tb.metad %>% select("Treatment        " = Group)

hp.anno.col.col <- c("#525252FF", "#8cc084FF", "#38a3a5FF", "#0466c8FF")
names(hp.anno.col.col) <- unique(hp.anno.col$`Treatment        `)

hp.anno.row <- 
  hp.top.tb %>% 
  select(`Level 4`, `PGPT Class` = cat, `Rel. abundance (%)` = mean) %>% 
  mutate(`Rel. abundance (%)` = `Rel. abundance (%)` * 100) %>% 
  mutate(`Level 4` = paste0(`Level 4`, "")) %>% 
  column_to_rownames("Level 4")

hp.anno.row.col <- colors.pb

hp.values <- 
  hp.top.tb %>% 
  select(`Level 4`, tb.metad$SAMPLE_ID) %>% 
  mutate(`Level 4` = paste0(`Level 4`, "                               ")) %>% 
  column_to_rownames("Level 4") %>%
  t() %>% 
  scale() %>% 
  t()

hp.breaks <- max(abs(min(hp.values, na.rm = T)), max(hp.values, na.rm = T))



col_fun1 = 
  colorRamp2::colorRamp2(c(seq(-hp.breaks, hp.breaks, .25)),
                         colorRampPalette(c("#364B9AFF", "#4A7BB7FF", "#6EA6CDFF", "#98CAE1FF", "#FFFFFFFF", "#FDB366FF", "#F67E4BFF", "#DD3D2DFF", "#A50026FF"))(length(seq(-hp.breaks, hp.breaks, .25))))

hp.max.rel <- max(hp.anno.row$`Rel. abundance (%)`, na.rm = T)

col_fun2 = 
  colorRamp2::colorRamp2(c(seq(1, hp.max.rel, 10)),
                         colorRampPalette(c("#FFFFFFFF",  "#000000FF"))(length(seq(1, hp.max.rel, 10))))


ha = ComplexHeatmap::rowAnnotation(`Rel. abundance (%)` = hp.anno.row$`Rel. abundance (%)`,
                                   `PGPT Class` = hp.anno.row$`PGPT Class`,
                                   col = list(`Rel. abundance (%)` = col_fun2,
                                              `PGPT Class` = hp.anno.row.col), 
                                   annotation_legend_param = list(`PGPT Class` = 
                                                                    list(
                                                                      labels = names(hp.anno.row.col),
                                                                      at = names(hp.anno.row.col))),
                                   
                                   show_annotation_name = FALSE)

hb = ComplexHeatmap::columnAnnotation(`Treatment        ` = hp.anno.col$`Treatment        `, 
                                      col = list(`Treatment        ` = hp.anno.col.col), show_annotation_name = FALSE)

# Plot
ComplexHeatmap::Heatmap(hp.values, 
                        name = "Z-Score        ",  
                        col = col_fun1,
                        cluster_rows = T, 
                        cluster_columns = T,
                        row_names_gp = gpar(fontsize = 8),
                        column_names_gp = gpar(fontsize = 8), 
                        left_annotation = ha,
                        top_annotation = hb) -> hp.plot

hp.plot

png(filename = paste0(path, "figures/Functional/Heatmap_PLABASE.png"), width = 5000 * 4, height = 5000 * 2, units = "px", res = 1200, type = "cairo")

draw(hp.plot, heatmap_legend_side = "left", annotation_legend_side = "left", padding = unit(c(1, 1, 1, 100), "mm"))

dev.off()

# 

hp.top.tb <- 
  tb.ann.count.ko.filt %>% 
  select(KO_ID, `Level 1`:`Level 2`, tb.metad$SAMPLE_ID) %>%
  mutate(across(.cols = c(`Level 1`:`Level 2`), ~ toupper(.)))  %>% 
  distinct() %>% 
  select(-KO_ID) %>% 
  group_by(across(`Level 1`:`Level 2`)) %>% 
  summarise(across(.cols = everything(), ~ sum(.x, na.rm = TRUE))) %>% 
  ungroup() %>% 
  mutate(across(.cols = where(is.numeric), ~ ./sum(.x, na.rm = TRUE))) %>% 
  mutate(`Level 1` = factor(`Level 1`, levels = c("METABOLISM", "ENVIRONMENTAL INFORMATION PROCESSING", "GENETIC INFORMATION PROCESSING", "CELLULAR PROCESSES"))) %>% 
  mutate(cat = paste0(.$`Level 1`, ": ", .$`Level 2`)) %>% 
  select(cat, `Level 1`:`Level 2`, tb.metad$SAMPLE_ID) %>% 
  mutate(mean = rowMeans(select(., where(is.numeric)))) 

hp.anno.col <- tb.metad %>% select("Treatment        " = Group)

hp.anno.col.col <- c("#525252FF", "#8cc084FF", "#38a3a5FF", "#0466c8FF")
names(hp.anno.col.col) <- unique(hp.anno.col$`Treatment        `)

hp.anno.row <- 
  hp.top.tb %>% 
  select(`Level 1` , `Level 2`, `KEGG Class` = cat, `Rel. abundance (%)` = mean) %>% 
  mutate(`Rel. abundance (%)` = `Rel. abundance (%)` * 100) %>% 
  mutate(`Level 2` = paste0(`Level 2`, "")) %>% 
  column_to_rownames("Level 2")

hp.anno.row.col <- c("#D73027FF", "#2166ACFF", "#1A9850FF", "#7D4F73FF")
names(hp.anno.row.col) <- levels(hp.top.tb$`Level 1`)

hp.values <- 
  hp.top.tb %>% 
  select(`Level 2`, tb.metad$SAMPLE_ID) %>% 
  mutate(`Level 2` = paste0(`Level 2`, "                               ")) %>% 
  column_to_rownames("Level 2") %>%
  t() %>% 
  scale() %>% 
  t()

hp.breaks <- max(abs(min(hp.values, na.rm = T)), max(hp.values, na.rm = T))

col_fun1 = 
  colorRamp2::colorRamp2(c(seq(-hp.breaks, hp.breaks, .25)),
                         colorRampPalette(c("#364B9AFF", "#4A7BB7FF", "#6EA6CDFF", "#98CAE1FF", "#FFFFFFFF", "#FDB366FF", "#F67E4BFF", "#DD3D2DFF", "#A50026FF"))(length(seq(-hp.breaks, hp.breaks, .25))))

hp.max.rel <- max(hp.anno.row$`Rel. abundance (%)`, na.rm = T)

col_fun2 = 
  colorRamp2::colorRamp2(c(seq(1, hp.max.rel, 10)),
                         colorRampPalette(c("#FFFFFFFF",  "#000000FF"))(length(seq(1, hp.max.rel, 10))))


ha = ComplexHeatmap::rowAnnotation(`Rel. abundance (%)` = hp.anno.row$`Rel. abundance (%)`,
                                   `KEGG Class` = hp.anno.row$`Level 1`,
                                   col = list(`Rel. abundance (%)` = col_fun2,
                                              `KEGG Class` = hp.anno.row.col), 
                                   annotation_legend_param = list(`KEGG Class` = 
                                                                    list(
                                                                      labels = names(hp.anno.row.col),
                                                                      at = names(hp.anno.row.col))),
                                   
                                   show_annotation_name = FALSE)

hb = ComplexHeatmap::columnAnnotation(`Treatment        ` = hp.anno.col$`Treatment        `, 
                                      col = list(`Treatment        ` = hp.anno.col.col), show_annotation_name = FALSE)

# Plot
ComplexHeatmap::Heatmap(hp.values, 
                        name = "Z-Score        ",  
                        col = col_fun1,
                        cluster_rows = T, 
                        cluster_columns = T,
                        row_names_gp = gpar(fontsize = 8),
                        column_names_gp = gpar(fontsize = 8), 
                        left_annotation = ha,
                        top_annotation = hb) -> hp.plot

hp.plot

png(filename = paste0(path, "figures/Functional/Heatmap_KO.png"), width = 5000 * 3.65, height = 5000 * 1.15, units = "px", res = 1200, type = "cairo")

draw(hp.plot, heatmap_legend_side = "left", annotation_legend_side = "left", padding = unit(c(1, 1, 1, 100), "mm"))

dev.off()

# "Phyloseqfying" data ---------------------------------------------------

library("phyloseq")
library("metagMisc")
library("microbiome")
library("ggvenn")


tb.metad.ps <- tb.metad %>% dplyr::rename(ID = SAMPLE_ID)

tb.count.ps <- 
  tb.count.ko %>% 
  rownames_to_column("ID") %>% 
  select(ID, tb.metad.ps$ID) %>% 
  distinct()

tb.taxon.ps <- 
  tb.count.ps %>% 
  select(ID) %>% 
  distinct() %>% 
  mutate("Kingdom" = ID, 
         "Phylum" = ID, 
         "Class" = ID, 
         "Order" = ID, 
         "Family" = ID, 
         "Genus" = ID, 
         "Species" = ID)

row.names(tb.count.ps) <- tb.count.ps$ID
row.names(tb.taxon.ps) <- tb.taxon.ps$ID
row.names(tb.metad.ps) <- tb.metad.ps$ID

tb.count.ps <- tb.count.ps %>% select(-ID)
tb.taxon.ps <- tb.taxon.ps %>% select(-ID)

tb.count.ps <- as.matrix(tb.count.ps)
tb.taxon.ps <- as.matrix(tb.taxon.ps)

ps.count <- otu_table(tb.count.ps, taxa_are_rows = TRUE)
ps.taxon <- tax_table(tb.taxon.ps)
ps.metad <- sample_data(tb.metad.ps)

ps.table.ko <- phyloseq(ps.count, ps.taxon, ps.metad)
ps.table.ko.tss <- phyloseq_standardize_otu_abundance(physeq = ps.table.ko, method = 'total')

# Alpha

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

gg.venn.a <- fun.venn(ps.table.ko, "KEGG Orthologs")

ggsave(filename = paste0(path, "figures/Functional/VennDiagram.png"), plot = gg.venn.a, bg = "white", dpi = 1200, width = 5000*1, height = 5000*1, units = "px")
ggsave(filename = paste0(path, "figures/Functional/VennDiagram.svg"), plot = gg.venn.a, bg = "white", dpi = 1200, width = 5000*1, height = 5000*1, units = "px")

# Alpha diversity ---------------------------------------------------------

adiv.table <- microbiome::alpha(ps.table.ko, index = c("observed", "shannon"))
names(adiv.table) <- c("Richness", "Shannon's diversity")

adiv.table <- merge(sample_data(ps.table.ko), adiv.table, by = "row.names", sort = F)[-1]

adiv.table.longer <- 
  adiv.table %>% 
  mutate(across(where(is.character), ~factor(., unique(.)))) %>% 
  pivot_longer(cols = Richness:`Shannon's diversity`, names_to = "index", values_to = "values") %>% 
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
                        ncol = 1)) +
  my_theme + 
  theme(legend.position = "right", 
        legend.box.background = element_rect(fill = "white", colour = "black", linewidth = .75),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

adiv.gg

ggsave(filename = paste0(path, "figures/Functional/AlphaDiv.png"), plot = adiv.gg, dpi = 1200, width = 5000*1.75, height = 5000*.75, units = "px")
ggsave(filename = paste0(path, "figures/Functional/AlphaDiv.svg"), plot = adiv.gg, dpi = 1200, width = 5000*1.75, height = 5000*.75, units = "px")

ls.adiv <- 
  
  list(
    
    "Alpha (Measures)" = 
      adiv.table %>% 
      select(ID, Group, Richness:`Shannon's diversity`) %>%
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

openxlsx::write.xlsx(x = ls.adiv, file = paste0(path, "tables/FuncAlphaDiv.xlsx"), rowNames = F)

# DAA ---------------------------------------------------------------------

daa.fun <- function(ps.table.subset) {
  
  tbDESeq2final = data.frame(check.names = F)
  
  tbDESeq2tmp = ps.table.subset
  
  g1 = sort(unique(sample_data(tbDESeq2tmp)$Group))[1]
  g2 = sort(unique(sample_data(tbDESeq2tmp)$Group))[2] 
  
  DESeq2tmp = phyloseq_to_deseq2(tbDESeq2tmp, ~ Group)
  DESeq2tmp = DESeq(DESeq2tmp, test = 'Wald', fitType="mean")
  
  DESeq2tmp.res = results(DESeq2tmp, cooksCutoff = FALSE)
  DESeq2tmp.sig = DESeq2tmp.res
  
  if (DESeq2tmp.sig@nrows == 0) { 
    
    next 
    
  } else {
    
    DESeq2tmp.sig2 = DESeq2tmp.sig %>% 
      as.data.frame(check.names = F) %>% 
      dplyr::mutate("Group" = if_else(log2FoldChange > 0, g2, g1)) %>% 
      dplyr::mutate("ID" = rownames(.))
    
    tbDESeq2tmp = as.data.frame(check.names = F, tax_table(tbDESeq2tmp)) %>% dplyr::mutate(ID = row.names(.))
    
    DESeq2tmp.sig2 = merge(DESeq2tmp.sig2, tbDESeq2tmp, by = "ID", all.x = T)
    
    DESeq2tmp.sig2 = DESeq2tmp.sig2 %>% 
      dplyr::rename(., "FunID" = 1) %>% 
      dplyr::mutate("comp" = paste0(g1, " × ", g2))
    
    rownames(DESeq2tmp.sig2) = NULL
    
    tbDESeq2final = rbind(tbDESeq2final, DESeq2tmp.sig2 %>% 
                            select(FunID, Group, comp, log2FoldChange:padj)) %>% 
      arrange(-log2FoldChange)
    
  }
  
  return(tbDESeq2final %>% mutate(padj = p.adjust(pvalue, "fdr")))
  
}

daa.ps.cont2.1 <- merge_phyloseq(subset_samples(ps.table.ko, Group == "CTRL"), subset_samples(ps.table.ko, Group == "SM05"))
daa.ps.cont2.2 <- merge_phyloseq(subset_samples(ps.table.ko, Group == "CTRL"), subset_samples(ps.table.ko, Group == "SM10"))
daa.ps.cont2.3 <- merge_phyloseq(subset_samples(ps.table.ko, Group == "CTRL"), subset_samples(ps.table.ko, Group == "SM20"))

sample_data(daa.ps.cont2.1) <- data.frame(check.names = F, sample_data(daa.ps.cont2.1)) %>% mutate(Group = ifelse(Group == "SM05", 2, 1)) %>% arrange(Group)
sample_data(daa.ps.cont2.2) <- data.frame(check.names = F, sample_data(daa.ps.cont2.2)) %>% mutate(Group = ifelse(Group == "SM10", 2, 1)) %>% arrange(Group)
sample_data(daa.ps.cont2.3) <- data.frame(check.names = F, sample_data(daa.ps.cont2.3)) %>% mutate(Group = ifelse(Group == "SM20", 2, 1)) %>% arrange(Group)

daa.cont2.1 <- daa.fun(daa.ps.cont2.1) %>% mutate(Group = ifelse(Group == 2, "Enriched", "Depleted"), comp = "SM05") %>% select(FunID, Group, comp, everything())
daa.cont2.2 <- daa.fun(daa.ps.cont2.2) %>% mutate(Group = ifelse(Group == 2, "Enriched", "Depleted"), comp = "SM10") %>% select(FunID, Group, comp, everything())
daa.cont2.3 <- daa.fun(daa.ps.cont2.3) %>% mutate(Group = ifelse(Group == 2, "Enriched", "Depleted"), comp = "SM20") %>% select(FunID, Group, comp, everything())

daa.cont <- 
  rbind(daa.cont2.1, daa.cont2.2, daa.cont2.3) %>% 
  mutate(psignif = 
           case_when(
             pvalue < .001 & abs(log2FoldChange) > log2(2) ~ "***", 
             pvalue < .01 & abs(log2FoldChange) > log2(2) ~ "**", 
             pvalue < .05 & abs(log2FoldChange) > log2(2) ~ "*", 
             pvalue >= .05 | abs(log2FoldChange) < log2(2) ~ NA_character_)) %>% 
  mutate(pgroup = 
           case_when(
             is.na(psignif) ~ "NS",
             !is.na(psignif) & log2FoldChange > 0 ~ comp,
             !is.na(psignif) & log2FoldChange < 0 ~ "CTRL"))

daa.glom.gp.fun <- function(my.ps.table, group) {
  
  return(
    
    merge(data.frame(check.names = F, tax_table(phyloseq_standardize_otu_abundance(merge_samples(my.ps.table, group, fun=sum))))["Species"],
          data.frame(check.names = F, t(otu_table(phyloseq_standardize_otu_abundance(merge_samples(my.ps.table, group, fun=sum))))), by = 0)[-1] %>% 
      dplyr::mutate(FunID = Species) %>% select(FunID, unique(tb.metad.ps[[group]]))
    
  )
  
}

daa.ko.glom0 <- daa.glom.gp.fun(ps.table.ko, "ID")
daa.ko.glom1 <- daa.glom.gp.fun(ps.table.ko, "Group")

daa.ko.glom.id <- merge(daa.cont, daa.ko.glom0) %>% arrange(comp, log2FoldChange)
daa.ko.glom.group <- merge(daa.cont, daa.ko.glom1) %>% arrange(comp, log2FoldChange)

# Graphs

# VP

tr.colors.1
vp.cols <- c("NS" = "grey", tr.colors.1) 


vp.labels <-
  daa.ko.glom.group %>% 
  filter(pgroup != "NS") %>% 
  select(comp, pgroup) %>% 
  table() %>% 
  data.frame() %>% 
  filter(Freq != 0) %>% 
  mutate(my.x = ifelse(pgroup == "CTRL", -3, 3))

gg.vp <-
  ggplot(data = daa.ko.glom.group %>% filter(!is.na(Group))) + 
  facet_grid(. ~ comp) + 
  geom_point(aes(x = ifelse(log2FoldChange > 0, sqrt(log2FoldChange), -sqrt(abs(log2FoldChange))), y = -log10(pvalue), color = pgroup), shape = 16, alpha = .75, size = 1.5) + 
  geom_hline(yintercept = -log10(0.05), color = "black", linetype = "dashed", size = 0.25) + 
  geom_vline(xintercept = c(-1, 1), color = "black", linetype = "dashed", size = 0.25) +
  geom_richtext(mapping = aes(label = str_glue("**{Freq}**"), x = my.x, y = 20, color = pgroup), 
                data = vp.labels, size = 3, vjust = 1, hjust = 0.5, 
                inherit.aes = T, label.padding = unit(0, "pt"), fill = NA, label.color = NA) + 
  scale_color_manual(values = vp.cols) + 
  scale_x_continuous(limits = ~ c(-1, 1) * max(abs(.x)), labels = function(x) x^2, breaks = scales::breaks_extended(n = 7)) +
  theme_bw(base_family = "Liberation Sans") + 
  ylab("-Log₁₀ (P-value)<br>") +
  xlab("<br>|Log₂ (Fold-Change)|") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), title = ""), alpha = "none", size = "none") +
  theme_bw() + 
  my_theme + 
  theme(legend.position = "none", 
        axis.title.x = element_markdown(size = 10, face = 'bold', colour = 'black'),
        axis.title.y = element_markdown(size = 10, face = 'bold', colour = 'black'))

gg.vp 

gg.vp + 
  ggpubr::get_legend(
    ggplot(data = tb.metad.ps, aes(x = 1,y = 1, color = FullGroup)) +
      geom_point() +   
      scale_color_manual(values = tr.colors.2, name = NULL) + 
      guides(color = 
               guide_legend(order = 1, 
                            override.aes = list(shape = 22, 
                                                fill = tr.colors, 
                                                color = "black", 
                                                size = 3,
                                                linetype = 0),
                            ncol = 2)) +
      theme_bw() + 
      my_theme + 
      theme(legend.justification = c(.6, .9), 
            legend.box.margin = margin(0, 0, 0, 0),
            legend.box.background = element_rect(fill = "white", colour = "black", linewidth = .75))) +
  plot_layout(ncol = 1, heights = c(1, .5)) -> gg.vp.leg

gg.vp.leg

ggsave(filename = paste0(path, "figures/Functional/VolcanoPlotDEF.png"), plot = gg.vp.leg, dpi = 1200, width = 5000*1.75, height = 5000*1, units = "px")
ggsave(filename = paste0(path, "figures/Functional/VolcanoPlotDEF.svg"), plot = gg.vp.leg, dpi = 1200, width = 5000*1.75, height = 5000*1, units = "px")

# Cat DAA -----------------------------------------------------------------

# PLABASE

ls.alluv.tb <- list()

for (i in c("SM05", "SM10", "SM20")) {
  
  tb.pb.alluv.all.tmp <- 
    tb.ann.count.pb %>%
    select(KO_ID, `Level 2`:`Level 3`, starts_with(i)) %>%
    distinct()
  
  tb.pb.alluv.all.tmp %>% 
    select(-KO_ID) %>% 
    mutate(across(.cols = where(is.numeric), ~ ./sum(.x, na.rm = TRUE)))  %>% 
    pivot_longer(cols = -c(`Level 2`, `Level 3`), names_to = "SAMPLE_ID", values_to = "Freq") %>%
    mutate(Group = str_remove_all(SAMPLE_ID, "\\.[:digit:]")) %>% 
    select(-SAMPLE_ID) %>% 
    group_by(`Level 2`, `Level 3`, Group) %>% 
    summarise(across(.cols = everything(), ~ sum(.x, na.rm = TRUE))) %>% 
    ungroup() %>% 
    mutate(`Level 2` = factor(`Level 2`, levels = unique(.$`Level 2`))) %>% 
    mutate(`Level 2` = fct_relevel(`Level 2`, "NON-PGPT", after = Inf)) %>% 
    mutate(cat = ifelse(
      `Level 2` %in% c("DIRECT EFFECTS", "INDIRECT EFFECTS"),
      paste0(.$`Level 2`, ": ", .$`Level 3`),
      paste0(.$`Level 2`))) %>% 
    arrange(Group, `Level 2`, -Freq) %>% 
    mutate(across(where(is.character), ~factor(., unique(.)))) %>% 
    mutate(`Level 3` = as.character(`Level 3`)) %>% 
    arrange(-Freq) %>% 
    mutate(`Level 3` = factor(`Level 3`, unique(.$`Level 3`))) %>% 
    arrange(Group, `Level 2`, -Freq) %>% 
    mutate(status = "All") -> tb.pb.alluv.all
  
  # DR
  
  my.dr.pb <- 
    daa.ko.glom.id %>% 
    filter(pgroup != "NS" & comp == i & Group == "Depleted") %>% 
    select(KO_ID = FunID, Group, comp, starts_with(i)) 
  
  tb.pb.alluv.all.tmp %>% 
    filter(KO_ID %in% my.dr.pb$KO_ID) %>% 
    select(-KO_ID) %>% 
    mutate(across(.cols = where(is.numeric), ~ ./sum(.x, na.rm = TRUE)))  %>% 
    pivot_longer(cols = -c(`Level 2`, `Level 3`), names_to = "SAMPLE_ID", values_to = "Freq") %>%
    mutate(Group = str_remove_all(SAMPLE_ID, "\\.[:digit:]")) %>% 
    select(-SAMPLE_ID) %>% 
    group_by(`Level 2`, `Level 3`, Group) %>% 
    summarise(across(.cols = everything(), ~ sum(.x, na.rm = TRUE))) %>% 
    ungroup() %>% 
    mutate(`Level 2` = factor(`Level 2`, levels = unique(.$`Level 2`))) %>% 
    mutate(`Level 2` = fct_relevel(`Level 2`, "NON-PGPT", after = Inf)) %>% 
    mutate(cat = ifelse(
      `Level 2` %in% c("DIRECT EFFECTS", "INDIRECT EFFECTS"),
      paste0(.$`Level 2`, ": ", .$`Level 3`),
      paste0(.$`Level 2`))) %>% 
    arrange(Group, `Level 2`, -Freq) %>% 
    mutate(across(where(is.character), ~factor(., unique(.)))) %>% 
    mutate(`Level 3` = as.character(`Level 3`)) %>% 
    arrange(-Freq) %>% 
    mutate(`Level 3` = factor(`Level 3`, unique(.$`Level 3`))) %>% 
    arrange(Group, `Level 2`, -Freq) %>% 
    mutate(status = "Depleted") -> tb.pb.alluv.dep
  
  # UP
  
  my.dr.ur <- 
    daa.ko.glom.id %>% 
    filter(pgroup != "NS" & comp == i & Group == "Enriched") %>% 
    select(KO_ID = FunID, Group, comp, starts_with(i)) 
  
  tb.pb.alluv.all.tmp %>% 
    filter(KO_ID %in% my.dr.ur$KO_ID) %>% 
    select(-KO_ID) %>% 
    mutate(across(.cols = where(is.numeric), ~ ./sum(.x, na.rm = TRUE)))  %>% 
    pivot_longer(cols = -c(`Level 2`, `Level 3`), names_to = "SAMPLE_ID", values_to = "Freq") %>%
    mutate(Group = str_remove_all(SAMPLE_ID, "\\.[:digit:]")) %>% 
    select(-SAMPLE_ID) %>% 
    group_by(`Level 2`, `Level 3`, Group) %>% 
    summarise(across(.cols = everything(), ~ sum(.x, na.rm = TRUE))) %>% 
    ungroup() %>% 
    mutate(`Level 2` = factor(`Level 2`, levels = unique(.$`Level 2`))) %>% 
    mutate(`Level 2` = fct_relevel(`Level 2`, "NON-PGPT", after = Inf)) %>% 
    mutate(cat = ifelse(
      `Level 2` %in% c("DIRECT EFFECTS", "INDIRECT EFFECTS"),
      paste0(.$`Level 2`, ": ", .$`Level 3`),
      paste0(.$`Level 2`))) %>% 
    arrange(Group, `Level 2`, -Freq) %>% 
    mutate(across(where(is.character), ~factor(., unique(.)))) %>% 
    mutate(`Level 3` = as.character(`Level 3`)) %>% 
    arrange(-Freq) %>% 
    mutate(`Level 3` = factor(`Level 3`, unique(.$`Level 3`))) %>% 
    arrange(Group, `Level 2`, -Freq) %>% 
    mutate(status = "Enriched") -> tb.pb.alluv.enr
  
  rbind(
    tb.pb.alluv.dep,
    tb.pb.alluv.all, 
    tb.pb.alluv.enr
  ) %>% 
    mutate(status = factor(status, levels = unique(.$status))) %>% 
    mutate(cat = factor(cat, levels = names(colors.pb))) -> tb.pb.alluv.daa
  
  
  ls.alluv.tb[[i]] <- tb.pb.alluv.daa
  
}

daa.pb.summary <- rbind(ls.alluv.tb[["SM05"]], ls.alluv.tb[["SM10"]], ls.alluv.tb[["SM20"]]) %>% mutate(Group = factor(Group, c("SM05", "SM10", "SM20")))

ggplot(daa.pb.summary,
       aes(x = status, stratum = `Level 3`, alluvium = `Level 3`,
           y = Freq/6, fill = cat)) +
  facet_wrap(.~Group, scales = "fixed", nrow = 1) + 
  geom_stratum(color = "#000000", alpha = 1, width = .333) +
  geom_flow() +
  scale_x_discrete(expand = c(0, 0), labels = c("<span style='color:#364B9AFF'>Depleted<br>Orthologs</span>", 
                                                "<span style='color:black'>All<br>Orthologs</span>", 
                                                "<span style='color:#A50026FF'>Enriched<br>Orthologs</span>")) +
  scale_y_continuous(expand = c(0, 0), position = "left", name = "Relative abundance (%)<br>",  
                     labels = scales::percent_format(accuracy = 1),
                     sec.axis = sec_axis(trans=~., 
                                         labels = scales::percent_format(accuracy = 1),
                                         name = "")) + 
  scale_fill_manual(values = colors.pb, name = NULL) +
  guides(fill = 
           guide_legend(order = 1, 
                        override.aes = list(shape = 16, 
                                            fill = colors.pb, 
                                            color = "black", 
                                            size = 1,
                                            linetype = 0),
                        ncol = 1)) +
  theme_bw() +
  my_theme + 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_markdown(size = 7, face = "bold", margin = unit(c(5, 0, 0, 0), "pt")), 
        panel.spacing.x = unit(25, "pt"),
        axis.ticks.x = element_blank()) -> gg.alluv.tmp.1

gg.alluv.tmp.1 <- 
  gg.alluv.tmp.1 +
  plot_grid(
    ggpubr::get_legend(addSmallLegend(gg.alluv.tmp.1 + theme(legend.position = "right")) + 
                         theme(legend.justification = c(0, 0.5), 
                               legend.box.margin = margin(0, 0, 0, 0),
                               legend.box.background = element_rect(fill = "white", colour = "black", linewidth = .75))),
    ncol = 2, align = "hv") + plot_layout(widths = c(1, 1))

# KOs

ls.alluv.tb <- list()

for (i in c("SM05", "SM10", "SM20")) {
  
  tb.pb.alluv.all.tmp <- 
    tb.ann.count.ko.filt %>%
    select(KO_ID, `Level 1`:`Level 2`, starts_with(i)) %>%
    distinct()
  
  tb.pb.alluv.all.tmp %>% 
    select(-KO_ID) %>% 
    mutate(across(.cols = where(is.numeric), ~ ./sum(.x, na.rm = TRUE)))  %>% 
    mutate(across(.cols = c(`Level 1`:`Level 2`), ~ toupper(.)))  %>% 
    pivot_longer(cols = -c(`Level 1`, `Level 2`), names_to = "SAMPLE_ID", values_to = "Freq") %>%
    mutate(Group = str_remove_all(SAMPLE_ID, "\\.[:digit:]")) %>% 
    select(-SAMPLE_ID) %>% 
    group_by(`Level 1`, `Level 2`, Group) %>% 
    summarise(across(.cols = everything(), ~ sum(.x, na.rm = TRUE))) %>% 
    ungroup() %>% 
    mutate(`Level 1` = factor(`Level 1`, levels = c("METABOLISM", "ENVIRONMENTAL INFORMATION PROCESSING", "GENETIC INFORMATION PROCESSING", "CELLULAR PROCESSES"))) %>% 
    mutate(cat = paste0(.$`Level 1`, ": ", .$`Level 2`)) %>% 
    arrange(Group, `Level 1`, -Freq) %>% 
    mutate(across(where(is.character), ~factor(., unique(.)))) %>% 
    mutate(`Level 2` = as.character(`Level 2`)) %>% 
    arrange(-Freq) %>% 
    mutate(`Level 2` = factor(`Level 2`, unique(.$`Level 2`))) %>% 
    arrange(Group, `Level 1`, -Freq) %>% 
    mutate(status = "All") -> tb.pb.alluv.all
  
  # DR
  
  my.dr.pb <- 
    daa.ko.glom.id %>% 
    filter(pgroup != "NS" & comp == i & Group == "Depleted") %>% 
    select(KO_ID = FunID, Group, comp, starts_with(i)) 
  
  tb.pb.alluv.all.tmp %>% 
    filter(KO_ID %in% my.dr.pb$KO_ID) %>% 
    select(-KO_ID) %>% 
    mutate(across(.cols = where(is.numeric), ~ ./sum(.x, na.rm = TRUE)))  %>% 
    mutate(across(.cols = c(`Level 1`:`Level 2`), ~ toupper(.)))  %>% 
    pivot_longer(cols = -c(`Level 1`, `Level 2`), names_to = "SAMPLE_ID", values_to = "Freq") %>%
    mutate(Group = str_remove_all(SAMPLE_ID, "\\.[:digit:]")) %>% 
    select(-SAMPLE_ID) %>% 
    group_by(`Level 1`, `Level 2`, Group) %>% 
    summarise(across(.cols = everything(), ~ sum(.x, na.rm = TRUE))) %>% 
    ungroup() %>% 
    mutate(`Level 1` = factor(`Level 1`, levels = c("METABOLISM", "ENVIRONMENTAL INFORMATION PROCESSING", "GENETIC INFORMATION PROCESSING", "CELLULAR PROCESSES"))) %>% 
    mutate(cat = paste0(.$`Level 1`, ": ", .$`Level 2`)) %>% 
    arrange(Group, `Level 1`, -Freq) %>% 
    mutate(across(where(is.character), ~factor(., unique(.)))) %>% 
    mutate(`Level 2` = as.character(`Level 2`)) %>% 
    arrange(-Freq) %>% 
    mutate(`Level 2` = factor(`Level 2`, unique(.$`Level 2`))) %>% 
    arrange(Group, `Level 1`, -Freq) %>% 
    mutate(status = "Depleted") -> tb.pb.alluv.dep
  
  # UP
  
  my.dr.ur <- 
    daa.ko.glom.id %>% 
    filter(pgroup != "NS" & comp == i & Group == "Enriched") %>% 
    select(KO_ID = FunID, Group, comp, starts_with(i)) 
  
  tb.pb.alluv.all.tmp %>% 
    filter(KO_ID %in% my.dr.ur$KO_ID) %>% 
    select(-KO_ID) %>% 
    mutate(across(.cols = where(is.numeric), ~ ./sum(.x, na.rm = TRUE)))  %>% 
    mutate(across(.cols = c(`Level 1`:`Level 2`), ~ toupper(.)))  %>% 
    pivot_longer(cols = -c(`Level 1`, `Level 2`), names_to = "SAMPLE_ID", values_to = "Freq") %>%
    mutate(Group = str_remove_all(SAMPLE_ID, "\\.[:digit:]")) %>% 
    select(-SAMPLE_ID) %>% 
    group_by(`Level 1`, `Level 2`, Group) %>% 
    summarise(across(.cols = everything(), ~ sum(.x, na.rm = TRUE))) %>% 
    ungroup() %>% 
    mutate(`Level 1` = factor(`Level 1`, levels = c("METABOLISM", "ENVIRONMENTAL INFORMATION PROCESSING", "GENETIC INFORMATION PROCESSING", "CELLULAR PROCESSES"))) %>% 
    mutate(cat = paste0(.$`Level 1`, ": ", .$`Level 2`)) %>% 
    arrange(Group, `Level 1`, -Freq) %>% 
    mutate(across(where(is.character), ~factor(., unique(.)))) %>% 
    mutate(`Level 2` = as.character(`Level 2`)) %>% 
    arrange(-Freq) %>% 
    mutate(`Level 2` = factor(`Level 2`, unique(.$`Level 2`))) %>% 
    arrange(Group, `Level 1`, -Freq) %>% 
    mutate(status = "Enriched") -> tb.pb.alluv.enr
  
  rbind(
    tb.pb.alluv.dep,
    tb.pb.alluv.all, 
    tb.pb.alluv.enr
  ) %>% 
    mutate(status = factor(status, levels = unique(.$status))) %>% 
    mutate(cat = factor(cat, levels = names(colors.ko))) -> tb.pb.alluv.daa
  
  
  ls.alluv.tb[[i]] <- tb.pb.alluv.daa
  
}

daa.ko.summary <- rbind(ls.alluv.tb[["SM05"]], ls.alluv.tb[["SM10"]], ls.alluv.tb[["SM20"]]) %>% mutate(Group = factor(Group, c("SM05", "SM10", "SM20")))

ggplot(daa.ko.summary,
       aes(x = status, stratum = `Level 2`, alluvium = `Level 2`,
           y = Freq/6, fill = cat)) +
  facet_wrap(.~Group, scales = "fixed", nrow = 1) + 
  geom_stratum(color = "#000000", alpha = 1, width = .333) +
  geom_flow() +
  scale_x_discrete(expand = c(0, 0), labels = c("<span style='color:#364B9AFF'>Depleted<br>Orthologs</span>", 
                                                "<span style='color:black'>All<br>Orthologs</span>", 
                                                "<span style='color:#A50026FF'>Enriched<br>Orthologs</span>")) +
  scale_y_continuous(expand = c(0, 0), position = "left", name = "Relative abundance (%)<br>",  
                     labels = scales::percent_format(accuracy = 1),
                     sec.axis = sec_axis(trans=~., 
                                         labels = scales::percent_format(accuracy = 1),
                                         name = "")) + 
  scale_fill_manual(values = colors.ko, name = NULL) +
  guides(fill = 
           guide_legend(order = 1, 
                        override.aes = list(shape = 16, 
                                            fill = colors.ko, 
                                            color = "black", 
                                            size = 1,
                                            linetype = 0),
                        ncol = 1)) +
  theme_bw() +
  my_theme + 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_markdown(size = 7, face = "bold", margin = unit(c(5, 0, 0, 0), "pt")), 
        panel.spacing.x = unit(25, "pt"),
        axis.ticks.x = element_blank()) -> gg.alluv.tmp.2

gg.alluv.tmp.2 <- 
  gg.alluv.tmp.2 +
  plot_grid(
    ggpubr::get_legend(addSmallLegend(gg.alluv.tmp.2 + theme(legend.position = "right")) + 
                         theme(legend.justification = c(0, 0.5), 
                               legend.box.margin = margin(0, 0, 0, 0),
                               legend.box.background = element_rect(fill = "white", colour = "black", linewidth = .75))),
    ncol = 2, align = "hv") + plot_layout(widths = c(1, 1))

gg.alluv.final <- gg.alluv.tmp.1 / gg.alluv.tmp.2

ggsave(filename = paste0(path, "figures/Functional/ProfileDAA.png"), plot = gg.alluv.final, width = 5000*3.5, height = 5000 * 2.5, units = 'px', dpi = 1200)
ggsave(filename = paste0(path, "figures/Functional/ProfileDAA.svg"), plot = gg.alluv.final, width = 5000*3.5, height = 5000 * 2.5, units = 'px', dpi = 1200)

# Tables ------------------------------------------------------------------

names(tb.ann.count.ko) <- names(tb.ann.count.ko) %>% str_replace_all(., "Level", "KO: Level")
names(tb.ann.count.pb) <- names(tb.ann.count.pb) %>% str_replace_all(., "Level", "PLaBAse: Level")

daa.ko.glom.id <- daa.ko.glom.id %>% dplyr::rename("KO_ID" = 1)
daa.ko.glom.group <- daa.ko.glom.group %>% dplyr::rename("KO_ID" = 1)

daa.pb.summary <- 
  daa.pb.summary %>% 
  data.frame(check.names = T) %>% 
  select(-cat) %>% 
  mutate(Freq = (Freq/6) * 100) %>% 
  arrange(status, Group) %>% 
  dplyr::rename("PLaBAse: Class" = 1, "PLaBAse: Subclass" = 2, "Status" = "status")

daa.ko.summary <- 
  daa.ko.summary %>% 
  data.frame(check.names = T) %>% 
  select(-cat) %>% 
  mutate(Freq = (Freq/6) * 100) %>% 
  arrange(status, Group) %>% 
  dplyr::rename("KO: Class" = 1, "KO: Subclass" = 2, "Status" = "status")

daa.tb.final <-
  
  list(
    
    "KO features" = tb.ann.count.ko,
    "PLaBAse features" = tb.ann.count.pb %>% as.data.frame(check.names = F),
    
    "DAA (Sample)" = daa.ko.glom.id,
    "DAA (Group)" = daa.ko.glom.group,
    
    "DAA KO (Summary)" = daa.ko.summary,
    "DAA PLaBAse (Summary)" = daa.pb.summary
    
  )

openxlsx::write.xlsx(x = daa.tb.final, file = paste0(path, "tables/FuncDiffAbundAnalysis.xlsx"), rowNames = F)

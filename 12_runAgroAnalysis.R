# Agro Analysis -----------------------------------------------------------

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

# Germination -------------------------------------------------------------

get_kruskal <- function(tb) {
  
  tb %>%
    group_by(Time) %>%
    nest() %>%
    mutate(kruskal = map(.x = data, ~kruskal.test(Value ~ Treatment, data = .x) %>% broom::tidy())) %>%
    unnest(kruskal) %>%
    mutate(across(where(is.numeric), ~round(., 4))) -> p
  
  return(p)
  
}

get_posthoc_lsd <- function(tb, tb.summary, tb.stat, ls.invert = NULL) {
  
  tb.compm <- data.frame(Time = "", Treatment = "", group = "")[-1,]
  
  for (i in unique(tb$Time)) {
    
    tmp.tb <- tb %>% filter(Time == i)
    
    if (stats::var(tmp.tb$Value) == 0) {
      
      next
      
    }
    
    if (i %in% ls.invert) {
      
      tmp.tb$Value <- -1 * tmp.tb$Value
      
    }
    
    tmp.lsd <- with(tmp.tb, agricolae::kruskal(tmp.tb$Value, tmp.tb$Treatment, alpha = 0.05))
    tmp.lsd <-
      cbind(Treatment = row.names(tmp.lsd[["groups"]]), tmp.lsd[["groups"]]) %>% 
      remove_rownames() %>% 
      mutate(Time = i, group = str_trim(groups)) %>% 
      select(Time, Treatment, group)
    
    if (tb.stat %>% filter(Time == i) %>% pull(p.value) > 0.05) {
      
      tmp.lsd <- mutate(tmp.lsd, group = NA_character_)
      
    }
    
    tb.compm <- bind_rows(tb.compm, tmp.lsd)
    
  }
  
  tb.compm <- merge(tb.summary, tb.compm, sort = F, all.x = T)
  
  return(tb.compm)
  
}

get_ggdensity <- function(tb) {
  
  tb %>% 
    ggplot(aes(Value)) +
    facet_wrap(~ Variable, scales = "free") + 
    geom_density(alpha = 0.8, fill = "gray30") + 
    xlab("Values") + ylab("Density") +
    theme_bw(base_family = "Arial") + 
    my_theme -> p
  
  return(p)
  
}

get_gghistogram <- function(tb) {
  
  tb %>% 
    ggplot(aes(Value)) +
    facet_wrap(~ Variable, scales = "free") + 
    geom_histogram(fill = "gray30") + 
    xlab("Values") + ylab("Density") +
    theme_bw(base_family = "Arial") + 
    my_theme -> p
  
  return(p)
  
}

get_qqplot <- function(tb) {
  
  tb %>% 
    nest(data = -Variable) %>% 
    mutate(test = map(.x=data, ~residuals(aov(Value ~ Treatment, data=.x)) %>% broom::tidy())) %>% 
    unnest(test) %>% 
    ggpubr::ggqqplot(., "x") +
    facet_wrap(~ Variable, scales = "free") + 
    theme_bw(base_family = "Arial") + 
    my_theme -> p
  
  return(p)
  
}

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

# Curves --------------------------------------------------------------------------------------

tb.curve.raw <-
  openxlsx::read.xlsx(xlsxFile = "misc/metadata2.xlsx", sheet = 1, check.names = F, sep.names = " ")

tb.curve.longer <- 
  tb.curve.raw %>% 
  mutate(across(Sample:Time, ~as.character(.))) %>% 
  mutate(across(Sample:Time, ~str_trim(.))) %>% 
  mutate(across(where(is.character), ~factor(., unique(.)))) %>% 
  pivot_longer(cols = where(is.numeric), names_to = "Variable", values_to = "Value") %>% 
  mutate(Value = (Value/7)) %>% 
  select(Variable, Treatment, Time, Value) %>% 
  arrange(Variable, Treatment, Time)

tb.curve.summary <-
  tb.curve.longer %>% 
  group_by(Variable, Treatment, Time) %>% 
  summarise(temp = list(c(psych::describe(Value), quantile(Value, na.rm = T), IQR = iqr(Value, na.rm = T)))) %>% 
  unnest_wider(temp) %>% 
  select(-c(vars, trimmed, mad)) %>% 
  mutate(across(where(is.numeric), ~case_when(is.nan(.) ~ NA_integer_, is.infinite(.) ~ NA_integer_,.default = .))) %>% 
  data.frame(check.names = F)

# ---------------------------------------------------------------------------------------------

tb.curve.assumptions <-
  tb.curve.longer %>% 
  na.omit(.) %>% 
  group_by(Variable) %>% 
  summarise(temp = list(c(shapiro.test(residuals(aov(Value ~ Treatment))), bartlett.test(Value ~ Treatment)))) %>% 
  unnest_wider(temp, names_repair = "universal") %>% 
  mutate(across(where(is.numeric), ~round(., 4))) %>% 
  select(Variable, 
         normality.stat = 2, normality.pvalue = 3, normality.method = 4, 
         homoscedasticity.stat = 6, homoscedasticity.pvalue = 8, homoscedasticity.method = 10)

get_ggdensity(tb.curve.longer)
get_gghistogram(tb.curve.longer)
get_qqplot(tb.curve.longer)

# Analysis of Variances -----------------------------------------------------------------------

tb.stat.curve <- 
  tb.curve.longer %>% 
  mutate(across(where(is.factor), ~droplevels(.))) %>% 
  filter(Variable == "Germination")

tb.stat.curve.summary <- 
  tb.curve.summary %>% 
  filter(Variable == "Germination") %>% 
  mutate(across(where(is.factor), ~droplevels(.)))

tb.kruskal.curve <- 
  get_kruskal(tb.stat.curve)

tb.stat.curve.compm.lsd <- 
  get_posthoc_lsd(tb.stat.curve, tb.stat.curve.summary, tb.kruskal.curve) %>% arrange(Time)

# Graphs --------------------------------------------------------------------------------------

tr.colors <- c("#525252FF", "#8cc084FF", "#38a3a5FF", "#0466c8FF")
names(tr.colors) <- unique(tb.curve.summary$Treatment)

tb.stat.curve.compm.lsd %>% 
  mutate(sig = ifelse(!is.na(group), "*", NA_character_)) %>% 
  group_by(Time) %>% 
  mutate(vpos = max(`100%`)) %>% 
  ungroup() %>% data.frame() -> tb.stat.curve.compm.lsd.2

# ---------------------------------------------------------------------------------------------

ggplot(tb.curve.summary, aes(x = Time, y = mean, group = Treatment)) +
  geom_line(aes(color = Treatment), size = 1, show.legend = F) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.05, position = position_dodge(.2)) + 
  geom_point(aes(fill = Treatment), size = 2, shape = 21, position = position_dodge(.2)) +
  geom_text(mapping = aes(label = sig, x = Time, y = vpos), 
            data = tb.stat.curve.compm.lsd.2,
            size = 4, vjust = -1.5, hjust = 0.5, color = "black", inherit.aes = T) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.02)), breaks = scales::extended_breaks(n = 6), 
                     labels = scales::percent_format(accuracy = 1)) +
  scale_x_discrete(expand = c(0.025, 0.025)) +
  scale_fill_manual(values = tr.colors) +
  scale_color_manual(values = tr.colors) +
  ylab('Seed germination (%)<br>') +
  xlab('<br>Time (hours)') +
  guides(fill = 
           guide_legend(order = 1, 
                        override.aes = list(shape = 22, 
                                            fill = tr.colors, 
                                            color = "black", 
                                            size = 3,
                                            linetype = 0),
                        ncol = 2)) +
  my_theme + 
  theme(legend.position = "top", 
        legend.box.background = element_rect(fill = "white", colour = "black", linewidth = .75)) -> g
g

ggsave(filename = paste0(path, "figures/Agro/Lineplot_Germ1.png"), plot = g, device = "png", width = 5000 * 1.5, height = 5000, units = "px", dpi = 1200, limitsize = F)
ggsave(filename = paste0(path, "figures/Agro/Lineplot_Germ1.svg"), plot = g, device = "svg", width = 5000 * 1.5, height = 5000, units = "px", dpi = 1200, limitsize = F)

ggplot(tb.curve.summary, aes(x = Time, y = mean, group = Treatment)) +
  geom_line(aes(color = Treatment), size = 1, show.legend = F) +
  #geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.05) + 
  geom_point(aes(fill = Treatment), size = 2, shape = 21) +
  geom_text(mapping = aes(label = sig, x = Time, y = vpos), 
            data = tb.stat.curve.compm.lsd.2,
            size = 4, vjust = -1.5, hjust = 0.5, color = "black", inherit.aes = T) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.02)), breaks = scales::extended_breaks(n = 6), 
                     labels = scales::percent_format(accuracy = 1)) +
  scale_x_discrete(expand = c(0.025, 0.025)) +
  scale_fill_manual(values = tr.colors) +
  scale_color_manual(values = tr.colors) +
  ylab('Seed germination (%)<br>') +
  xlab('<br>Time (hours)') +
  guides(fill = 
           guide_legend(order = 1, 
                        override.aes = list(shape = 22, 
                                            fill = tr.colors, 
                                            color = "black", 
                                            size = 3,
                                            linetype = 0),
                        ncol = 2)) +
  my_theme + 
  theme(legend.position = "top", 
        legend.box.background = element_rect(fill = "white", colour = "black", linewidth = .75)) -> g
g

ggsave(filename = paste0(path, "figures/Agro/Lineplot_Germ2.png"), plot = g, device = "png", width = 5000 * 1.5, height = 5000, units = "px", dpi = 1200, limitsize = F)
ggsave(filename = paste0(path, "figures/Agro/Lineplot_Germ2.svg"), plot = g, device = "svg", width = 5000 * 1.5, height = 5000, units = "px", dpi = 1200, limitsize = F)

# Agro --------------------------------------------------------------------

get_kruskal <- function(tb) {
  
  tb %>%
    group_by(Variable) %>%
    nest() %>%
    mutate(kruskal = map(.x = data, ~kruskal.test(Value ~ Treatment, data = .x) %>% broom::tidy())) %>%
    unnest(kruskal) %>%
    mutate(across(where(is.numeric), ~round(., 4))) -> p
  
  return(p)
  
}

get_posthoc_lsd <- function(tb, tb.summary, tb.stat, ls.invert = NULL) {
  
  tb.compm <- data.frame(Variable = "", Treatment = "", group = "")[-1,]
  
  for (i in unique(tb$Variable)) {
    
    tmp.tb <- tb %>% filter(Variable == i)
    
    if (stats::var(tmp.tb$Value) == 0) {
      
      next
      
    }
    
    if (i %in% ls.invert) {
      
      tmp.tb$Value <- -1 * tmp.tb$Value
      
    }
    
    tmp.lsd <- with(tmp.tb, agricolae::kruskal(tmp.tb$Value, tmp.tb$Treatment, alpha = 0.05))
    tmp.lsd <-
      cbind(Treatment = row.names(tmp.lsd[["groups"]]), tmp.lsd[["groups"]]) %>% 
      remove_rownames() %>% 
      mutate(Variable = i, group = str_trim(groups)) %>% 
      select(Variable, Treatment, group)
    
    if (tb.stat %>% filter(Variable == i) %>% pull(p.value) > 0.05) {
      
      tmp.lsd <- mutate(tmp.lsd, group = NA_character_)
      
    }
    
    tb.compm <- bind_rows(tb.compm, tmp.lsd)
    
  }
  
  tb.compm <- merge(tb.summary, tb.compm, sort = F, all.x = T)
  
  return(tb.compm)
  
}

# -------------------------------------------------------------------------

tb.agro.raw <-
  openxlsx::read.xlsx(xlsxFile = "misc/metadata2.xlsx", sheet = 2, check.names = F, sep.names = " ") %>% 
  select(-c(Treat1, Treat2, Pot)) %>% 
  group_by(Sample, Group, Treatment) %>% 
  summarise(across(where(is.numeric), mean)) %>% 
  data.frame(check.names = F)

tb.agro.longer <- 
  tb.agro.raw %>% 
  mutate(across(Sample:Treatment, ~as.character(.))) %>% 
  mutate(across(Sample:Treatment, ~str_trim(.))) %>% 
  mutate(across(where(is.character), ~factor(., unique(.)))) %>% 
  pivot_longer(cols = where(is.numeric), names_to = "Variable", values_to = "Value") %>% 
  mutate(Abbr = case_when(
    Variable == "Shoot fresh weight (g)" ~ "ShootFW",
    Variable == "Shoot dry weight (g)" ~ "ShootDW",
    Variable == "Root fresh weight (g)" ~ "RootFW",
    Variable == "Root dry weight (g)" ~ "RootDW",
    Variable == "Leaves (n)" ~ "Leaves",
    Variable == "Shoot length (cm)" ~ "ShootL",
    Variable == "Root length (cm)" ~ "RootL",
    Variable == "Chlorophyll A (mg g<sup>−1</sup>)" ~ "ChlA",
    Variable == "Chlorophyll B (mg g<sup>−1</sup>)" ~ "ChlB",
    Variable == "Carotenoids (mg g<sup>−1</sup>)" ~ "Car")) %>% 
  select(Sample:Variable, Abbr, Value) %>% 
  mutate(across(where(is.character), ~factor(., unique(.)))) %>% 
  arrange(Variable, Treatment)

tb.agro.summary <-
  tb.agro.longer %>% 
  group_by(Variable, Treatment) %>% 
  summarise(temp = list(c(psych::describe(Value), quantile(Value, na.rm = T), IQR = iqr(Value, na.rm = T)))) %>% 
  unnest_wider(temp) %>% 
  select(-c(vars, trimmed, mad)) %>% 
  mutate(across(where(is.numeric), ~case_when(is.nan(.) ~ NA_integer_, is.infinite(.) ~ NA_integer_,.default = .))) %>% 
  data.frame(check.names = F)

# ---------------------------------------------------------------------------------------------

tb.agro.assumptions <-
  tb.agro.longer %>% 
  na.omit(.) %>% 
  group_by(Variable) %>% 
  summarise(temp = list(c(shapiro.test(residuals(aov(Value ~ Treatment))), bartlett.test(Value ~ Treatment)))) %>% 
  unnest_wider(temp, names_repair = "universal") %>% 
  mutate(across(where(is.numeric), ~round(., 4))) %>% 
  select(Variable, 
         normality.stat = 2, normality.pvalue = 3, normality.method = 4, 
         homoscedasticity.stat = 6, homoscedasticity.pvalue = 8, homoscedasticity.method = 10)

get_ggdensity(tb.agro.longer)
get_gghistogram(tb.agro.longer)
get_qqplot(tb.agro.longer)

# -------------------------------------------------------------------------

tb.stat.agro <- 
  tb.agro.longer %>% 
  mutate(across(where(is.factor), ~droplevels(.))) 

tb.stat.agro.summary <- 
  tb.agro.summary %>% 
  mutate(Group = str_remove_all(Treatment, ":.*") %>% factor(., unique(.)))

tb.kruskal.agro <- 
  get_kruskal(tb.stat.agro) %>% 
  mutate(p.full = paste0("Kruskal-Wallis, *p*: ", ifelse(p.value < 0.001, "<0.001", sprintf("%1.3f", round(p.value, 3)))))

tb.stat.agro.compm.lsd <- 
  get_posthoc_lsd(tb.stat.agro, tb.stat.agro.summary, tb.kruskal.agro) %>% 
  mutate(Group = str_remove_all(Treatment, ":.*") %>% factor(., unique(.)))


ggplot() +
  facet_wrap( ~ Variable, scales = "free") +
  # geom_boxplot(data = tb.agro.longer, aes(x = Group, y = Value, fill = Treatment, color = Treatment), 
  #              width = .5, outliers = F, show.legend = F, alpha = .5) + 
  geom_point(data = tb.agro.longer, aes(x = Group, y = Value, color = Treatment),
             position = position_jitterdodge(dodge.width=0.9), size = -1, shape = 16, alpha = 1) +
  geom_bar(data = tb.stat.agro.summary, aes(x = Group, y = mean, fill = Treatment), 
           color = "black", width = .75, stat = 'identity', show.legend = F) + 
  geom_errorbar(data = tb.stat.agro.summary, aes(x = Group, ymin=mean-sd, ymax=mean+sd), width=.05) +
  geom_richtext(mapping = aes(label = p.full, x = 2.5, y = Inf, vjust = 2.5),
                data = tb.kruskal.agro, position = position_dodge(width = 0),
                size = 2.5, hjust = 0.5, color = "black", inherit.aes = T, 
                label.padding = unit(0, "pt"), fill = NA, label.color = NA) +
  geom_richtext(mapping = aes(label = group, x = as.numeric(Group), y = max), 
                data = tb.stat.agro.compm.lsd,
                size = 2.5, vjust = -1, hjust = 0.5, color = "black", inherit.aes = T,
                label.padding = unit(0, "pt"), fill = NA, label.color = NA) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.5)), breaks = scales::extended_breaks(n = 6)) +
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
  theme(legend.position = "inside", 
        legend.position.inside = c(.75, .15),
        legend.box.background = element_rect(fill = "white", colour = "black", linewidth = .75),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) -> g3

g3

ggsave(filename = paste0(path, "figures/Agro/Barplot_Agro.png"), plot = g3, device = "png", width = 5000 * 2.15, height = 5000 * 2, units = "px", dpi = 1200, limitsize = F)
ggsave(filename = paste0(path, "figures/Agro/Barplot_Agro.svg"), plot = g3, device = "svg", width = 5000 * 2.15, height = 5000 * 2, units = "px", dpi = 1200, limitsize = F)

# -------------------------------------------------------------------------

get_corr <- function(tb) {
  
  tb %>%
    group_by(Variable) %>%
    nest() %>%
    mutate(corr = map(.x = data, ~cor.test(.$`*Spirulina maxima* (g L<sup>−1</sup>)`, .$Value, data = .x, method = "spearman") %>% broom::tidy())) %>%
    unnest(corr) %>%
    select(Variable, `r` = estimate, statistic, p.value, method) %>% 
    data.frame(check.names = T) -> p
  
  return(p)
  
}

tb.agro.cor <- 
  tb.agro.longer %>% 
  mutate(`*Spirulina maxima* (g L<sup>−1</sup>)` = case_when(
    Group == "CTRL" ~ 0.0,
    Group == "SM05" ~ 0.5,
    Group == "SM10" ~ 1.0,
    Group == "SM20" ~ 2.0
  ))

tb.agro.cor.stat <- get_corr(tb.agro.cor)

ggplot(tb.agro.cor, aes(x = `*Spirulina maxima* (g L<sup>−1</sup>)`, y = Value, fill = Treatment)) +
  facet_wrap( ~ Variable, scales = "free") +
  geom_smooth(method = lm, color = "#343a40", fill = "#ced4da", se = TRUE, alpha = .75) +
  geom_point(size = 2, shape = 21) +
  geom_richtext(data = tb.agro.cor.stat,
                mapping = aes(label = paste0("*R* = ", round(r, 3), "<br>*p* = ", ifelse(p.value < 0.001, "<0.001", sprintf("%1.3f", round(p.value, 3)))), x = 0, y = Inf),
                fill = NA, label.color = NA, position = position_dodge(width = 0), 
                size = 2.5, vjust = 1.25, hjust = 0, color = "black", inherit.aes = T) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)), breaks = scales::extended_breaks(n = 6)) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05)), breaks = c(0, 0.5, 1.0, 2.0), labels = c("0.0", "0.5", "1.0", "2.0")) +
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
                        ncol = 1)) +
  my_theme + 
  theme(legend.position = "inside", 
        legend.position.inside = c(.75, .15),
        legend.box.background = element_rect(fill = "white", colour = "black", linewidth = .75),
        axis.ticks.x = element_blank())  -> g4

g4

ggsave(filename = paste0(path, "figures/Agro/CorrLinear_Agro.png"), plot = g4, device = "png", width = 5000 * 2.5, height = 5000 * 1.75, units = "px", dpi = 1200, limitsize = F)
ggsave(filename = paste0(path, "figures/Agro/CorrLinear_Agro.svg"), plot = g4, device = "svg", width = 5000 * 2.5, height = 5000 * 1.75, units = "px", dpi = 1200, limitsize = F)

# MV: PCA and Heatmap -----------------------------------------------------

library("FactoMineR")
library("factoextra")

tb.agro.raw
tb.agro.longer

mv.factor <-
  tb.agro.longer %>%
  select(Sample, Group, Treatment) %>% 
  distinct() %>% 
  data.frame(check.names = F)

mv.matrix <-
  tb.agro.longer %>% 
  select(Sample, Abbr, Value) %>% 
  pivot_wider(names_from = Abbr, values_from = Value) %>% 
  data.frame(check.names = F) %>% 
  column_to_rownames("Sample")

tr.colors.1 <- c("#525252FF", "#8cc084FF", "#38a3a5FF", "#0466c8FF")
names(tr.colors.1) <- unique(mv.factor$Group)

tr.colors.2 <- c("#525252FF", "#8cc084FF", "#38a3a5FF", "#0466c8FF")
names(tr.colors.2) <- unique(mv.factor$Treatment)

# PCA using FactoMineR
res.pca <- PCA(mv.matrix, graph = F, scale.unit = T)

# Screeplot of explained variances
fviz_screeplot(res.pca, addlabels = T, linecolor ="red", barfill="#495057", barcolor ="#000000") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), breaks = scales::extended_breaks(n = 6)) +
  scale_x_discrete(expand = expansion(mult = c(0.1, 0.1))) +
  my_theme -> pca.scree

pca.scree

fviz_contrib(res.pca, choice="var", axes = c(1,2), fill="#495057", color ="#000000") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), breaks = scales::extended_breaks(n = 6)) +
  scale_x_discrete(expand = expansion(mult = c(0.1, 0.1))) +
  my_theme -> pca.contr

pca.contr

# PCA - IND

fviz_pca_biplot(res.pca, 
                geom.ind = "point", fill.ind = mv.factor$Treatment, col.ind = "black",
                pointshape = 16, pointsize = 2, labelsize = 2, addEllipses = F,
                invisible = "quali", repel = T, col.var = "contrib",
                gradient.cols = rev(paletteer::paletteer_c("grDevices::SunsetDark", 30)),
                legend.title = list(color = "Contribution (%)")) +
  geom_point(size = 2.25, shape = 21, aes(fill = mv.factor$Treatment), color = "#00000000", inherit.aes = T) +
  xlab(paste0("<br>PC1 (",round(res.pca$eig[1,3], 2), "%)")) +
  ylab(paste0("PC2 (",round(res.pca$eig[2,3]-res.pca$eig[1,3], 2),"%)<br>")) +
  scale_fill_manual(values = tr.colors.2) +
  labs(title = "") +
  guides(
    fill = "none",
    color = 
      guide_colourbar(barheight = .5, barwidth = 10, 
                      label.position = "bottom", title.position = "top", 
                      direction = "horizontal", title.hjust = .5, 
                      draw.ulim = T, draw.llim = T)) +
  theme(legend.position = "top", 
        legend.title = element_text(color = "black", size = 10, face = "bold"),
        
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
        axis.text.y = element_markdown(size = 7, colour = 'black')) -> gg.pca

gg.pca +
  plot_grid(
    ggpubr::get_legend(ggplot(mv.factor, aes(x = 1, y = 1, fill = Treatment)) +
                         geom_point(shape = 22) +
                         scale_fill_manual(values = tr.colors.2) +
                         guides(fill = 
                                  guide_legend(order = 1, 
                                               override.aes = list(shape = 22, 
                                                                   fill = tr.colors, 
                                                                   color = "black", 
                                                                   size = 3,
                                                                   linetype = 0, 
                                                                   alpha = 1),
                                               ncol = 1)) + 
                         my_theme +
                         theme(legend.position = "right", 
                               legend.box.background = element_rect(fill = "white", colour = "black", linewidth = .75),
                               legend.justification = c(0.1, 0.5))),
    ncol = 1, align = 'hv') + plot_layout(widths = c(1, 1)) -> gg.pca

gg.pca

ggsave(filename = paste0(path, "figures/Agro/PCA_Agro.png"), plot = gg.pca, device = "png", width = 5000 * 2, height = 5000 * 1.25, units = "px", dpi = 1200, limitsize = F)
ggsave(filename = paste0(path, "figures/Agro/PCA_Agro.svg"), plot = gg.pca, device = "svg", width = 5000 * 2, height = 5000 * 1.25, units = "px", dpi = 1200, limitsize = F)

# HPs

# Anno. Col
hp.anno.col <- mv.factor %>% select("Treatment        " = Group)

hp.anno.col.col <- c("#525252FF", "#8cc084FF", "#38a3a5FF", "#0466c8FF")
names(hp.anno.col.col) <- unique(hp.anno.col$`Treatment        `)

hp.anno.col.col.ls <- 
  list(
    `Treatment        ` = hp.anno.col.col
  )

# Mat. prep.
hp.values <- mv.matrix %>% scale() %>% t()

hp.breaks <- max(abs(min(hp.values, na.rm = T)), max(hp.values, na.rm = T))

# Plot
ComplexHeatmap::pheatmap(hp.values, name = "Z-Score        ",
                         scale = "none",
                         breaks = c(seq(-ceiling(hp.breaks), ceiling(hp.breaks), .2)),
                         cluster_row = T, 
                         cluster_cols = T, 
                         fontsize_col = 10, 
                         fontsize_row = 10,
                         border_color = NA,
                         annotation_names_col = F,
                         annotation_names_row = F,
                         annotation_col = hp.anno.col, 
                         annotation_colors = hp.anno.col.col.ls,
                         color = colorRampPalette(c("#364B9AFF", "#4A7BB7FF", "#6EA6CDFF", "#98CAE1FF", "#FFFFFFFF", "#FDB366FF", "#F67E4BFF", "#DD3D2DFF", "#A50026FF"))(length(seq(-ceiling(hp.breaks), ceiling(hp.breaks), .2))-1),
                         clustering_method = "ward.D2") -> hp.plot

hp.plot

png(filename = paste0(path, "figures/Agro/Heatmap_Agro.png"), width = 5000 * 2.5, height = 5000 * 1, units = "px", res = 1200, type = "cairo")

hp.plot

dev.off()

# Correlations ------------------------------------------------------------

library("psych")
library("circlize")

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    from = rownames(cormat)[row(cormat)[ut]],
    to = rownames(cormat)[col(cormat)[ut]],
    r  = cormat[ut],
    p = pmat[ut]
  )
}

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

# Get top tables

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

hp.top.tb.p <- get.top.fun.table(ps.rarefied.tss, "Phylum", 10)
hp.top.tb.c <- get.top.fun.table(ps.rarefied.tss, "Class", 20)
hp.top.tb.o <- get.top.fun.table(ps.rarefied.tss, "Order", 30)
hp.top.tb.f <- get.top.fun.table(ps.rarefied.tss, "Family", 50)

hp.top.tb.c <- hp.top.tb.c %>% mutate(Phylum = ifelse(Phylum %in% hp.top.names.p, Phylum, "Others"))
hp.top.tb.o <- hp.top.tb.o %>% mutate(Phylum = ifelse(Phylum %in% hp.top.names.p, Phylum, "Others"))
hp.top.tb.f <- hp.top.tb.f %>% mutate(Phylum = ifelse(Phylum %in% hp.top.names.p, Phylum, "Others"))

hp.top.tb.p.corr <- hp.top.tb.p %>% select(Taxa, tb.metad$ID) %>% column_to_rownames("Taxa") %>% t()
hp.top.tb.c.corr <- hp.top.tb.c %>% select(Taxa, tb.metad$ID) %>% column_to_rownames("Taxa") %>% t()
hp.top.tb.o.corr <- hp.top.tb.o %>% select(Taxa, tb.metad$ID) %>% column_to_rownames("Taxa") %>% t()
hp.top.tb.f.corr <- hp.top.tb.f %>% select(Taxa, tb.metad$ID) %>% column_to_rownames("Taxa") %>% t()

tb.agro.p <- merge(hp.top.tb.p.corr, mv.matrix, by = 0)
tb.agro.c <- merge(hp.top.tb.c.corr, mv.matrix, by = 0)
tb.agro.o <- merge(hp.top.tb.o.corr, mv.matrix, by = 0)
tb.agro.f <- merge(hp.top.tb.f.corr, mv.matrix, by = 0)

tb.agro.p.corr <- corr.test(tb.agro.p[-1], method = "pearson", use = "pairwise", adjust = "fdr", alpha = 0.05)
tb.agro.c.corr <- corr.test(tb.agro.c[-1], method = "pearson", use = "pairwise", adjust = "fdr", alpha = 0.05)
tb.agro.o.corr <- corr.test(tb.agro.o[-1], method = "pearson", use = "pairwise", adjust = "fdr", alpha = 0.05)
tb.agro.f.corr <- corr.test(tb.agro.f[-1], method = "pearson", use = "pairwise", adjust = "fdr", alpha = 0.05)

tb.agro.p.corr.flat <- 
  flattenCorrMatrix(tb.agro.p.corr$r, tb.agro.p.corr$p) %>% 
  filter(from %in% colnames(hp.top.tb.p.corr)) %>% 
  filter(to %in% colnames(mv.matrix)) %>% 
  dplyr::rename(pearson = r, p.value = p) %>% 
  mutate(sig = if_else(abs(pearson) > 0.5 & p.value < 0.05, "*", ""))

tb.agro.c.corr.flat <- 
  flattenCorrMatrix(tb.agro.c.corr$r, tb.agro.c.corr$p) %>% 
  filter(from %in% colnames(hp.top.tb.c.corr)) %>% 
  filter(to %in% colnames(mv.matrix)) %>% 
  dplyr::rename(pearson = r, p.value = p) %>% 
  mutate(sig = if_else(abs(pearson) > 0.5 & p.value < 0.05, "*", ""))

tb.agro.o.corr.flat <- 
  flattenCorrMatrix(tb.agro.o.corr$r, tb.agro.o.corr$p) %>% 
  filter(from %in% colnames(hp.top.tb.o.corr)) %>% 
  filter(to %in% colnames(mv.matrix)) %>% 
  dplyr::rename(pearson = r, p.value = p) %>% 
  mutate(sig = if_else(abs(pearson) > 0.5 & p.value < 0.05, "*", ""))

tb.agro.f.corr.flat <- 
  flattenCorrMatrix(tb.agro.f.corr$r, tb.agro.f.corr$p) %>% 
  filter(from %in% colnames(hp.top.tb.f.corr)) %>% 
  filter(to %in% colnames(mv.matrix)) %>% 
  dplyr::rename(pearson = r, p.value = p) %>% 
  mutate(sig = if_else(abs(pearson) > 0.5 & p.value < 0.05, "*", ""))

# HPs

# Phylum

hp.anno.row <- 
  hp.top.tb.p %>% 
  select(Taxa, Phylum) %>% 
  mutate(Phylum = factor(Phylum, levels = c(hp.top.names.p, "Others"))) %>% 
  mutate(Taxa = paste0(Taxa, "        ")) %>% 
  column_to_rownames("Taxa")

hp.anno.row.col.p <- colorRampPalette(c(
  
  paletteer::paletteer_d("futurevisions::pegasi")
  
))(length(unique(hp.top.names.p)))

names(hp.anno.row.col.p) <- unique(hp.top.names.p)

hp.anno.row.col.p["Others"] <- "#E6E6E6FF"

hp.values <- 
  tb.agro.p.corr.flat %>% 
  select(from, to, pearson) %>% 
  pivot_wider(names_from = to, values_from = pearson) %>% 
  mutate(from = paste0(from, "        ")) %>% 
  column_to_rownames("from")

hp.labels <- 
  tb.agro.p.corr.flat %>% 
  select(from, to, p.value) %>% 
  pivot_wider(names_from = to, values_from = p.value) %>% 
  column_to_rownames("from")

col_fun = colorRamp2(c(seq(-1, 1, .25)), c("#364B9AFF", "#4A7BB7FF", "#6EA6CDFF", "#98CAE1FF", "#FFFFFFFF", "#FDB366FF", "#F67E4BFF", "#DD3D2DFF", "#A50026FF"))

ha = ComplexHeatmap::rowAnnotation(Phylum = hp.anno.row$Phylum, col = list(Phylum = hp.anno.row.col.p), show_annotation_name = FALSE)

# Plot
ComplexHeatmap::Heatmap(hp.values, 
                        name = "Pearson's correlation (r)        ", 
                        col = col_fun,
                        cluster_rows = F, 
                        cluster_columns = F,
                        row_names_gp = gpar(fontsize = 10),
                        column_names_gp = gpar(fontsize = 10), 
                        left_annotation = ha,
                        cell_fun = function(j, i, x, y, w, h, fill){
                          if(hp.labels[i, j] < 0.05) {
                            grid.text('*', x, y)
                          }
                        }) -> hp.plot

hp.plot

png(filename = paste0(path, "figures/Agro/CorrHeatmap_1-Phylum.png"), width = 5000 * 1.65, height = 5000 * .9, units = "px", res = 1200, type = "cairo")

hp.plot

dev.off()

# Class

hp.anno.row <- 
  hp.top.tb.c %>% 
  select(Taxa, Phylum) %>% 
  mutate(Phylum = factor(Phylum, levels = c(hp.top.names.p, "Others"))) %>% 
  mutate(Taxa = paste0(Taxa, "        ")) %>% 
  column_to_rownames("Taxa")

hp.anno.row.col.p <- colorRampPalette(c(
  
  paletteer::paletteer_d("futurevisions::pegasi")
  
))(length(unique(hp.top.names.p)))

names(hp.anno.row.col.p) <- unique(hp.top.names.p)

hp.anno.row.col.p["Others"] <- "#E6E6E6FF"

hp.values <- 
  tb.agro.c.corr.flat %>% 
  select(from, to, pearson) %>% 
  pivot_wider(names_from = to, values_from = pearson) %>% 
  mutate(from = paste0(from, "        ")) %>% 
  column_to_rownames("from")

hp.labels <- 
  tb.agro.c.corr.flat %>% 
  select(from, to, p.value) %>% 
  pivot_wider(names_from = to, values_from = p.value) %>% 
  column_to_rownames("from")

col_fun = colorRamp2(c(seq(-1, 1, .25)), c("#364B9AFF", "#4A7BB7FF", "#6EA6CDFF", "#98CAE1FF", "#FFFFFFFF", "#FDB366FF", "#F67E4BFF", "#DD3D2DFF", "#A50026FF"))

ha = ComplexHeatmap::rowAnnotation(Phylum = hp.anno.row$Phylum, col = list(Phylum = hp.anno.row.col.p), show_annotation_name = FALSE)

# Plot
ComplexHeatmap::Heatmap(hp.values, 
                        name = "Pearson's correlation (r)        ", 
                        col = col_fun,
                        cluster_rows = F, 
                        cluster_columns = F,
                        row_names_gp = gpar(fontsize = 10),
                        column_names_gp = gpar(fontsize = 10), 
                        left_annotation = ha,
                        cell_fun = function(j, i, x, y, w, h, fill){
                          if(hp.labels[i, j] < 0.05) {
                            grid.text('*', x, y)
                          }
                        }) -> hp.plot

hp.plot

png(filename = paste0(path, "figures/Agro/CorrHeatmap_2-Class.png"), width = 5000 * 1.65, height = 5000 * 1.35, units = "px", res = 1200, type = "cairo")

hp.plot

dev.off()

# Order

hp.anno.row <- 
  hp.top.tb.o %>% 
  select(Taxa, Phylum) %>% 
  mutate(Phylum = factor(Phylum, levels = c(hp.top.names.p, "Others"))) %>% 
  mutate(Taxa = paste0(Taxa, "        ")) %>% 
  column_to_rownames("Taxa")

hp.anno.row.col.p <- colorRampPalette(c(
  
  paletteer::paletteer_d("futurevisions::pegasi")
  
))(length(unique(hp.top.names.p)))

names(hp.anno.row.col.p) <- unique(hp.top.names.p)

hp.anno.row.col.p["Others"] <- "#E6E6E6FF"

hp.values <- 
  tb.agro.o.corr.flat %>% 
  select(from, to, pearson) %>% 
  pivot_wider(names_from = to, values_from = pearson) %>% 
  mutate(from = paste0(from, "        ")) %>% 
  column_to_rownames("from")

hp.labels <- 
  tb.agro.o.corr.flat %>% 
  select(from, to, p.value) %>% 
  pivot_wider(names_from = to, values_from = p.value) %>% 
  column_to_rownames("from")

col_fun = colorRamp2(c(seq(-1, 1, .25)), c("#364B9AFF", "#4A7BB7FF", "#6EA6CDFF", "#98CAE1FF", "#FFFFFFFF", "#FDB366FF", "#F67E4BFF", "#DD3D2DFF", "#A50026FF"))

ha = ComplexHeatmap::rowAnnotation(Phylum = hp.anno.row$Phylum, col = list(Phylum = hp.anno.row.col.p), show_annotation_name = FALSE)

# Plot
ComplexHeatmap::Heatmap(hp.values, 
                        name = "Pearson's correlation (r)        ", 
                        col = col_fun,
                        cluster_rows = F, 
                        cluster_columns = F,
                        row_names_gp = gpar(fontsize = 10),
                        column_names_gp = gpar(fontsize = 10), 
                        left_annotation = ha,
                        cell_fun = function(j, i, x, y, w, h, fill){
                          if(hp.labels[i, j] < 0.05) {
                            grid.text('*', x, y)
                          }
                        }) -> hp.plot

hp.plot

png(filename = paste0(path, "figures/Agro/CorrHeatmap_3-Order.png"), width = 5000 * 1.65, height = 5000 * 1.75, units = "px", res = 1200, type = "cairo")

hp.plot

dev.off()

# Family

hp.anno.row <- 
  hp.top.tb.f %>% 
  select(Taxa, Phylum) %>% 
  mutate(Phylum = factor(Phylum, levels = c(hp.top.names.p, "Others"))) %>% 
  mutate(Taxa = paste0(Taxa, "        ")) %>% 
  column_to_rownames("Taxa")

hp.anno.row.col.p <- colorRampPalette(c(
  
  paletteer::paletteer_d("futurevisions::pegasi")
  
))(length(unique(hp.top.names.p)))

names(hp.anno.row.col.p) <- unique(hp.top.names.p)

hp.anno.row.col.p["Others"] <- "#E6E6E6FF"

hp.values <- 
  tb.agro.f.corr.flat %>% 
  select(from, to, pearson) %>% 
  pivot_wider(names_from = to, values_from = pearson) %>% 
  mutate(from = paste0(from, "        ")) %>% 
  column_to_rownames("from")

hp.labels <- 
  tb.agro.f.corr.flat %>% 
  select(from, to, p.value) %>% 
  pivot_wider(names_from = to, values_from = p.value) %>% 
  column_to_rownames("from")

col_fun = colorRamp2(c(seq(-1, 1, .25)), c("#364B9AFF", "#4A7BB7FF", "#6EA6CDFF", "#98CAE1FF", "#FFFFFFFF", "#FDB366FF", "#F67E4BFF", "#DD3D2DFF", "#A50026FF"))

ha = ComplexHeatmap::rowAnnotation(Phylum = hp.anno.row$Phylum, col = list(Phylum = hp.anno.row.col.p), show_annotation_name = FALSE)

# Plot
ComplexHeatmap::Heatmap(hp.values, 
                        name = "Pearson's correlation (r)        ", 
                        col = col_fun,
                        cluster_rows = F, 
                        cluster_columns = F,
                        row_names_gp = gpar(fontsize = 10),
                        column_names_gp = gpar(fontsize = 10), 
                        left_annotation = ha,
                        cell_fun = function(j, i, x, y, w, h, fill){
                          if(hp.labels[i, j] < 0.05) {
                            grid.text('*', x, y)
                          }
                        }) -> hp.plot

hp.plot

png(filename = paste0(path, "figures/Agro/CorrHeatmap_4-Family.png"), width = 5000 * 1.65, height = 5000 * 2.5, units = "px", res = 1200, type = "cairo")

hp.plot

dev.off()

# Tables ------------------------------------------------------------------

ls.final.tables.curve <-
  
  list(
    
    "Raw" = 
      tb.curve.longer %>% select(Variable, Treatment, Time, Value) %>% 
      mutate(across(where(is.factor), ~as.character(.))) %>% 
      mutate(across(where(is.character), ~str_remove_all(., "</?sup>"))) %>% 
      mutate(across(where(is.character), ~str_remove_all(., "\\*"))),
    
    "Kruskal" = 
      tb.kruskal.curve %>% 
      select(-data) %>% 
      data.frame(check.names = F) %>% 
      mutate(across(where(is.factor), ~as.character(.))) %>% 
      mutate(across(where(is.character), ~str_remove_all(., "</?sup>"))) %>% 
      mutate(across(where(is.character), ~str_remove_all(., "\\*"))),
    
    "Summary" = 
      tb.stat.curve.compm.lsd %>% 
      select(-c(`0%`, `50%`, `100%`)) %>%
      mutate(across(where(is.factor), ~as.character(.))) %>% 
      mutate(across(where(is.character), ~str_remove_all(., "</?sup>"))) %>% 
      mutate(across(where(is.character), ~str_remove_all(., "\\*"))) %>% 
      dplyr::rename("Q1" = "25%", "Q3" = "75%") %>% 
      mutate(across(where(is.numeric), ~sprintf("%1.2f", round(., 2)))) %>% 
      mutate(full.mean = str_trim(str_remove(paste0(mean, " (± ", sd, ") ", group), "NA"))),
    
    "CompM" = 
      tb.stat.curve.compm.lsd %>% 
      mutate(across(where(is.factor), ~as.character(.))) %>% 
      mutate(across(where(is.character), ~str_remove_all(., "</?sup>"))) %>% 
      mutate(across(where(is.character), ~str_remove_all(., "\\*"))) %>% 
      mutate(across(where(is.numeric), ~sprintf("%1.2f", round(., 2)))) %>% 
      mutate(full.mean = str_trim(str_remove(paste0(mean, " (± ", sd, ") ", group), "NA"))) %>% 
      select(Time, Treatment, full.mean) %>% 
      pivot_wider(names_from = Time, values_from = full.mean) %>% 
      data.frame(check.names = F)
    
  )

openxlsx::write.xlsx(x = ls.final.tables.curve, file = paste0(path, "tables/AgroAnalysis_Germination.xlsx"), rowNames = F)

ls.final.tables.agro <-
  
  list(
    
    "Raw" = 
      tb.agro.longer %>% select(Variable, Abbr, Treatment, Value) %>% 
      mutate(across(where(is.factor), ~as.character(.))) %>% 
      mutate(across(where(is.character), ~str_remove_all(., "</?sup>"))) %>% 
      mutate(across(where(is.character), ~str_remove_all(., "\\*"))),
    
    "Kruskal" = 
      tb.kruskal.agro %>% 
      select(-data) %>% 
      data.frame(check.names = F) %>% 
      mutate(across(where(is.factor), ~as.character(.))) %>% 
      mutate(across(where(is.character), ~str_remove_all(., "</?sup>"))) %>% 
      mutate(across(where(is.character), ~str_remove_all(., "\\*"))),
    
    "Summary" = 
      tb.stat.agro.compm.lsd %>% 
      select(-c(`0%`, `50%`, `100%`)) %>%
      mutate(across(where(is.factor), ~as.character(.))) %>% 
      mutate(across(where(is.character), ~str_remove_all(., "</?sup>"))) %>% 
      mutate(across(where(is.character), ~str_remove_all(., "\\*"))) %>% 
      dplyr::rename("Q1" = "25%", "Q3" = "75%") %>% 
      mutate(across(where(is.numeric), ~sprintf("%1.2f", round(., 2)))) %>% 
      mutate(full.mean = str_trim(str_remove(paste0(mean, " (± ", sd, ") ", group), "NA"))),
    
    "CompM" = 
      tb.stat.agro.compm.lsd %>% 
      mutate(across(where(is.factor), ~as.character(.))) %>% 
      mutate(across(where(is.character), ~str_remove_all(., "</?sup>"))) %>% 
      mutate(across(where(is.character), ~str_remove_all(., "\\*"))) %>% 
      mutate(across(where(is.numeric), ~sprintf("%1.2f", round(., 2)))) %>% 
      mutate(full.mean = str_trim(str_remove(paste0(mean, " (± ", sd, ") ", group), "NA"))) %>% 
      select(Variable, Treatment, full.mean) %>% 
      pivot_wider(names_from = Variable, values_from = full.mean) %>% 
      data.frame(check.names = F)
    
  )

openxlsx::write.xlsx(x = ls.final.tables.agro, file = paste0(path, "tables/AgroAnalysis_PlantTest.xlsx"), rowNames = F)

ls.final.tables.corr  <-
  
  list(
    
    "Phylum" = tb.agro.p.corr.flat %>% 
      dplyr::rename("Taxa" = "from", "Variable" = "to"),
    
    "Class" = tb.agro.c.corr.flat %>% 
      dplyr::rename("Taxa" = "from", "Variable" = "to"),
    
    "Order" = tb.agro.o.corr.flat %>% 
      dplyr::rename("Taxa" = "from", "Variable" = "to"),
    
    "Family" = tb.agro.f.corr.flat %>% 
      dplyr::rename("Taxa" = "from", "Variable" = "to")
      
  )

openxlsx::write.xlsx(x = ls.final.tables.corr, file = paste0(path, "tables/AgroAnalysis_Correlations.xlsx"), rowNames = F)

# End ---------------------------------------------------------------------

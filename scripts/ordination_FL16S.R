### Setup env
library(tidyverse)
library(vegan)
library(ape)
library(patchwork)

source('scripts/MFD_colors.R')

mfdo2.palette <- readRDS("data/palette_mfd_hab2.rds") %>%
  unlist()

setwd("/mfd_diversity")

### Load data
data <- data.table::fread('data/2025-02-13_MFD_arcbac_genus_rarefaction_rel.csv', na.strings = "")

metadata <- readxl::read_excel("data/2025-04-14_mfd_db.xlsx")

### Beta diversity
## Calculate Hellinger-transformed Bray-Curtis dissimilarity in parallel
dist <- data %>%
  # mutate(across(where(is.numeric), ~./sum(.)*100)) %>%
  filter(!Genus == "Unclassified") %>%
  select(where(is.numeric), Genus) %>%
  column_to_rownames(var = "Genus") %>%
  t() %>% 
  decostand(., method = "hellinger") %>%
  parallelDist::parDist(., method = "bray", threads = 100) %>%
  as.matrix() %>%
  data.frame()

## Save Bray-Curtis matrix
data.table::fwrite(dist, sep = ",", row.names = TRUE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_BC_distance_full.csv"))


### Grid subset
## 10km grid subset
filter <- data.table::fread("output/2025-04-15_MFD_FL16S_combined_diversity_MFDO1.tsv") %>%
  select(complex) %>%
  distinct() %>%
  pull(complex)

## Grid metadata
metadata.grid <- data.table::fread('data/2025-02-19_MFD_samples_grid_10km.tsv', na.strings = "") %>%
  select(project_id:sampling_comment) %>%
  select(fieldsample_barcode, everything()) %>%
  mutate(complex = str_c(mfd_sampletype, mfd_areatype, mfd_hab1, sep = ", "))

metadata.grid.filt <- metadata.grid %>%
  mutate(across(mfd_hab1, ~str_replace(., "Sclerophyllous scrub", "Temperate heath and scrub")),
         across(mfd_hab2, ~str_replace(., "Scrub", "Sclerophyllous scrub"))) %>%
  mutate(complex = str_c(mfd_sampletype, mfd_areatype, mfd_hab1, sep = ", ")) %>%
  filter(complex %in% filter)

samples.grid <- metadata.grid %>%
  pull(fieldsample_barcode)

samples.grid.filt <- metadata.grid.filt %>%
  pull(fieldsample_barcode)

## Subset distance matrix
dist.grid <- dist %>%
  select(all_of(samples.grid)) %>%
  filter(rownames(.) %in% samples.grid)

dist.grid.filt <- dist %>%
  select(all_of(samples.grid.filt)) %>%
  filter(rownames(.) %in% samples.grid.filt)

## Write to output
data.table::fwrite(metadata.grid, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_metadata_BC_grid.tsv"))

data.table::fwrite(dist.grid, sep = ",", row.names = TRUE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_BC_distance_grid.csv"))

data.table::fwrite(metadata.grid.filt, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_metadata_BC_grid_filt.tsv"))

data.table::fwrite(dist.grid.filt, sep = ",", row.names = TRUE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_BC_distance_grid_filt.csv"))


### Create subset for PCOA
## Filter metadata and create "complex" corresponding to full MFDO1 string
metadata.select <- metadata %>%
  filter(fieldsample_barcode %in% colnames(dist),
    !is.na(mfd_hab1)) %>%
  select(fieldsample_barcode, mfd_sampletype:mfd_hab3) %>%
  mutate(complex = str_c(mfd_sampletype, mfd_areatype, mfd_hab1, sep = ", ")) %>%
  # filter(complex %in% filter) %>%
  mutate(across(complex, ~factor(., levels = names(mfdo1.palette)))) %>%
  group_by(complex) %>%
  mutate(complex_size = n()) %>%
  ungroup()

## Create group summary
groups.summary <- metadata.select %>%
  group_by(complex) %>%
  summarise(size = n()) %>%
  arrange(size)

## Subset metadata to remove groups with low sample counts
metadata.subset <- metadata.select %>%
  filter(!complex %in% c("Other, Urban, Landfill",
                         "Soil, Urban, Roadside",
                         "Water, Subterranean, Freshwater",
                         "Other, Urban, Drinking water",
                         "Water, Urban, Drinking water",
                         "Soil, Subterranean, Urban",
                         "Sediment, Subterranean, Saltwater")) %>%
  mutate(across(complex, ~droplevels(.)))

## Create new group summary of the subset
groups.summary.subset <- metadata.subset %>%
  group_by(complex) %>%
  summarise(size = n()) %>%
  arrange(size)

## Extract sample IDs from samples with MFDO1 information
samples.subset <- metadata.subset %>%
  filter(!is.na(mfd_hab1)) %>%
  pull(fieldsample_barcode)

## Extract MFDO1 levels
levels.subset <- metadata.subset %>%
  pull(complex) %>%
  droplevels() %>%
  levels()

## Subset color-palette based on MFDO1 levels
mfdo1.palette.subset <- mfdo1.palette[levels.subset]

## Subset distance matrix
dist.subset <- dist %>%
  select(all_of(samples.subset)) %>%
  filter(rownames(.) %in% samples.subset)

## Write to output
data.table::fwrite(metadata.subset, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_metadata_BC_subset.tsv"))

data.table::fwrite(dist.subset, sep = ",", row.names = TRUE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_BC_distance_subset.csv"))

## Subset distance matrix
dist.grid.subset <- dist.grid %>%
  select(any_of(samples.subset)) %>%
  filter(rownames(.) %in% samples.subset)

metadata.grid.subset <- metadata.grid %>%
  filter(fieldsample_barcode %in% samples.subset)

## Write to output
data.table::fwrite(metadata.grid.subset, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_metadata_BC_grid_subset.tsv"))

data.table::fwrite(dist.grid.subset, sep = ",", row.names = TRUE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_BC_distance_grid_subset.csv"))



### PCoA - Principal Coordinate Analysis
## Import Bray-Curtis matrix (if uncommented)
# dist <- data.table::fread("output/2024-03-19_BC-distance-subset.csv") %>% column_to_rownames(var = "V1")

## Perform the pcoa
PCOA <- pcoa(dist.subset)

## plot the eigenvalues and interpret 
barplot(PCOA$values$Relative_eig[1:10])

## Get percentage of variance explained by the first 3 principal coordinates
PCO1 <- round(sum(as.vector(PCOA$value$Relative_eig)[1])*100, 1)
PCO2 <- round(sum(as.vector(PCOA$value$Relative_eig)[2])*100, 1)
PCO3 <- round(sum(as.vector(PCOA$value$Relative_eig)[3])*100, 1)

## PCO1 and PCO2 explains 30% of the variation in the communities
sum(as.vector(PCOA$value$Relative_eig)[1:2])

## Extract the scores for plotting with ggplot
PCOA.scores <- PCOA$vectors %>%
  as.data.frame() %>%
  select(1:6) %>%
  rename_with(., ~str_replace(., "Axis.", "PCO")) %>%
  cbind(fieldsample_barcode = colnames(dist.subset)) %>%
  relocate(fieldsample_barcode, .before = "PCO1") %>%
  left_join(metadata.subset)

### Stats based on spatial subset
## Betadisper
betadisper <- betadisper(as.dist(dist.subset), metadata.subset$complex)

betadisper.aov <- anova(betadisper)

betadisper.aov

TukeyHSD(betadisper)

## Anosim
set.seed(123)
anosim <- anosim(as.dist(dist.subset), metadata.subset$complex, permutations = 999, parallel = 10)

anosim

## PERMANOVA - import
df.permanova.total <- data.table::fread("output/2025-03-19_MFD_PERMANOVA.tsv")

df.permanova.contrasts <- data.table::fread("output/2025-03-19_MFD_PERMANOVA_contrasts.tsv")


## Add labels to scores-object
PCOA.scores <- PCOA.scores %>%
  mutate(across(complex, ~factor(., levels = names(mfdo1.palette.subset)))) %>%
  filter(!is.na(complex)) %>%
  arrange(complex) %>%
  left_join(df.permanova.contrasts, by = c("complex" = "Model")) %>%
  mutate(label = str_c(mfd_sampletype, ", ", mfd_areatype, "<br>", mfd_hab1, ",<br>*n* = ", complex_size, sep = "",
                       "<br>", "__PERMANOVA:__ ", "<br>*R^2^* = ", 
                       round(R2, 2), "%, ", "*p* = ", pval)) %>%
  select(-R2, -pval) %>%
  cbind(df.permanova.total %>% filter(!is.na(pval))) %>%
  mutate(label_all = str_c("MFDO1 subset, ", "<br>*n* = ", nrow(.), "<br>", "__ANOSIM:__ ", "*R* = ", 
                           round(anosim$statistic, 2), ", ", "*p* = ", anosim$signif,
                           "<br>", "__PERMANOVA:__ ", "*R^2^* = ", 
                           round(R2, 2), "%, ", "*p* = ", pval)) %>%
  select(-Model, -R2, -pval) %>%
  mutate(across(label, ~str_replace(., "reclaimed lowland", "lowland")),
         across(label, ~factor(.)))

## Write object to output
data.table::fwrite(PCOA.scores, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_pcoa_scores.tsv"))


### Visualise results
## Plot of all samples in the subset - colored by MFDO1
p.pcoa.1v2.all <- ggplot() +
  geom_point(data = PCOA.scores, aes(x = PCO1, y = PCO2, fill = complex), 
             size = 4, alpha = 1, color = "black", pch = 21) +
  theme_minimal(base_size = 19) +
  scale_fill_manual(values = mfdo1.palette.subset) +
  guides(fill = guide_legend(override.aes = list(shape = 22, size = 5, alpha = 1))) +
  theme(legend.position = "none",
        aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(colour = "black"),
        strip.text.x = element_text(size = 16),
        plot.margin = margin(0,0,0,0),
        strip.text = ggtext::element_markdown()) +
  xlab(str_c("PCO1 - ", PCO1, "%", sep = "")) +
  ylab(str_c("PCO2 - ", PCO2, "%", sep = "")) +
  facet_wrap(vars(label_all), ncol = 1)

## Render plot
p.pcoa.1v2.all

## Plot of all samples in the subset - add contour to show density
p.pcoa.1v2.all.contour <- ggplot() +
  geom_point(data = PCOA.scores, aes(x = PCO1, y = PCO2, fill = complex),
             size = 5, alpha = 1, color = "black", pch = 21)  +
  scale_x_continuous(limits = c(-0.625, 0.4)) +
  scale_y_continuous(limits = c(-0.3, 0.525)) +
  scale_fill_manual(values = mfdo1.palette.subset) +
  ggnewscale::new_scale_fill() +
  stat_density_2d(data = PCOA.scores, aes(x = PCO1, y = PCO2, fill = after_stat(level)), 
                  geom = "polygon", colour = "black", adjust = 2.25, alpha = 0.2) +
  scale_fill_distiller(palette = "Greys", direction = 1) +
  theme_minimal(base_size = 19) +
  guides(fill = guide_legend(override.aes = list(shape = 22, size = 5, alpha = 1))) +
  theme(legend.position = "none",
        aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(colour = "black"),
        strip.text.x = element_text(size = 16),
        plot.margin = margin(0,0,0,0),
        strip.text = ggtext::element_markdown()) +
  xlab(str_c("PCO1 - ", PCO1, "%", sep = "")) +
  ylab(str_c("PCO2 - ", PCO2, "%", sep = "")) +
  facet_wrap(vars(label_all), ncol = 1)

## Render plot
p.pcoa.1v2.all.contour

## Plot of all samples in the subset - faceted panel
p.pcoa.1v2.facet <- PCOA.scores %>%
  ggplot() +
  geom_point(data = PCOA.scores[-15], aes(x = PCO1, y = PCO2), 
             size = 2, alpha = 0.5, color = "black", fill = "white", pch = 21) +
  scale_x_continuous(limits = c(-0.625, 0.4)) +
  scale_y_continuous(limits = c(-0.3, 0.525)) +
  geom_point(data = PCOA.scores, aes(x = PCO1, y = PCO2, fill = complex),
             size = 2, alpha = 1, color = "black", pch = 21) +
  theme_minimal(base_size = 19) +
  scale_fill_manual(values = mfdo1.palette.subset) +
  guides(fill = guide_legend(override.aes = list(shape = 22, size = 5, alpha = 1))) +
  theme(legend.position = "none",
        aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(colour = "black"),
        strip.text.x = element_text(size = 16),
        plot.margin = margin(0,0,0,0),
        strip.text = ggtext::element_markdown()) +
  xlab(str_c("PCO1 - ", PCO1, "%", sep = "")) +
  ylab(str_c("PCO2 - ", PCO2, "%", sep = "")) +
  facet_wrap(~factor(label, levels = levels(label)), ncol = 6)

## Render plot
p.pcoa.1v2.facet

## Create subset for Bogs, mires and fens
PCOA.scores.filt.bmf <- PCOA.scores %>%
  filter(complex == levels.subset[8])  %>%
  mutate(across(mfd_hab2, ~str_remove(., "[ \\s]*\\(non-habitat type\\)")),
         across(mfd_hab2, ~na_if(., "Mire")),
         across(mfd_hab2, ~replace_na(., "Other"))) %>%
  mutate(across(mfd_hab2, ~factor(., levels = c("Calcareous fens", "Sphagnum acid bogs", "Wet thicket", "Other")))) %>%
  mutate(label = str_c(mfd_sampletype, ", ", mfd_areatype, "<br>", mfd_hab1, ",<br>*n* = ", complex_size, sep = "")) %>%
  mutate(across(label, ~factor(.)))

#PCOA.scores.filt.bmf.tmp <- PCOA.scores %>%
#  filter(!complex == levels.subset[7]) %>%
#  select(-label) %>%
#  mutate(label = levels(PCOA.scores.filt.bmf$label))

## Extract MFDO2 levels of Bogs, mires and fens
levels.bmf <- PCOA.scores.filt.bmf %>%
  pull(mfd_hab2) %>%
  unique()

## Create new color palette for MFDO2 
mfdo2.palette.bmf <- c("#cdad00", "#32849f", "#7301A8FF", "grey60")
names(mfdo2.palette.bmf) <- c("Calcareous fens", "Sphagnum acid bogs", "Wet thicket", "Other")

## Plot of all Bogs, mires and fens samples - colored by MFDO2
p.pcoa.1v2.bmf <- ggplot() +
  geom_point(data = PCOA.scores[-15], aes(x = PCO1, y = PCO2),
             size = 5, alpha = 0.5, color = "black", fill = "white", pch = 21) +
  geom_point(data = PCOA.scores.filt.bmf, aes(x = PCO1, y = PCO2, fill = mfd_hab2),
             size = 5, alpha = 1, color = "black", pch = 21) +
  scale_fill_manual(values = mfdo2.palette.bmf, labels = names(mfdo2.palette.bmf)) +
  scale_x_continuous(limits = c(-0.625, 0.4)) +
  scale_y_continuous(limits = c(-0.3, 0.525)) +
  theme_minimal(base_size = 19) +
  guides(fill = guide_legend(title = "MFDO2", nrow = 2)) +
  theme(legend.position = "bottom",
        aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(colour = "black"),
        strip.text.x = element_text(size = 16),
        strip.text = ggtext::element_markdown(),
        plot.margin = margin(0,0,0,0, unit = "cm")) +
  xlab(str_c("PCO1 - ", PCO1, "%", sep = "")) +
  ylab(str_c("PCO2 - ", PCO2, "%", sep = "")) +
  facet_wrap(vars(label), ncol = 1)

## Render plot
p.pcoa.1v2.bmf


### Arrange plots 
ggarranged1 <- p.pcoa.1v2.all.contour / p.pcoa.1v2.bmf

## Render plot
ggarranged1

## Write arranged plots to output
png(file = 'output/ordination_left_panel.png',
    width = 800,
    height = 1200) 
ggarranged1
dev.off()

pdf(file = 'output/ordination_left_panel.pdf',
    width = 8,
    height = 12) 
ggarranged1
dev.off()

tiff(file = 'output/ordination_left_panel.tiff',
     width = 800,
     height = 1200) 
ggarranged1
dev.off()

ggsave("output/ordination_left_panel.svg", 
       plot = ggarranged1, width = 8, height = 12, 
       units = "in", dpi = "retina")

## Write arranged plots to output
png(file = 'output/ordination_righ_panel.png',
    width = 1900,
    height = 1200) 
p.pcoa.1v2.facet
dev.off()

pdf(file = 'output/ordination_right_panel.pdf',
    width = 19,
    height = 12) 
p.pcoa.1v2.facet
dev.off()

tiff(file = 'output/ordination_right_panel.tiff',
     width = 1900,
     height = 1200) 
p.pcoa.1v2.facet
dev.off()

ggsave("output/ordination_right_panel.svg", 
       plot = p.pcoa.1v2.facet, width = 19, height = 12, 
       units = "in", dpi = "retina")


### Create individual MFDO2 plots
## Extract labels
levels.label <- PCOA.scores %>%
  pull(label) %>%
  levels()

## Create empty data list
var.list = list()

## Change data to list format based on the MFDO1 categories
for (i in 1:length(levels.label)) {
  
  var.list[[i]] = PCOA.scores[PCOA.scores$label == levels.label[i],]
  
}

## Create empty plot list
plot.list = list()

## Create plot for each MFDO1 category
for (i in 1:length(var.list)) {
  data.filt <- PCOA.scores[PCOA.scores$label == levels.label[i],]
  
  p.1v2 <- ggplot() +
    geom_point(data = PCOA.scores[-15], aes(x = PCO1, y = PCO2), 
               size = 2, alpha = 0.5, color = "black", fill = "white", pch = 21) +
    geom_point(data = data.filt, aes(x = PCO1, y = PCO2, fill = mfd_hab2),
               size = 2, alpha = 1, color = "black", pch = 21) +
    scale_x_continuous(limits = c(-0.625, 0.4)) +
    scale_y_continuous(limits = c(-0.3, 0.525)) +
    theme_minimal(base_size = 12) +
    guides(fill = guide_legend(title = "MFDO2")) +
    theme(legend.position = "right",
          aspect.ratio = 1,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_line(colour = "black"),
          strip.text.x = element_text(size = 12),
          plot.margin = margin(0,0,0,0, unit = "cm")) +
    xlab(str_c("PCO1 - ", PCO1, "%", sep = "")) +
    ylab(str_c("PCO2 - ", PCO2, "%", sep = "")) +
    facet_wrap(vars(label), ncol = 1)
  
  p.1v3 <- ggplot() +
    geom_point(data = PCOA.scores[-15], aes(x = PCO1, y = PCO3), 
               size = 2, alpha = 0.5, color = "black", fill = "white", pch = 21) +
    geom_point(data = data.filt, aes(x = PCO1, y = PCO3, fill = mfd_hab2),
               size = 2, alpha = 1, color = "black", pch = 21) +
    scale_x_continuous(limits = c(-0.55, 0.3)) +
    scale_y_continuous(limits = c(-0.325, 0.375)) +
    theme_minimal(base_size = 12) +
    guides(fill = guide_legend(title = "MFDO2")) +
    theme(legend.position = "right",
          aspect.ratio = 1,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_line(colour = "black"),
          strip.text.x = element_text(size = 12),
          plot.margin = margin(0,0,0,0, unit = "cm")) +
    xlab(str_c("PCO1 - ", PCO1, "%", sep = "")) +
    ylab(str_c("PCO3 - ", PCO3, "%", sep = "")) +
    facet_wrap(vars(label), ncol = 1)
  
  p <- p.1v2 + p.1v3 + plot_layout(guides = "collect")
  
  plot.list[[i]] = p
}

### Save as multi-page pdf, with each page representing af MFDO1 category
## landscape format
pdf("output/ordinations_MFDO1_MFDO2_groups.pdf", paper = "a4r", width = 19, height = 12)
for (i in 1:length(var.list)) {
  print(plot.list[[i]])
}
dev.off()


### Save image
save.image(paste0('output/', format(Sys.time(), "%Y-%m-%d"), "_MFD_ordinations.RData"))


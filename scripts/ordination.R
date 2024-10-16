### Setup env
library(tidyverse)
library(vegan)
library(ape)
library(ampvis2)
library(patchwork)

source('MFD_colors.R')

mfdo2.palette <- readRDS("data/palette_mfd_hab2.rds") %>%
  unlist()

setwd("/mfd_diversity")

### Load data
load('data/2024-03-07_mfd-ampvis-arcbac-data.RData')

## Extract sample IDs from samples with MFDO1 information
samples.subset <- mfd.ampvis.arcbac.ra$metadata %>%
  filter(!is.na(mfd_hab1)) %>%
  pull(fieldsample_barcode)

## Filter metadata and create "complex" correponding to full MFDO1 string
metadata.filt <- mfd.ampvis.arcbac.ra$metadata %>%
  filter(fieldsample_barcode %in% samples.subset) %>%
  select(fieldsample_barcode, mfd_sampletype:mfd_hab3) %>%
  mutate(across(mfd_hab1, ~str_replace(., "Sclerophyllous scrub", "Temperate heath and scrub")),
         across(mfd_hab2, ~str_replace(., "Scrub", "Sclerophyllous scrub"))) %>%
  mutate(complex = str_c(mfd_sampletype, mfd_areatype, mfd_hab1, sep = ", ")) %>%
  mutate(across(complex, ~factor(., levels = names(mfdo1.palette)))) %>%
  group_by(complex) %>%
  mutate(complex_size = n()) %>%
  ungroup()

## Create group summary
groups.summary <- metadata.filt %>%
  group_by(complex) %>%
  summarise(size = n())

## Subset metadata to remove groups with low sample counts
metadata.subset <- metadata.filt %>%
  filter(!complex %in% c("Water, Subterranean, Freshwater",
                         "Water, Urban, Sandfilter",
                         "Soil, Urban, Other",
                         "Sediment, Urban, Other",
                         "Soil, Urban, Roadside",
                         "Soil, Subterranean, Urban"))

## Create new group summary of the subset
groups.summary.sub <- metadata.subset %>%
  group_by(complex) %>%
  summarise(size = n())

## Extract MFDO1 levels
levels.subset <- metadata.subset %>%
  pull(complex) %>%
  droplevels() %>%
  levels()

## Subset color-palette based on MFDO1 levels
mfdo1.palette.subset <- mfdo1.palette[levels.subset]


### Beta diversity
## Subset the observational table
data <- mfd.ampvis.arcbac.ra$abund %>%
  cbind(mfd.ampvis.arcbac.ra$tax) %>%
  select(all_of(samples.subset), Kingdom:Genus) %>%
  remove_rownames()

## Calculate Hellinger-transformed Bray-Curtis dissimilarity in parallel
dist <- data %>%
  mutate(across(where(is.numeric), ~./sum(.)*100)) %>%
  filter(!Genus == "Unclassified") %>%
  select(where(is.numeric), Genus) %>%
  column_to_rownames(var = "Genus") %>%
  t() %>% 
  decostand(., method = "hellinger") %>%
  parallelDist::parDist(., method = "bray", threads = 100) %>%
  as.matrix() %>%
  data.frame()

## Save Bray-Curtis matrix
data.table::fwrite(dist, "output/2024-03-19_BC-distance-subset.csv", row.names = TRUE)


### PCoA - Principal Coordinate Analysis
## Import Bray-Curtis matrix (if uncommented)
# dist <- data.table::fread("output/2024-03-19_BC-distance-subset.csv") %>% column_to_rownames(var = "V1")

## Perform the pcoa
PCOA <- pcoa(dist)

## plot the eigenvalues and interpret 
barplot(PCOA$values$Relative_eig[1:10])

## Get percentage of variance explained by the first 3 principal coordinates
sum(as.vector(PCOA$value$Relative_eig)[1])
sum(as.vector(PCOA$value$Relative_eig)[2])
sum(as.vector(PCOA$value$Relative_eig)[3])

## PCO1 and PCO2 explains 30% of the variation in the communities
sum(as.vector(PCOA$value$Relative_eig)[1:2])

## Extract the scores for plotting with ggplot
PCOA.scores <- PCOA$vectors %>%
  as.data.frame() %>%
  select(1:6) %>%
  rename_with(., ~str_replace(., "Axis.", "PCO")) %>%
  cbind(fieldsample_barcode = colnames(dist)) %>%
  relocate(fieldsample_barcode, .before = "PCO1") %>%
  left_join(metadata.filt)

## Add "complex" correponding to full MFDO1 string
PCOA.scores.subset <- PCOA.scores %>%
  mutate(across(complex, ~factor(., levels = names(mfdo1.palette.subset)))) %>%
  filter(!is.na(complex)) %>%
  arrange(complex) %>%
  mutate(label = str_c(mfd_sampletype, ", ", mfd_areatype, "\n", mfd_hab1, ",\nn = ", complex_size, sep = ""),
         label_all = str_c("MFDO1 subset, ", "\nn = ", nrow(.))) %>%
  mutate(across(label, ~factor(.)))

## Write object to output
data.table::fwrite(PCOA.scores.subset, "output/2024-03-19_pcoa-scores.csv")

### Visualise results
## Plot of all samples in the subset - colored by MFDO1
p.pcoa.1v2.all <- ggplot() +
  geom_point(data = PCOA.scores.subset, aes(x = PCO1, y = PCO2, fill = complex), 
             size = 4, alpha = 1, color = "black", pch = 21) +
  theme_minimal(base_size = 19) +
  scale_fill_manual(values = mfdo1.palette) +
  guides(fill = guide_legend(override.aes = list(shape = 22, size = 5, alpha = 1))) +
  theme(legend.position = "none",
        aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(colour = "black"),
        strip.text.x = element_text(size = 19),
        plot.margin = margin(0,0,0,0)) +
  xlab("PCO1 - 19%") +
  ylab("PCO2 - 11%") +
  facet_wrap(vars(label_all), ncol = 1)

## Render plot
p.pcoa.1v2.all

## Plot of all samples in the subset - add contour to show density
p.pcoa.1v2.all.contour <- ggplot() +
  geom_point(data = PCOA.scores.subset, aes(x = PCO1, y = PCO2),
             size = 5, alpha = 0.5, color = "black", pch = 21, fill = "white")  +
  scale_x_continuous(limits = c(-0.625, 0.4)) +
  scale_y_continuous(limits = c(-0.3, 0.525)) +
  ggnewscale::new_scale_fill() +
  stat_density_2d(data = PCOA.scores.subset, aes(x = PCO1, y = PCO2, fill = after_stat(level)), 
                  geom = "polygon", colour = "black", adjust = 2.25, alpha = 0.2) +
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
        strip.text.x = element_text(size = 19),
        plot.margin = margin(0,0,0,0)) +
  xlab("PCO1 - 19%") +
  ylab("PCO2 - 11%") +
  facet_wrap(vars(label_all), ncol = 1)

## Render plot
p.pcoa.1v2.all.contour

## Plot of all samples in the subset - faceted panel
p.pcoa.1v2.facet <- PCOA.scores.subset %>%
  ggplot() +
  geom_point(data = PCOA.scores.subset[-15], aes(x = PCO1, y = PCO2), 
             size = 2, alpha = 0.5, color = "black", fill = "white", pch = 21) +
  scale_x_continuous(limits = c(-0.625, 0.4)) +
  scale_y_continuous(limits = c(-0.3, 0.525)) +
  geom_point(data = PCOA.scores.subset, aes(x = PCO1, y = PCO2, fill = complex),
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
        strip.text.x = element_text(size = 19),
        plot.margin = margin(0,0,0,0)) +
  xlab("PCO1 - 19%") +
  ylab("PCO2 - 11%") +
  facet_wrap(~factor(label, levels = rev(levels(label))), ncol = 6)

## Render plot
p.pcoa.1v2.facet

## Create subset for Bogs, mires and fens
PCOA.scores.filt.bmf <- PCOA.scores.subset %>%
  filter(complex == levels.subset[7])  %>%
  mutate(across(mfd_hab2, ~str_remove(., "[ \\s]*\\(non-habitat type\\)")),
         across(mfd_hab2, ~na_if(., "Mire")),
         across(mfd_hab2, ~replace_na(., "Other"))) %>%
  mutate(across(mfd_hab2, ~factor(., levels = c("Calcareous fens", "Sphagnum acid bogs", "Wet thicket", "Other")))) %>%
  mutate(across(label, ~factor(.)))

#PCOA.scores.filt.bmf.tmp <- PCOA.scores.subset %>%
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
  geom_point(data = PCOA.scores.subset[-15], aes(x = PCO1, y = PCO2),
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
        strip.text.x = element_text(size = 19),
        plot.margin = margin(0,0,0,0, unit = "cm")) +
  xlab("PCO1 - 19%") +
  ylab("PCO3 - 11%") +
  facet_wrap(vars(label), ncol = 1)

## Render plot
p.pcoa.1v2.bmf


### Arrange plots 
ggarranged <- p.pcoa.1v2.facet + (p.pcoa.1v2.all.contour / p.pcoa.1v2.bmf) + plot_layout(widths = c(3.8,1))

## Render plot
ggarranged

## Write arranged plots to output
png(file = 'output/ordination-mfdo1-frag-combined.png',
    width = 1900,
    height = 1200) 
ggarranged
dev.off()

tiff(file = 'output/ordination-mfdo1-frag-combined.tiff',
    width = 1900,
    height = 1200) 
ggarranged
dev.off()

ggsave("output/ordination-mfdo1-frag-combined.svg", 
       plot = ggarranged, width = 14, height = 8, 
       units = "in", dpi = "retina")

## Write individual plots to output
# tiff(file = 'output/ordination-mfdo1-frag-all.tiff',
#      width = 1900,
#      height = 1200) 
# p.pcoa.1v2.all
# dev.off()
# 
# tiff(file = 'output/ordination-mfdo1-frag-MFDO1.tiff',
#     width = 1900,
#     height = 1200) 
# p.pcoa.1v2.facet
# dev.off()
# 
# tiff(file = 'output/ordination-mfdo1-frag-contour.tiff',
#      width = 600,
#      height = 600) 
# p.pcoa.1v2.all.contour
# dev.off()
# 
# tiff(file = 'output/ordination-mfdo1-frag-blowout.tiff',
#      width = 600,
#      height = 600) 
# p.pcoa.1v2.bmf
# dev.off()


### Create individual MFDO2 plots
## Extract labels
levels.label <- PCOA.scores.subset %>%
  pull(label) %>%
  levels()

## Create empty data list
var.list = list()

## Change data to list format based on the MFDO1 categories
for (i in 1:length(levels.label)) {
  
  var.list[[i]] = PCOA.scores.subset[PCOA.scores.subset$label == levels.label[i],]
  
}

## Create empty plot list
plot.list = list()

## Create plot for each MFDO1 category
for (i in 1:length(var_list)) {
  data.filt <- PCOA.scores.subset[PCOA.scores.subset$label == levels.label[i],]
  
  p.1v2 <- ggplot() +
    geom_point(data = PCOA.scores.subset[-15], aes(x = PCO1, y = PCO2), 
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
    xlab("PCO1 - 19%") +
    ylab("PCO2 - 11%") +
    facet_wrap(vars(label), ncol = 1)
  
  p.1v3 <- ggplot() +
    geom_point(data = PCOA.scores.subset[-15], aes(x = PCO1, y = PCO3), 
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
    xlab("PCO1 - 19%") +
    ylab("PCO3 - 4%") +
    facet_wrap(vars(label), ncol = 1)
  
  p <- p.1v2 + p.1v3 + plot_layout(guides = "collect")
  
  plot_list[[i]] = p
}

### Save as multi-page pdf, with each page representing af MFDO1 category
## landscape format
pdf("output/ordinations_MFDO1_MFDO2_groups.pdf", paper = "a4r", width = 19, height = 12)
for (i in 1:length(var_list)) {
  print(plot_list[[i]])
}
dev.off()


### Save image
save.image('output/2024-03-19_ordinations.RData')

### Setup env
library(tidyverse)
library(ampvis2)

source('scripts/MFD_colors.R')

setwd("/mfd_diversity")


### Import data
filter <- data.table::fread("output/2025-04-15_MFD_FL16S_combined_diversity_MFDO1.tsv") %>%
  select(complex) %>%
  distinct() %>%
  pull(complex)

## Filter and create "complex" correponding to full MFDO1 string for 10 km reference grid samples
metadata <- data.table::fread('data/2025-02-19_MFD_samples_grid_10km.tsv', na.strings = "") %>%
  mutate(across(mfd_hab1, ~str_replace(., "Sclerophyllous scrub", "Temperate heath and scrub")),
         across(mfd_hab2, ~str_replace(., "Scrub", "Sclerophyllous scrub"))) %>%
  filter(complex %in% filter) %>%
  mutate(across(complex, ~factor(., levels = names(mfdo1.palette)))) %>%
  select(fieldsample_barcode, everything())

## Create group summary
groups.summary <- metadata %>%
  group_by(complex) %>%
  summarise(size = n())

## Pull sample IDs
samples.subset <- metadata %>%
  pull(fieldsample_barcode)

## Import count data of the genus-aggregated observational table of 16S dervied from metagenomes
data <- data.table::fread('data/2025-02-13_MFD_arcbac_genus_rarefaction_rel.csv', na.strings = "")

## Subset data to 10 km reference samples
data.subset <- data %>%
  select(all_of(samples.subset), Kingdom:Genus) %>%
  filter(rowSums(across(where(is.numeric)))!=0)

rm(data)

## Create ampvis2 object with data subset
ampvis.grid <- amp_load(otutable = data.subset,
                        metadata = metadata)

## Normalise and remove unclassified fraction
ampvis.grid.subset <- ampvis.grid %>%
  amp_subset_taxa(normalise = FALSE, tax_vector = c("Unclassified"), remove = TRUE)

## Create phylum level aggregated wide data frame using ampvis2
heatmap.phylum <- ampvis.grid.subset %>%
  amp_heatmap(group_by = "complex",
              tax_aggregate = "Phylum",
              tax_show = 20,
              normalise = FALSE,
              showRemainingTaxa = FALSE, 
              plot_values_size = 3,
              plot_na = TRUE,
              plot_values = FALSE,
              textmap = TRUE) %>%
  rownames_to_column(var = "Phylum")

## Pivot longer
# data.heatmap.phylum <- heatmap.phylum %>%
#   rownames_to_column(var = "Phylum") %>%
#   pivot_longer(!Phylum, names_to = "complex", values_to = "abund")

## Write table to output
data.table::fwrite(heatmap.phylum, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_phylum_heatmap_grid_MFDO1.tsv"))


### Save iamge
save.image(paste0('output/', format(Sys.time(), "%Y-%m-%d"), "_heatmap.RData"))


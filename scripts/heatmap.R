### Setup env
library(tidyverse)
library(ampvis2)

source('MFD_colors.R')

setwd("/mfd_diversity")


### Import data
## Filter and create "complex" correponding to full MFDO1 string for 10 km reference grid samples
metadata <- data.table::fread('data/2024-05-10_samples-grid-10km.csv', na.strings = "") %>%
  mutate(complex = str_c(mfd_sampletype, mfd_areatype, mfd_hab1, sep = ", ")) %>%
  group_by(complex) %>%
  mutate(complex_size = n()) %>%
  mutate(label = str_c(complex, ", n = ", complex_size, sep = "")) %>%
  ungroup() %>%
  filter(!complex %in% c("Water, Subterranean, Freshwater",
                         "Water, Urban, Sandfilter",
                         "Soil, Urban, Other",
                         "Sediment, Urban, Other",
                         "Soil, Urban, Roadside",
                         "Soil, Subterranean, Urban")) %>%
  mutate(across(complex, ~factor(., levels = names(mfdo1.palette)))) %>%
  select(-project_id)

## Create group summary
groups.summary <- metadata %>%
  group_by(complex) %>%
  summarise(size = n())

## Pull sample IDs
samples.subset <- metadata %>%
  pull(fieldsample_barcode)

## Import count data of the genus-aggregated observational table of 16S dervied from metagenomes
data <- data.table::fread('data/2024-03-07_arcbac-rarefaction-count.csv', na.strings = "")

## Subset data to 10 km reference samples
data.subset <- data %>%
  select(all_of(samples.subset), Kingdom:Genus) %>%
  filter(rowSums(across(where(is.numeric)))!=0)

## Create ampvis2 object with data subset
ampvis.grid <- amp_load(otutable = data.subset,
                        metadata = metadata)

## Normalise and remove unclassified fraction
ampvis.grid.subset <- ampvis.grid %>%
  amp_subset_taxa(normalise = TRUE, tax_vector = c("Unclassified"), remove = TRUE)

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
              textmap = TRUE)

## Pivot longer
data.heatmap.phylum <- heatmap.phylum %>%
  rownames_to_column(var = "Phylum") %>%
  pivot_longer(!Phylum, names_to = "complex", values_to = "abund")

## Write table to output
data.table::fwrite(data.heatmap.phylum, "output/2024-03-19_phylum-heatmap-grid.csv")


### Save iamge
save.image(file = 'output/2024-03-19_heatmap.RData')

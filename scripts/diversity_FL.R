### Setup env
library(tidyverse)
library(vegan)
library(ampvis2)
library(iNEXT)
library(patchwork)
library(rstatix)

### Import color palettes
source('MFD_colors.R')

setwd("/mfd_diversity")


### Import data
## Load metadata
metadata <- data.table::fread("data/2024-05-10_MFD-FL-metadata-subset.csv", sep = ",", header = TRUE) %>%

## Load observational table submitted to random subsampling without replacement
fl.asv.subset.ra <- data.table::fread("data/2024-05-10_MFD_FL16S_OTU_subset-ra.csv", sep = ",", header = TRUE) %>%

### Alpha diversity
fl.alpha.diversity <- fl.asv.subset.ra %>%
  mutate(across(everything(), ~+as.logical(.x))) %>%
  filter(rowSums(across(where(is.numeric)))!=0) %>%
  colSums() %>%
  data.frame('observed_species' = .) %>%
  rownames_to_column(var = "fieldsample_barcode") %>%
  left_join(metadata %>% select(fieldsample_barcode, complex)) %>%
  select(-fieldsample_barcode) %>%
  filter(!complex %in% setdiff(metadata$complex, levels.subset)) %>%
  mutate(across(complex, ~factor(., levels = levels)))

## Visualise disribution of alpha diversity across MFDO1
plot.alpha <- fl.alpha.diversity %>%
  ggplot(aes(x = complex, y = observed_species, fill = complex)) +
  geom_boxplot() +
  scale_fill_manual(values = mfdo1.palette.subset,
                    breaks = levels,
                    name = "MFDO1") +
  xlab('MFDO1') +
  ylab('Observed species') +
  coord_flip() +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "none")

## Render plot
plot.alpha

## Write alpha table to output
data.table::fwrite(fl.alpha.diversity, "output/2024-03-19_alpha-diversity-fl.csv")


### Statistical analysis of differences in alpha diversity
## Get summary
fl.alpha.diversity %>%
  filter(complex %in% c("Soil, Natural, Forests", "Soil, Urban, Greenspaces")) %>%
  group_by(complex) %>%
  summarise(median = median(observed_species), 
            mean = mean(observed_species),
            sd = sd(observed_species))

## Identify outliers
fl.alpha.diversity %>%
  filter(complex %in% c("Soil, Natural, Forests", "Soil, Urban, Greenspaces")) %>%
  group_by(complex) %>%
  identify_outliers(observed_species)

## qqplot
fl.alpha.diversity %>%
  filter(complex %in% c("Soil, Natural, Forests", "Soil, Urban, Greenspaces")) %>%
  ggqqplot(., x = "observed_species", facet.by = "complex")

## Test for normality (Shapiro-Wilk)
fl.alpha.diversity %>%
  filter(complex %in% c("Soil, Natural, Forests", "Soil, Urban, Greenspaces")) %>%
  group_by(complex) %>%
  shapiro_test(observed_species)

## Test for homoscedasticity
fl.alpha.diversity %>%
  filter(complex %in% c("Soil, Natural, Forests", "Soil, Urban, Greenspaces")) %>%
  levene_test(observed_species ~ complex)

## Welch t-test 
fl.alpha.diversity %>%
  filter(complex %in% c("Soil, Urban, Greenspaces", "Soil, Natural, Forests")) %>%
  mutate(across(complex, ~factor(., levels = c("Soil, Urban, Greenspaces", "Soil, Natural, Forests")))) %>%
  t_test(observed_species ~ complex, alternative = "greater", conf.level = 0.95, detailed = TRUE)

## Cohen’s d for Welch t-test
fl.alpha.diversity %>%
  filter(complex %in% c("Soil, Urban, Greenspaces", "Soil, Natural, Forests")) %>%
  mutate(across(complex, ~factor(., levels = c("Soil, Urban, Greenspaces", "Soil, Natural, Forests")))) %>%
  cohens_d(weight ~ group, var.equal = FALSE)

## Student t-test
fl.alpha.diversity %>%
  filter(complex %in% c("Soil, Urban, Greenspaces", "Soil, Natural, Forests")) %>%
  mutate(across(complex, ~factor(., levels = c("Soil, Urban, Greenspaces", "Soil, Natural, Forests")))) %>%
  t_test(observed_species ~ complex, alternative = "greater", conf.level = 0.95, detailed = TRUE, var.equal = TRUE)

## Cohen’s d for Welch t-test
fl.alpha.diversity %>%
  filter(complex %in% c("Soil, Urban, Greenspaces", "Soil, Natural, Forests")) %>%
  mutate(across(complex, ~factor(., levels = c("Soil, Urban, Greenspaces", "Soil, Natural, Forests")))) %>%
  cohens_d(weight ~ group, var.equal = TRUE)

## Wilcoxon signed-rank test
fl.alpha.diversity %>%
  filter(complex %in% c("Soil, Urban, Greenspaces", "Soil, Natural, Forests")) %>%
  mutate(across(complex, ~factor(., levels = c("Soil, Urban, Greenspaces", "Soil, Natural, Forests")))) %>%
  wilcox_test(observed_species ~ complex, alternative = "greater", conf.level = 0.95, detailed = TRUE)



### Gamma
## Load gamma table
data.gamma <- data.table::fread("data/2024-05-10_MFD-FL-iNEXT-richness-MFDO1.csv", sep = ",", header = TRUE)

## Create summary on MFDO1 after filtering
groups.summary.subset <- metadata %>%
  filter(reads >= min.reads) %>%
  group_by(complex) %>%
  summarise(size = n())

## Combine into data frame
fl.gamma.diversity <- data.gamma %>%
  mutate(across(Diversity, ~factor(., levels = c("Species richness", "Shannon diversity", "Simpson diversity")))) %>%
  rename(complex = Assemblage) %>% 
  left_join(groups.summary.subset %>% select(-size)) %>%
  mutate(across(Diversity, ~replace_na(., "Species richness"))) %>%
  complete(complex, Diversity) %>%
  mutate(across(complex, ~factor(., levels = levels.subset))) %>%
  filter(!is.na(complex)) %>%
  left_join(groups.summary.subset)

## Plot of total richness (q=0)
plot.gamma.total <- fl.gamma.diversity %>%
  ggplot(aes(x = complex, y = Estimator, group = Diversity, fill = complex)) +
  geom_bar(position=position_dodge(), aes(y=Estimator), stat="identity") +
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0.5, position=position_dodge(width=0.9)) +
  scale_fill_manual(values = mfdo1.palette.subset,
                    breaks = levels,
                    name = "MFDO1") +
  facet_grid(rows = vars(Diversity), scales = "free_y") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

## Plot of Shannon diversity (q=1)
plot.gamma.shannon <- fl.gamma.diversity %>%
  filter(Diversity == "Shannon diversity") %>%
  ggplot(aes(x = complex, y = Estimator, group = Diversity, fill = complex)) +
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0.5, position=position_dodge(width=0.9)) +
  geom_point(shape = 21) +
  geom_text(aes(y = -2000, label = str_c("N = ", size)), size = 4,
            fontface = "bold") +
  scale_fill_manual(values = mfdo1.palette.subset,
                    breaks = levels,
                    name = "MFDO1") +
  scale_y_continuous(limits = c(-3000,50000)) +
  xlab('') +
  ylab('Shannon diversity') +
  coord_flip() +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_blank())

## Render plots
plot.gamma.total + plot.gamma.shannon

## Write gamma table to output
data.table::fwrite(fl.gamma.diversity, "output/2024-03-19_gamma-diversity-fl.csv")

### Combine alpha and gamma table
table.diversity <- fl.alpha.diversity %>%
  group_by(complex) %>%
  summarise(median_diversity = round(median(observed_species), 0),
            mean_diversity = round(mean(observed_species), 0),
            sd_diversity = round(sd(observed_species), 0)) %>%
  left_join(fl.gamma.diversity)

## Write gamma table to output
data.table::fwrite(table.diversity, "output/2024-03-19_combined-diversity-fl.csv")


### Save image
save.image(file = 'output/2024-03-19_FL-diversity.RData')
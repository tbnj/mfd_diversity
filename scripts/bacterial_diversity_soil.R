### Setup env
library(tidyverse)
library(patchwork)
library(rstatix)
library(multcompView)

source('MFD_colors.R')

setwd("/mfd_diversity")

### Import data
## Alpha diversity from FL16S amplicons
fl.alpha.diversity <- data.table::fread("output/2024-03-19_alpha-diversity-fl.csv")

## Gamma diversity from FL16S amplicons
fl.gamma.diversity <- data.table::fread("output/2024-03-19_gamma-diversity-fl.csv")

## Extract soil MFDO1 categories
levels.subset <- fl.alpha.diversity %>%
  filter(str_detect(complex, "Soil"),
         !str_detect(complex, "Subterranean|Coastal")) %>%
  group_by(complex) %>%
  reframe(mean = mean(observed_species),
          median = median(observed_species)) %>%
  arrange(desc(median)) %>%
  pull(complex) %>%
  unique()

## Subset the color palette
mfdo1.palette.subset <- mfdo1.palette[levels.subset]

## Subset alpha diversity
fl.alpha.diversity.subset <- fl.alpha.diversity %>%
  filter(complex %in% levels.subset) %>%
  mutate(across(complex, ~factor(., levels = levels.subset)))

## Subset gamma diversity
fl.gamma.diversity.subset <- fl.gamma.diversity %>%
  filter(complex %in% levels.subset) %>%
  mutate(across(complex, ~factor(., levels = levels.subset)))

## Summary of soil MFDO1 categories
fl.alpha.diversity.subset %>%
  group_by(complex) %>%
  summarise(group_size = n(),
            median = median(observed_species),
            mean = mean(observed_species),
            sd = sd(observed_species))

### Statitical test
## Test for normality (Shapiro-Wilk)
fl.alpha.diversity.subset %>%
  group_by(complex) %>%
  shapiro_test(observed_species)

## Test for homoscedasticity
fl.alpha.diversity.subset %>%
  levene_test(observed_species ~ complex)

## analysis of variance
anova <- aov(observed_species ~ complex, data = fl.alpha.diversity.subset)
summary(anova)

## Tukey's tukey (Tukey HSD)
tukey_hsd(anova)

## Create data frame of Tukey HSD results
tukey <- TukeyHSD(anova)

## Compact letter display
cld <- multcompLetters4(anova, tukey)

## Attanged table of alpha diversity
dt <- fl.alpha.diversity.subset %>%
  group_by(complex) %>%
  summarise(mean=mean(observed_species), sd = sd(observed_species)) %>%
  arrange(desc(mean))

## Extract the compact letter display and adding to the table
cld <- as.data.frame.list(cld$complex)
dt$cld <- cld$Letters

## Investigate
print(dt)

## Write to output
data.table::fwrite(dt, "output/2024-03-19_soil_diveristy_analysis.csv")


### Visualisation
## Boxplot of alpha diversity
plot.alpha <- fl.alpha.diversity.subset %>%
  ggplot(aes(x = complex, y = observed_species, fill = complex)) +
  geom_boxplot() +
  geom_text(data = dt, aes(x = complex, y = -100, label = cld),
            size = 4, fontface = "bold") +
  scale_fill_manual(values = mfdo1.palette.subset,
                    breaks = levels,
                    name = "MFDO1") +
  xlab('MFDO1') +
  ylab('Observed species') +
  # coord_flip() +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "none")

## Render plot
plot.alpha

### Barplot of gamma diversity (Shannon)
plot.gamma <- fl.gamma.diversity.subset %>%
  filter(Diversity == "Shannon diversity") %>%
  mutate(across(complex, ~factor(., levels = levels.subset))) %>%
  ggplot(aes(x = complex, y = Estimator, group = Diversity, fill = complex)) +
  geom_bar(position=position_dodge(), aes(y=Estimator), stat="identity", color = "black") +
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0, position=position_dodge(width=0.9), size = 1) +
  # geom_point(shape = 21) +
  geom_text(aes(y = -2000, label = str_c("N = ", size)), size = 4,
            fontface = "bold") +
  scale_fill_manual(values = mfdo1.palette.subset,
                    breaks = levels,
                    name = "MFDO1") +
  scale_y_continuous(limits = c(-3000,50000)) +
  xlab('') +
  ylab('Shannon diversity') +
  # coord_flip() +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_blank())

## Combined plot
plot.gamma / plot.alpha

### Setup env
library(tidyverse)
library(knitr)
library(openxlsx)
library(ggraph)
library(igraph)
library(gtools)
library(treedataverse)
library(dendroextras)
library(ggvenn)
library(ggtreeExtra)
library(ggnewscale)
library(wesanderson)
library(ape)
library(gridExtra)
library(pvclust)
library(vegan)
library(parallelDist)

set.seed(123)

options(stringsAsFactors = F, gsubfn.engine = "R")

Sys.setenv("LANGUAGE"="En")

source("MFD_colors.R")

setwd("/mfd_diversity")

wd <- getwd()

data.path <- paste0(wd, '/data')
output.path <- paste0(wd, '/output')

### Load data
## MFD Ontology
mfd_ontology <- read.xlsx(paste0(data.path, "/2023-12-21_mfd-habitat-ontology.xlsx"), sheet = 1)

## Bray-Curtis matix
BC.df <- data.table::fread(paste0(data.path, "/2024-03-19_BC-distance-subset.csv"))

## Set rownames for BC matrix
BC.df <- BC.df %>%
  `rownames<-`(colnames(BC.df))

## MFD db
mfd_db <- read.xlsx(paste0(data.path, "/2024-02-13_mfd_db.xlsx"), sheet = 1) %>%
  filter(fieldsample_barcode %in% colnames(BC.df))

## 10 km grid representatives
rep.10km.df <- read_csv(paste0(data.path, "/2024-05-10_samples-grid-10km.csv")) %>%
  filter(fieldsample_barcode %in% colnames(BC.df))

## Heatmap data frame
phylum.df <- read_csv(paste0(data.path, "/2024-03-19_phylum-heatmap-grid.csv"))

## Alpha diversity data frame
diversity.df <- read_csv(paste0(data.path, "/2024-03-19_alpha-diversity-fl.csv"))

## Gamma diversity data frame
gamma.df <- read_csv(paste0(data.path, "/2024-03-19_gamma-diversity-fl.csv"))

## Genus aggregated observational table
reduced.genus <-  data.table::fread(paste0(data.path, "/2024-03-07_arcbac-rarefaction-rel.csv")) %>%
  select(any_of(rep.10km.df$fieldsample_barcode), Kingdom:Genus) %>%
  filter(rowSums(across(where(is.numeric)))!=0)

## Perform Hellinger transformation
reduced.genus.hell <- reduced.genus %>%
  filter(!Genus == "Unclassified") %>%
  select(where(is.numeric), Genus) %>%
  column_to_rownames(var = "Genus") %>%
  t() %>%
  decostand(., method = "hellinger") %>%
  as.data.frame()


### Habitat dissimilarities
## Habitat parsing functions
mfd_db_ontology <- mfd_db %>%
  select(fieldsample_barcode, mfd_sampletype, mfd_areatype, mfd_hab1, mfd_hab2, mfd_hab3) %>%
  mutate(`sample type cumulative` = mfd_sampletype,
         `area type cumulative` = paste0(mfd_sampletype, ", ", mfd_areatype),
         `MFDO1 cumulative` = paste0(mfd_sampletype, ", ", mfd_areatype, ", ", mfd_hab1),
         `MFDO2 cumulative` = paste0(mfd_sampletype, ", ", mfd_areatype, ", ", mfd_hab1, ", ", mfd_hab2),
         `MFDO3 cumulative` = paste0(mfd_sampletype, ", ", mfd_areatype, ", ", mfd_hab1, ", ", mfd_hab2, ", ", mfd_hab3))

rep_10k_ontology <- rep.10km.df %>%
  select(fieldsample_barcode, mfd_sampletype, mfd_areatype, mfd_hab1, mfd_hab2, mfd_hab3) %>%
  mutate(`sample type cumulative` = mfd_sampletype,
         `area type cumulative` = paste0(mfd_sampletype, ", ", mfd_areatype),
         `MFDO1 cumulative` = paste0(mfd_sampletype, ", ", mfd_areatype, ", ", mfd_hab1),
         `MFDO2 cumulative` = paste0(mfd_sampletype, ", ", mfd_areatype, ", ", mfd_hab1, ", ", mfd_hab2),
         `MFDO3 cumulative` = paste0(mfd_sampletype, ", ", mfd_areatype, ", ", mfd_hab1, ", ", mfd_hab2, ", ", mfd_hab3))

get_av_dist <- function(mat, lab.df){
  labs <- unique(lab.df$lab)
  mat.to_return <- matrix(-1, nrow = length(labs), ncol = length(labs)) %>%
    `colnames<-`(labs) %>%
    `rownames<-`(labs)
  for(i in labs){
    samples.i <- lab.df %>% filter(lab == i) %>% pull(fieldsample_barcode)
    for(j in labs){
      samples.j <- lab.df %>% filter(lab == j) %>% pull(fieldsample_barcode)
      #print(mat[samples.i, samples.j])
      mat.to_return[i, j] <- mean(as.matrix(mat[samples.i, samples.j]))
    }
  }
  return(mat.to_return)
}

## Dissimilarity functions
BC.wrap <- function(counts, lab){
  mat <- parallelDist(counts,
                      method = "bray",
                      binary = F,
                      diag = T,
                      upper = T,
                      threads = 10) %>%
  as.matrix()
  mat.to_return <- get_av_dist(mat, lab) %>%
    as.dist()
  return(mat.to_return)
}

### Hab profiles
lab.df <- rep_10k_ontology %>%
             mutate(lab = `MFDO1 cumulative`) %>%
  mutate(across(lab,
                ~str_replace(.,
                             "Soil, Natural, Sclerophyllous scrub",
                             "Soil, Natural, Temperate heath and scrub"))) %>%
  filter(!lab %in% c("Water, Subterranean, Freshwater",
                     "Water, Urban, Sandfilter",
                     "Soil, Urban, Other",
                     "Sediment, Urban, Other",
                     "Soil, Urban, Roadside",
                     "Soil, Subterranean, Urban",
                     "Soil, Natural, NA")) %>%
  select(fieldsample_barcode, lab)

## MFDO1 summary (mean dissimilarity)
mfdo1.profiles.hell <- reduced.genus %>%
  filter(!Genus == "Unclassified") %>%
  select(where(is.numeric), Genus) %>%
  column_to_rownames(var = "Genus") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("fieldsample_barcode") %>%
  left_join(lab.df,
            by = "fieldsample_barcode") %>%
  filter(!is.na(lab)) %>%
  dplyr::select(-fieldsample_barcode) %>%
  pivot_longer(names_to = "Genera", values_to = "value", -lab) %>%
  group_by(lab, Genera) %>%
  summarise(value.mean = mean(value)) %>%
  pivot_wider(names_from = "Genera", values_from = "value.mean") %>%
  column_to_rownames("lab") %>%
  decostand(., method = "hellinger") %>%
  as.data.frame()

## Bootstrap SBHC large
ab.df <- rep_10k_ontology %>%
             mutate(lab = `MFDO1 cumulative`) %>%
  mutate(across(lab,
                ~str_replace(.,
                             "Soil, Natural, Sclerophyllous scrub",
                             "Soil, Natural, Temperate heath and scrub"))) %>%
  filter(!lab %in% c("Water, Subterranean, Freshwater",
                     "Water, Urban, Sandfilter",
                     "Soil, Urban, Other",
                     "Sediment, Urban, Other",
                     "Soil, Urban, Roadside",
                     "Soil, Subterranean, Urban",
                     "Soil, Natural, NA")) %>%
  select(fieldsample_barcode, lab)

## hclust with BC-wrapper
createHclustObject <- function(x)hclust(BC.wrap(counts = as.matrix(x),
                                                lab = (lab.df %>%
                                                         filter(fieldsample_barcode%in%rownames(reduced.genus.hell)))), "ave")

## bootstrap
set.seed(123)
b <- bootstrap::bootstrap(reduced.genus.hell,
                          fun = createHclustObject,
                          n = 10L,
                          mc.cores = 15)

## Write to output
write.csv(paste0(output.path, "/bootstrap_weigths.csv"))

## Plot results
hc <- createHclustObject(reduced.genus.hell)
plot(hc)

## Draw bootstrap values to corresponding node
bootlabels.hclust(hc, b, col="blue")


### Hab summarisation
## Average distance for full dataset
mfdo1.dist.full <- get_av_dist((BC.df %>% as.data.frame() %>% `rownames<-`(colnames(.))), (mfd_db_ontology %>%
                                    mutate(lab = `MFDO1 cumulative`) %>%
                                    mutate(across(lab,
                                                  ~str_replace(.,
                                                               "Soil, Natural, Sclerophyllous scrub",
                                                               "Soil, Natural, Temperate heath and scrub"))) %>%
                                    filter(!lab %in% c("Water, Subterranean, Freshwater",
                                                       "Water, Urban, Sandfilter",
                                                       "Soil, Urban, Other",
                                                       "Sediment, Urban, Other",
                                                       "Soil, Urban, Roadside",
                                                       "Soil, Subterranean, Urban",
                                                       "Soil, Natural, NA")) %>%
                                    select(fieldsample_barcode, lab)))

## Average distance for 10 km reference subset
mfdo1.dist.rep <- get_av_dist((BC.df %>% as.data.frame() %>% `rownames<-`(colnames(.))), (rep_10k_ontology %>%
                                    mutate(lab = `MFDO1 cumulative`) %>%
                                    mutate(across(lab,
                                                  ~str_replace(.,
                                                               "Soil, Natural, Sclerophyllous scrub",
                                                               "Soil, Natural, Temperate heath and scrub"))) %>%
                                    filter(!lab %in% c("Water, Subterranean, Freshwater",
                                                       "Water, Urban, Sandfilter",
                                                       "Soil, Urban, Other",
                                                       "Sediment, Urban, Other",
                                                       "Soil, Urban, Roadside",
                                                       "Soil, Subterranean, Urban",
                                                       "Soil, Natural, NA")) %>%
                                    select(fieldsample_barcode, lab)))


### Visualisation - Beta and fingerprint using geom_facet

## Transform phylum level 16S metagenome profiles
phylum.long <- (phylum.df %>%
                  pivot_longer(names_to = "hab", values_to = "rel.ab", -Phylum) %>%
                  mutate(rel.ab = as.numeric(rel.ab),
                         Phylum = factor(Phylum, levels = phylum.df$Phylum)))[c("hab", "rel.ab", "Phylum")]
phylum.long <- (phylum.df %>%
  mutate(hab = complex,
         rel.ab = abund))[c("hab", "rel.ab", "Phylum")]

## Sample count summary of 16S metagenome profiles
dd.MG.n <- rep_10k_ontology %>%
                                    mutate(lab = `MFDO1 cumulative`) %>%
                                    mutate(across(lab,
                                                  ~str_replace(.,
                                                               "Soil, Natural, Sclerophyllous scrub",
                                                               "Soil, Natural, Temperate heath and scrub"))) %>%
                                    filter(!lab %in% c("Water, Subterranean, Freshwater",
                                                       "Water, Urban, Sandfilter",
                                                       "Soil, Urban, Other",
                                                       "Sediment, Urban, Other",
                                                       "Soil, Urban, Roadside",
                                                       "Soil, Subterranean, Urban",
                                                       "Soil, Natural, NA")) %>%
  group_by(`MFDO1 cumulative`) %>%
  summarise(n=n())

## Dendogram of 16S metagenome profiles
p.rep <- hc %>%
  ggtree(ladderize = T) +
  theme_tree2() +
  scale_x_continuous(labels = abs)

b.df <- data.frame(x = -(hc$height),
                   bootstrap = b) %>%
  left_join((p.rep$data %>% select(x, node) %>% as.data.frame()),
            by = "x")

p.colors.rep <- p.rep %<+% data.frame(label = colnames(mfdo1.dist.rep),
                              beta.within = diag(mfdo1.dist.rep)) +
  geom_tiplab(aes(label = label), align = T, linesize = .5, as_ylab = T) +
  geom_tippoint(aes(fill = label, size = beta.within), shape = 21) +
  scale_fill_manual(values = mfdo1.palette) +
  guides(fill = "none") +
  new_scale_fill()

p.bootstrap.rep <- p.colors.rep %<+% b.df +
  geom_point(data = (. %>% filter(!isTip)), aes(fill=bootstrap), shape = 23, size = 3, color = "black") +
  scale_fill_distiller(type = "seq", palette = "Greys", direction = 1,  na.value = NA) +
  scale_fill_stepsn(colors = rev(grey.colors(8))[c(1, 8)], na.value = NA,
                       limits = c(0.7, 1),
                       breaks = c(0.7, 0.8, 0.9, 1.0),
                       labels = c("0.7", "0.8", "0.9", "1.0"),
                       name = "Bootstrap") +
  scale_color_continuous(na.value = NA) +
  new_scale_fill() +
  new_scale_color()

## Phylum level 16S metagenome profiles
p.fingerprint <- p.bootstrap.rep +
  geom_facet(panel = "Fingerprint",
             data = phylum.long,
             geom = geom_tile,
             aes(x = Phylum, fill = rel.ab),
             color = "black",
            lwd = .25,
            linetype = 1) +
  scale_x_discrete() +
  scale_fill_viridis_c() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1)) +
  guides(size = "none") +
  theme(legend.position="bottom")

## Sample count summary of 16S metagenome profiles
p.n.MG <- p.bootstrap.rep +
  geom_facet(panel = "# MG",
             data = dd.MG.n,
             geom = geom_text,
             aes(x = 1, label = paste0("n = ", n)))

## Alpha diversity based on FL16S amplicons
p.a_div <- p.n.MG +
  new_scale_fill() +
  geom_facet(panel = "Alpha diversity",
             data = diversity.df[c("complex", "observed_species")],
             geom = geom_jitter,
             orientation = "y",
             aes(x = observed_species, color = label),
             size = 1.0) +
  scale_fill_manual(values = mfdo1.palette) +
  scale_color_manual(values = mfdo1.palette) +
  geom_facet(panel = "Alpha diversity",
             data = diversity.df[c("complex", "observed_species")],
             geom = geom_boxplot,
             orientation = "y",
             aes(x = observed_species, fill = label),
             color = "black",
             width = 0.8,
             alpha = 0.6,
             outlier.size = -1) +
  scale_fill_manual(values = mfdo1.palette) +
  guides(fill = "none",
         color = "none")

## Gamma diversity based on FL16S amplicons
p.c_div <- p.a_div +
  geom_facet(pane = "Gamma diversity",
             data = (gamma.df %>% filter(Diversity == "Shannon diversity")),
             geom = geom_errorbar,
             aes(x= Estimator, xmin = LCL, xmax = UCL),
             color = "black",
             width = 0.5,
             orientation = "y") +
  geom_facet(panel = "Gamma diversity",
             data = (gamma.df %>% filter(Diversity == "Shannon diversity")),
             geom = geom_point,
             aes(x = Estimator, fill = label),
             pch = 21,
             color = "black",
             width = 1,
             alpha = 1,
             orientation = "y")

## Sample count summary of FL16S amplicons
p.n.16S <- p.c_div +
  geom_facet(panel = "# FL16S",
             data = (gamma.df %>% filter(Diversity == "Shannon diversity")),
             geom = geom_text,
             aes(x = 1, label = paste0("n = ", size))) +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))

## Render and save plots in svg format
p.n.16S

ggsave(p.n.16S, filename = paste0(output.path, "/tree_cont_scale.svg"), width = 12, height = 4, dpi = "retina")

p.fingerprint

ggsave(p.fingerprint, filename = paste0(output.path, "/tree_disc_scale.svg"), width = 8, height = 6, dpi = "retina")

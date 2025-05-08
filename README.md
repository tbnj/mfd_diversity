# mfd_diversity
The scripts in this repository are part of the [Microflora Danica project](https://github.com/cmc-aau/mfd_wiki/wiki). 

Be advised, that for the metagenomic-derived data, the term "OTU" is only used due to format requirements by ampvis2, and they do not represent classical OTUs. 
The generated profiles can be thought of as taxonomic bins. 

All scritps are needed to make the combined figure panel at the end. 

## Scripts
### Amplicon 16S data 
`scripts/diversity_FL16S.R` calculates alpha diversity using the near full-length 16S OTUs from the UMI Nanopore data across the MFDO1 ontology level. It performs the statistical analysis on alpha diveristy and combines the alpha with the gamma diversity calculated for the same data.

### Metagenomic 16S data 
`scripts/heatmap_16S.R` recreates the prokaryotic fingerprint profile across the MFDO1 ontology level. The script also renders the figures used in the manuscript. 


`scripts/permanova_contrasts_16S.R` performs the PERMANOVA and contrasts abalysis across the MFDO1 ontology level.


`scripts/ordination_16S.R` calculates the Hellinger-transformed Bray-Curtis dissimlarity matrix and recreates the prokaryotic ordination panels with addition of results from ANOSIM and PERMANOVA. 


`scripts/combine_diversity_FL16S.R` combines the different aspects of diveristy and calcualtes the between-group Hellinger-transformed Bray-Curtis dissimilarites across the MFDO1 ontology level. The script also renders the figures used in the manuscript. 

## Data
The scripts rely on data files available from the MFD [github](https://github.com/cmc-aau/mfd_metadata) and the MFD Zenodo [repo](https://zenodo.org/records/12605769). 

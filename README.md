# LiverwortAssociatedDiazotrophs
This repository host the metadata, microbial community matrix, and the scripts used for the analyses that support the findings of the study:

Adriel M Sierra, Dennis Alejandro Escolástico-Ortiz, Charles E Zartman, Nicolas Derome, Connie Lovejoy, Juan Carlos Villarreal A, Assembly and co-occurrence networks of nitrogen-fixing bacteria associated with epiphyllous liverworts in fragmented tropical forests, ISME Communications, Volume 5, Issue 1, January 2025, ycaf173, https://doi.org/10.1093/ismeco/ycaf173

### Abstract
Understanding the spatial dynamics of plant-associated microbial communities is increasingly urgent in the context of habitat loss and the biodiversity crisis. However, the influence of reduced habitat size and connectivity on the assembly mechanisms underlying microbial associations is fundamental to advancing microbial ecology and conservation. In the Brazilian Amazon, we investigated nitrogen-fixing (diazotrophic) bacterial communities associated with two epiphyllous liverworts, Cololejeunea surinamensis and Radula flaccida, across 11 forest sites within the Biological Dynamics of Forest Fragments Project landscape. Using amplicon sequencing targeting the nitrogenase gene (nifH), we characterized diazotroph community diversity, inferred assembly mechanisms through null models, and analyzed co-occurrence network structure. Host-specific associations were evident: C. surinamensis predominantly hosted Hassallia, while R. flaccida was primarily associated with Fischerella. Despite habitat fragmentation, diazotrophic richness and composition remained similar across habitats of different sizes, consistent with strong homogenizing dispersal. Network analyses revealed that smaller fragments harbored more modular communities with fewer module hubs, pronounced shifts in key species relative abundance, and reduced network robustness. Our findings underscore the influence of habitat size on the stability of liverwort-associated diazotrophs, with smaller fragments exhibiting lower key species specificity and disruption of microbe-microbe interactions. Our results emphasize the importance of conserving large, connected forest habitats to maintain the functional integrity of phyllosphere N-fixing microbiota.

### Data availability
Raw sequence data were deposited in the NCBI Sequence Read Archive (SRA) with their respective accession numbers under the BioProject: PRJNA1169898.

### Data
Phyloseq objects with sample metadata, microbial community matrix of amplicon sequence variants (ASVs), ASVs taxonomic classification assigned using the adapted nifH ARB database v1.0.3 accessed in 2022 (Heller et al. 2014; Moynihan 2020), ASVs phylogenetic tree and full refseq alignment.  

- [Complete diazotroph community matrix](https://github.com/adrielmsierra/LiverwortAssociatedDiazotrophs/blob/main/01_Datasets/nifH_ps_complete_dataset_tree_abr23.rds). This is the complete dataset, including samples with uneven sequencing depths and without prior rarefaction, retaining all detected ASVs, including singletons.
- [Rarefied dataset](https://github.com/adrielmsierra/LiverwortAssociatedDiazotrophs/blob/main/01_Datasets/nifH_ps_rarefied_dataset_tree_jun24.rds) to an even sequencing depth of 514 reads for downstream analysis. The dataset comprised 146 samples, with a mean sequencing depth of 7003.9 reads (min = 571; max = 41 496) and a total of 4897 ASVs (3267 ASVs from 70 samples of C. surinamensis, and 2033 ASVs from 76 samples of R. flaccida).
- [Host-specific diazotrophs](https://github.com/adrielmsierra/LiverwortAssociatedDiazotrophs/blob/main/01_Datasets/Bryophyte_nifH_coremicrobiota_Jun24.rds). Filtered dataset corresponding to 75 ASVs. 

### Scripts
Scripts in this repository are organized in the following directories:
- Relative abundance of diazotrophic microbiota
- Host-specific diazotrophic microbiota
- Analyses of diazotrophic community diversity
  * Alpha Diversity
  * Beta diversity and ordination
- Metacommunity assembly processes
- Network analyses and detection of module hubs

### References
- Heller   P, Tripp   HJ, Turk-Kubo   K. et al.  ARBitrator: a software pipeline for on-demand retrieval of auto-curated nifH sequences from GenBank. Bioinformatics  2014;30:2883–90. 10.1093/bioinformatics/btu417
- Moynihan   MA. nifHdada2 GitHub repository. Zenodo  2020. 10.5281/zenodo.3958370

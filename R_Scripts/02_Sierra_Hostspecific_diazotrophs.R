###################################################################################
#                           Bryophyte NifH (Diazotrophs)
#                       Host-specific diazotrophic microbiota
#
#
#                                  Adriel M. Sierra
#
#
#                                  Jun 2024
#
###################################################################################

# Dataset
nifH.ps.com.dataset = readRDS("./02_Phyloseq_data_analysis/nifH_ps_complete_dataset_tree_abr23.rds")
# rarefied dataset 
nifh.dataset.rarefied = readRDS("./02_Phyloseq_data_analysis/nifH_ps_rarefied_dataset_tree_jun24.rds")

###core microbiome from complete dataset. script from https://bocasbiome.github.io/wf1.html
#Step one is to generate an ASV table without row names.

data.IndVal.host <- data.frame(otu_table(nifh.dataset.rarefied))

data.IndVal.host <- data.frame(otu_table(nifH.ps.com.dataset))
data.IndVal.host
data.IndVal.ASV <- tibble::remove_rownames(data.IndVal.host)

##Indicator value group file
library(stringr)
data.IndVal.group <- data.frame(sample_data(nifH.ps.com.dataset)) %>%
  select(Species)
data.IndVal.group$Species = as.factor(c("Radula flaccida", "Cololejeunea surinamensis"))
#data.IndVal.group$Status <- data.IndVal.group$Species
data.IndVal.group <- data.IndVal.group[, c(2,1)]
data.IndVal.group$Status <- str_replace(data.IndVal.group$Status, "Radula flaccida", "1")
data.IndVal.group$Status <- str_replace(data.IndVal.group$Status, "Cololejeunea surinamensis", "2")
data.IndVal.group <- tibble::rownames_to_column(data.IndVal.group, "Label")
data.IndVal.group$Status <- as.integer(c(rep("1", 76), rep("2", 70)))
data.IndVal.group$Status <-  as.integer(data.IndVal.group$Status)
data.IndVal.group$Species <-  as.character(data.IndVal.group$Species)

# calculate the indicator values. We set a seed for reproducibility and saved a table of results.
install.packages("labdsv")
library(labdsv)
detach("package:nlme", unload=TRUE, character.only=TRUE)
#unloadNamespace("nlme")
#unloadNamespace("ape")
#unloadNamespace("phangorn")

#install.packages("labdsv")
set.seed(1280)
iva <- indval(data.IndVal.ASV, data.IndVal.group$Status)
#Table of the significant indicator species at p= 0.01
gr <- iva$maxcls[iva$pval <= 0.01]
iv <- iva$indcls[iva$pval <= 0.01]
pv <- iva$pval[iva$pval <= 0.01]
fr <- apply(data.IndVal.ASV > 0, 2, sum)[iva$pval <= 0.01]
indval.out <- data.frame(group = gr, indval = iv, pval = pv, freq = fr)
indval.out <- indval.out[order(indval.out$group, -indval.out$indval),]
indval.out
dim(indval.out)


indval.out$prob.corrected = p.adjust(indval.out$pval, "bonferroni")
write.csv(indval.out,
          file = "./Tables_pub/IndVal_microbiome_host_fromRarefiedData_bonf_jun24.csv")
# Export phyloseq of core microbiota
ind_asvs_list <- row.names(indval.out)
ps.indv.core <- subset_taxa(nifH.ps.com.dataset,
                            rownames(tax_table(nifH.ps.com.dataset)) %in%
                              ind_asvs_list)

sample_sums(ps.indv.core)
taxa_sums(ps.indv.core)
table(sample_sums(ps.indv.core)==0)
#core microbiota absent in 58 samples
nifH.core.ps.ex <- subset_samples(ps.indv.core, sample_sums(ps.indv.core)!=0)
saveRDS(nifH.core.ps.ex, "Bryophyte_nifH_coremicrobiota_Jun24.rds")
nifH.core.ps.ex <- readRDS("Bryophyte_nifH_coremicrobiota_Jun24.rds")
table(nifH.core.ps.ex@sam_data$Species,nifH.core.ps.ex@sam_data$Size)
# remove samples without taxa
nifH.core.ps.ex <- subset_samples(nifH.core.ps, sample_sums(nifH.core.ps)!=0)
taxa_sums(nifH.core.ps)
#host specific core microbiome
rad_ind_asvs <- indval.out[indval.out$group == c(1,2), ]
rad_ind_asvs_list <- row.names(rad_ind_asvs)
ps.indv.rad.core <- subset_taxa(nifH.ps.com.dataset,
                                rownames(tax_table(nifH.ps.com.dataset)) %in%
                                  rad_ind_asvs_list)
ps.indv.rad.core
tax_table(ps.indv.rad.core)
#ps.indv.rad.core@phy_tree = NULL
col_ind_asvs <- indval.out[indval.out$group == 2, ]
col_ind_asvs_list <- row.names(col_ind_asvs)
ps.indv.col.core <- subset_taxa(nifH.ps.com.dataset,
                                rownames(tax_table(nifH.ps.com.dataset)) %in%
                                  col_ind_asvs_list)
ps.indv.col.core
tax_table(ps.indv.col.core)
#ps.indv.col.core@phy_tree = NULL
#merge core microbiome phyloseq
nifH.core.ps = merge_phyloseq(ps.indv.rad.core, ps.indv.col.core)
sample_sums(nifH.core.ps)
taxa_sums(nifH.core.ps)
table(sample_sums(nifH.core.ps)==0)
# remove samples without taxa
nifH.core.ps.ex <- subset_samples(nifH.core.ps, sample_sums(nifH.core.ps)!=0)
taxa_sums(nifH.core.ps)

#extract sequences to compared on BLAST

remotes::install_github("alexpiper/seqateurs")
library(seqateurs)
core_taxa_names = taxa_names(nifH.core.ps) 
ps_to_fasta(nifH.core.ps, seqnames = "core_taxa_names", width = 1000)


meco_nifHdataset.rare <- phyloseq2meco(nifh.dataset.rarefied)
meco_nifHdataset.rare
class(meco_nifHdataset.rare)
rownames(meco_nifHdataset.rare$tax_table) -> meco_nifHdataset.rare$tax_table$Species

#LefSe
# use default parameters to calculate reltaive abundance
meco_nifHdataset.rare$sample_sums() %>% range
meco_nifHdataset.rare$cal_abund(rel = TRUE) 



t1 <- trans_diff$new(dataset = meco_nifHdataset.rare, method = "lefse", group ="Host", alpha = 0.01, lefse_subgroup = NULL, taxa_level = "Genus")
t1 <- trans_diff$new(dataset = meco_nifHdataset.rare, method = "lefse", group ="Class.size", alpha = 0.01, lefse_subgroup = NULL, taxa_level = "Species")
t1 <- trans_diff$new(dataset = meco_nifHdataset.rare, method = "lefse", group = "Size", taxa_level = "Genus")
# see t1$res_diff for the result
# From v0.8.0, threshold is used for the LDA score selection.
t1$plot_diff_bar(threshold = 2)
# we show 20 taxa with the highest LDA (log10)
t1$plot_diff_bar(use_number = 1:30, width = 0.8)

t1$res_diff

# Subseting dataset per host species
radu.rare = subset_samples(nifH.ps.com.dataset, Host %in% c("Radula flaccida"))
colo.rare = subset_samples(nifH.ps.com.dataset, Host %in% c("Cololejeunea surinamensis"))
radu.rare = prune_taxa(taxa_sums(radu.rare) > 1, radu.rare)
colo.rare = prune_taxa(taxa_sums(colo.rare) > 1, colo.rare)

meco_radu.rare <- phyloseq2meco(radu.rare)
meco_radu.rare

meco_colo.rare <- phyloseq2meco(colo.rare)
meco_colo.rare
# use default parameters to calculate relative abundance
meco_radu.rare$sample_sums() %>% range
meco_radu.rare$cal_abund(rel = TRUE) 

meco_colo.rare$sample_sums() %>% range
meco_colo.rare$cal_abund(rel = TRUE) 

radu.rare.t1 <- trans_diff$new(dataset = meco_radu.rare, method = "lefse", group ="Size", alpha = 0.05, lefse_subgroup = NULL, taxa_level = "Species")
colo.rare.t1 <- trans_diff$new(dataset = meco_colo.rare, method = "lefse", group ="Size", alpha = 0.05, lefse_subgroup = NULL, taxa_level = "Species")

# see t1$res_diff for the result
# From v0.8.0, threshold is used for the LDA score selection.
radu.rare.t1$plot_diff_bar(threshold = 2)
colo.rare.t1$plot_diff_bar(threshold = 2)
# we show 20 taxa with the highest LDA (log10)
radu.rare.t1$plot_diff_bar(use_number = 1:30, width = 0.8)
colo.rare.t1$plot_diff_bar(use_number = 1:30, width = 0.8)

radu.rare.t1$res_diff
colo.rare.t1$res_diff

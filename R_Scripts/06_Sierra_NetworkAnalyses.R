###################################################################################
#                           Bryophyte NifH (Diazotrophs)
#                     Network analyses and detection of module hubs
#
#                                  Adriel M. Sierra
#
#
#                                  Jun 2024
#
###################################################################################
# Dataset
# rarefied dataset 
nifh.dataset.rarefied = readRDS("./02_Phyloseq_data_analysis/nifH_ps_rarefied_dataset_tree_jun24.rds")

# Subset datasets
rarefied.dt.small = subset_samples(nifh.dataset.rarefied, Class.size %in% c("Small"))
rarefied.dt.small = prune_taxa(taxa_sums(rarefied.dt.small) > 1, rarefied.dt.small)
rarefied.dt.large = subset_samples(nifh.dataset.rarefied, Class.size %in% c("Large"))
rarefied.dt.large = prune_taxa(taxa_sums(rarefied.dt.large) > 1, rarefied.dt.large)



# to work with the package microeco we need to convert the datafrom phyloseq to microtable object
nifH.ps.dataset # data
meco_rarefied.dt.small <- phyloseq2meco(rarefied.dt.small)
meco_rarefied.dt.small
meco_rarefied.dt.large <- phyloseq2meco(rarefied.dt.large)
meco_rarefied.dt.large
#small network
network.rare.small <- trans_network$new(dataset = meco_rarefied.dt.small, cor_method = "spearman", filter_thres = 0.001)
# construct network; require igraph package
network.rare.small$cal_network(COR_p_thres = 0.01, COR_optimization = TRUE)
# use arbitrary coefficient threshold to contruct network
#network.rare.small$cal_network(COR_p_thres = 0.01, COR_cut = 0.7)
# return t1$res_network
network.rare.small$cal_module(method = "cluster_fast_greedy")
#plot network
network.rare.small$plot_network(method = "ggraph", node_color = "Genus")
network.rare.small$plot_network(method = "networkD3", node_color = "module")
network.rare.small$plot_network(method = "networkD3", node_color = "Phylum")

# calculate network attributes
network.rare.small$cal_network_attr()
network.rare.small$res_network_attr

# get node properties
network.rare.small$get_node_table(node_roles = TRUE)
# return t1$res_node_table
network.rare.small$res_node_table
# get edge properties
network.rare.small$get_edge_table()
# return t1$res_edge_table 
network.rare.small$res_edge_table

network.rare.small$get_adjacency_matrix()
# return t1$res_adjacency_matrix
network.rare.small$res_adjacency_matrix
# add_label = TRUE can be used to directly add text label for points
network.rare.small$plot_taxa_roles(use_type = 1)
network.rare.small$plot_taxa_roles()
# plot node roles with phylum information
network.rare.small$plot_taxa_roles(use_type = 2)

network.rare.small$cal_sum_links(taxa_level = "Genus")
#install.packages("chorddiag")
library(chorddiag)
# interactive visualization; require chorddiag package; see https://github.com/mattflor/chorddiag
network.rare.small$plot_sum_links(plot_pos = TRUE, plot_num = 10, color_values = RColorBrewer::brewer.pal(10, "Paired"))
# From v1.2.0, method = "circlize" is available for conveniently saving the static plot
# If circlize package is not installed, first run: install.packages("circlize")
network.rare.small$plot_sum_links(method = "circlize", transparency = 0.2, annotationTrackHeight = circlize::mm_h(c(5, 5)))



#large network
network.rare.large <- trans_network$new(dataset = meco_rarefied.dt.large, cor_method = "spearman", filter_thres = 0.001)
# construct network; require igraph package
network.rare.large$cal_network(COR_p_thres = 0.01, COR_optimization = TRUE)
# use arbitrary coefficient threshold to construct network
network.rare.large$cal_network(COR_p_thres = 0.01, COR_cut = 0.7)
# return t1$res_network
network.rare.large$res_network
network.rare.large$cal_module(method = "cluster_fast_greedy")
#plot network
network.rare.large$plot_network(method = "ggraph", node_color = "module")
network.rare.large$plot_network(method = "networkD3", node_color = "module")
network.rare.large$plot_network(method = "networkD3", node_color = "Phylum")

# calculate network attributes
network.rare.large$cal_network_attr()
network.rare.large$res_network_attr
network.rare.small$res_network_attr

# get node properties
network.rare.large$get_node_table(node_roles = TRUE)
# return t1$res_node_table

# get edge properties
network.rare.large$get_edge_table()
# return t1$res_edge_table 
network.rare.large$get_adjacency_matrix()
# return t1$res_adjacency_matrix

# add_label = TRUE can be used to directly add text label for points
network.rare.large$plot_taxa_roles(use_type = 1)

# plot node roles with phylum information
network.rare.large$plot_taxa_roles(use_type = 2)

network.rare.large$cal_sum_links(taxa_level = "Genus")
#install.packages("chorddiag")
# interactive visualization; require chorddiag package; see https://github.com/mattflor/chorddiag
network.rare.large$plot_sum_links(plot_pos = TRUE, plot_num = 10, color_values = RColorBrewer::brewer.pal(10, "Paired"))
# From v1.2.0, method = "circlize" is available for conveniently saving the static plot
# If circlize package is not installed, first run: install.packages("circlize")
network.rare.large$plot_sum_links(method = "circlize", transparency = 0.2, annotationTrackHeight = circlize::mm_h(c(5, 5)))




#install.packages("rgexf")
library(rgexf)

#### comparing networks
# install the required packages
# aplot: one dependency of the trans_venn class of microeco package
# agricolae: for Duncan's new multiple range test
packages <- c("meconetcomp", "rgexf", "pheatmap", "aplot", "agricolae")
# Now check or install
for(x in packages){
  if(!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
  }
}

library(microeco)
library(meconetcomp)
# use pipe operator in magrittr package
library(magrittr)
library(igraph)
library(ggplot2)
library(devtools)
library(remotes)
#remotes::install_github("zdk123/SpiecEasi")
library(SpiecEasi)
# dataset with amplicon seq
meco_nifHdataset.rare


# first create a list
bryo_nifh_network <- list()
# select samples of "Radula flaccida small" group
# use clone to get a deep copy of our microtable (R6 object)
tmp <- clone(meco_nifHdataset.rare)
# change sample_table directly
tmp$sample_table %<>% subset(Host=="Radula flaccida")%<>%subset(Class.size=="Small")
# trim all files in the object
tmp$tidy_dataset()
# use filter_thres parameter to filter the feature with low relative abundance
tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.0005)
#tmp <- trans_network$new(dataset = tmp, cor_method = "sparcc", use_sparcc_method = "SpiecEasi", filter_thres = 0.001)

# COR_p_thres represents the p value threshold
# COR_cut denotes the correlation coefficient thresholdCOR_cut = 0.6,
tmp$cal_network(COR_p_thres = 0.01, COR_optimization = TRUE,
                COR_optimization_low_high = c(0.01, 0.8))#p-value adjustment for default method
tmp$res_network
# put the network into the list
bryo_nifh_network$Rflaccida_Small <- tmp

# select samples of "Radula flaccida large" group
# use clone to get a deep copy of our microtable (R6 object)
tmp <- clone(meco_nifHdataset.rare)
# change sample_table directly
tmp$sample_table %<>% subset(Host=="Radula flaccida") %<>%subset(Class.size=="Large")
# trim all files in the object
tmp$tidy_dataset()
# use filter_thres parameter to filter the feature with low relative abundance
tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.0005)
# COR_p_thres represents the p value threshold
# COR_cut denotes the correlation coefficient threshold
tmp$cal_network(COR_p_thres = 0.01,  COR_optimization = TRUE,
                COR_optimization_low_high = c(0.01, 0.8))
tmp$res_network
# put the network into the list
bryo_nifh_network$Rflaccida_Large <- tmp

# select samples of "Radula flaccida small" group
# use clone to get a deep copy of our microtable (R6 object)
tmp <- clone(meco_nifHdataset.rare)
# change sample_table directly
tmp$sample_table %<>% subset(Host=="Cololejeunea surinamensis") %<>%subset(Class.size=="Small")
# trim all files in the object
tmp$tidy_dataset()
# use filter_thres parameter to filter the feature with low relative abundance
tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.0005)
# COR_p_thres represents the p value threshold
# COR_cut denotes the correlation coefficient threshold
tmp$cal_network(COR_p_thres = 0.01,  COR_optimization = TRUE,
                COR_optimization_low_high = c(0.01, 0.8))
tmp$res_network
# put the network into the list
bryo_nifh_network$Csurinamensis_Small <- tmp

# select samples of "Radula flaccida large" group
# use clone to get a deep copy of our microtable (R6 object)
tmp <- clone(meco_nifHdataset.rare)
# change sample_table directly
tmp$sample_table %<>% subset(Host=="Cololejeunea surinamensis") %<>%subset(Class.size=="Large")
# trim all files in the object
tmp$tidy_dataset()
# use filter_thres parameter to filter the feature with low relative abundance
tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.0005)
# COR_p_thres represents the p value threshold
# COR_cut denotes the correlation coefficient threshold
tmp$cal_network(COR_p_thres = 0.01,   COR_optimization = TRUE,
                COR_optimization_low_high = c(0.01, 0.8))
tmp$res_network
# put the network into the list
bryo_nifh_network$Csurinamensis_Large <- tmp
bryo_nifh_network$Rflaccida_Small
str(bryo_nifh_network)
# network modularity
# partition modules for all the networks in the list.
bryo_nifh_network %<>% cal_module(undirected_method = "cluster_fast_greedy")
bryo_nifh_network$Rflaccida_Small$res_network_attr
# network topological attributes
tmp <- cal_network_attr(bryo_nifh_network)
bryo_nifh_network$Rflaccida_Small$res_network_attr
#View(tmp)
# tmp is a data.frame object

#node and edge properties extraction for all networks
bryo_nifh_network %<>% get_node_table(node_roles = TRUE) %>% get_edge_table
bryo_nifh_network$Rflaccida_Small$
  # compare nodes across networks
  # obtain the node distributions by searching the res_node_table in the object
  tmp <- node_comp(bryo_nifh_network, property = "name")

# obtain nodes intersection
tmp1 <- trans_venn$new(tmp, ratio = "numratio")
g1.n <- tmp1$plot_venn(fill_color = FALSE)
ggsave("bryo_amp_node_overlap.pdf", g1.n, width = 7, height = 6, dpi = 350)
# calculate jaccard distance to reflect the overall differences of networks
tmp$cal_betadiv(method = "jaccard")
tmp$beta_diversity$jaccard

# get the edge distributions across networks
tmp <- edge_comp(bryo_nifh_network)
# obtain edges intersection
tmp1 <- trans_venn$new(tmp, ratio = "numratio")
g1.e <- tmp1$plot_venn(fill_color = FALSE)
ggsave("bryo_amp_edge_overlap.pdf", g1.e, width = 7, height = 6)
# calculate jaccard distance
tmp$cal_betadiv(method = "jaccard")
tmp$beta_diversity$jaccard

#11.6 Extract overlapped edges of networks to a new network
# first obtain edges distribution and intersection
tmp <- edge_comp(bryo_nifh_network)
tmp1 <- trans_venn$new(tmp)
# convert intersection result to a microtable object
tmp2 <- tmp1$trans_comm()
tmp2$sample_table
# extract the intersection of all the four networks ("IW", "TW" and "CW")
# please use colnames(tmp2$otu_table) to find the required name
tmp2$otu_table
Intersec_rfla <- subset_network(bryo_nifh_network, venn = tmp2, name = "Rflaccida_Small&Rflaccida_Large")

Intersec_csur <- subset_network(bryo_nifh_network, venn = tmp2, name = "Csurinamensis_Small&Csurinamensis_Large")
# Intersec_all is a trans_network object
# for example, save Intersec_all as gexf format
Intersec_rfla$plot_network(method = "ggraph", node_color = "module")
Intersec_rfla$plot_network(method = "networkD3", node_color = "module")

Intersec_rfla$save_network("Intersec_rfla.gexf")
Intersec_csur$plot_network(method = "ggraph", node_color = "module")
Intersec_csur$save_network("Intersec_csur")

ggarrange(Intersec_rfla$plot_network(method = "ggraph", node_color = "phylum"), Intersec_csur$plot_network(method = "ggraph", node_color = "module"),
          ncol = 1, nrow = 2, labels = c("R. flaccida", "C. surinamensis"))

# safe networks for gephi
bryo_nifh_network$Rflaccida_Small$save_network(filepath = "Rflaccida_small_network.gexf")
bryo_nifh_network$Rflaccida_Large$save_network(filepath = "Rflaccida_large_network.gexf")
bryo_nifh_network$Csurinamensis_Small$save_network(filepath = "Csurinamensis_small_network.gexf")
bryo_nifh_network$Csurinamensis_Large$save_network(filepath = "Csurinamensis_large_network.gexf")


#11.7 Compare phylogenetic distances of paired nodes in edges
# I decide not to include this part
# filter useless features to speed up the calculation
node_names <- unique(unlist(lapply(bryo_nifh_network, function(x){colnames(x$data_abund)})))
filter_bryo_amp <- microeco::clone(meco_nifHdataset.rare)
filter_bryo_amp$otu_table <- filter_bryo_amp$otu_table[node_names, ]
filter_bryo_amp$tidy_dataset()
# obtain phylogenetic distance matrix
phylogenetic_distance_bryo <- as.matrix(cophenetic(filter_bryo_amp$phylo_tree))
# use both the positive and negative labels
tmp <- edge_node_distance$new(network_list = bryo_nifh_network, dis_matrix = phylogenetic_distance_bryo, label = c("+", "-"))
tmp$cal_diff(method = "anova")
tmp$res_diff
tmp$data_table

# visualization
g1 <- tmp$plot(violinplot_add = "none", add_sig = TRUE, add_sig_text_size = 5) + ylab("Phylogenetic distance")
ggsave("bryo_amp_phylo_distance.pdf", g1, width = 7, height = 6)

# show different modules with at least 10 nodes and positive edges
tmp <- edge_node_distance$new(network_list = bryo_nifh_network, dis_matrix = phylogenetic_distance_bryo, 
                              label = "+", with_module = TRUE, module_thres = 10)
tmp$cal_diff(method = "anova")
g1 <- tmp$plot(boxplot_add = "none", add_sig = TRUE, add_sig_text_size = 5) + ylab("Phylogenetic distance")
ggsave("bryo_amp_phylo_distance_modules.pdf", g1, width = 8, height = 6)

#11.8 Compare node sources of edges across networks
bryo_amp_network_edgetax <- edge_tax_comp(bryo_nifh_network, taxrank = "Phylum", label = "+", rel = TRUE)
# filter the features with small number
bryo_amp_network_edgetax <- bryo_amp_network_edgetax[apply(bryo_amp_network_edgetax, 1, mean) > 0.01, ]
colnames(bryo_amp_network_edgetax) = c("R. flaccida Small", "R. flaccida Large", "C. surinamensis Small", "C. surinamensis Large")
#col_annot = data.frame(c("R. flaccida Small", "R. flaccida Large", "C. surinamensis Small", "C. surinamensis Large"))
# visualization
g1 <- pheatmap::pheatmap(bryo_amp_network_edgetax,   display_numbers = TRUE)
ggsave("bryo_amp_edge_tax_comp_phylum.pdf", g1, width = 7, height = 7, dpi = 350)

#Compare topological properties of sub-networks
# calculate global properties of all sub-networks
tmp <- subnet_property(bryo_nifh_network)
# then prepare the data for the correlation analysis
# use sample names (second column) as rownames
rownames(tmp) <- tmp[, 2]
# delete first two columns (network name and sample name)
tmp <- tmp[, -c(1:2)]
# load ready-made abiotic factor and diversity table
#data
class(meco_nifHdataset.rare$sample_table$Size)

tmp1 <- trans_env$new(dataset = meco_nifHdataset.rare, add_data = meco_nifHdataset.rare$sample_table)
tmp1$cal_cor(use_data = meco_nifHdataset.rare, by_group = "Size", add_abund_table = tmp, cor_method = "spearman")
tmp1$res_cor
# generate correlation heatmap
g1 <- tmp1$plot_cor()
ggsave("bryo_amp_subnet_property.pdf", g1, width = 11, height = 5)

#Robustness of network

tmp <- robustness$new(bryo_nifh_network, remove_strategy = c("edge_rand", "edge_strong", "node_rand", "node_degree_high",'node_hub'), 
                      remove_ratio = seq(0, 0.99, 0.1), measure = c("Eff", "Eigen", "Pcr"), run = 10)
View(tmp$res_table)
View(tmp$res_summary)
tmp$plot(linewidth = 1)+theme_classic2()+
  theme(text = element_text(size = 14))

theme_pubr()
theme_minimal()


tmp1 <- tmp$res_table %>% .[.$remove_strategy == "node_rand" & .$measure == "Eigen", ]
library(microeco)
library(ggplot2)
theme_set(theme_bw())

# generate correlation heatmap
g1 <- tmp1$plot_cor()
t1 <- trans_env$new(dataset = NULL, add_data = tmp1)
t1$dataset$sample_table <- t1$data_env
t1$plot_scatterfit()
t1$plot_scatterfit(x = "remove_ratio", y = "value", type = "cor", group = "Network") + 
  xlab("Ratio of randomly removed nodes") + ylab("Network connectivity") +theme_classic()+ stat_cor(label.y = c(0.022,0.023,0.024,0.025), size = 6) +theme(axis.title = element_text(size = 15))

View(bryo_nifh_network$Rflaccida_Large$res_node_table) 

# Vulnerability of nodes
vul_table <- vulnerability(bryo_nifh_network)
head(sort(vul_table$vulnerability))
write.csv(vul_table, "network_vulnerability_aug24.csv")
ggplot(vul_table,aes(vulnerability, fill=Network))+
  geom_density()
#Cohesion
t1 <- cohesionclass$new(bryo_nifh_network)
View(t1$res_list$sample)
View(t1$res_list$feature)
summary(t1$cal_diff(method = "anova"))
t1$res_diff
summary(aov(as.numeric(t1$res_list$sample$c_pos) ~ as.factor(t1$res_list$sample$network)))
TukeyHSD(aov(as.numeric(t1$res_list$sample$c_pos) ~ as.factor(t1$res_list$sample$network)))

t1$plot(measure = "r_pos")
t1$plot(measure = "r_pos", boxplot_add = "none")
write.csv(t1$res_list$sample, "Cohesion_results_Aug24.csv")
t1$res_list$sample > 0
table(t1$res_list$feature$r_neg > 0)
#Plot visualization

#plot network
ggarrange(bryo_nifh_network$Rflaccida_Small$plot_network(method = "ggraph", node_color = "Phylum"),
          bryo_nifh_network$Rflaccida_Large$plot_network(method = "ggraph", node_color = "Phylum"),
          bryo_nifh_network$Csurinamensis_Small$plot_network(method = "ggraph", node_color = "Phylum"),
          bryo_nifh_network$Csurinamensis_Large$plot_network(method = "ggraph", node_color = "Phylum"), 
          nrow = 2, ncol = 2, labels = c("R. flaccida - Small", "R. flaccida - Large", "C. surinamensis - Small", "C. surinamensis - Large"),common.legend = TRUE)
network.rare.large$plot_network(method = "networkD3", node_color = "module")
network.rare.large$plot_network(method = "networkD3", node_color = "Phylum")

# calculate network attributes
network.rare.large$cal_network_attr()
network.rare.large$res_network_attr
network.rare.small$res_network_attr

# get node properties
network.rare.large$get_node_table(node_roles = TRUE)
# return t1$res_node_table

# get edge properties
network.rare.large$get_edge_table()
# return t1$res_edge_table 
network.rare.large$get_adjacency_matrix()
# return t1$res_adjacency_matrix

# add_label = TRUE can be used to directly add text label for points
bryo_nifh_network$Rflaccida_Small$plot_taxa_roles(use_type = 1)
bryo_nifh_network$Rflaccida_Large$plot_taxa_roles(use_type = 1)

bryo_nifh_network$Csurinamensis_Small$plot_taxa_roles(use_type = 1)
bryo_nifh_network$Csurinamensis_Large$plot_taxa_roles(use_type = 1)
# plot node roles with phylum information
ggarrange(bryo_nifh_network$Csurinamensis_Small$plot_taxa_roles(use_type = 2, use_level = "Phylum")+theme_classic(),
          bryo_nifh_network$Csurinamensis_Large$plot_taxa_roles(use_type = 2, use_level = "Phylum")+theme_classic(),
          bryo_nifh_network$Rflaccida_Small$plot_taxa_roles(use_type = 2, use_level = "Phylum", pch=17)+theme_classic(),
          bryo_nifh_network$Rflaccida_Large$plot_taxa_roles(use_type = 2, use_level = "Phylum")+theme_classic(),
          nrow =2, ncol = 2, legend = "right", common.legend = TRUE,
          labels = c("C. surinamensis - Small", "C. surinamensis - Large", "R. flaccida - Small", "R. flaccida - Large"))

zipi = ggarrange(bryo_nifh_network$Csurinamensis_Small$plot_taxa_roles(use_type = 1, cex=3, alpha=0.7, roles_color_values = c("#00798c", "#D16103"))+theme_classic()+
                   theme(text = element_text(size = 16)),
                 bryo_nifh_network$Csurinamensis_Large$plot_taxa_roles(use_type = 1, cex=3, alpha=0.7, roles_color_values = c("#00798c", "#D16103"))+theme_classic()+
                   theme(text = element_text(size = 16)), 
                 bryo_nifh_network$Rflaccida_Small$plot_taxa_roles(use_type = 1, cex=3, alpha=0.7, roles_color_values = c("#D16103", "white"))+theme_classic()+
                   theme(text = element_text(size = 16)),
                 bryo_nifh_network$Rflaccida_Large$plot_taxa_roles(use_type = 1, cex=3, alpha=0.7, roles_color_values = c("#00798c", "#D16103"))+theme_classic()+
                   theme(text = element_text(size = 16)), 
                 nrow =2, ncol = 2, legend = "right", common.legend = TRUE,
                 labels = c("C. surinamensis - Small", "C. surinamensis - Large","R. flaccida - Small", "R. flaccida - Large"))
zipi
table(bryo_nifh_network$Csurinamensis_Large$res_node_table$taxa_roles)



bryo_nifh_network$Rflaccida_Small$cal_sum_links(taxa_level = "Genus")
dim(bryo_nifh_network$Rflaccida_Small$res_edge_table[1])#total number of interactions +
bryo_nifh_network$Rflaccida_Small$res_sum_links_pos = (bryo_nifh_network$Rflaccida_Small$res_sum_links_pos/3984)*100#transform to percentage

bryo_nifh_network$Rflaccida_Large$cal_sum_links(taxa_level = "Genus")
dim(bryo_nifh_network$Rflaccida_Large$res_edge_table[1])#total number of interactions +
bryo_nifh_network$Rflaccida_Large$res_sum_links_pos = (bryo_nifh_network$Rflaccida_Large$res_sum_links_pos/5681)*100#transform to percentage

bryo_nifh_network$Csurinamensis_Small$cal_sum_links(taxa_level = "Genus")
dim(bryo_nifh_network$Csurinamensis_Small$res_edge_table[1])#total number of interactions +
bryo_nifh_network$Csurinamensis_Small$res_sum_links_pos = (bryo_nifh_network$Csurinamensis_Small$res_sum_links_pos/5084)*100#transform to percentage

bryo_nifh_network$Csurinamensis_Large$cal_sum_links(taxa_level = "Genus")
dim(bryo_nifh_network$Csurinamensis_Large$res_edge_table[1])#total number of interactions +
bryo_nifh_network$Csurinamensis_Large$res_sum_links_pos = (bryo_nifh_network$Csurinamensis_Large$res_sum_links_pos/6405)*100
#devtools::install_github("mattflor/chorddiag")
library(chorddiag)
# interactive visualization; require chorddiag package; see https://github.com/mattflor/chorddiag
bryo_nifh_network$Rflaccida_Small$plot_sum_links(plot_pos = TRUE, plot_num = 10, color_values = RColorBrewer::brewer.pal(10, "Paired"))
bryo_nifh_network$Rflaccida_Large$plot_sum_links(plot_pos = TRUE, plot_num = 10, color_values = RColorBrewer::brewer.pal(10, "Paired"))
bryo_nifh_network$Csurinamensis_Small$plot_sum_links(plot_pos = TRUE, plot_num = 10, color_values = RColorBrewer::brewer.pal(10, "Paired"))
bryo_nifh_network$Csurinamensis_Large$plot_sum_links(plot_pos = TRUE, plot_num = 10, color_values = RColorBrewer::brewer.pal(10, "Paired"))
# From v1.2.0, method = "circlize" is available for conveniently saving the static plot
# If circlize package is not installed, first run: install.packages("circlize")
bryo_nifh_network$Rflaccida_Small$plot_sum_links(method = "circlize", transparency = 0.2, annotationTrackHeight = circlize::mm_h(c(5, 5)))
bryo_nifh_network$Rflaccida_Large$plot_sum_links(method = "circlize", transparency = 0.2, annotationTrackHeight = circlize::mm_h(c(5, 5)))

bryo_nifh_network$Csurinamensis_Small$plot_sum_links(method = "circlize", transparency = 0.2, annotationTrackHeight = circlize::mm_h(c(5, 5)))
bryo_nifh_network$Csurinamensis_Large$plot_sum_links(method = "circlize", transparency = 0.2, annotationTrackHeight = circlize::mm_h(c(5, 5)))

library(circlize)
chordDiagram(bryo_nifh_network$Csurinamensis_Large$res_sum_links_pos)+ theme(axis.text.y = element_text(angle = 180))




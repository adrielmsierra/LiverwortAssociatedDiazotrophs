###################################################################################
#                           Bryophyte NifH (Diazotrophs)
#                               Diversity metrics
#                               Betadiversity
#                        Community composition (Ordination)
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
# Host-specific diazotrophic microbiota
nifH.core.ps.ex <- readRDS("Bryophyte_nifH_coremicrobiota_Jun24.rds")

# Calculate relative abundance
nifH.ps.relative.rare <- transform_sample_counts(nifh.dataset.rarefied, function(OTU) OTU/sum(OTU))
nifH.ps.relative.core <- transform_sample_counts(nifH.core.ps.ex, function(OTU) OTU/sum(OTU))


# Calculate distance matrix
# Bray-Curtis Non phylogenetic and Weigthed UNifrac for phylogetically


#Ordination
######
ordination.PCoA.bray.rare <- ordinate(nifh.dataset.rarefied, method="PCoA", distance="bray", na.rm=TRUE)
ordination.PCoA.wunifrac.rare <- ordinate(nifh.dataset.rarefied, method="PCoA", distance="unifrac", weighted=TRUE)
#core
ordination.PCoA.bray.core <- ordinate(nifH.ps.relative.core , method="PCoA", distance="bray", na.rm=TRUE)
ordination.PCoA.wunifrac.core <- ordinate(nifH.ps.relative.core , method="PCoA", distance="unifrac", weighted=TRUE)

#extract distance matrix
# Bray-Curtis and WUnifrac

set.seed(1911)
#Bray-Curtis

distance.matrix.bray.rare <- phyloseq::distance(nifH.ps.relative.rare, method = "bray")
distance.matrix.bray.core <- phyloseq::distance(nifH.ps.relative.core, method = "bray")
#WUnifrac
distance.matrix.wunifrac.rare <- UniFrac(nifH.ps.relative.rare, weighted = TRUE)
distance.matrix.wunifrac.core <- UniFrac(nifH.ps.relative.core, weighted = TRUE)
wunifrac_dist <- UniFrac(nifH.ps.relative.rare, weighted = TRUE)

sample_data(nifH.ps.relative.rare)$Group = as.factor(paste(sampledf$Reserve, sampledf$Size))

dt.size = data.frame(sampledf.pops[1:76,], ordination.PCoA.wunifrac.rare)
str(dt.size)

library(dplyr)

dt.sizepop = dt.size[1:76,] %>% 
  group_by(sampledf.pops) %>%
  summarise_all("mean")

View(dt.size[1:76,])#radula
View(dt.size[77:146,])#cololejeunea
head(dt.sizepop)
dt.sizepop = data.frame(dt.sizepop)

rownames(dt.sizepop) = dt.sizepop$sampledf.pops
dist.bray = vegdist(dt.sizepop[-1], dist="bray")
write.csv(dist.bray, "Braycurtis_rarefied_dist.csv")
sampledf <- data.frame(sample_data(nifH.ps.relative.rare))
sampledf.core <- data.frame(sample_data(nifH.ps.relative.core))


# Define your groups
groups <- sample_data(nifH.ps.relative.rare)$Group  # Replace 'Group' with your grouping variable

# Function to calculate the average distance within each group
average_within_group <- function(physeq, group) {
  # Subset phyloseq object by group
  physeq_group <- prune_samples(sample_data(physeq)$Group == group, physeq)
  
  # Calculate weighted UniFrac distance for the group
  dist_group <- UniFrac(physeq_group, weighted=TRUE)
  
  # Return the mean distance
  return(mean(as.matrix(dist_group)))
}

# Apply the function to each group
group_distances <- sapply(unique(groups), function(g) average_within_group(nifH.ps.relative.rare, g))

# Print results
group_distances

nifH.ps.relative.rare.subs = subset_samples(nifH.ps.relative.rare, Species=="Radula flaccida")
nifH.ps.relative.rare.subs = subset_samples(nifH.ps.relative.rare, Species=="Cololejeunea surinamensis")
# Define your groups
groups <- sample_data(nifH.ps.relative.rare.subs)$Group  # Replace 'Group' with your grouping variable
# Create a matrix to store the average distances between each pair of groups
group_combinations <- combn(unique(groups), 2)
between_group_distances <- matrix(NA, nrow=ncol(group_combinations), ncol=1)

# Function to calculate the average distance between two groups
average_between_groups <- function(group1, group2, dist_matrix, groups) {
  # Get sample names for each group
  samples_group1 <- rownames(sample_data(nifH.ps.relative.rare.subs))[groups == group1]
  samples_group2 <- rownames(sample_data(nifH.ps.relative.rare.subs))[groups == group2]
  
  # Subset the distance matrix for these samples
  dist_subset <- dist_matrix[samples_group1, samples_group2]
  
  # Return the mean distance between the two groups
  return(mean(dist_subset))
}

# Apply the function to each pair of groups
for (i in 1:ncol(group_combinations)) {
  group1 <- group_combinations[1, i]
  group2 <- group_combinations[2, i]
  between_group_distances[i] <- average_between_groups(group1, group2, as.matrix(wunifrac_dist), groups)
}

# Convert to a more interpretable format
between_group_df <- data.frame(Group1 = group_combinations[1, ],
                               Group2 = group_combinations[2, ],
                               Average_Distance = between_group_distances)

# Print results
View(between_group_df)
#l'analyse multivariée permutationnelle de la variance (PERMANOVA) 
#colnames(sampledf)=c("Voucher", "Species", "Site", "Size", "Host")

sampledf$Voucher = as.factor(sampledf$Voucher)
sampledf$Host = as.factor(sampledf$Host)
sampledf$Class.size = as.factor(sampledf$Class.size)
str(sampledf)

sampledf.core$Voucher = as.factor(sampledf.core$Voucher)
sampledf.core$Host = as.factor(sampledf.core$Host)
sampledf.core$Class.size = as.factor(sampledf.core$Class.size)
str(sampledf.core)

# Permutational Multivariate Analysis of Variance Using Distance Matrices

library("GUniFrac")
# bray-curtis
mad.bact.bray.rare = adonis3(distance.matrix.bray.rare~sampledf$Host+sampledf$Size, method ="bray", perm=999)
mad.bact.bray.rare$aov.tab

mad.bact.bray.core = adonis3(distance.matrix.bray.core~sampledf.core$Host+sampledf.core$Size, method ="bray", perm=999)
mad.bact.bray.core$aov.tab
#permanova_pairwise(distance.matrix.bray.rare, sampledf$Size,permutations = 999, method = "bray",padj = "bonferroni")
list(unique(colnames(sampledf)))
pairwise.adonis(distance.matrix.bray.rare, sampledf$Size, sim.function = "vegdist", 
                sim.method = "bray", p.adjust.m = "bonferroni")
pairwise.adonis(distance.matrix.bray.core, sampledf.core$Size, sim.function = "vegdist", 
                sim.method = "bray", p.adjust.m = "bonferroni")


dt.sizepop = distance.matrix.bray.rare[77:146,] %>% 
  group_by(sampledf.pops) %>%
  summarise_all("sum")

# Weighted Unifrac
mad.bact.wu.rare = adonis3(distance.matrix.wunifrac.rare~sampledf$Host+sampledf$Size, perm=999)
mad.bact.wu.rare = adonis3(distance.matrix.wunifrac.rare~sampledf$sampledf.pops, perm=999)
mad.bact.wu.rare$aov.tab

mad.bact.wu.core = adonis3(distance.matrix.wunifrac.core~sampledf.core$Host+sampledf.core$Size, perm=999)
mad.bact.wu.core$aov.tab

pairwise.adonis(distance.matrix.wunifrac.rare, sampledf$Size, p.adjust.m = "bonferroni")
pairwise.adonis(distance.matrix.wunifrac.core, sampledf.core$Size, sim.function = "vegdist", 
                sim.method = "bray", p.adjust.m = "bonferroni")


# Multivariate homogeneity of groups dispersions (variances)
# complete data
# Multivariate homogeneity of groups dispersions (variances)

# rarefied data
# Multivariate homogeneity of groups dispersions (variances)
#betapcoa.bray.rare.h<-betadisper(distance.matrix.bray.rare, sampledf$Host)
#anova(betapcoa.bray.rare.h)
betapcoa.bray.rare.s<-betadisper(distance.matrix.bray.rare, sampledf$Size, bias.adjust = TRUE)
anova(betapcoa.bray.rare.s)
write.csv(betapcoa.bray.rare.s$vectors, "Pcoa_vectors_rarefied_braycurtis.csv")
#betapcoa.bray.rare.cs<-betadisper(distance.matrix.bray.rare, sampledf$Class.size)
#anova(betapcoa.bray.rare.cs)
#betapcoa.wu.rare.h<-betadisper(distance.matrix.wunifrac.rare, sampledf$Host)
#anova(betapcoa.wu.rare.h)
betapcoa.wu.rare.s<-betadisper(distance.matrix.wunifrac.rare, sampledf$Size, bias.adjust = TRUE)
write.csv(betapcoa.wu.rare.s$vectors, "Pcoa_vectors_rarefied_wu.csv")

anova(betapcoa.wu.rare.s)
#betapcoa.wu.rare.cs<-betadisper(distance.matrix.wunifrac.rare, sampledf$Class.size)
#anova(betapcoa.wu.rare.cs)

# core data
#betapcoa.bray.core.h<-betadisper(distance.matrix.bray.core, sampledf.core$Host)
#anova(betapcoa.bray.core.h)
betapcoa.bray.core.s<-betadisper(distance.matrix.bray.core, sampledf.core$Size, bias.adjust = TRUE)
#anova(betapcoa.bray.core.s)
#betapcoa.bray.core.cs<-betadisper(distance.matrix.bray.core, sampledf.core$Class.size)
#anova(betapcoa.bray.core.cs)
#betapcoa.wu.core.h<-betadisper(distance.matrix.wunifrac.core, sampledf.core$Host)
#anova(betapcoa.wu.core.h)
betapcoa.wu.core.s<-betadisper(distance.matrix.wunifrac.core, sampledf.core$Size, bias.adjust = TRUE)
#anova(betapcoa.wu.core.s)
#betapcoa.wu.core.cs<-betadisper(distance.matrix.wunifrac.core, sampledf.core$Class.size)
#anova(betapcoa.wu.core.cs)
#
betaperm = permutest(betapcoa.bray.core.s, pairwise = TRUE)
(mod.HSD.betaperm <- TukeyHSD(betapcoa.bray.core.s))



mad.bact.rc = adonis2(distance.matrix.bray.rare~sampledf$Size, method ="raup", perm=999)
mad.bact.rc


# Extract distances

fragdist.bray.rare = data.frame(betapcoa.bray.rare.s$distances, sampledf[,c(2:7)])
colnames(fragdist.bray.rare) = c("Distance", "Voucher", "Species", "Reserve", "Class.size", "Size", "Pop")

fragdist.wunifrac.rare = data.frame(betapcoa.wu.rare.s$distances, sampledf[,c(2:7)])
colnames(fragdist.wunifrac.rare) = c("Distance", "Voucher", "Species", "Reserve", "Class.size", "Size", "Pop")

fragdist.bray.core = data.frame(betapcoa.bray.core.s$distances, sampledf.core[,c(2:7)])
colnames(fragdist.bray.core) = c("Distance", "Voucher", "Species", "Reserve", "Class.size", "Size", "Pop")

fragdist.wunifrac.core = data.frame(betapcoa.wu.core.s$distances, sampledf.core[,c(2:7)])
colnames(fragdist.wunifrac.core) = c("Distance", "Voucher", "Species", "Reserve", "Class.size", "Size", "Pop")

#

teste.sta.b.rare = aov(Distance~Size*Species, data=fragdist.bray.rare)
summary(teste.sta.b.rare)

teste.sta.wu.rare = aov(Distance~Size*Species, data=fragdist.wunifrac.rare)
summary(teste.sta.wu.rare)

teste.sta.b.core = aov(Distance~Size*Species, data=fragdist.bray.core)
summary(teste.sta.b.core)

teste.sta.wu.core = aov(Distance~Size*Species, data=fragdist.wunifrac.core)
summary(teste.sta.wu.core)




#rarefied
# matrices
#distance matrices with Bray-Curtis and Weighted Unifrac
rare.pcoa<-pcoa(distance.matrix.bray.rare)
rare.pcoa<-pcoa(distance.matrix.wunifrac.rare)
#
rare.pcoa.vectors<-as.data.frame(rare.pcoa$vectors)
rare.pcoa.vectors$Size<-as.factor(sampledf$Size)
rare.pcoa.vectors$host<-as.factor(sampledf$Host)
rare.pcoa.vectors$class.size<-as.factor(sampledf$Class.size)
# percentage explain variability
explainvar1rel.rare <- rare.pcoa$values[1,1]/sum(rare.pcoa$values[1]) * 100
explainvar2rel.rare <- rare.pcoa$values[2,1]/sum(rare.pcoa$values[1]) * 100

head(rare.pcoa.vectors)
xlabr.rare = paste("PCoA Axis 1 (", format(round(explainvar1rel.rare, 1), nsmall=1)," %)", sep="")
ylabr.rare = paste("PCoA Axis 2 (", format(round(explainvar2rel.rare, 1), nsmall=1)," %)", sep="")

# Plot ordination
p.bray.pcoa.rare = ggplot(data = rare.pcoa.vectors, aes(x = Axis.1, y = Axis.2, fill = host))+#braycurtis
  ggtitle("Rarefied", subtitle = "Bray-Curtis")+
  #p.wu.pcoa.rare = ggplot(data = rare.pcoa.vectors, aes(x = Axis.1, y = Axis.2, fill = host))+# weigthed unifrac 
  #ggtitle("", subtitle = "Weighted Unifrac")+
  # points
  geom_point(aes(color = class.size, shape= host), size = 3, stat="identity", pos="identity", show.legend = TRUE)+
  # stat ellipse
  stat_ellipse(aes(color=class.size, linetype=host), level=0.95, size=0.75, show.legend = TRUE, alpha=0.25)+
  scale_color_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")+
  labs(x = xlabr.rare, y = ylabr.rare)+
  theme_bw() +
  theme(text = element_text(size = 14))+
  theme(panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())# custom theme

#final separate plots
p.bray.pcoa.rare
p.wu.pcoa.rare
#supp figure ordination aug24
supp.fig.ord = ggarrange(p.bray.pcoa, p.wu.pcoa, p.bray.pcoa.rare,
                         p.wu.pcoa.rare, ncol = 2, nrow = 2, align = c("hv"), 
                         common.legend = TRUE, legend="right")
###################################
#core microbiota
# matrices
#core.pcoa<-pcoa(distance.matrix.bray.core)
core.pcoa<-pcoa(distance.matrix.wunifrac.core)
head(core.pcoa)
core.pcoa.vectors<-as.data.frame(core.pcoa$vectors)
core.pcoa.vectors$Size<-as.factor(sampledf.core$Size)
core.pcoa.vectors$host<-as.factor(sampledf.core$Host)
core.pcoa.vectors$class.size<-as.factor(sampledf.core$Class.size)
# perct explain variability
explainvar1rel.core <- core.pcoa$values[1,1]/sum(core.pcoa$values[1]) * 100
explainvar2rel.core <- core.pcoa$values[2,1]/sum(core.pcoa$values[1]) * 100

head(core.pcoa.vectors)
xlabr.core = paste("PCoA Axis 1 (", format(round(explainvar1rel.core, 1), nsmall=1)," %)", sep="")
ylabr.core = paste("PCoA Axis 2 (", format(round(explainvar2rel.core, 1), nsmall=1)," %)", sep="")


#p.bray.pcoa.core = ggplot(data = core.pcoa.vectors, aes(x = Axis.1, y = Axis.2, fill = host))+#braycurtis
#ggtitle("Core", subtitle = "Bray-Curtis")+
p.uw.core.pcoa = ggplot(data = core.pcoa.vectors, aes(x = Axis.1, y = Axis.2, fill = host))+# weigthed unifrac 
  ggtitle("", subtitle = "Weighted Unifrac")+
  # points
  geom_point(aes(color = class.size, shape= host), size = 3, stat="identity", pos="identity", show.legend = TRUE)+
  # stat ellipse
  stat_ellipse(aes(color=class.size, linetype=class.size), level=0.95, size=0.75, show.legend = TRUE, alpha=0.25)+
  scale_color_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")+
  labs(x = xlabr.core, y = ylabr.core)+
  theme_bw() +
  theme(text = element_text(size = 16))+
  theme(panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())# custom theme
#Final plots separate
p.bray.pcoa.core
p.uw.core.pcoa
#

#Compiling figure with all PCoAs
FigureS_pcoas = ggarrange(p.bray.pcoa, p.wu.pcoa, p.bray.pcoa.rare, 
                          p.wu.pcoa.rare, p.bray.pcoa.core, p.uw.core.pcoa, 
                          ncol = 2, nrow = 3, align = c("hv"), 
                          common.legend = TRUE, legend="top")
FigureS_pcoas
#######
#distance to centroid regression

library(introdataviz)

compare_means(Distance ~ Size, data = fragdist.bray.rare[1:76,])
levels(fragdist.bray.rare$Size)=as.factor(c("1 ha","10 ha", "100 ha","Continuous"))
my_compar <- list(c("1 ha","10 ha"), c("1 ha","100 ha"), c("1 ha","Continuous"), c("10 ha","100 ha"), c("10 ha","Continuous"), c("100 ha","Continuous"))
pdiazo.bray.rare.rad = ggplot(data=fragdist.bray.rare[1:76,], aes(x=Size, y=Distance,  group=Size)) + 
  #geom_point(aes(shape=Species), position=pd, size=3, color="darkgrey", fill="black")+ labs(x="Forest size", y = "Distance to centroid")+
  #geom_smooth(aes(linetype= Species), method=lm, color="black", fill="lightgrey", se=TRUE)+
  ggtitle("R. flaccida", subtitle = "Bray-Curtis")+
  #geom_split_violin(trim=FALSE)+
  geom_violin(trim=FALSE, fill="#e6ab02", color="black")+
  #scale_color_brewer(palette="Dark2")+
  #scale_fill_brewer(palette="Dark2")+
  geom_boxplot(width = 0.04, )+
  #stat_cor(label.y = c(0.8,0.82), size = 8)+ 
  #stat_regline_equation(label.y = c(0.85,0.87), size = 8)+
  theme_classic()+
  theme(text = element_text(size = 14))+
  #facet_grid(. ~ variable, space = "free_x", scales = "free_x") +
  stat_compare_means(comparisons = my_compar,  label = "p.signif", 
                     size = 6, label.y = c(0.81, 0.79, 0.77)) # Increased label size
#theme(legend.justification=c(),
#      legend.position=c())
pdiazo.bray.rare.rad

compare_means(Distance ~ Size, data = fragdist.bray.rare[77:146,])

pdiazo.bray.rare.col = ggplot(data=fragdist.bray.rare[77:146,], aes(x=Size, y=Distance, group=Size)) + 
  #geom_point(aes(shape=Species), position=pd, size=3, color="darkgrey", fill="black")+ labs(x="Forest size", y = "Distance to centroid")+
  #geom_smooth(aes(linetype= Species), method=lm, color="black", fill="lightgrey", se=TRUE)+
  ggtitle("C. surinamensis", subtitle = "Bray-Curtis")+
  geom_violin(trim=FALSE, fill="#e6ab02", color="black")+
  scale_color_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")+
  geom_boxplot(width = 0.04)+
  #stat_cor(label.y = c(0.8,0.82), size = 8)+ 
  #stat_regline_equation(label.y = c(0.85,0.87), size = 8)+
  theme_classic()+
  theme(text = element_text(size = 14))+
  stat_compare_means(comparisons = my_compar,  label = "p.signif", 
                     size = 6, label.y = c(0.81, 0.79, 0.77, 0.75, 0.73, 0.71)) # Increased label size
#theme(legend.justification=c(),
#legend.position=c())
pdiazo.bray.rare.col
levels(fragdist.wunifrac.rare$Size)=as.factor(c("1 ha","10 ha", "100 ha","Continuous"))

compare_means(Distance ~ Size, data = fragdist.wunifrac.rare[1:76,])

pdiazo.wu.rare.rad = ggplot(data=fragdist.wunifrac.rare[1:76,], aes(x=Size, y=Distance, group=Size)) + 
  ggtitle("R. flaccida", subtitle = "Weighted Unifrac")+
  #geom_point(aes(shape=Species), position=pd, size=3, color="darkgrey", fill="black")+ labs(x="Forest size", y = "Distance to centroid")+
  #geom_smooth(aes(linetype= Species), method=lm, color="black", fill="lightgrey", se=TRUE)+
  geom_violin(trim=FALSE, fill="#e6ab02", color="black")+
  #scale_color_brewer(palette="Dark2")+
  #scale_fill_brewer(palette="Dark2")+
  geom_boxplot(width = 0.04)+
  #geom_boxplot(aes(fill=Host), width = 0.8)+
  #stat_cor(label.y = c(0.8,0.82), size = 8)+ 
  #stat_regline_equation(label.y = c(0.85,0.87), size = 8)+
  theme_classic()+
  theme(text = element_text(size = 14))+
  stat_compare_means(comparisons = my_compar,  label = "p.signif", 
                     size = 6, label.y = c(0.94, 0.9, 0.86, 0.82)) # Increased label size
#theme(legend.justification=c(),
#     legend.position=c())
pdiazo.wu.rare.rad
compare_means(Distance ~ Size, data = fragdist.wunifrac.rare[77:146,])
pdiazo.wu.rare.col = ggplot(data=fragdist.wunifrac.rare[77:146,], aes(x=Size, y=Distance, group=Size)) + 
  ggtitle("C. surinamensis", subtitle = "Weighted Unifrac")+
  #geom_point(aes(shape=Species), position=pd, size=3, color="darkgrey", fill="black")+ labs(x="Forest size", y = "Distance to centroid")+
  #geom_smooth(aes(linetype= Species), method=lm, color="black", fill="lightgrey", se=TRUE)+
  geom_violin(trim=FALSE, fill="#e6ab02", color="black")+
  #scale_color_brewer(palette="Dark2")+
  #scale_fill_brewer(palette="Dark2")+
  geom_boxplot(width = 0.04)+
  #geom_boxplot(aes(fill=Host), width = 0.8)+
  #stat_cor(label.y = c(0.8,0.82), size = 8)+ 
  #stat_regline_equation(label.y = c(0.85,0.87), size = 8)+
  theme_classic()+
  theme(text = element_text(size = 14))
#stat_compare_means(comparisons = my_compar,  label = "p.signif", 
#                   size = 6, label.y = c(0.94, 0.9, 0.86, 0.82)) # Increased label size 
#theme(legend.justification=c(),
#     legend.position=c())
pdiazo.wu.rare.col



pdiazo.bray.core = ggplot(data=fragdist.bray.core, aes(x=Size, y=Distance, fill=Species, linetype= Species, group=Size)) + 
  #geom_point(aes(shape=Species), position=pd, size=3, color="darkgrey", fill="black")+ labs(x="Forest size", y = "Distance to centroid")+
  #geom_smooth(aes(linetype= Species), method=lm, color="black", fill="lightgrey", se=TRUE)+
  geom_violin(trim=FALSE, fill="#e6ab02", color="black")+
  scale_color_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")+
  geom_boxplot(width = 0.07)+
  stat_cor(label.y = c(0.8,0.82), size = 8)+ 
  stat_regline_equation(label.y = c(0.85,0.87), size = 8)+
  theme_classic()+
  theme(text = element_text(size = 14))+
  theme(legend.justification=c(),
        legend.position=c())
pdiazo.bray.core

pdiazo.uw.core = ggplot(data=fragdist.wunifrac.core, aes(x=Size, y=Distance, fill=Species, linetype= Species, group=Size)) + 
  #geom_point(aes(shape=Species), position=pd, size=3, color="darkgrey", fill="black")+ labs(x="Forest size", y = "Distance to centroid")+
  #geom_smooth(aes(linetype= Species), method=lm, color="black", fill="lightgrey", se=TRUE)+
  geom_violin(trim=FALSE, fill="#e6ab02", color="black")+
  scale_color_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")+
  geom_boxplot(width = 0.07)+
  stat_cor(label.y = c(0.8,0.82), size = 8)+ 
  stat_regline_equation(label.y = c(0.85,0.87), size = 8)+
  theme_classic()+
  theme(text = element_text(size = 14))+
  theme(legend.justification=c(),
        legend.position=c())
pdiazo.uw.core

Figure2_pcoas = ggarrange(ggarrange(p.bray.pcoa.,p.bray.pcoa., 
                                    p.bray.pcoa.rare, p.wu.pcoa.rare,  
                                    ncol = 1, nrow = 2, align = c("hv"), 
                                    common.legend = TRUE, legend="right"), ggarrange(pdiazo.bray.rare,pdiazo.wu.rare,ncol = 1, nrow = 2), ncol=2)
Figure2_pcoas

#Core microbiome
pdiazo.bray.core <- ggplot(data=fragdist.bray.core, aes(x=betapcoa.bray.core.group, y=betapcoa.bray.core.distances)) + 
  geom_violin(trim=FALSE, fill="#e6ab02", color="black")+
  scale_color_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")+
  geom_boxplot(width = 0.07)+
  #geom_boxplot(aes(fill=Host), width = 0.8)+
  theme_classic()
#pdiazo.bray = pdiazo.bray + labs(x="Forest size", y = "Distance to centroid")
pdiazo.bray.core = pdiazo.bray.core + labs(x="Taille de la forêt", y = "Distance au centroïde")
pdiazo.bray.core + theme(text = element_text(size = 18)) 

pdiazo.unifrac.core <- ggplot(data=fragdist.unifrac.core, aes(x=betapcoa.unifrac.core.group, y=betapcoa.unifrac.core.distances)) + 
  geom_violin(trim=FALSE, fill="#e6ab02", color="black")+
  scale_color_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")+
  geom_boxplot(width = 0.07)+
  #geom_boxplot(aes(fill=Host), width = 0.8)+
  theme_classic()
#pdiazo.unifrac = pdiazo.unifrac + labs(x="Forest size", y = "Distance to centroid")
pdiazo.unifrac.core = pdiazo.unifrac.core + labs(x="Taille de la forêt", y = "Distance au centroïde")
pdiazo.unifrac.core + theme(text = element_text(size = 18)) 

fig.dist.pub = ggarrange(pdiazo.bray, pdiazo.unifrac,
                         pdiazo.bray.core, pdiazo.unifrac.core,
                         ncol = 2, nrow = 2, align = c("hv"), common.legend = TRUE, legend="right")
fig.dist.pub

pdiazo.bray.rare.rad

fig.dist.pub.aug24 = ggarrange(pdiazo.bray.rare.rad, pdiazo.bray.rare.col,
                               pdiazo.wu.rare.rad, pdiazo.wu.rare.col,
                               ncol = 1, nrow = 4, align = c("hv"), common.legend = TRUE, legend="right")
fig.dist.pub.aug24

ggarrange(ggarrange(p.bray.pcoa.rare,
                    p.wu.pcoa.rare, ncol = 1, nrow = 2, align = c("hv"), 
                    common.legend = TRUE, legend="right"), 
          fig.dist.pub.aug24,
          ncol = 2, nrow = 1, align = c("hv"), common.legend = TRUE, legend="right")
ggarrange(ggarrange(p.bray.pcoa.rare,
                    p.wu.pcoa.rare, ncol = 1, nrow = 2, align = c("hv"), 
                    common.legend = TRUE, legend="top"), 
          ggarrange(pdiazo.bray.rare.rad, pdiazo.bray.rare.col,
                    pdiazo.wu.rare.rad, pdiazo.wu.rare.col,
                    ncol = 2, nrow = 2, align = c("hv")), ncol = 2, nrow = 1)

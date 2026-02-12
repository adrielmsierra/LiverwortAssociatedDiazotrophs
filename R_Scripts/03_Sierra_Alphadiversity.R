###################################################################################
#                           Bryophyte NifH (Diazotrophs)
#                               Diversity metrics
#                               Alpha diversity
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

#### sequencing depth diversity bias to consider rarefying
####
####
#Alpha data estimated for the whole dataset
###alpha diversity
#observe species hill numbers
# complete microbiote community
hillq0 <- estimate_richness(nifH.ps.com.dataset, measures = "Observed")

nifH.ps.com.dataset.norm <- microbiome::transform(nifH.ps.com.dataset, transform = "compositional")

shannon.hill <- exp(vegan::diversity(otu_table(nifH.ps.com.dataset.norm), index = "shannon"))
shannon.hill.df <- as.data.frame(shannon.hill)

simpson.hill <- 1/(1-(vegan::diversity(otu_table(nifH.ps.com.dataset.norm), index = "simpson")))
simpson.hill.df <- as.data.frame(simpson.hill)

pielou.hill = vegan::diversity(otu_table(nifH.ps.com.dataset.norm), index = "shannon")/log(specnumber(otu_table(nifH.ps.com.dataset.norm)))
pielou.hill.df = as.data.frame(pielou.hill)

calcDiv(otu_table(nifH.ps.com.dataset.norm))

library(picante)
faith_phydiv = pd(otu_table(nifH.ps.com.dataset.norm), phy_tree(nifH.ps.com.dataset.norm), include.root=FALSE)
faith_phydiv.df = as.data.frame(faith_phydiv)
alpha_data<- cbind(hillq0, shannon.hill.df, simpson.hill.df, pielou.hill.df, faith_phydiv.df)

##### sequencing depth diversity bias correlation
summary(lm(alpha_data$Observed~sample_sums(nifH.ps.com.dataset)))
seqdepth.div = data.frame(alpha_data$Observed, sample_sums(nifH.ps.com.dataset))
#plot
ggplot(seqdepth.div, aes(y=alpha_data.Observed, x=sample_sums.nifH.ps.com.dataset.))+ 
  geom_point(size=2, stat='identity')+ labs(x="Sequencing depth", y ="Observed number of species")+
  geom_smooth(aes(x=sample_sums.nifH.ps.com.dataset., y=alpha_data.Observed), method= "lm", color="Black", fill="lightgrey", se=TRUE)+
  theme_classic()+
  theme(axis.text.x=element_text(size = 16), axis.text.y=element_text(size = 16))+
  stat_cor(label.y = 210, size = 8)+ 
  stat_regline_equation(label.y = 190, size = 8)+
  theme(text = element_text(size = 14)) 

# rarefied to the sequencing depth of 90% of the minimum sample depth in the dataset this is 513.9 reads
nifh.dataset.rarefied = rarefy_even_depth(nifH.ps.com.dataset, rngseed=1223, sample.size=0.9*min(sample_sums(nifH.ps.com.dataset)), replace=F)
saveRDS(nifh.dataset.rarefied, "./02_Phyloseq_data_analysis/nifH_ps_rarefied_dataset_tree_jun24.rds")
nifh.dataset.rarefied = readRDS("./02_Phyloseq_data_analysis/nifH_ps_rarefied_dataset_tree_jun24.rds")

#
# Diversity estimates with hill numbers
hillq0.rare <- estimate_richness(nifh.dataset.rarefied, measures = "Observed")

nifH.ps.rare.dataset.norm <- microbiome::transform(nifh.dataset.rarefied, transform = "compositional")

shannon.rare.hill <- exp(vegan::diversity(otu_table(nifH.ps.rare.dataset.norm), index = "shannon"))
shannon.rare.hill.df <- as.data.frame(shannon.rare.hill)

simpson.rare.hill <- 1/(1-(vegan::diversity(otu_table(nifH.ps.rare.dataset.norm), index = "simpson")))
simpson.rare.hill.df <- as.data.frame(simpson.rare.hill)

pielou.rare.hill <- vegan::diversity(otu_table(nifH.ps.rare.dataset.norm), index = "shannon")/log(specnumber(otu_table(nifH.ps.rare.dataset.norm)))
pielou.rare.hill.df = as.data.frame(pielou.rare.hill)

faith_phydiv.rare = pd(otu_table(nifH.ps.rare.dataset.norm), phy_tree(nifH.ps.rare.dataset.norm), include.root=FALSE)
faith_phydiv.rare.df = as.data.frame(faith_phydiv.rare)

faith_phydiv.rare.df
alpha.rare<- cbind(hillq0.rare, shannon.rare.hill.df, simpson.rare.hill.df, pielou.rare.hill.df, faith_phydiv.rare.df, alpha_data[,8:14])
write.csv(alpha.rare, "./Tables_pub/alpha_rarefiedDiversity_Jun24.csv")
alpha.rare = read.csv("./Tables_pub/alpha_rarefiedDiversity_Jun24.csv")

#Shapiro test
head(alpha.rare)
str(alpha.rare)
test_names <- c( "Observed", "shannon.rare.hill", "simpson.rare.hill", "pielou.rare.hill", "PD")
shapiro.test(alpha.rare$Observed)
shapiro.test(alpha.rare$shannon.rare.hill)
shapiro.test(alpha.rare$simpson.rare.hill)
shapiro.test(alpha.rare$PD)
shapiro.test(alpha.rare$pielou.rare.hill)
library(ggpubr)
ggqqplot(alpha.rare$Observed)
for (test in test_names) {
  test_formula <- paste(test,"~","Host")
  print(paste("Mesure d'alpha diversité:", test))
  print(summary(glm(as.formula(test_formula), data=alpha.rare, family = poisson)))
}
#Radula
for (test in test_names) {
  test_formula <- paste(test,"~","Size")
  print(paste("Mesure d'alpha diversité:", test))
  print(summary(glm(as.formula(test_formula), data=alpha.rare[c(1:77),])))
}

summary(glm(Observed~Size,data= alpha.rare[c(1:77),]))
#Cololejeunea
for (test in test_names) {
  test_formula <- paste(test,"~","Size")
  print(paste("Mesure d'alpha diversité:", test))
  print(summary(glm(as.formula(test_formula), data=alpha.rare[c(77:146),])))
}

pd
Fig.div.rare.ms = 
  ggplot(alpha.rare, aes(x=Size, y=Observed, fill=Size, color=Size, linetype= Host, group=Host)) + 
  #ggplot(alpha_data, aes(x=Size, y=Observed)) + 
  ggtitle("Alpha Diversity", subtitle="Rarefied diazotrophic microbiota")+
  geom_point(position=pd, size=5)+ labs(x=" ", y ="Number of species (q = 0)")+
  geom_smooth(method=lm, color="black", fill="lightgrey", se=TRUE)+
  scale_fill_brewer(palette="Dark2")+
  scale_color_brewer(palette="Dark2")+
  theme_classic()+
  theme(text = element_text(size = 14), axis.text.x=element_text(angle=20, hjust=1))+
  theme(legend.justification=c(),
        legend.position=c()) 
Fig.div.rare.ms

Fig.sha.rare.ms = 
  ggplot(alpha.rare, aes(x=Size, y=shannon.rare.hill, fill=Size, color=Size, linetype= Host, group=Host)) + 
  #ggplot(alpha_data, aes(x=Size, y=shannon.hill)) + 
  geom_point(position=pd, size=5)+ labs(x=" ", y ="Shannon entropy (q = 1)")+
  geom_smooth(method=lm, color="black", fill="lightgrey", se=TRUE)+
  scale_fill_brewer(palette="Dark2")+
  scale_color_brewer(palette="Dark2")+
  theme_classic()+
  theme(text = element_text(size = 14), axis.text.x=element_text(angle=20, hjust=1))+
  theme(legend.justification=c(),
        legend.position=c()) 
Fig.sha.rare.ms

Fig.sim.rare.ms = ggplot(alpha.rare, aes(x=Size, y=simpson.rare.hill, fill=Size, color=Size, linetype= Host, group=Host)) + 
  geom_point(position=pd, size=5)+ labs(x="Forest fragment size", y ="inverse Simpson index (q = 2)")+
  geom_smooth(method=lm, color="black", fill="lightgrey", se=TRUE)+
  scale_fill_brewer(palette="Dark2")+
  scale_color_brewer(palette="Dark2")+
  theme_classic()+
  theme(text = element_text(size = 14), axis.text.x=element_text(angle=20, hjust=1))+
  theme(legend.justification=c(),
        legend.position=c())
Fig.sim.rare.ms

Fig.pie.rare.ms = ggplot(alpha.rare, aes(x=Size, y=pielou.rare.hill, fill=Size, color=Size, linetype= Host, group=Host)) + 
  geom_point(position=pd, size=5)+ labs(x=" ", y ="Pielou Index")+
  geom_smooth(method=lm, color="black", fill="lightgrey", se=TRUE)+
  scale_fill_brewer(palette="Dark2")+
  scale_color_brewer(palette="Dark2")+
  theme_classic()+
  theme(text = element_text(size = 14), axis.text.x=element_text(angle=20, hjust=1))+
  theme(legend.justification=c(),
        legend.position=c())
Fig.pie.rare.ms

Fig.pd.rare.ms = ggplot(alpha.rare, aes(x=Size, y=PD, fill=Size, color=Size, linetype= Host, group=Host)) + 
  geom_point(position=pd, size=5)+ labs(x=" ", y ="Faith's PD")+
  geom_smooth(method=lm, color="black", fill="lightgrey", se=TRUE)+
  scale_fill_brewer(palette="Dark2")+
  scale_color_brewer(palette="Dark2")+
  theme_classic()+
  theme(text = element_text(size = 14), axis.text.x=element_text(angle=20, hjust=1))+
  theme(legend.justification=c(),
        legend.position=c())
Fig.pd.rare.ms

#######################
###alpha core microbiota
###
###
######################
nifH.core.ps.ex


hillq0.core <- estimate_richness(nifH.core.ps.ex, measures = "Observed")

nifH.ps.core.dataset.norm <- microbiome::transform(nifH.core.ps.ex, transform = "compositional")

shannon.core.hill <- exp(vegan::diversity(otu_table(nifH.ps.core.dataset.norm), index = "shannon"))
shannon.core.hill.df <- as.data.frame(shannon.core.hill)

simpson.core.hill <- 1/(1-(vegan::diversity(otu_table(nifH.ps.core.dataset.norm), index = "simpson")))
simpson.core.hill.df <- as.data.frame(simpson.core.hill)

pielou.core.hill <- vegan::diversity(otu_table(nifH.ps.core.dataset.norm), index = "shannon")/log(specnumber(otu_table(nifH.ps.core.dataset.norm)))
pielou.core.hill.df = as.data.frame(pielou.core.hill)

library(picante)
faith_phydiv.core = pd(otu_table(nifH.ps.core.dataset.norm), phy_tree(nifH.ps.core.dataset.norm), include.root=FALSE)
faith_phydiv.core.df = as.data.frame(faith_phydiv.core)

faith_phydiv.core.df
alpha.core<- cbind(hillq0.core, shannon.core.hill.df, simpson.core.hill.df, pielou.core.hill.df, faith_phydiv.core.df)
write.csv(alpha.core, "./Tables_pub/alpha_corediversity_Abr23.csv")
alpha.core = read.csv("./Tables_pub/alpha_corediversity_Abr23.csv")

prich.core.size = plot_richness(nifH.core.ps.ex, x="Size", measures=c("Observed")) +
  geom_point(pch = 21, fill = "White", col="white")+
  geom_boxplot(aes(fill = Host))+xlab("Size")+
  scale_fill_brewer(palette = "Dark2")+
  theme_classic() 
prich.core.size + theme(text = element_text(size = 18)) 

#alpha_data<- estimate_richness(nifH.merge.ps, measures = c("Shannon", "Simpson", "InvSimpson","Chao1"))
alpha.core$Host = nifH.core.ps.ex@sam_data$Host
alpha.core$Size = nifH.core.ps.ex@sam_data$Size
alpha.core$Reserve = nifH.core.ps.ex@sam_data$Reserve
alpha.core$Class.size = nifH.core.ps.ex@sam_data$Class.size
alpha.core$Pop = nifH.core.ps.ex@sam_data$Pop
alpha.core$sites = nifH.core.ps.ex@sam_data$Reserve
head(alpha.core)
str(alpha.core)
#test_names <- c( "Shannon", "Simpson", "InvSimpson","Chao1")
test_names <- c( "Observed", "shannon.core.hill", "simpson.core.hill", "pielou.core.hill", "PD")
shapiro.test(alpha.core$Observed)
shapiro.test(alpha.core$shannon.core.hill)
shapiro.test(alpha.core$simpson.core.hill)
shapiro.test(alpha.core$PD)
shapiro.test(alpha.core$pielou.core.hill)

par(mfrow=c(1,5))
hobsc = hist(alpha.core$Observed)
hshac = hist(alpha.core$shannon.core.hill)
hsimc = hist(alpha.core$simpson.core.hill)
hpiec = hist(alpha.core$pielou.core.hill)
hphyc = hist(alpha.core$PD)
head(alpha.core)
core.test_names <- c("Observed", "shannon.core.hill", "simpson.core.hill", "pielou.core.hill", "PD")
for (test in core.test_names) {
  test_formula <- paste(test,"~","Host")
  print(paste("Testing homogeneity of variance:", test))
  print(leveneTest(as.formula(test_formula), data=alpha.core))
}

head(alpha.core)
meandiv.core = data.frame(dfrich.obs.plot <- data_summary(alpha.core, varname="Observed", groupnames=c("Host", "Size", "sites")),
                          dfrich.sha.plot <- data_summary(alpha.core, varname="shannon.core.hill", groupnames=c("Host", "Size", "sites")),
                          dfrich.simp.plot <- data_summary(alpha.core, varname="simpson.core.hill", groupnames=c("Host", "Size", "sites")),
                          dfrich.pie.plot <- data_summary(alpha.core, varname="pielou.core.hill", groupnames=c("Host", "Size", "sites")),
                          dfrich.pd.plot <- data_summary(alpha.core, varname="PD", groupnames=c("Host", "Size", "sites")))

write.csv(meandiv.core, "Tables_pub/MeanRichness_CoreMicrobiota_Jun24.csv")
meandiv.core.Size = data.frame(dfrich.obs.plot <- data_summary(alpha.core, varname="Observed", groupnames=c("Host", "Size")),
                               dfrich.sha.plot <- data_summary(alpha.core, varname="shannon.core.hill", groupnames=c("Host", "Size")),
                               dfrich.simp.plot <- data_summary(alpha.core, varname="simpson.core.hill", groupnames=c("Host", "Size")),
                               dfrich.pie.plot <- data_summary(alpha.core, varname="pielou.core.hill", groupnames=c("Host", "Size")),
                               dfrich.pd.plot <- data_summary(alpha.core, varname="PD", groupnames=c("Host", "Size")))
table(alpha.core$Host, alpha.core$Size)
write.csv(meandiv.core.Size, "Tables_pub/MeanRichness_CoreMicrobiota_Size_Jun24.csv")

core.test_names <- c("Observed", "shannon.core.hill", "simpson.core.hill", "pielou.core.hill", "PD")

for (test in core.test_names) {
  test_formula <- paste(test,"~","Host")
  print(paste("Mesure d'alpha diversité:", test))
  print(summary(glm(as.formula(test_formula), data=alpha.core)))
}

#Fragment size stats

for (test in core.test_names) {
  test_formula <- paste(test,"~","Size")
  print(paste("Mesure d'alpha diversité:", test))
  print(summary(glm(as.formula(test_formula), data=alpha.core[c(1:19),])))
}
View(alpha.core)
for (test in core.test_names) {
  test_formula <- paste(test,"~","Size")
  print(paste("Mesure d'alpha diversité:", test))
  print(summary(glm(as.formula(test_formula), data=alpha.core[c(30:88),])))
}

#Multiple tests non paramétriques rangs de Wilcox
for (test in test_names) {
  test_formula <- paste(test,"~","Host")
  print(paste("Mesure d'alpha diversité:", test))
  print(pairwise.wilcox.test(alpha.core[,test], alpha.core$Host, p.adjust.method = "bonferroni"))
}

#Fragment size stats

for (test in test_names) {
  test_formula <- paste(test,"~","Size")
  print(paste("Mesure d'alpha diversité:", test))
  print(kruskal.test(as.formula(test_formula), data=alpha.core))
}

#Multiple tests non paramétriques rangs de Wilcox
for (test in test_names) {
  test_formula <- paste(test,"~","Size")
  print(paste("Mesure d'alpha diversité:", test))
  print(pairwise.wilcox.test(alpha.core[,test], alpha.core$Size, p.adjust.method = "bonferroni"))
}

######
# core Diversity figures
#
####
head(alpha.core)
Fig.div.core.ms = 
  ggplot(alpha.core, aes(x=Size, y=Observed, fill=Size, color=Size, linetype= Host, group=Host)) + 
  #ggplot(alpha_data, aes(x=Size, y=Observed)) + 
  ggtitle("Alpha Diversity", subtitle="Core diazotrophic microbiota")+
  geom_point(position=pd, size=5)+ labs(x=" ", y ="Number of species (q = 0)")+
  geom_smooth(method=lm, color="black", fill="lightgrey", se=TRUE)+
  scale_fill_brewer(palette="Dark2")+
  scale_color_brewer(palette="Dark2")+
  theme_classic()+
  theme(text = element_text(size = 14), axis.text.x=element_text(angle=20, hjust=1))+
  theme(legend.justification=c(),
        legend.position=c()) 
Fig.div.core.ms

Fig.sha.core.ms = 
  ggplot(alpha.core, aes(x=Size, y=shannon.core.hill, fill=Size, color=Size, linetype= Host, group=Host)) + 
  #ggplot(alpha_data, aes(x=Size, y=shannon.hill)) + 
  geom_point(position=pd, size=5)+ labs(x=" ", y ="Shannon entropy (q = 1)")+
  geom_smooth(method=lm, color="black", fill="lightgrey", se=TRUE)+
  scale_fill_brewer(palette="Dark2")+
  scale_color_brewer(palette="Dark2")+
  theme_classic()+
  theme(text = element_text(size = 14), axis.text.x=element_text(angle=20, hjust=1))+
  theme(legend.justification=c(),
        legend.position=c()) 
Fig.sha.core.ms

Fig.sim.core.ms = ggplot(alpha.core, aes(x=Size, y=simpson.core.hill, fill=Size, color=Size, linetype= Host, group=Host)) + 
  geom_point(position=pd, size=5)+ labs(x="Forest fragment size", y ="inverse Simpson index (q = 2)")+
  geom_smooth(method=lm, color="black", fill="lightgrey", se=TRUE)+
  scale_fill_brewer(palette="Dark2")+
  scale_color_brewer(palette="Dark2")+
  theme_classic()+
  theme(text = element_text(size = 14), axis.text.x=element_text(angle=20, hjust=1))+
  theme(legend.justification=c(),
        legend.position=c())
Fig.sim.core.ms

Fig.pie.core.ms = ggplot(alpha.core, aes(x=Size, y=pielou.core.hill, fill=Size, color=Size, linetype= Host, group=Host)) + 
  geom_point(position=pd, size=5)+ labs(x=" ", y ="Pielou Index")+
  geom_smooth(method=lm, color="black", fill="lightgrey", se=TRUE)+
  scale_fill_brewer(palette="Dark2")+
  scale_color_brewer(palette="Dark2")+
  theme_classic()+
  theme(text = element_text(size = 14), axis.text.x=element_text(angle=20, hjust=1))+
  theme(legend.justification=c(),
        legend.position=c())
Fig.pie.core.ms

Fig.pd.core.ms = ggplot(alpha.core, aes(x=Size, y=PD, fill=Size, color=Size, linetype= Host, group=Host)) + 
  geom_point(position=pd, size=5)+ labs(x=" ", y ="Faith's PD")+
  geom_smooth(method=lm, color="black", fill="lightgrey", se=TRUE)+
  scale_fill_brewer(palette="Dark2")+
  scale_color_brewer(palette="Dark2")+
  theme_classic()+
  theme(text = element_text(size = 14), axis.text.x=element_text(angle=20, hjust=1))+
  theme(legend.justification=c(),
        legend.position=c())
Fig.pd.core.ms

##Arrange manuscript figure
fig.div.pub = ggarrange(Fig.div.ms, Fig.sha.ms, Fig.sim.ms, 
                        align = c("hv"), common.legend = TRUE, legend="right")
fig.div.pub

fig.div.pub = ggarrange(Fig.div.ms, Fig.sha.ms, Fig.sim.ms, Fig.pie.ms, Fig.pd.ms,
                        Fig.div.core.ms, Fig.sha.core.ms, Fig.sim.core.ms, Fig.pie.core.ms, Fig.pd.core.ms,
                        ncol = 5, nrow = 2, align = c("hv"), common.legend = TRUE, legend="bottom")
fig.div.pub

fig.div.pub.core.rare = ggarrange(Fig.div.rare.ms, Fig.sha.rare.ms, Fig.sim.rare.ms, Fig.pie.rare.ms, Fig.pd.rare.ms,
                                  Fig.div.core.ms, Fig.sha.core.ms, Fig.sim.core.ms, Fig.pie.core.ms, Fig.pd.core.ms,
                                  ncol = 5, nrow = 2, align = c("hv"), common.legend = TRUE, legend="bottom")
fig.div.pub.core.rare

# linear regression
colnames(alpha_data)
colnames(alpha.rare) = c("SpecimenID", "Observed", "shannon.hill", "simpson.hill", "pielou.hill", "PD", "SR", "Host", "Size",
                         "Reserve", "Class.size", "Pop", "sites")
str(alpha.rare)
alpha.rare$Class.size = as.factor(alpha.rare$Class.size)
alpha.rare$Observed = as.numeric(alpha.rare$Observed)
colnames(alpha.core) = c("SpecimenID", "Observed", "shannon.hill", "simpson.hill", "pielou.hill", "PD", "SR", "Host", "Size",
                         "Reserve", "Class.size", "Pop", "sites")
str(alpha.core)
alpha.core$Observed = as.numeric(alpha.core$Observed)
alpha.core$Class.size = as.factor(alpha.core$Class.size)
str(alpha_data)

#Observedglm(Observed~Host*Class.size, data=alpha_data),
obs.cof = plot_summs(glm(Observed~Class.size, data=alpha.rare[c(77:146),]),glm(Observed~Class.size, data=alpha.rare[c(1:76),]), glm(Observed~Class.size, data=alpha.core[c(30:88),]), glm(Observed~Class.size, data=alpha.core[c(1:29),]),
                     plot.distributions = FALSE, inner_ci_level = .9, point.size=6, point.shape = c(21,22,21,22), line.size=c(1,2.5), 
                     colors=c("#E69F00", "#E69F00", "#56B4E9", "#56B4E9"), legend.title="Microbial dataset", model.names =c("C surinamensis Rarefied", "R flaccida Rarefied",  "C surinamensis Core", "R flaccida Core"))+
  ggtitle("Alpha Diversity", subtitle="Number of species (q = 0)")+
  labs(y = "1- & 10-ha Fragments")+
  theme_classic()+
  theme(axis.text.y.left = element_text(color = "white"))+
  theme(text = element_text(size = 14), axis.text.y=element_text(angle=90, vjust=1, size = 0.1), axis.text.x=element_text(hjust=1))+
  theme(legend.justification=c(),
        legend.position=c()) 
obs.cof

sha.cof = plot_summs(glm(shannon.hill~Class.size, data=alpha.rare[c(77:146),]),glm(shannon.hill~Class.size, data=alpha.rare[c(1:76),]), glm(shannon.hill~Class.size, data=alpha.core[c(30:88),]), glm(shannon.hill~Class.size, data=alpha.core[c(1:29),]),
                     plot.distributions = FALSE, inner_ci_level = .9, point.size=6, point.shape = c(21,22,21,22), line.size=c(1,2.5), 
                     colors=c("#E69F00", "#E69F00", "#56B4E9", "#56B4E9"), legend.title="Microbial dataset", model.names =c("C surinamensis Rarefied", "R flaccida Rarefied",  "C surinamensis Core", "R flaccida Core"))+
  ggtitle("", subtitle="Shannon entropy (q = 1)")+
  labs(y = "1- & 10-ha Fragments")+
  theme_classic()+
  theme(axis.text.y.left = element_text(color = "white"))+
  theme(text = element_text(size = 14), axis.text.y=element_text(angle=90, vjust=1, size = 0.1), axis.text.x=element_text(hjust=1))+
  theme(legend.justification=c(),
        legend.position=c()) 
sha.cof

simp.cof = plot_summs(glm(simpson.hill~Class.size, data=alpha.rare[c(77:146),]),glm(simpson.hill~Class.size, data=alpha.rare[c(1:76),]), glm(simpson.hill~Class.size, data=alpha.core[c(30:88),]), glm(simpson.hill~Class.size, data=alpha.core[c(1:29),]),
                      plot.distributions = FALSE, inner_ci_level = .9, point.size=6, point.shape = c(21,22,21,22), line.size=c(1,2.5), 
                      colors=c("#E69F00", "#E69F00", "#56B4E9", "#56B4E9"), legend.title="Microbial dataset", model.names =c("C surinamensis Rarefied", "R flaccida Rarefied",  "C surinamensis Core", "R flaccida Core"))+
  ggtitle("", subtitle="Inverse Simpson index (q = 2)")+
  labs(y = "1- & 10-ha Fragments")+
  theme_classic()+
  theme(axis.text.y.left = element_text(color = "white"))+
  theme(text = element_text(size = 14), axis.text.y=element_text(angle=90, vjust=1, size = 0.1), axis.text.x=element_text(hjust=1))+
  theme(legend.justification=c(),
        legend.position=c()) 
simp.cof

pd.cof = plot_summs(glm(PD~Class.size, data=alpha.rare[c(77:146),]),glm(PD~Class.size, data=alpha.rare[c(1:76),]), glm(PD~Class.size, data=alpha.core[c(30:88),]), glm(PD~Class.size, data=alpha.core[c(1:29),]),
                    plot.distributions = FALSE, inner_ci_level = .9, point.size=6, point.shape = c(21,22,21,22), line.size=c(1,2.5), 
                    colors=c("#E69F00", "#E69F00", "#56B4E9", "#56B4E9"), legend.title="Microbial dataset", model.names =c("C surinamensis Rarefied", "R flaccida Rarefied",  "C surinamensis Core", "R flaccida Core"))+
  ggtitle("", subtitle="Faith's PD")+
  labs(y = "1- & 10-ha Fragments")+
  theme_classic()+
  theme(axis.text.y.left = element_text(color = "white"))+
  theme(text = element_text(size = 14), axis.text.y=element_text(angle=90, vjust=1, size = 0.1), axis.text.x=element_text(hjust=1))+
  theme(legend.justification=c(),
        legend.position=c()) 
pd.cof

install.packages("openxlsx", dependencies = TRUE)

export_summs(glm(Observed~Class.size*Host, data=alpha.rare),#glm(PD~Class.size, data=alpha.rare[c(1:76),]), 
             glm(Observed~Class.size*Host, data=alpha.core), #glm(PD~Class.size, data=alpha.core[c(1:29),]), 
             scale = TRUE,
             error_format = "[{conf.low}, {conf.high}]",statistics = "all",
             to.file = "xlsx", file.name = "Observed_glm_sept24.xlsx")
export_summs(glm(shannon.hill~Class.size*Host, data=alpha.rare),#glm(PD~Class.size, data=alpha.rare[c(1:76),]), 
             glm(shannon.hill~Class.size*Host, data=alpha.core), #glm(PD~Class.size, data=alpha.core[c(1:29),]), 
             scale = TRUE,
             error_format = "[{conf.low}, {conf.high}]", statistics = "all",
             to.file = "xlsx", file.name = "Shannon_glm_sept24.xlsx")
export_summs(glm(simpson.hill~Class.size*Host, data=alpha.rare),#glm(PD~Class.size, data=alpha.rare[c(1:76),]), 
             glm(simpson.hill~Class.size*Host, data=alpha.core), #glm(PD~Class.size, data=alpha.core[c(1:29),]), 
             scale = TRUE, statistics = "all",
             error_format = "[{conf.low}, {conf.high}]", to.file = "xlsx", file.name = "Simpson_glm_sept24.xlsx")
export_summs(glm(PD~Class.size*Host, data=alpha.rare),#glm(PD~Class.size, data=alpha.rare[c(1:76),]), 
             glm(PD~Class.size*Host, data=alpha.core), #glm(PD~Class.size, data=alpha.core[c(1:29),]), 
             scale = TRUE, statistics = "all",
             error_format = "[{conf.low}, {conf.high}]", to.file = "xlsx", file.name = "faith_glm_sept24.xlsx")
###Final figure
#Figure 1 coeficient
Figure1B = ggarrange(obs.cof, sha.cof, simp.cof, pd.cof,zipi,
                     ncol = 2, nrow = 2, align = c("hv"), 
                     common.legend = TRUE, legend="right")
Figure1B

# view all plots individually (not shown)
plot_list1
plot_list2
# Combine all plots
#install.packages("cowplot")
library(cowplot)
library(purrr)
plot_grid(plotlist = plot_list1,
          ncol = 2)

plot_grid(plotlist = plot_list2,
          ncol = 2)

head(alpha_data)
summary(glm(PD~Class.size, data=alpha_data))
fit.core = glm(shannon.core.hill~Class.size, data=alpha.core)
fit.rare = glm(shannon.rare.hill~Class.size, data=alpha.rare)


summary(fit.com)
summary(fit.core)
plot_summs(fit.com, plot.distributions = TRUE, inner_ci_level = .9)
plot_summs(fit.com, fit.rare, fit.core, plot.distributions = FALSE, inner_ci_level = .9)
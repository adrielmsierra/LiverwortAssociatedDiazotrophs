###################################################################################
#                           Bryophyte NifH (Diazotrophs)
#                 Relative abundance of diazotrophic microbiota
#
#                                  Adriel M. Sierra
#
#
#                                  Jun 2024
#
###################################################################################
# Dataset
nifH.ps.com.dataset = readRDS("./02_Phyloseq_data_analysis/nifH_ps_complete_dataset_tree_abr23.rds")

# Review number of taxa by different ranks 
#number of taxa in the dataset
unique(tax_table(nifH.ps.com.dataset)[, "Phylum"])
unique(tax_table(nifH.ps.com.dataset)[, "Class"])
unique(tax_table(nifH.ps.com.dataset)[, "Order"])
unique(tax_table(nifH.ps.com.dataset)[, "Family"])
unique(tax_table(nifH.ps.com.dataset)[, "Genus"])
#which are the most abundant taxa
View(sort(table(tax_table(nifH.ps.com.dataset)[, "Phylum"], exclude = NULL)))
View(sort(table(tax_table(nifH.ps.com.dataset)[, "Class"], exclude = NULL)))
View(sort(table(tax_table(nifH.ps.com.dataset)[, "Order"], exclude = NULL)))
View(sort(table(tax_table(nifH.ps.com.dataset)[, "Family"], exclude = NULL)))
View(sort(table(tax_table(nifH.ps.com.dataset)[, "Genus"], exclude = NULL)))
View(sort(taxa_sums(nifH.ps.com.dataset), decreasing = TRUE))
View(tax_table(nifH.ps.com.dataset@tax_table@.Data))

# Explore the relative abundance of the Diazotrophs taxa
#abundance community composition
nifH.ps.relative <- transform_sample_counts(nifH.ps.com.dataset, function(OTU) OTU/sum(OTU))

phylum.sp.abun =  plot_bar(nifH.ps.relative, fill="Phylum") + facet_wrap(~Host, scales="free_x")+
  theme(text = element_text(size = 6))+
  #scale_fill_brewer(palette = "Set3")+
  ylab("Relative abundance")+xlab("Samples")+
  theme_classic()
phylum.sp.abun + theme(axis.text.x=element_text(angle=90, hjust=1))

##Visualize Bacterial Orders
dat.ord <- nifH.ps.relative %>% tax_glom(taxrank = "Order") %>% psmelt()
dat.ord
str(dat.ord)
dat.ord$Phylum = as.factor(dat.ord$Phylum)
dat.ord$Order = as.factor(dat.ord$Order)
dat.ord$Sample = as.factor(dat.ord$Sample)
#dat.ord$Family = as.factor(dat.ord$Family)

dat.ord.cyano = subset(dat.ord, Phylum=="Cyanobacteria")
dat.ord.cyano.cs = subset(dat.ord.cyano, Species=="Cololejeunea surinamensis")
fit <- lm(Abundance~Size, data = dat.ord.cyano.cs)

plot_coefs(fit)
efs.phy <- hedg_g(dat.ord, Abundance~Size+Species+Phylum, keep_d = TRUE)
str(efs.phy)
var=data.frame(efs.phy$hedg_g, as.factor(efs.phy$Size_ref))
str(var)
colnames(var) <- c("hedg_g", "Size")

str(var)
binned_plot(dat.ord, Abundance~Class.size+Species+Phylum)+theme_classic() 
ecdf_plot(dat.ord, Abundance~Size+Species+Phylum)+theme_classic() 
auc(dat.ord, Abundance~Size+Species+Phylum)
avg.abun.sp.ord = data.frame(aggregate(dat.ord$Abundance, by=list(Category=dat.ord$Phylum, dat.ord$Order, dat.ord$Host, dat.ord$Sample, dat.ord$Size), FUN=mean))
colnames(avg.abun.sp.ord) = c("Phylum", "Order", "Host_Species", "Sample", "Size", "Abundance")
str(avg.abun.sp.ord)
avg.abun.sp.ord$Abundance

#plot
Fig.Abun.ord.com <- ggplot(avg.abun.sp.ord, aes(x=Phylum, y=Abundance, fill=Host_Species)) + 
  scale_fill_brewer(palette="Dark2")+
  #geom_violin()+
  geom_boxplot(width = 0.7)+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.2)+
  theme_classic()

Fig.Abun.ord.com + labs(x="Bacterial Order", y = "Relative Abundance")+ 
  theme(axis.text.x=element_text(angle=90, size = 8))
levels(avg.abun.sp.ord$Order)

Fig.Abun.ord.com.size <- ggplot(avg.abun.sp.ord, aes(x=Phylum, y=Abundance, fill=Size)) + 
  scale_fill_brewer(palette="Dark2")+
  #geom_violin()+
  geom_boxplot(width = 0.7)+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.2)+
  theme_classic()
Fig.Abun.ord.com.size + labs(x="Bacterial Order", y = "Relative Abundance")+ 
  theme(axis.text.x=element_text(angle=90, size = 8))

levels(avg.abun.sp.ord$Order)

abundant.orders = sort(table(tax_table(nifH.ps.relative)[, "Order"], exclude = NULL), decreasing = TRUE)>45
orders.abunt <- abundant.orders == TRUE
table(orders.abunt)
filt.ord = orders.abunt[1:6]
filt.ord = data.frame(filt.ord [-2])
class(filt.ord)
ordenes = row.names(filt.ord)

filt.ord.others = orders.abunt[7:17]
class(filt.ord.others)
ordenes.less = as.factor(row.names(filt.ord.others))
levels(ordenes.less)
print(ordenes.less)
avg.abun.sp.ord$Classorder = avg.abun.sp.ord$Order
table(avg.abun.sp.ord$Order == ordenes)
str(avg.abun.sp.ord)

####
addorder = function(x){
  # If population label x is present function will output y
  if(x=="Nostocales") y = "Nostocales" 
  if(x=="Clostridiales") y = "Clostridiales" 
  if(x=="<NA>") y = "<NA>"
  if(x=="Rhizobiales") y = "Rhizobiales" 
  if(x=="Bacillales") y = "Bacillales" 
  if(x=="Oscillatoriales") y = "Oscillatoriales" 
  if(x=="Rhodobacterales"|x=="Enterobacterales"|x=="Selenomonadales"|x=="Burkholderiales"|x=="Rhodospirillales" |x=="Nevskiales" |x=="Chroococcidiopsidales"|x=="Pseudomonadales"|x=="Chromatiales"|x=="Desulfovibrionales"|x=="Hydrogenophilales") y = "Others" 
  return(y)
}

avg.abun.sp.ord$Classorder =  sapply(avg.abun.sp.ord$Order, addorder) 
avg.abun.sp.ord$Classorder = as.factor(avg.abun.sp.ord$Classorder)
str(avg.abun.sp.ord)


Fig.Abun.ms.ord <- ggplot(avg.abun.sp.ord, aes(x=Classorder, y=Abundance, fill=Host_Species)) + 
  scale_fill_brewer(palette="Dark2")+
  #geom_violin()+
  #geom_bar(stat = "identity", position=position_dodge())+
  #geom_errorbar(aes(ymin=Abundance-se, ymax=Abundance+se), width=.2,
  #              position=position_dodge(.9))+
  geom_boxplot(width = 0.4)+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1)+
  theme_classic()
Fig.Abun.ms.ord = Fig.Abun.ms.ord + labs(x="Bacterial orders", y = "Relative Abundance")
Fig.Abun.ms.ord

Fig.Abun.ms.ord.size <- ggplot(avg.abun.sp.ord, aes(x=Classorder, y=Abundance, fill=Size)) + 
  scale_fill_brewer(palette="Dark2")+
  #geom_violin()+
  #geom_bar(stat = "identity", position=position_dodge())+
  #geom_errorbar(aes(ymin=Abundance-se, ymax=Abundance+se), width=.2,
  #              position=position_dodge(.9))+
  geom_boxplot(width = 0.4)+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.05)+
  theme_classic()
Fig.Abun.ms.ord.size = Fig.Abun.ms.ord.size + labs(x="Bacterial orders", y = "Relative Abundance")
Fig.Abun.ms.ord.size
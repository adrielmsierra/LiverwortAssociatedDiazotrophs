###################################################################################
#                           Bryophyte NifH (Diazotrophs)
#                         Metacommunity assembly processes
#
#                                  Adriel M. Sierra
#
#
#                                  Jun 2024
#
###################################################################################

#### ecological processes dominating the community assembly 

#https://chiliubio.github.io/microeco_tutorial/model-based-class.html#trans_nullmodel-class
(meco_nifHdataset.rare$sample_table$Size)
(meco_nifHdataset.rare$beta_diversity)

var.data = data.frame(meco_nifHdataset.rare$sample_table)
str(var.data)
# generate trans_nullmodel object
# as an example, we only use high abundance OTU with mean relative abundance > 0.0005
t1 <- trans_nullmodel$new(meco_nifHdataset.rare, filter_thres = 0.0005, add_data = var.data)

# use pH as the test variable
t1$cal_mantel_corr(use_env = "Size")
# return t1$res_mantel_corr
# plot the mantel correlogram
t1$plot_mantel_corr()

# see null.model parameter for other null models
# null model run 500 times for the example
t1$cal_ses_betampd(runs = 500, abundance.weighted = TRUE)
# return t1$res_ses_betampd
meco_nifHdataset.rare$sample_table$Size
# add betaNRI matrix to beta_diversity list
meco_nifHdataset.rare$beta_diversity[["betaNRI"]] <- t1$res_ses_betampd
# create trans_beta class, use measure "betaNRI"
trans_beta$new(meco_nifHdataset.rare, measure = "Bray")
meco_nifHdataset.rare$sample_table = meco_nifHdataset.rare$sample_table
t2 <- trans_beta$new(meco_nifHdataset.rare, group = "Size", measure = "betaNRI")
# transform the distance for each group
t2$cal_group_distance()
# see the help document for more methods, e.g. "anova" and "KW_dunn"
t2$cal_group_distance_diff(method = "wilcox")
# plot the results
g1 <- t2$plot_group_distance(add = "mean")
g1 + geom_hline(yintercept = -2, linetype = 2) + geom_hline(yintercept = 2, linetype = 2)


# null model run 500 times
t1$cal_ses_betamntd(runs = 500, abundance.weighted = TRUE, null.model = "taxa.labels")
# return t1$res_ses_betamntd


# result stored in t1$res_rcbray
t1$cal_rcbray(runs = 1000)
# return t1$res_rcbray

# use betaNTI and rcbray to evaluate processes
t1$cal_process(use_betamntd = TRUE, group = "Size")

t1$cal_process(use_betamntd = TRUE)

View(t1$res_process)
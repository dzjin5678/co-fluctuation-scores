# MODEL FITTING: whole-brain co-fluc amplitude DEPENDENT CHANGES IN regional co-fluctuation score.
library(dplyr)
library(R.matlab)
library(ggsegGlasser)
library(ggsegSchaefer)
library(ggseg)
library(ggplot2)
library(ggseg3d)
library(cifti)
library(stringr)
library(factoextra)
library(matrixStats)
library(scales)
library(Hmisc)
library(tidyr)
library(cocor)
source("./GAMs_fixed.R")


cs_name = "1_region_level_ratioOfMeans"
# cs_name = "2_region_level_meanOfRatios"
# cs_name = "3_region_level_leaveOneRegionOut"
# cs_name = "4_region_level_regressOneRegionOut"
# cs_name = "5_region_level_norm_regress_beta"
# cs_name = "6_region_level_norm_zscore"
# cs_name = "7_region_level_means"

# default settings.
high_fd = 0; # remove high head motion frames.
num_of_bins = 20; # 20 bins.
atlas_i=2; # schaefer400x17.
regress_cofounds=0; # regress fd, gs, hr and br.
gsr=0; # no global signal regression.

# extract atlas name and number of regions.
atlas_names=c("schaefer200", "schaefer400", "hcpmmp");
atlas_rois=c(200, 400, 360);
atlas = atlas_names[atlas_i];
num_of_rois = atlas_rois[atlas_i];
# extract pipeline name
if (gsr==1){
  pipeline = 'proc_regress_fixWglob_ts_mat';
}else if (gsr==0){
  pipeline = 'proc_regress_fix_ts_mat';
}

dataset <- "HCP"
gam_data_filename=sprintf("bins_gam_data_zscored_bins_rss_high_fd-%s_num_of_bins-%s_atlas-%s_regress_cofounds-%s_gsr-%s.csv", 
                      high_fd, num_of_bins, atlas, regress_cofounds, pipeline);
print(gam_data_filename)
gam_amp_filename=sprintf("bins_gam_amp_zscored_bins_rss_high_fd-%s_num_of_bins-%s_atlas-%s_regress_cofounds-%s_gsr-%s.csv", 
                          high_fd, num_of_bins, atlas, regress_cofounds, pipeline);
gam_smooths_filename=sprintf("bins_gam_smooths_zscored_bins_rss_high_fd-%s_num_of_bins-%s_atlas-%s_regress_cofounds-%s_gsr-%s.csv", 
                         high_fd, num_of_bins, atlas, regress_cofounds, pipeline);
gam_estimated_smooths_filename=sprintf("bins_gam_estimated_smooths_zscored_bins_rss_high_fd-%s_num_of_bins-%s_atlas-%s_regress_cofounds-%s_gsr-%s.csv", 
                             high_fd, num_of_bins, atlas, regress_cofounds, pipeline);
smooth_var <- "zscored_bins_rss"
covariates <- "sex + age + tSNR + mean_fd + mean_gs + mean_br + mean_hr"
# begin gam.
gam.amp.schaefer400x17 <- matrix(data=NA, nrow=num_of_rois, ncol=12) # matrix to save gam.fit output to
for(row in c(1:num_of_rois)){
  region <- sprintf("V%s", row)
  # read data.
  gam.data <- read.csv(sprintf("../results_bins/DRIVER/%s/%s/%s", 
                               dataset, cs_name, gam_data_filename), header = FALSE)
  print(gam.data)
  colnames(gam.data)[401] <- smooth_var
  colnames(gam.data)[402] <- "sex"
  colnames(gam.data)[403] <- "subj_id"
  colnames(gam.data)[404] <- "mean_fd"
  colnames(gam.data)[405] <- "mean_gs"
  colnames(gam.data)[406] <- "mean_br"
  colnames(gam.data)[407] <- "mean_hr"
  colnames(gam.data)[408] <- "age"
  colnames(gam.data)[409] <- "tSNR"
  #run the gam.fit.smooth function
  GAM.RESULTS <- gam.fit.smooth(region = region, smooth_var = smooth_var, covariates = covariates,
                                knots = 3, set_fx = TRUE, stats_only = FALSE)
  gam.amp.schaefer400x17[row,] <- GAM.RESULTS
}
gam.amp.schaefer400x17 <- as.data.frame(gam.amp.schaefer400x17)
colnames(gam.amp.schaefer400x17) <- c("label", #region name
                               "GAM.age.Fvalue", #GAM F-value for the age smooth term
                               "GAM.age.pvalue", #GAM p-value for the age smooth term
                               "GAM.age.partialR2", #partial Rsq from age and age-null models
                               "Anova.age.pvalue", #Anova p-value comparing age and age-null models
                               "age.onsetchange", #age at which fluctuation amplitude starts significantly changing (first significant derivative)
                               "age.peakchange", #age at which fluctuation amplitude exhibits maximal change (largest significant derivative)
                               "minage.decrease", #age at which fluctuation amplitude starts significantly decreasing (first significant negative derivative)
                               "maxage.increase", #age at which fluctuation amplitude stops significantly increasing (last significant positive derivative)
                               "age.maturation", #age at which fluctuation amplitude stops changing (last significant derivative)
                               "mean.curvature", 
                               "mean.derivative_2l")
cols = c(2:12)
gam.amp.schaefer400x17[,cols] = apply(gam.amp.schaefer400x17[,cols], 2, function(x) as.numeric(as.character(x))) #format as numeric
write.csv(gam.amp.schaefer400x17,
          sprintf("../results_bins/DRIVER/%s/%s/%s", dataset, cs_name, gam_amp_filename),
          row.names = F, quote = F)


#### PREDICT GAM SMOOTH FITTED VALUES ####
##Function to predict fitted values of a measure based on a
# fitted GAM smooth (measure ~ s(smooth_var, k = knots, fx = set_fx) + covariates)) and a prediction df
# number of predictions to make
np <- 200
gam.smooths.schaefer400x17 <- matrix(data=NA, ncol=7)
colnames(gam.smooths.schaefer400x17) <- c(smooth_var,"fitted","se","lower","upper","index","label")

for(row in c(1:400)){
  region <- sprintf("V%s", row)
  GAM.SMOOTH <- gam.smooth.predict(region = region, smooth_var = smooth_var, covariates = covariates,
                                   knots = 3, set_fx = TRUE, increments = np) #run the gam.smooth.predict function
  preddata <- as.data.frame(GAM.SMOOTH[2]) #get predicted.smooth df from function output
  preddata$index <- rep(x=row, np) #region index
  preddata$label <- rep(x=GAM.SMOOTH[1], np) #label
  gam.smooths.schaefer400x17 <- rbind(gam.smooths.schaefer400x17, preddata)
}
gam.smooths.schaefer400x17 <- gam.smooths.schaefer400x17[-1,] #remove empty initialization row
gam.smooths.schaefer400x17$label <- as.character(gam.smooths.schaefer400x17$label)

write.csv(gam.smooths.schaefer400x17,
          sprintf("../results_bins/DRIVER/%s/%s/%s", dataset, cs_name, gam_smooths_filename),
          row.names = F, quote = F)


#### Region-wise GAM Smooth Estimates ####
np <- 200 # number of ages to evaluate the smooth at
gam.estimated.smooths.schaefer400x17 <- matrix(data=NA, ncol=4)
colnames(gam.estimated.smooths.schaefer400x17) <- c(smooth_var,"est","index","label")
for(row in c(1:400)){
  region <- sprintf("V%s", row)
  GAM.ESTIMATES <- gam.estimate.smooth(region = region, smooth_var = smooth_var, covariates = covariates,
                                       knots = 3, set_fx = TRUE, increments = np) # run the gam.estimate.smooth function
  GAM.ESTIMATES$index <- rep(x=row, np) # region index
  GAM.ESTIMATES$label <- rep(x=region, np) # label
  gam.estimated.smooths.schaefer400x17 <- rbind(gam.estimated.smooths.schaefer400x17, GAM.ESTIMATES)
}
gam.estimated.smooths.schaefer400x17 <- gam.estimated.smooths.schaefer400x17[-1,] #remove empty initialization row

write.csv(gam.estimated.smooths.schaefer400x17,
          sprintf("../results_bins/DRIVER/%s/%s/%s", dataset, cs_name, gam_estimated_smooths_filename),
          row.names = F, quote = F)
gc()






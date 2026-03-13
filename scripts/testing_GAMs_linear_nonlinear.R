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
library(mgcv)
library(gratia)
library(tidyverse)
library(numDeriv)


#### FIT GAM SMOOTH ####
gam.fit.smooth <- function(region, smooth_var, covariates, knots, set_fx = FALSE, stats_only = FALSE){
  
  parcel <- region
  
  #Fit the gam nonlinear
  modelformula <- as.formula(sprintf("%s ~ s(%s, k = %s, fx = %s) + %s", region, smooth_var, knots, set_fx, covariates))
  gam.model_nonlinear <- gam(modelformula, method = "REML", data = gam.data)
  
  #Fit the gam linear
  modelformula_linear <- as.formula(sprintf("%s ~ %s + %s", region, smooth_var, covariates))
  gam.model_linear <- gam(modelformula_linear, method = "REML", data = gam.data)
  
  res <- anova(gam.model_linear, gam.model_nonlinear, test="Chisq")
  return(res)
}


# default settings.
high_fd = 0; # remove high head motion frames.
num_of_bins = 20; # 20 bins.
atlas_i=1; # schaefer200x17.
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
gam_amp_filename=sprintf("bins_gam_amp_zscored_bins_rss_high_fd-%s_num_of_bins-%s_atlas-%s_regress_cofounds-%s_gsr-%s.csv", 
                         high_fd, num_of_bins, atlas, regress_cofounds, pipeline);
gam_smooths_filename=sprintf("bins_gam_smooths_zscored_bins_rss_high_fd-%s_num_of_bins-%s_atlas-%s_regress_cofounds-%s_gsr-%s.csv", 
                             high_fd, num_of_bins, atlas, regress_cofounds, pipeline);
gam_estimated_smooths_filename=sprintf("bins_gam_estimated_smooths_zscored_bins_rss_high_fd-%s_num_of_bins-%s_atlas-%s_regress_cofounds-%s_gsr-%s.csv", 
                                       high_fd, num_of_bins, atlas, regress_cofounds, pipeline);

# ratioOfMeans, meanOfRatios, norm_zscore, norm_regress_beta, norm_regress_bAndBeta, leaveOneRegionOut
cs_name <- "1_region_level_ratioOfMeans"
smooth_var <- "zscored_bins_rss"
covariates <- "sex + age + tSNR + mean_fd + mean_gs + mean_br + mean_hr"
# begin gam.
gam.nonlinear_linear_compare_res <- matrix(data=NA, nrow=200, ncol=2)
for(row in c(1:200)){
  region <- sprintf("V%s", row)
  # read data.
  gam.data <- read.csv(sprintf("../results_bins/DRIVER/%s/%s/%s", 
                               dataset, cs_name, gam_data_filename), header = FALSE)
  # print(gam.data)
  colnames(gam.data)[201] <- smooth_var
  colnames(gam.data)[202] <- "sex"
  colnames(gam.data)[203] <- "subj_id"
  colnames(gam.data)[204] <- "mean_fd"
  colnames(gam.data)[205] <- "mean_gs"
  colnames(gam.data)[206] <- "mean_br"
  colnames(gam.data)[207] <- "mean_hr"
  colnames(gam.data)[208] <- "age"
  colnames(gam.data)[209] <- "tSNR"
  #run the gam.fit.smooth function
  res <- gam.fit.smooth(region = region, smooth_var = smooth_var, covariates = covariates,
                                knots = 3, set_fx = TRUE, stats_only = FALSE)
  dev <- res$Deviance[2]
  p   <- res$`Pr(>Chi)`[2]
  # print(dev)
  # print(p)
  gam.nonlinear_linear_compare_res[row, 1] <- dev
  gam.nonlinear_linear_compare_res[row, 2] <- p
}

gam_nonlinear_linear_compare_res_filename <- sprintf("gam_nonlinear_linear_compare_high_fd-%s_num_of_bins-%s_atlas-%s_regress_cofounds-%s_gsr-%s.csv", 
                          high_fd, num_of_bins, atlas, regress_cofounds, pipeline);

write.csv(gam.nonlinear_linear_compare_res,
          sprintf("../results_bins/DRIVER/%s/%s/%s", dataset, cs_name, gam_nonlinear_linear_compare_res_filename),
          row.names = F, quote = F)


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

for (cs_name in c("1_region_level_ratioOfMeans")){
  
# for (cs_name in c("1_region_level_ratioOfMeans", "2_region_level_meanOfRatios", "3_region_level_leaveOneRegionOut", "4_region_level_regressOneRegionOut", 
#                   "5_region_level_norm_regress_beta", "6_region_level_norm_zscore", "7_region_level_means")){
  # default settings.
  high_fd = 0; # remove high head motion frames.
  num_of_bins = 20; # 20 bins.
  atlas_i=1; # schaefer200x17.
  regress_cofounds=1; # regress fd, gs, hr and br.
  gsr=0; # no global signal regression.
  
  # gam settings.
  for (k in c(5)){
  
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
    gam_amp_filename=sprintf("bins_gam_k_%s_amp_zscored_bins_rss_high_fd-%s_num_of_bins-%s_atlas-%s_regress_cofounds-%s_gsr-%s.csv", 
                              k, high_fd, num_of_bins, atlas, regress_cofounds, pipeline);
    gam_smooths_filename=sprintf("bins_gam_k_%s_smooths_zscored_bins_rss_high_fd-%s_num_of_bins-%s_atlas-%s_regress_cofounds-%s_gsr-%s.csv", 
                              k, high_fd, num_of_bins, atlas, regress_cofounds, pipeline);
    gam_estimated_smooths_filename=sprintf("bins_gam_k_%s_estimated_smooths_zscored_bins_rss_high_fd-%s_num_of_bins-%s_atlas-%s_regress_cofounds-%s_gsr-%s.csv", 
                              k, high_fd, num_of_bins, atlas, regress_cofounds, pipeline);
    smooth_var <- "zscored_bins_rss"
    covariates <- "sex + age + tSNR + mean_fd + mean_gs + mean_br + mean_hr"
    # covariates <- "sex + age + tSNR"
    # begin gam.
    gam.amp.schaefer200x17 <- matrix(data=NA, nrow=200, ncol=12) # matrix to save gam.fit output to
    gam.amp.cov.t.schaefer200x17 <- matrix(data=NA, nrow=200, ncol=8) # matrix to save gam.fit output to
    gam.amp.cov.p.schaefer200x17 <- matrix(data=NA, nrow=200, ncol=8) # matrix to save gam.fit output to
    gam.amp.k_check.schaefer200x17 <- matrix(data=NA, nrow=200, ncol=5) # matrix to save gam.fit output to
    
    # read data.
    gam.data <- read.csv(sprintf("../results_bins/DRIVER/%s/%s/%s", 
                                 dataset, cs_name, gam_data_filename), header = FALSE)
    colnames(gam.data)[201] <- smooth_var
    colnames(gam.data)[202] <- "sex"
    colnames(gam.data)[203] <- "subj_id"
    colnames(gam.data)[204] <- "mean_fd"
    colnames(gam.data)[205] <- "mean_gs"
    colnames(gam.data)[206] <- "mean_br"
    colnames(gam.data)[207] <- "mean_hr"
    colnames(gam.data)[208] <- "age"
    colnames(gam.data)[209] <- "tSNR"
    
    for(row in c(1:200)){
      region <- sprintf("V%s", row)
      print(region)
      #run the gam.fit.smooth function
      GAM.RESULTS <- gam.fit.smooth(region = region, smooth_var = smooth_var, covariates = covariates,
                                    knots = k, set_fx = TRUE, stats_only = FALSE)
      gam.amp.schaefer200x17[row,] <- GAM.RESULTS$results
      gam.amp.cov.t.schaefer200x17[row,] <- GAM.RESULTS$p.t
      gam.amp.cov.p.schaefer200x17[row,] <- GAM.RESULTS$p.pv
      gam.amp.k_check.schaefer200x17[row, 1:4] <- GAM.RESULTS$gam_check_k
      gam.amp.k_check.schaefer200x17[row, 5] <- GAM.RESULTS$aic_values
    }
    gam.amp.schaefer200x17 <- as.data.frame(gam.amp.schaefer200x17)
    colnames(gam.amp.schaefer200x17) <- c("label", #region name
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
    gam.amp.schaefer200x17[,cols] = apply(gam.amp.schaefer200x17[,cols], 2, function(x) as.numeric(as.character(x))) #format as numeric
    write.csv(gam.amp.schaefer200x17,
              sprintf("../results_bins/DRIVER/%s/%s/%s", dataset, cs_name, gam_amp_filename),
              row.names = F, quote = F)
    
    gam.amp.cov.t.schaefer200x17 <- as.data.frame(gam.amp.cov.t.schaefer200x17)
    colnames(gam.amp.cov.t.schaefer200x17) <- c("intercept","sex","age", "tSNR","mean_fd","mean_gs","mean_br","mean_hr")
    write.csv(gam.amp.cov.t.schaefer200x17,
              sprintf("../results_bins/DRIVER/%s/%s/cov_t_%s", dataset, cs_name, gam_amp_filename),
              row.names = F, quote = F)
    gam.amp.cov.p.schaefer200x17 <- as.data.frame(gam.amp.cov.p.schaefer200x17)
    colnames(gam.amp.cov.p.schaefer200x17) <- c("intercept","sex","age", "tSNR","mean_fd","mean_gs","mean_br","mean_hr")
    write.csv(gam.amp.cov.p.schaefer200x17,
              sprintf("../results_bins/DRIVER/%s/%s/cov_p_%s", dataset, cs_name, gam_amp_filename),
              row.names = F, quote = F)
    gam.amp.k_check.schaefer200x17 <- as.data.frame(gam.amp.k_check.schaefer200x17)
    colnames(gam.amp.k_check.schaefer200x17) <- c("k","edf","k_index", "p_value")
    write.csv(gam.amp.k_check.schaefer200x17,
              sprintf("../results_bins/DRIVER/%s/%s/k_check_aic_%s", dataset, cs_name, gam_amp_filename),
              row.names = F, quote = F)
    
    #### PREDICT GAM SMOOTH FITTED VALUES ####
    ##Function to predict fitted values of a measure based on a
    # fitted GAM smooth (measure ~ s(smooth_var, k = knots, fx = set_fx) + covariates)) and a prediction df
    # number of predictions to make
    np <- 200
    gam.smooths.schaefer200x17 <- matrix(data=NA, ncol=7)
    colnames(gam.smooths.schaefer200x17) <- c(smooth_var,"fitted","se","lower","upper","index","label")
    
    for(row in c(1:200)){
      region <- sprintf("V%s", row)
      GAM.SMOOTH <- gam.smooth.predict(region = region, smooth_var = smooth_var, covariates = covariates,
                                       knots = k, set_fx = TRUE, increments = np) #run the gam.smooth.predict function
      preddata <- as.data.frame(GAM.SMOOTH[2]) #get predicted.smooth df from function output
      preddata$index <- rep(x=row, np) #region index
      preddata$label <- rep(x=GAM.SMOOTH[1], np) #label
      gam.smooths.schaefer200x17 <- rbind(gam.smooths.schaefer200x17, preddata)
    }
    gam.smooths.schaefer200x17 <- gam.smooths.schaefer200x17[-1,] #remove empty initialization row
    gam.smooths.schaefer200x17$label <- as.character(gam.smooths.schaefer200x17$label)
    
    write.csv(gam.smooths.schaefer200x17,
              sprintf("../results_bins/DRIVER/%s/%s/%s", dataset, cs_name, gam_smooths_filename),
              row.names = F, quote = F)
    
    
    #### Region-wise GAM Smooth Estimates ####
    np <- 200 # number of ages to evaluate the smooth at
    gam.estimated.smooths.schaefer200x17 <- matrix(data=NA, ncol=4)
    colnames(gam.estimated.smooths.schaefer200x17) <- c(smooth_var,"est","index","label")
    for(row in c(1:200)){
      region <- sprintf("V%s", row)
      GAM.ESTIMATES <- gam.estimate.smooth(region = region, smooth_var = smooth_var, covariates = covariates,
                                           knots = k, set_fx = TRUE, increments = np) # run the gam.estimate.smooth function
      GAM.ESTIMATES$index <- rep(x=row, np) # region index
      GAM.ESTIMATES$label <- rep(x=region, np) # label
      gam.estimated.smooths.schaefer200x17 <- rbind(gam.estimated.smooths.schaefer200x17, GAM.ESTIMATES)
    }
    gam.estimated.smooths.schaefer200x17 <- gam.estimated.smooths.schaefer200x17[-1,] #remove empty initialization row
    
    write.csv(gam.estimated.smooths.schaefer200x17,
              sprintf("../results_bins/DRIVER/%s/%s/%s", dataset, cs_name, gam_estimated_smooths_filename),
              row.names = F, quote = F)
    gc()
  }
}


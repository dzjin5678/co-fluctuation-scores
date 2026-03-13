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

# default settings.
high_fd = 0; # remove high head motion frames.
num_of_bins = 20; # 20 bins.
atlas_i=1; # schaefer200x17.
regress_cofounds=0; # regress fd, gs, hr and br.
gsr=0; # no global signal regression.

# extract atlas name and number of regions.
atlas_names=c("schaefer200");
atlas_rois=c(190);
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
gam.amp.schaefer200x17 <- matrix(data=NA, nrow=200, ncol=12) # matrix to save gam.fit output to
for(row in c(1:190)){
  region <- sprintf("V%s", row)
  # read data.
  gam.data <- read.csv(sprintf("../results_bins/tSNR/%s/%s", 
                               cs_name, gam_data_filename), header = FALSE)
  # print(gam.data)
  colnames(gam.data)[201-10] <- smooth_var
  colnames(gam.data)[202-10] <- "sex"
  colnames(gam.data)[203-10] <- "subj_id"
  colnames(gam.data)[204-10] <- "mean_fd"
  colnames(gam.data)[205-10] <- "mean_gs"
  colnames(gam.data)[206-10] <- "mean_br"
  colnames(gam.data)[207-10] <- "mean_hr"
  colnames(gam.data)[208-10] <- "age"
  colnames(gam.data)[209-10] <- "tSNR"
  #run the gam.fit.smooth function
  GAM.RESULTS <- gam.fit.smooth(region = region, smooth_var = smooth_var, covariates = covariates,
                                knots = 3, set_fx = TRUE, stats_only = FALSE)
  gam.amp.schaefer200x17[row,] <- GAM.RESULTS
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
          sprintf("../results_bins/tSNR/%s/%s", cs_name, gam_amp_filename),
          row.names = F, quote = F)


#### PREDICT GAM SMOOTH FITTED VALUES ####
##Function to predict fitted values of a measure based on a
# fitted GAM smooth (measure ~ s(smooth_var, k = knots, fx = set_fx) + covariates)) and a prediction df
# number of predictions to make
np <- 200
gam.smooths.schaefer200x17 <- matrix(data=NA, ncol=7)
colnames(gam.smooths.schaefer200x17) <- c(smooth_var,"fitted","se","lower","upper","index","label")

for(row in c(1:190)){
  region <- sprintf("V%s", row)
  GAM.SMOOTH <- gam.smooth.predict(region = region, smooth_var = smooth_var, covariates = covariates,
                                   knots = 3, set_fx = TRUE, increments = np) #run the gam.smooth.predict function
  preddata <- as.data.frame(GAM.SMOOTH[2]) #get predicted.smooth df from function output
  preddata$index <- rep(x=row, np) #region index
  preddata$label <- rep(x=GAM.SMOOTH[1], np) #label
  gam.smooths.schaefer200x17 <- rbind(gam.smooths.schaefer200x17, preddata)
}
gam.smooths.schaefer200x17 <- gam.smooths.schaefer200x17[-1,] #remove empty initialization row
gam.smooths.schaefer200x17$label <- as.character(gam.smooths.schaefer200x17$label)

write.csv(gam.smooths.schaefer200x17,
          sprintf("../results_bins/tSNR/%s/%s", cs_name, gam_smooths_filename),
          row.names = F, quote = F)


#### Region-wise GAM Smooth Estimates ####
np <- 200 # number of ages to evaluate the smooth at
gam.estimated.smooths.schaefer200x17 <- matrix(data=NA, ncol=4)
colnames(gam.estimated.smooths.schaefer200x17) <- c(smooth_var,"est","index","label")
for(row in c(1:190)){
  region <- sprintf("V%s", row)
  GAM.ESTIMATES <- gam.estimate.smooth(region = region, smooth_var = smooth_var, covariates = covariates,
                                       knots = 3, set_fx = TRUE, increments = np) # run the gam.estimate.smooth function
  GAM.ESTIMATES$index <- rep(x=row, np) # region index
  GAM.ESTIMATES$label <- rep(x=region, np) # label
  gam.estimated.smooths.schaefer200x17 <- rbind(gam.estimated.smooths.schaefer200x17, GAM.ESTIMATES)
}
gam.estimated.smooths.schaefer200x17 <- gam.estimated.smooths.schaefer200x17[-1,] #remove empty initialization row

write.csv(gam.estimated.smooths.schaefer200x17,
          sprintf("../results_bins/tSNR/%s/%s", cs_name, gam_estimated_smooths_filename),
          row.names = F, quote = F)
gc()


###################################### plot ###################################### 
library(dplyr)
library(R.matlab)
library(ggsegGlasser)
library(ggsegSchaefer)
library(ggseg)
library(ggplot2)
library(ggdensity)
library(ggseg3d)
library(cifti)
library(stringr)
library(factoextra)
library(matrixStats)
library(scales)
library(Hmisc)
library(tidyr)
library(cocor)
library(devEMF)

hcp3t_unrelated100_tSNR_mean <- readMat("../outputs_bins/HCP/tSNR/hcp3t_unrelated100_tSNR_mean.mat")
hcp3t_unrelated100_tSNR_mean$low.tSNR.roi.idx

schaefer200x17_atlas <- read.csv(sprintf("../../atlas/schaefer200x17_atlas.CSV"))
schaefer200x17_atlas <- schaefer200x17_atlas[-hcp3t_unrelated100_tSNR_mean$low.tSNR.roi.idx, ]
# schaefer200x17_label <- read.csv("G:/datasets/atlas/Schaefer/Schaefer2018_200Parcels_17Networks_order_FSLMNI152_2mm.Centroid_RAS.csv")
# schaefer200x17_roi_names <- schaefer200x17_label$ROI.Name

scheafer200x17_mat <- readMat("../../atlas/Schaefer200x17.mat");
yeo8 <- scheafer200x17_mat$lab17to8
yeo8 <- yeo8[-hcp3t_unrelated100_tSNR_mean$low.tSNR.roi.idx, ]
# color settings.
color_arr <- array(data = NA, dim = 8, dimnames = NULL)
color_arr[1] <- rgb(162/255, 81/255, 172/255) # VIS
color_arr[2] <- rgb(120/255, 154/255, 192/255) # SMN
color_arr[3] <- rgb(64/255, 152/255, 50/255) # DAN
color_arr[4] <- rgb(224/255, 102/255, 254/255) # SVAN
color_arr[5] <- rgb(169/255, 169/255, 169/255) # LIM
color_arr[6] <- rgb(238/255, 185/255, 67/255) # CONT
color_arr[7] <- rgb(217/255, 113/255, 125/255) # DMN
color_arr[8] <- rgb(0/255, 0/255, 128/255) # TP

gam_fig_filename=sprintf("bins_gam_data_zscored_bins_rss_high_fd-%s_num_of_bins-%s_atlas-%s_regress_cofounds-%s_gsr-%s", 
                         high_fd, num_of_bins, atlas, regress_cofounds, pipeline);
dir.create(sprintf("../results_bins/tSNR/%s/%s", cs_name, gam_fig_filename))
gam.smoothestimates.schaefer200x17 <- gam.estimated.smooths.schaefer200x17
gam.smoothestimates.schaefer200x17 <- merge(gam.smoothestimates.schaefer200x17, gam.amp.schaefer200x17, by="label", sort = F)
gam.smoothestimates.schaefer200x17 <- merge(gam.smoothestimates.schaefer200x17, schaefer200x17_atlas, by="label", sort = F)

min_limit <- min(gam.smoothestimates.schaefer200x17$est);
max_limit <- max(gam.smoothestimates.schaefer200x17$est);

for(subnet_id in c(1:8)){
  subnet_data <- gam.smoothestimates.schaefer200x17[gam.smoothestimates.schaefer200x17$yeo8_label == subnet_id, ]
  file_path_name = sprintf("../results_bins/tSNR/%s/%s/%s.jpeg", cs_name, gam_fig_filename, as.character(subnet_id))
  jpeg(file = file_path_name,
       width = 100*15, height = 100*15, units = "px", res = 7*15);
  a <- ggplot(subnet_data, aes(zscored_bins_rss, est, group=index, color=color_arr[subnet_id])) +
    geom_line(size=4.0, alpha = 0.7) +
    theme_classic() +
    ylab(NULL) +
    xlab(NULL) +
    theme(legend.position = "none", plot.margin = unit(c(1,1,1,1), "lines")) +
    theme(axis.text = element_text(size=105, family = "Arial", color = c("black")),
          axis.title = element_text(size=95, family = "Arial", color = c("black")),
          axis.line = element_line(size=.22), axis.ticks = element_line(size=.22)) +
    scale_y_continuous(limits = c(min_limit, max_limit), breaks = c(-0.02, 0, 0.02)) +
    scale_color_manual(values=color_arr[subnet_id])
  print(a)
  dev.off()
}

for(subnet_id in c(1:8)){
  subnet_data <- gam.smoothestimates.schaefer200x17[gam.smoothestimates.schaefer200x17$yeo8_label == subnet_id, ]
  file_path_name = sprintf("../results_bins/tSNR/%s/%s/%s.emf", cs_name, gam_fig_filename, as.character(subnet_id))
  a <- ggplot(subnet_data, aes(zscored_bins_rss, est, group=index, color=color_arr[subnet_id])) +
    geom_line(size=0.5, alpha = 0.7) +
    theme_classic() +
    ylab(NULL) +
    xlab(NULL) +
    theme(legend.position = "none", plot.margin = unit(c(1,1,1,1), "lines")) +
    theme(axis.text = element_text(size=17, family = "Arial", color = c("black")),
          axis.line = element_line(size=.22), axis.ticks = element_line(size=.22)) +
    scale_y_continuous(limits = c(min_limit, max_limit), breaks = c(-0.02, 0, 0.02)) +
    scale_color_manual(values=color_arr[subnet_id])
  emf(file_path_name, width=2.3, height=2, emfPlus=FALSE)
  print(a)
  dev.off()
}



# MODEL FITTING: whole-brain co-fluc amplitude DEPENDENT CHANGES IN regional co-fluctuation score.
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


cs_name = "1_region_level_ratioOfMeans"
y_breaks = c(-0.02, 0, 0.02);
smooth_var <- "zscored_bins_rss"
covariates <- "sex + age + tSNR + mean_fd + mean_gs + mean_br + mean_hr"

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

schaefer200x17_atlas <- read.csv(sprintf("../../atlas/schaefer200x17_atlas.CSV"))
schaefer200x17_label <- read.csv("G:/datasets/atlas/Schaefer/Schaefer2018_200Parcels_17Networks_order_FSLMNI152_2mm.Centroid_RAS.csv")
schaefer200x17_roi_names <- schaefer200x17_label$ROI.Name

scheafer200x17_mat <- readMat("../../atlas/Schaefer200x17.mat");
yeo8 <- scheafer200x17_mat$lab17to8

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

gam_fig_filename=sprintf("bins_gam_six_regions_data_zscored_bins_rss_high_fd-%s_num_of_bins-%s_atlas-%s_regress_cofounds-%s_gsr-%s",
                         high_fd, num_of_bins, atlas, regress_cofounds, pipeline);
save_dir <- sprintf("../results_bins/DRIVER/%s/%s/%s", dataset, cs_name, gam_fig_filename)
dir.create(save_dir)

k=3;
gam_amp_filename=sprintf("../results_bins/DRIVER/%s/%s/bins_gam_k_%s_amp_zscored_bins_rss_high_fd-%s_num_of_bins-%s_atlas-%s_regress_cofounds-%s_gsr-%s.csv", 
                         dataset, cs_name, k, high_fd, num_of_bins, atlas, regress_cofounds, pipeline);
k_check_aic_3_filename=sprintf("../results_bins/DRIVER/%s/%s/k_check_aic_bins_gam_k_%s_amp_zscored_bins_rss_high_fd-%s_num_of_bins-%s_atlas-%s_regress_cofounds-%s_gsr-%s.csv", 
                         dataset, cs_name, k, high_fd, num_of_bins, atlas, regress_cofounds, pipeline);
gam_estimated_smooths_filename=sprintf("../results_bins/DRIVER/%s/%s/bins_gam_k_%s_estimated_smooths_zscored_bins_rss_high_fd-%s_num_of_bins-%s_atlas-%s_regress_cofounds-%s_gsr-%s.csv", 
                                       dataset, cs_name, k, high_fd, num_of_bins, atlas, regress_cofounds, pipeline);
# GAM fit
gam.k.3.amp.schaefer200x17 <- read.csv(gam_amp_filename)
# GAM smooth estimates
gam.k.3.smoothestimates.schaefer200x17 <- read.csv(gam_estimated_smooths_filename)
gam.k.3.smoothestimates.schaefer200x17 <- merge(gam.k.3.smoothestimates.schaefer200x17, gam.k.3.amp.schaefer200x17, by="label", sort = F)
gam.k.3.smoothestimates.schaefer200x17 <- merge(gam.k.3.smoothestimates.schaefer200x17, schaefer200x17_atlas, by="label", sort = F)
k_check_aic_3 <- read.csv(k_check_aic_3_filename)

k=5;
gam_amp_filename=sprintf("../results_bins/DRIVER/%s/%s/bins_gam_k_%s_amp_zscored_bins_rss_high_fd-%s_num_of_bins-%s_atlas-%s_regress_cofounds-%s_gsr-%s.csv", 
                         dataset, cs_name, k, high_fd, num_of_bins, atlas, regress_cofounds, pipeline);
k_check_aic_5_filename=sprintf("../results_bins/DRIVER/%s/%s/k_check_aic_bins_gam_k_%s_amp_zscored_bins_rss_high_fd-%s_num_of_bins-%s_atlas-%s_regress_cofounds-%s_gsr-%s.csv", 
                      dataset, cs_name, k, high_fd, num_of_bins, atlas, regress_cofounds, pipeline);
gam_estimated_smooths_filename=sprintf("../results_bins/DRIVER/%s/%s/bins_gam_k_%s_estimated_smooths_zscored_bins_rss_high_fd-%s_num_of_bins-%s_atlas-%s_regress_cofounds-%s_gsr-%s.csv", 
                                       dataset, cs_name, k, high_fd, num_of_bins, atlas, regress_cofounds, pipeline);
# GAM fit
gam.k.5.amp.schaefer200x17 <- read.csv(gam_amp_filename)
# GAM smooth estimates
gam.k.5.smoothestimates.schaefer200x17 <- read.csv(gam_estimated_smooths_filename)
gam.k.5.smoothestimates.schaefer200x17 <- merge(gam.k.5.smoothestimates.schaefer200x17, gam.k.5.amp.schaefer200x17, by="label", sort = F)
gam.k.5.smoothestimates.schaefer200x17 <- merge(gam.k.5.smoothestimates.schaefer200x17, schaefer200x17_atlas, by="label", sort = F)
k_check_aic_5 <- read.csv(k_check_aic_5_filename)


k=10;
gam_amp_filename=sprintf("../results_bins/DRIVER/%s/%s/bins_gam_k_%s_amp_zscored_bins_rss_high_fd-%s_num_of_bins-%s_atlas-%s_regress_cofounds-%s_gsr-%s.csv", 
                         dataset, cs_name, k, high_fd, num_of_bins, atlas, regress_cofounds, pipeline);
k_check_aic_10_filename=sprintf("../results_bins/DRIVER/%s/%s/k_check_aic_bins_gam_k_%s_amp_zscored_bins_rss_high_fd-%s_num_of_bins-%s_atlas-%s_regress_cofounds-%s_gsr-%s.csv", 
                      dataset, cs_name, k, high_fd, num_of_bins, atlas, regress_cofounds, pipeline);
gam_estimated_smooths_filename=sprintf("../results_bins/DRIVER/%s/%s/bins_gam_k_%s_estimated_smooths_zscored_bins_rss_high_fd-%s_num_of_bins-%s_atlas-%s_regress_cofounds-%s_gsr-%s.csv", 
                                       dataset, cs_name, k, high_fd, num_of_bins, atlas, regress_cofounds, pipeline);
# GAM fit
gam.k.10.amp.schaefer200x17 <- read.csv(gam_amp_filename)
# GAM smooth estimates
gam.k.10.smoothestimates.schaefer200x17 <- read.csv(gam_estimated_smooths_filename)
gam.k.10.smoothestimates.schaefer200x17 <- merge(gam.k.10.smoothestimates.schaefer200x17, gam.k.10.amp.schaefer200x17, by="label", sort = F)
gam.k.10.smoothestimates.schaefer200x17 <- merge(gam.k.10.smoothestimates.schaefer200x17, schaefer200x17_atlas, by="label", sort = F)
k_check_aic_10 <- read.csv(k_check_aic_10_filename)



################################################################################
################################################################################
######################### 以下进入到画图环节 ###################################
################################################################################                                                                            
################################################################################

plot_data <- k_check_aic_3$k_index
min_lim <- min(plot_data)
max_lim <- max(plot_data)
plot_data[k_check_aic_3$p_value > 0.05/200] <- NA
some_data = tibble(
  region = schaefer200x17_roi_names,
  p = as.vector(plot_data)
)
for (lr in c("left", "right")){
  for (lm in c("lateral", "medial")){
    file_path_name = sprintf("../results_bins/DRIVER/%s/%s/%s/all_k_3_%s_%s.jpeg", dataset, cs_name, gam_fig_filename, lr, lm)
    jpeg(file = file_path_name, width=1500, height=1000, units = "px",res =800)
    p <- ggseg(.data = some_data, atlas = schaefer17_200, mapping=aes(fill=p, colour=I("white")), position = c("stacked"), hemisphere=lr, view=lm) +
      theme_void() +
      labs(fill="") +
      paletteer::scale_fill_paletteer_c("ggthemes::Red-Gold", direction = -1, limits = c(round(min_lim, 2), round(max_lim, 2)),breaks = c(round(min_lim, 2), round(max_lim, 2)), oob = squish) +
      theme(legend.text = element_text(family = "Arial", size = 20, color = c("black")), legend.position = "none")
    p
    print(p)
    dev.off()
  }
}
# single maps.
file_path_name = sprintf("../results_bins/DRIVER/%s/%s/%s/all_k_3.jpeg", dataset, cs_name, gam_fig_filename)
jpeg(file = file_path_name, width=1500, height=1500, units = "px",res =300)
p <- ggseg(.data = some_data, atlas = schaefer17_200, mapping=aes(fill=p, colour=I("white")), position = c("stacked")) +
  theme_void() +
  labs(fill="") +
  paletteer::scale_fill_paletteer_c("ggthemes::Red-Gold", direction = -1,
                                    limits = c(round(min_lim, 1), round(max_lim, 1)),
                                    breaks = c(round(min_lim, 1), round(max_lim, 1)),
                                    oob = squish) +
  theme(legend.text = element_text(family = "Arial", size = 15, color = c("black")))
p
print(p)
dev.off()



# ################################ plot six representative regions ##############################
# min_limit <- min(gam.k.10.smoothestimates.schaefer200x17$est);
# max_limit <- max(gam.k.10.smoothestimates.schaefer200x17$est);
# # for (rep_resion_i in c(12, 59, 55, 144, 169, 12)){
# # for (rep_resion_i in c(11, 12, 27, 28, 32, 35, 40, 41, 51, 53, 57, 58, 80, 88, 99, 100)){
# for (rep_resion_i in c(1:100)){
#     
#   max_lim <- 1
#   min_lim <- 0
#   plot_data <- rep(NA, 200)
#   plot_data[rep_resion_i] <- 1
#   some_data = tibble(
#     region = schaefer200x17_roi_names,
#     p = as.vector(plot_data)
#   )
#   for (lr in c("left", "right")){
#     for (lm in c("lateral", "medial")){
#       file_path_name = sprintf("../results_bins/DRIVER/%s/%s/%s/six_regions_%s_%s_%s.jpeg", dataset, cs_name, gam_fig_filename, rep_resion_i, lr, lm)
#       jpeg(file = file_path_name, width=1500, height=1000, units = "px",res =800)
#       p <- ggseg(.data = some_data, atlas = schaefer17_200, mapping=aes(fill=p, colour=I("white")), position = c("stacked"), hemisphere=lr, view=lm) +
#         theme_void() +
#         labs(fill="") +
#         paletteer::scale_fill_paletteer_c("ggthemes::Orange-Blue Diverging", direction = -1, limits = c(round(min_lim, 2), round(max_lim, 2)),breaks = c(round(min_lim, 2), round(max_lim, 2)), oob = squish) +
#         theme(legend.text = element_text(family = "Arial", size = 20, color = c("black")), legend.position = "none")
#       p
#       print(p)
#       dev.off()
#     }
#   }
# 
#   region_i_data_k_3 <- gam.k.3.smoothestimates.schaefer200x17[gam.k.3.smoothestimates.schaefer200x17$label == sprintf("V%s", rep_resion_i), ]
#   region_i_data_k_5 <- gam.k.5.smoothestimates.schaefer200x17[gam.k.5.smoothestimates.schaefer200x17$label == sprintf("V%s", rep_resion_i), ]
#   region_i_data_k_10 <- gam.k.10.smoothestimates.schaefer200x17[gam.k.10.smoothestimates.schaefer200x17$label == sprintf("V%s", rep_resion_i), ]
#   
#   region_i_data_k_3$index <- 1
#   region_i_data_k_5$index <- 2
#   region_i_data_k_10$index <- 3
#   
#   region_i_data <- rbind(region_i_data_k_3, region_i_data_k_5, region_i_data_k_10)
# 
#   # browser()
#   file_path_name = sprintf("../results_bins/DRIVER/%s/%s/%s/%s.jpeg", dataset, cs_name, gam_fig_filename, as.character(rep_resion_i))
#   jpeg(file = file_path_name,
#        width = 100*15, height = 100*15, units = "px", res = 7*15);
#   a <- ggplot(region_i_data, aes(zscored_bins_rss, est, group=index, color=region_i_data$index)) +
#     geom_line(size=4.0, alpha = 0.7) +
#     theme_classic() +
#     ylab(NULL) +
#     xlab(NULL) +
#     paletteer::scale_colour_paletteer_c("ggthemes::Sunset-Sunrise Diverging") + 
#     theme(legend.position = "none", plot.margin = unit(c(1,1,1,1), "lines")) +
#     theme(axis.text = element_text(size=105, family = "Arial", color = c("black")),
#           axis.title = element_text(size=95, family = "Arial", color = c("black")),
#           axis.line = element_line(size=.22), axis.ticks = element_line(size=.22)) +
#     scale_y_continuous(limits = c(min_limit, max_limit), breaks = y_breaks)
#   print(a)
#   dev.off()
# }


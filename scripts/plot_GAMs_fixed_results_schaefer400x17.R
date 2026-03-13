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
# cs_name = "2_region_level_meanOfRatios"
# cs_name = "3_region_level_leaveOneRegionOut"
# cs_name = "4_region_level_regressOneRegionOut"
# cs_name = "5_region_level_norm_regress_beta"
# cs_name = "6_region_level_norm_zscore"
# cs_name = "7_region_level_means"

y_breaks = c(-0.02, 0, 0.02);
# y_breaks = c(-0.4, 0, 0.4);
# y_breaks = c(-0.05, 0, 0.1);
# y_breaks = c(-12, 0, 30);

smooth_var <- "zscored_bins_rss"
covariates <- "sex + age + tSNR + mean_fd + mean_gs + mean_br + mean_hr"

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
gam_data_filename=sprintf("../results_bins/DRIVER/%s/%s/bins_gam_data_zscored_bins_rss_high_fd-%s_num_of_bins-%s_atlas-%s_regress_cofounds-%s_gsr-%s.csv", 
                          dataset, cs_name, high_fd, num_of_bins, atlas, regress_cofounds, pipeline);
gam_amp_filename=sprintf("../results_bins/DRIVER/%s/%s/bins_gam_amp_zscored_bins_rss_high_fd-%s_num_of_bins-%s_atlas-%s_regress_cofounds-%s_gsr-%s.csv", 
                         dataset, cs_name, high_fd, num_of_bins, atlas, regress_cofounds, pipeline);
gam_smooths_filename=sprintf("../results_bins/DRIVER/%s/%s/bins_gam_smooths_zscored_bins_rss_high_fd-%s_num_of_bins-%s_atlas-%s_regress_cofounds-%s_gsr-%s.csv", 
                             dataset, cs_name, high_fd, num_of_bins, atlas, regress_cofounds, pipeline);
gam_estimated_smooths_filename=sprintf("../results_bins/DRIVER/%s/%s/bins_gam_estimated_smooths_zscored_bins_rss_high_fd-%s_num_of_bins-%s_atlas-%s_regress_cofounds-%s_gsr-%s.csv", 
                                       dataset, cs_name, high_fd, num_of_bins, atlas, regress_cofounds, pipeline);
gam_fig_filename=sprintf("bins_gam_data_zscored_bins_rss_high_fd-%s_num_of_bins-%s_atlas-%s_regress_cofounds-%s_gsr-%s", 
                          high_fd, num_of_bins, atlas, regress_cofounds, pipeline);
save_dir <- sprintf("../results_bins/DRIVER/%s/%s/%s", dataset, cs_name, gam_fig_filename)
dir.create(save_dir)

schaefer400x17_label <- read.csv("G:/datasets/atlas/Schaefer/Schaefer2018_400Parcels_17Networks_order_FSLMNI152_2mm.Centroid_RAS.csv")
schaefer400x17_roi_names <- schaefer400x17_label$ROI.Name
schaefer400x17_roi_names[1:200]  <- paste0("lh_", schaefer400x17_roi_names[1:200])
schaefer400x17_roi_names[201:400] <- paste0("rh_", schaefer400x17_roi_names[201:400])

############################## load plot data  #################################
# GAM fit
gam.amp.schaefer400x17 <- read.csv(gam_amp_filename)
# GAM smooth estimates
gam.smoothestimates.schaefer400x17 <- read.csv(gam_estimated_smooths_filename)
gam.smoothestimates.schaefer400x17 <- merge(gam.smoothestimates.schaefer400x17, gam.amp.schaefer400x17, by="label", sort = F)
# GAM smooth.
gam.smooths.schaefer400x17 <- read.csv(gam_smooths_filename)
gam.smooths.schaefer400x17 <- merge(gam.smooths.schaefer400x17, gam.amp.schaefer400x17, by="label", sort = F)
# GAM data
gam.data <- read.csv(gam_data_filename, header = FALSE)
colnames(gam.data)[401] <- smooth_var
colnames(gam.data)[402] <- "sex"
colnames(gam.data)[403] <- "subj_id"
colnames(gam.data)[404] <- "mean_fd"
colnames(gam.data)[405] <- "mean_gs"
colnames(gam.data)[406] <- "mean_br"
colnames(gam.data)[407] <- "mean_hr"
colnames(gam.data)[408] <- "age"
colnames(gam.data)[409] <- "tSNR"



################################################################################
################################################################################
######################### 以下进入到画图环节 ###################################
################################################################################                                                                            
################################################################################

min_limit <- min(gam.smoothestimates.schaefer400x17$est);
max_limit <- max(gam.smoothestimates.schaefer400x17$est);

# # display all trajectories colored by partial R2.
# q_r2 <- quantile(gam.smoothestimates.schaefer400x17$GAM.age.partialR2, probs = c(0.25, 0.75), na.rm = FALSE)
# df <- gam.smoothestimates.schaefer400x17
# df_filtered <- df %>% filter(GAM.age.partialR2 <= q_r2[1] | GAM.age.partialR2 >= q_r2[2])
# 
# color_full <- gam.smoothestimates.schaefer400x17$GAM.age.partialR2
# color <- df_filtered$GAM.age.partialR2
# 
# file_path_name = sprintf("../results_bins/DRIVER/%s/%s/%s/all_primary_effect.jpeg", dataset, cs_name, gam_fig_filename)
# 
# jpeg(file = file_path_name, width=1600, height=1500, units = "px",res =105)
# ggplot(df_filtered, aes(zscored_bins_rss, est, group=index, color=color)) +
#   geom_line(size=2.0, alpha = 0.8) +
#   paletteer::scale_color_paletteer_c("ggthemes::Orange-Blue Diverging", direction = -1, limits = c(min(color_full), max(color_full))) +
#   theme_classic() +
#   ylab(NULL) +
#   xlab(NULL) +
#   theme(legend.position = "none") +
#   theme(axis.text = element_text(size=70, family = "Arial", color = c("black")),
#         axis.title = element_text(size=65, family = "Arial", color = c("black")),
#         axis.line = element_line(size=.22), axis.ticks = element_line(size=.22)) +
#   scale_y_continuous(limits = c(min_limit, max_limit), breaks = y_breaks)
# dev.off()
# 
# 
# q_r2 <- quantile(gam.smoothestimates.schaefer400x17$GAM.age.partialR2, probs = c(0.25, 0.75), na.rm = FALSE)
# df <- gam.smoothestimates.schaefer400x17
# df_filtered <- df %>% filter(GAM.age.partialR2 <= q_r2[1] | GAM.age.partialR2 >= q_r2[2])
# color_full <- gam.smoothestimates.schaefer400x17$GAM.age.partialR2
# color <- df_filtered$GAM.age.partialR2
# file_path_name = sprintf("../results_bins/DRIVER/%s/%s/%s/all_primary_effect.emf", dataset, cs_name, gam_fig_filename)
# p <- ggplot(df_filtered, aes(zscored_bins_rss, est, group=index, color=color)) +
#   geom_line(size=0.5, alpha = 0.8) +
#   paletteer::scale_color_paletteer_c("ggthemes::Orange-Blue Diverging", direction = -1, limits = c(min(color_full), max(color_full))) +
#   theme_classic() +
#   ylab(NULL) +
#   xlab(NULL) +
#   theme(legend.position = "none") +
#   theme(axis.text = element_text(size=17, family = "Arial", color = c("black")),
#         axis.line = element_line(size=.22), axis.ticks = element_line(size=.22)) +
#   scale_y_continuous(limits = c(min_limit, max_limit), breaks = y_breaks)
# emf(file_path_name, width=3, height=3, emfPlus=FALSE)
# print(p)
# dev.off()
# 
# 
# 
# # display all trajectories colored by second derivations
# q_r2 <- quantile(gam.smoothestimates.schaefer400x17$mean.derivative_2l, probs = c(0.25, 0.75), na.rm = FALSE)
# df <- gam.smoothestimates.schaefer400x17
# df_filtered <- df %>% filter(mean.derivative_2l <= q_r2[1] | mean.derivative_2l >= q_r2[2])
# 
# color_full <- gam.smoothestimates.schaefer400x17$mean.derivative_2l
# color <- df_filtered$mean.derivative_2l
# file_path_name = sprintf("../results_bins/DRIVER/%s/%s/%s/all_second_effect.jpeg", dataset, cs_name, gam_fig_filename)
# jpeg(file = file_path_name, width=1600, height=1500, units = "px",res =105)
# ggplot(df_filtered, aes(zscored_bins_rss, est, group=index, color=color)) +
#   geom_line(size=2, alpha = 0.8) +
#   # paletteer::scale_color_paletteer_c("ggthemes::Orange-Blue-White Diverging", direction = 1, limits = c(min(color_full), max(color_full))) +
#   paletteer::scale_color_paletteer_c("grDevices::Geyser", direction = 1, limits = c(min(color), max(color))) +
#   theme_classic() +
#   ylab(NULL) +
#   xlab(NULL) +
#   theme(legend.position = "none") +
#   theme(axis.text = element_text(size=70, family = "Arial", color = c("black")),
#         axis.title = element_text(size=65, family = "Arial", color = c("black")),
#         axis.line = element_line(size=.22), axis.ticks = element_line(size=.22)) +
#   scale_y_continuous(limits = c(min_limit, max_limit), breaks = y_breaks)
# dev.off()
# 
# 
# q_r2 <- quantile(gam.smoothestimates.schaefer400x17$mean.derivative_2l, probs = c(0.25, 0.75), na.rm = FALSE)
# df <- gam.smoothestimates.schaefer400x17
# df_filtered <- df %>% filter(mean.derivative_2l <= q_r2[1] | mean.derivative_2l >= q_r2[2])
# color_full <- gam.smoothestimates.schaefer400x17$mean.derivative_2l
# color <- df_filtered$mean.derivative_2l
# file_path_name = sprintf("../results_bins/DRIVER/%s/%s/%s/all_second_effect.emf", dataset, cs_name, gam_fig_filename)
# p <- ggplot(df_filtered, aes(zscored_bins_rss, est, group=index, color=color)) +
#   geom_line(size=0.5, alpha = 0.8) +
#   paletteer::scale_color_paletteer_c("grDevices::Geyser", direction = 1, limits = c(min(color), max(color))) +
#   theme_classic() +
#   ylab(NULL) +
#   xlab(NULL) +
#   theme(legend.position = "none") +
#   theme(axis.text = element_text(size=17, family = "Arial", color = c("black")),
#         axis.line = element_line(size=.22), axis.ticks = element_line(size=.22)) +
#   scale_y_continuous(limits = c(min_limit, max_limit), breaks = y_breaks)
# emf(file_path_name, width=3, height=3, emfPlus=FALSE)
# print(p)
# dev.off()


################################ amplitude effect ##############################
effect_data <- gam.amp.schaefer400x17
some_data = tibble(
  label = schaefer400x17_roi_names,
  p = as.vector(effect_data$GAM.age.partialR2)
)
max_lim <- max(some_data$p)
min_lim <- min(some_data$p)

# four map.
for (lr in c("left", "right")){
  for (lm in c("lateral", "medial")){
    file_path_name = sprintf("../results_bins/DRIVER/%s/%s/%s/all_primary_effect_%s_%s.jpeg", dataset, cs_name, gam_fig_filename, lr, lm)
    jpeg(file = file_path_name, width=1500, height=1000, units = "px",res =800)
    p <- ggseg(.data = some_data, atlas = schaefer17_400, mapping=aes(fill=p, colour=I("white")), position = c("stacked"), hemisphere=lr, view=lm) +
      theme_void() +
      labs(fill="") +
      paletteer::scale_fill_paletteer_c("ggthemes::Orange-Blue Diverging", direction = -1, limits = c(round(min_lim, 2), round(max_lim, 2)),breaks = c(round(min_lim, 2), round(max_lim, 2)), oob = squish) +
      theme(legend.text = element_text(family = "Arial", size = 20, color = c("black")), legend.position = "none")
    p
    print(p)
    dev.off()
  }
}

# single maps.
file_path_name = sprintf("../results_bins/DRIVER/%s/%s/%s/all_primary_effect_brain_map.jpeg", dataset, cs_name, gam_fig_filename)
jpeg(file = file_path_name, width=1500, height=1500, units = "px",res =300)
p <- ggseg(.data = some_data, atlas = schaefer17_400, mapping=aes(fill=p, colour=I("white")), position = c("stacked")) +
  theme_void() +
  labs(fill="") +
  paletteer::scale_fill_paletteer_c("ggthemes::Orange-Blue Diverging", direction = -1, 
                                    limits = c(round(min_lim, 2), round(max_lim, 2)), 
                                    breaks = c(round(min_lim, 2), round(max_lim, 2)), 
                                    oob = squish) +
  theme(legend.text = element_text(family = "Arial", size = 20, color = c("black")), legend.position = "bottom")
p
print(p)
dev.off()


############################# second derivations ###############################
second_derivatives_data <- gam.amp.schaefer400x17
some_data = tibble(
  label = schaefer400x17_roi_names,
  p = as.vector(second_derivatives_data$mean.derivative_2l)
)
max_lim <- max(some_data$p)
min_lim <- min(some_data$p)
# single map.
for (lr in c("left", "right")){
  for (lm in c("lateral", "medial")){
    file_path_name = sprintf("../results_bins/DRIVER/%s/%s/%s/all_second_effect_%s_%s.jpeg", dataset, cs_name, gam_fig_filename, lr, lm)
    jpeg(file = file_path_name, width = 1500, height = 1000, units = "px", res = 800)
    p <- ggseg(.data = some_data, atlas = schaefer17_400, mapping=aes(fill=p, colour=I("white")), position = c("stacked"), hemisphere=lr, view=lm) +
      theme_void() +
      labs(fill="") +
      paletteer::scale_fill_paletteer_c("grDevices::Geyser", direction = 1,
                                        limits = c(round(min_lim, 5), round(max_lim, 5)),
                                        breaks = c(round(min_lim, 5), round(max_lim, 5)), oob = squish) +
      theme(legend.text = element_text(family = "Arial", size = 20, color = c("black")), legend.position = "none")
    p
    print(p)
    dev.off()
  }
}
# four map.
file_path_name = sprintf("../results_bins/DRIVER/%s/%s/%s/all_second_effect_brain_map.jpeg", dataset, cs_name, gam_fig_filename)
jpeg(file = file_path_name, width=1500, height=1500, units = "px",res =300)
p <- ggseg(.data = some_data, atlas = schaefer17_400, mapping=aes(fill=p, colour=I("white")), position = c("stacked")) +
  theme_void() +
  labs(fill="") +
  paletteer::scale_fill_paletteer_c("grDevices::Geyser", direction = 1, 
                                    limits = c(round(min_lim, 5), round(max_lim, 5)), 
                                    breaks = c(round(min_lim, 5), round(max_lim, 5)), oob = squish) +
  theme(legend.text = element_text(family = "Arial", size = 20, color = c("black")), legend.position = "bottom")
p
print(p)
dev.off()


############################## mean curvature ######################################
effect_data <- gam.amp.schaefer400x17
some_data = tibble(
  label = schaefer400x17_roi_names,
  p = as.vector(effect_data$mean.curvature)
)
max_lim <- max(some_data$p)
min_lim <- min(some_data$p)

# four map.
for (lr in c("left", "right")){
  for (lm in c("lateral", "medial")){
    file_path_name = sprintf("../results_bins/DRIVER/%s/%s/%s/all_mean_curv_%s_%s.jpeg", dataset, cs_name, gam_fig_filename, lr, lm)
    jpeg(file = file_path_name, width=1500, height=1000, units = "px",res =800)
    p <- ggseg(.data = some_data, atlas = schaefer17_400, mapping=aes(fill=p, colour=I("white")), position = c("stacked"), hemisphere=lr, view=lm) +
      theme_void() +
      labs(fill="") +
      paletteer::scale_fill_paletteer_c("grDevices::Geyser", direction = -1, limits = c(round(min_lim, 2), round(max_lim, 2)),breaks = c(round(min_lim, 2), round(max_lim, 2)), oob = squish) +
      theme(legend.text = element_text(family = "Arial", size = 20, color = c("black")), legend.position = "none")
    p
    print(p)
    dev.off()
  }
}

# single maps.
file_path_name = sprintf("../results_bins/DRIVER/%s/%s/%s/all_mean_curv_brain_map.jpeg", dataset, cs_name, gam_fig_filename)
jpeg(file = file_path_name, width=1500, height=1500, units = "px",res =300)
p <- ggseg(.data = some_data, atlas = schaefer17_400, mapping=aes(fill=p, colour=I("white")), position = c("stacked")) +
  theme_void() +
  labs(fill="") +
  paletteer::scale_fill_paletteer_c("grDevices::Geyser", direction = -1, 
                                    limits = c(round(min_lim, 2), round(max_lim, 2)), 
                                    breaks = c(round(min_lim, 2), round(max_lim, 2)), 
                                    oob = squish) +
  theme(legend.text = element_text(family = "Arial", size = 20, color = c("black")), legend.position = "bottom")
p
print(p)
dev.off()



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
# y_breaks = c(-0.4, 0, 0.3);
# y_breaks = c(-0.05, 0, 0.1);
# y_breaks = c(-1, 0, 2);
# y_breaks = c(-2, 0, 1);
# y_breaks = c(-10, 0, 30);

smooth_var <- "zscored_bins_rss"
covariates <- "sex + age + tSNR + mean_fd + mean_gs + mean_br + mean_hr"

# default settings.
high_fd = 0; # remove high head motion frames.
num_of_bins = 20; # 20 bins.
atlas_i=1; # schaefer200x17.
regress_cofounds=1; # regress fd, gs, hr and br.
gsr=0; # no global signal regression.

# k="best";
k=5;

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
gam_amp_filename=sprintf("../results_bins/DRIVER/%s/%s/bins_gam_k_%s_amp_zscored_bins_rss_high_fd-%s_num_of_bins-%s_atlas-%s_regress_cofounds-%s_gsr-%s.csv", 
                         dataset, cs_name, k, high_fd, num_of_bins, atlas, regress_cofounds, pipeline);
gam_smooths_filename=sprintf("../results_bins/DRIVER/%s/%s/bins_gam_k_%s_smooths_zscored_bins_rss_high_fd-%s_num_of_bins-%s_atlas-%s_regress_cofounds-%s_gsr-%s.csv", 
                             dataset, cs_name, k, high_fd, num_of_bins, atlas, regress_cofounds, pipeline);
gam_estimated_smooths_filename=sprintf("../results_bins/DRIVER/%s/%s/bins_gam_k_%s_estimated_smooths_zscored_bins_rss_high_fd-%s_num_of_bins-%s_atlas-%s_regress_cofounds-%s_gsr-%s.csv", 
                                       dataset, cs_name, k, high_fd, num_of_bins, atlas, regress_cofounds, pipeline);
gam_fig_filename=sprintf("bins_gam_k_%s_data_zscored_bins_rss_high_fd-%s_num_of_bins-%s_atlas-%s_regress_cofounds-%s_gsr-%s",
                         k, high_fd, num_of_bins, atlas, regress_cofounds, pipeline);
save_dir <- sprintf("../results_bins/DRIVER/%s/%s/%s", dataset, cs_name, gam_fig_filename)
dir.create(save_dir)
dir.create(sprintf("../results_bins/DRIVER/%s/%s/k_%s_cov_t_p_corrected", dataset, cs_name, k))

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


############################## load plot data  #################################
# GAM fit
gam.amp.schaefer200x17 <- read.csv(gam_amp_filename)
# GAM smooth estimates
gam.smoothestimates.schaefer200x17 <- read.csv(gam_estimated_smooths_filename)
gam.smoothestimates.schaefer200x17 <- merge(gam.smoothestimates.schaefer200x17, gam.amp.schaefer200x17, by="label", sort = F)
gam.smoothestimates.schaefer200x17 <- merge(gam.smoothestimates.schaefer200x17, schaefer200x17_atlas, by="label", sort = F)
# GAM smooth.
gam.smooths.schaefer200x17 <- read.csv(gam_smooths_filename)
gam.smooths.schaefer200x17 <- merge(gam.smooths.schaefer200x17, gam.amp.schaefer200x17, by="label", sort = F)
gam.smooths.schaefer200x17 <- merge(gam.smooths.schaefer200x17, schaefer200x17_atlas, by="label", sort = F)
# GAM data
gam.data <- read.csv(gam_data_filename, header = FALSE)
colnames(gam.data)[201] <- smooth_var
colnames(gam.data)[202] <- "sex"
colnames(gam.data)[203] <- "subj_id"
colnames(gam.data)[204] <- "mean_fd"
colnames(gam.data)[205] <- "mean_gs"
colnames(gam.data)[206] <- "mean_br"
colnames(gam.data)[207] <- "mean_hr"
colnames(gam.data)[208] <- "age"
colnames(gam.data)[209] <- "tSNR"



################################################################################
################################################################################
######################### 以下进入到画图环节 ###################################
################################################################################                                                                            
################################################################################


# ################################ covariates effect ##############################
# gam_cov_p_amp_filename=sprintf("../results_bins/DRIVER/%s/%s/cov_p_bins_gam_k_%s_amp_zscored_bins_rss_high_fd-%s_num_of_bins-%s_atlas-%s_regress_cofounds-%s_gsr-%s.csv",
#                                 dataset, cs_name, k, high_fd, num_of_bins, atlas, regress_cofounds, pipeline);
# gam_cov_t_amp_filename=sprintf("../results_bins/DRIVER/%s/%s/cov_t_bins_gam_k_%s_amp_zscored_bins_rss_high_fd-%s_num_of_bins-%s_atlas-%s_regress_cofounds-%s_gsr-%s.csv",
#                                 dataset, cs_name, k, high_fd, num_of_bins, atlas, regress_cofounds, pipeline);
# gam.amp.cov.p.schaefer200x17 <- read.csv(gam_cov_p_amp_filename)
# gam.amp.cov.t.schaefer200x17 <- read.csv(gam_cov_t_amp_filename)
# for (cov_i in c("sex","age", "tSNR","mean_fd","mean_gs","mean_br","mean_hr")){
#   plot_p <- gam.amp.cov.p.schaefer200x17[[cov_i]]
#   plot_p <- plot_p * 200
#   plot_t <- gam.amp.cov.t.schaefer200x17[[cov_i]]
#   plot_t[plot_p >= 0.01] <- NA
#   some_data = tibble(
#     region = schaefer200x17_roi_names,
#     p = as.vector(plot_t)
#   )
#   max_lim <- max(some_data$p, na.rm = TRUE)
#   min_lim <- min(some_data$p, na.rm = TRUE)
#   upper_lim <- max(abs(max_lim), abs(min_lim))
#   max_lim <- upper_lim
#   min_lim <- -upper_lim
#   if (all(is.na(plot_t))){
#     max_lim <- 1
#     min_lim <- -1
#   }
#   # four map.
#   for (lr in c("left", "right")){
#     for (lm in c("lateral", "medial")){
#       file_path_name = sprintf("../results_bins/DRIVER/%s/%s/k_%s_cov_t_p_corrected/cov_t_%s_all_primary_effect_%s_%s.jpeg", dataset, cs_name, k, cov_i, lr, lm)
#       jpeg(file = file_path_name, width=1500, height=1000, units = "px",res =800)
#       p <- ggseg(.data = some_data, atlas = schaefer17_200, mapping=aes(fill=p, colour=I("white")), position = c("stacked"), hemisphere=lr, view=lm) +
#         theme_void() +
#         labs(fill="") +
#         paletteer::scale_fill_paletteer_c("ggthemes::Orange-Blue Diverging", direction = -1,
#                                           limits = c(round(min_lim, 1), round(max_lim, 1)),
#                                           breaks = c(round(min_lim, 1), round(max_lim, 1)),
#                                           oob = squish) +
#         theme(legend.text = element_text(family = "Arial", size = 20, color = c("black")), legend.position = "none")
#       p
#       print(p)
#       dev.off()
#     }
#   }
#   # single maps.
#   file_path_name = sprintf("../results_bins/DRIVER/%s/%s/k_%s_cov_t_p_corrected/cov_t_%s_all_primary_effect_brain_map.jpeg", dataset, cs_name, k, cov_i)
#   jpeg(file = file_path_name, width=1500, height=1500, units = "px",res =300)
#   p <- ggseg(.data = some_data, atlas = schaefer17_200, mapping=aes(fill=p, colour=I("white")), position = c("stacked")) +
#     theme_void() +
#     labs(fill="") +
#     paletteer::scale_fill_paletteer_c("ggthemes::Orange-Blue Diverging", direction = -1,
#                                       limits = c(round(min_lim, 1), round(max_lim, 1)),
#                                       breaks = c(round(min_lim, 1), round(max_lim, 1)),
#                                       oob = squish) +
#     theme(legend.text = element_text(family = "Arial", size = 20, color = c("black")), legend.position = "bottom")
#   p
#   print(p)
#   dev.off()
# }
# browser()


# ########################### amplitude effect, bh corrected p values ##############################
# GAM.age.pvalue <- gam.amp.schaefer200x17$GAM.age.pvalue
# GAM.age.pvalue <- GAM.age.pvalue * 200
# GAM.age.pvalue.corrected <- GAM.age.pvalue < 0.05
# GAM.age.partialR2 <- gam.amp.schaefer200x17$GAM.age.partialR2
# GAM.age.partialR2[GAM.age.pvalue >= 0.05] <- NA
# some_data = tibble(
#   region = schaefer200x17_roi_names,
#   p = as.vector(GAM.age.partialR2)
# )
# max_lim <- max(some_data$p)
# min_lim <- min(some_data$p)
# # four map.
# for (lr in c("left", "right")){
#   for (lm in c("lateral", "medial")){
#     file_path_name = sprintf("../results_bins/DRIVER/%s/%s/%s/corrected_primary_effect_%s_%s.jpeg", dataset, cs_name, gam_fig_filename, lr, lm)
#     jpeg(file = file_path_name, width=1500, height=1000, units = "px",res =800)
#     p <- ggseg(.data = some_data, atlas = schaefer17_200, mapping=aes(fill=p, colour=I("white")), position = c("stacked"), hemisphere=lr, view=lm) +
#       theme_void() +
#       labs(fill="") +
#       paletteer::scale_fill_paletteer_c("ggthemes::Orange-Blue Diverging", direction = -1, limits = c(round(min_lim, 2), round(max_lim, 2)),breaks = c(round(min_lim, 2), round(max_lim, 2)), oob = squish) +
#       theme(legend.text = element_text(family = "Arial", size = 20, color = c("black")), legend.position = "none")
#     p
#     print(p)
#     dev.off()
#   }
# }
# print(schaefer200x17_roi_names[GAM.age.pvalue >= 0.05])
# browser()


############################# fitted trajectories ##############################
# display trajectory by yeo8.
min_limit <- min(gam.smoothestimates.schaefer200x17$est);
max_limit <- max(gam.smoothestimates.schaefer200x17$est);
for(subnet_id in c(1:8)){
  subnet_data <- gam.smoothestimates.schaefer200x17[gam.smoothestimates.schaefer200x17$yeo8_label == subnet_id, ]
  file_path_name = sprintf("../results_bins/DRIVER/%s/%s/%s/%s.jpeg", dataset, cs_name, gam_fig_filename, as.character(subnet_id))
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
    scale_y_continuous(limits = c(min_limit, max_limit), breaks = y_breaks) +
    scale_color_manual(values=color_arr[subnet_id])
  print(a)
  dev.off()
}
for(subnet_id in c(1:8)){
  subnet_data <- gam.smoothestimates.schaefer200x17[gam.smoothestimates.schaefer200x17$yeo8_label == subnet_id, ]
  file_path_name = sprintf("../results_bins/DRIVER/%s/%s/%s/%s.emf", dataset, cs_name, gam_fig_filename, as.character(subnet_id))
  a <- ggplot(subnet_data, aes(zscored_bins_rss, est, group=index, color=color_arr[subnet_id])) +
    geom_line(size=0.5, alpha = 0.7) +
    theme_classic() +
    ylab(NULL) +
    xlab(NULL) +
    theme(legend.position = "none", plot.margin = unit(c(1,1,1,1), "lines")) +
    theme(axis.text = element_text(size=17, family = "Arial", color = c("black")),
          axis.line = element_line(size=.22), axis.ticks = element_line(size=.22)) +
    scale_y_continuous(limits = c(min_limit, max_limit), breaks = y_breaks) +
    scale_color_manual(values=color_arr[subnet_id])
  emf(file_path_name, width=2.3, height=2, emfPlus=FALSE)
  print(a)
  dev.off()
}

# display all trajectories colored by partial R2.
q_r2 <- quantile(gam.smoothestimates.schaefer200x17$GAM.age.partialR2, probs = c(0.25, 0.75), na.rm = FALSE)
df <- gam.smoothestimates.schaefer200x17
df_filtered <- df %>% filter(GAM.age.partialR2 <= q_r2[1] | GAM.age.partialR2 >= q_r2[2])

# color_full <- gam.smoothestimates.schaefer200x17$GAM.age.partialR2
# min_lim <- min(color_full)
# max_lim <- max(color_full)

max_lim <- max(gam.smoothestimates.schaefer200x17$GAM.age.partialR2)
min_lim <- min(gam.smoothestimates.schaefer200x17$GAM.age.partialR2)
abs_max <- max(c(abs(max_lim), abs(min_lim)))
max_lim <- abs_max
min_lim <- -abs_max

color <- df_filtered$GAM.age.partialR2
file_path_name = sprintf("../results_bins/DRIVER/%s/%s/%s/all_primary_effect.jpeg", dataset, cs_name, gam_fig_filename)
jpeg(file = file_path_name, width=1600, height=1500, units = "px",res =105)
ggplot(df_filtered, aes(zscored_bins_rss, est, group=index, color=color)) +
  geom_line(size=2.0, alpha = 0.8) +
  paletteer::scale_color_paletteer_c("ggthemes::Orange-Blue Diverging", direction = -1, limits = c(min_lim, max_lim)) +
  theme_classic() +
  ylab(NULL) +
  xlab(NULL) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size=70, family = "Arial", color = c("black")),
        axis.title = element_text(size=65, family = "Arial", color = c("black")),
        axis.line = element_line(size=.22), axis.ticks = element_line(size=.22)) +
  scale_y_continuous(limits = c(min_limit, max_limit), breaks = y_breaks)
dev.off()

file_path_name = sprintf("../results_bins/DRIVER/%s/%s/%s/all_primary_effect.emf", dataset, cs_name, gam_fig_filename)
p <- ggplot(df_filtered, aes(zscored_bins_rss, est, group=index, color=color)) +
  geom_line(size=0.5, alpha = 0.8) +
  paletteer::scale_color_paletteer_c("ggthemes::Orange-Blue Diverging", direction = -1, limits = c(min_lim, max_lim)) +
  theme_classic() +
  ylab(NULL) +
  xlab(NULL) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size=17, family = "Arial", color = c("black")),
        axis.line = element_line(size=.22), axis.ticks = element_line(size=.22)) +
  scale_y_continuous(limits = c(min_limit, max_limit), breaks = y_breaks)
emf(file_path_name, width=3, height=3, emfPlus=FALSE)
print(p)
dev.off()
# browser()


# display all trajectories colored by second derivations
q_r2 <- quantile(gam.smoothestimates.schaefer200x17$mean.derivative_2l, probs = c(0.25, 0.75), na.rm = FALSE)
df <- gam.smoothestimates.schaefer200x17
df_filtered <- df %>% filter(mean.derivative_2l <= q_r2[1] | mean.derivative_2l >= q_r2[2])

# color_full <- gam.smoothestimates.schaefer200x17$mean.derivative_2l
# min_lim <- min(color_full)
# max_lim <- max(color_full)

max_lim <- max(gam.smoothestimates.schaefer200x17$mean.derivative_2l)
min_lim <- min(gam.smoothestimates.schaefer200x17$mean.derivative_2l)
abs_max <- max(c(abs(max_lim), abs(min_lim)))
max_lim <- abs_max
min_lim <- -abs_max

color <- df_filtered$mean.derivative_2l
file_path_name = sprintf("../results_bins/DRIVER/%s/%s/%s/all_second_effect.jpeg", dataset, cs_name, gam_fig_filename)
jpeg(file = file_path_name, width=1600, height=1500, units = "px",res =105)
ggplot(df_filtered, aes(zscored_bins_rss, est, group=index, color=color)) +
  geom_line(size=2, alpha = 0.8) +
  paletteer::scale_color_paletteer_c("grDevices::Geyser", direction = 1, limits = c(min_lim, max_lim)) +
  theme_classic() +
  ylab(NULL) +
  xlab(NULL) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size=70, family = "Arial", color = c("black")),
        axis.title = element_text(size=65, family = "Arial", color = c("black")),
        axis.line = element_line(size=.22), axis.ticks = element_line(size=.22)) +
  scale_y_continuous(limits = c(min_limit, max_limit), breaks = y_breaks)
dev.off()

file_path_name = sprintf("../results_bins/DRIVER/%s/%s/%s/all_second_effect.emf", dataset, cs_name, gam_fig_filename)
p <- ggplot(df_filtered, aes(zscored_bins_rss, est, group=index, color=color)) +
  geom_line(size=0.5, alpha = 0.8) +
  paletteer::scale_color_paletteer_c("grDevices::Geyser", direction = 1, limits = c(min_lim, max_lim)) +
  theme_classic() +
  ylab(NULL) +
  xlab(NULL) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size=17, family = "Arial", color = c("black")),
        axis.line = element_line(size=.22), axis.ticks = element_line(size=.22)) +
  scale_y_continuous(limits = c(min_limit, max_limit), breaks = y_breaks)
emf(file_path_name, width=3, height=3, emfPlus=FALSE)
print(p)
dev.off()


################################ amplitude effect ##############################
effect_data <- gam.amp.schaefer200x17
some_data = tibble(
  region = schaefer200x17_roi_names,
  p = as.vector(effect_data$GAM.age.partialR2)
)

max_lim <- max(some_data$p)
min_lim <- min(some_data$p)
abs_max <- max(c(abs(max_lim), abs(min_lim)))
max_lim <- abs_max
min_lim <- -abs_max

# four map.
for (lr in c("left", "right")){
  for (lm in c("lateral", "medial")){
    file_path_name = sprintf("../results_bins/DRIVER/%s/%s/%s/all_primary_effect_%s_%s.jpeg", dataset, cs_name, gam_fig_filename, lr, lm)
    jpeg(file = file_path_name, width=1500, height=1000, units = "px", res =800)
    p <- ggseg(.data = some_data, atlas = schaefer17_200, mapping=aes(fill=p, colour=I("white")), position = c("stacked"), hemisphere=lr, view=lm) +
      theme_void() +
      labs(fill="") +
      paletteer::scale_fill_paletteer_c("ggthemes::Orange-Blue Diverging", direction = -1, limits = c(round(min_lim, 2), round(max_lim, 2)),breaks = c(round(min_lim, 2), round(max_lim, 2)), oob = squish) +
      theme(legend.text = element_text(family = "Arial", size = 20, color = c("black")), legend.position = "none")
    p
    print(p)
    dev.off()
  }
}
# for (lr in c("right")){
#   for (lm in c("lateral", "medial")){
#     file_path_name = sprintf("../results_bins/DRIVER/%s/%s/%s/all_primary_effect_aa_%s.jpeg", dataset, cs_name, gam_fig_filename, lm)
#     jpeg(file = file_path_name, width=1500, height=1000, units = "px", res =800)
#     p <- ggseg(.data = some_data, atlas = schaefer17_200, mapping=aes(fill=p, colour=I("white")), position = c("stacked"), hemisphere=lr, view=lm) +
#       theme_void() +
#       labs(fill="") +
#       paletteer::scale_fill_paletteer_c("ggthemes::Orange-Blue Diverging", direction = -1, limits = c(round(min_lim, 2), round(max_lim, 2)),breaks = c(round(min_lim, 2), round(max_lim, 2)), oob = squish) +
#       theme(legend.text = element_text(family = "Arial", size = 20, color = c("black")), legend.position = "none")
#     p
#     print(p)
#     dev.off()
#   }
# }
# single maps.
file_path_name = sprintf("../results_bins/DRIVER/%s/%s/%s/all_primary_effect_brain_map.jpeg", dataset, cs_name, gam_fig_filename)
jpeg(file = file_path_name, width=1500, height=1500, units = "px",res =300)
p <- ggseg(.data = some_data, atlas = schaefer17_200, mapping=aes(fill=p, colour=I("white")), position = c("stacked")) +
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
second_derivatives_data <- gam.amp.schaefer200x17
some_data = tibble(
  region = schaefer200x17_roi_names,
  p = as.vector(second_derivatives_data$mean.derivative_2l)
)

max_lim <- max(some_data$p)
min_lim <- min(some_data$p)
abs_max <- max(c(abs(max_lim), abs(min_lim)))
max_lim <- abs_max
min_lim <- -abs_max

# single map.
for (lr in c("left", "right")){
  for (lm in c("lateral", "medial")){
    file_path_name = sprintf("../results_bins/DRIVER/%s/%s/%s/all_second_effect_%s_%s.jpeg", dataset, cs_name, gam_fig_filename, lr, lm)
    jpeg(file = file_path_name, width = 1500, height = 1000, units = "px", res = 800)
    p <- ggseg(.data = some_data, atlas = schaefer17_200, mapping=aes(fill=p, colour=I("white")), position = c("stacked"), hemisphere=lr, view=lm) +
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
p <- ggseg(.data = some_data, atlas = schaefer17_200, mapping=aes(fill=p, colour=I("white")), position = c("stacked")) +
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
effect_data <- gam.amp.schaefer200x17
some_data = tibble(
  region = schaefer200x17_roi_names,
  p = as.vector(effect_data$mean.curvature)
)

max_lim <- max(some_data$p)
min_lim <- min(some_data$p)
abs_max <- max(c(abs(max_lim), abs(min_lim)))
max_lim <- abs_max
min_lim <- -abs_max

# four map.
for (lr in c("left", "right")){
  for (lm in c("lateral", "medial")){
    file_path_name = sprintf("../results_bins/DRIVER/%s/%s/%s/all_mean_curv_%s_%s.jpeg", dataset, cs_name, gam_fig_filename, lr, lm)
    jpeg(file = file_path_name, width=1500, height=1000, units = "px",res =800)
    p <- ggseg(.data = some_data, atlas = schaefer17_200, mapping=aes(fill=p, colour=I("white")), position = c("stacked"), hemisphere=lr, view=lm) +
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
p <- ggseg(.data = some_data, atlas = schaefer17_200, mapping=aes(fill=p, colour=I("white")), position = c("stacked")) +
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


########################## partial R2 and second derivative ####################
# plot_data <- gam.amp.schaefer200x17 %>% select(GAM.age.partialR2, mean.derivative_2l)
# plot_data <- scale(plot_data)
# plot_data <- as.data.frame(plot_data)
# plot_data[3] <- yeo8
# colnames(plot_data)[1] <- "partialR2"
# colnames(plot_data)[2] <- "secondDerivative"
# colnames(plot_data)[3] <- "Group"
# plot_data$Group = factor(plot_data$Group,levels = c(1,2,3,4,5,6,7,8))
# scatterplot <- ggplot(data=plot_data, aes(x = partialR2, y = secondDerivative , color = Group)) +
#   geom_point(size = 0.2) +
#   # geom_hline(yintercept = 0, color = "gray", linetype = "solid", size = 0.3, alpha = 0.5) +
#   # geom_vline(xintercept = 0, color = "gray", linetype = "solid", size = 0.3, alpha = 0.5) +
#   scale_color_manual(values = color_arr)+  #设置填充的颜色
#   labs(x=element_blank(), y=element_blank()) +
#   theme_minimal() +
#   theme(legend.position = "none",
#         panel.grid = element_blank(),
#         # panel.border = element_line(color = "black", size = 0.2),
#         axis.line = element_line(color = "black", size = 0.1),
#         axis.text = element_text(size=8, family = "Arial", color = c("black")))
# p7 <- scatterplot + geom_hdr_lines(linewidth=0.3, probs = c(0.41, 0.4),
#                                    method = "kde", xlim=c(-4, 3), ylim=c(-3, 4))
# file_path_name = sprintf("../results_bins/DRIVER/%s/%s/%s/scatter.jpeg", dataset, cs_name, gam_fig_filename)
# jpeg(file = file_path_name, width = 300, height = 300, units = "px", res = 1500)
# print(p7)
# dev.off()


plot_data <- gam.amp.schaefer200x17 %>% select(GAM.age.partialR2, mean.derivative_2l)
plot_data <- scale(plot_data)
plot_data <- as.data.frame(plot_data)
plot_data[3] <- yeo8
colnames(plot_data)[1] <- "partialR2"
colnames(plot_data)[2] <- "secondDerivative"
colnames(plot_data)[3] <- "Group"
plot_data$Group = factor(plot_data$Group,levels = c(1,2,3,4,5,6,7,8))
scatterplot <- ggplot(data=plot_data, aes(x = partialR2, y = secondDerivative , color = Group)) +
  geom_point(size = 2.2) +
  # geom_hline(yintercept = 0, color = "gray", linetype = "solid", size = 0.3, alpha = 0.5) +
  # geom_vline(xintercept = 0, color = "gray", linetype = "solid", size = 0.3, alpha = 0.5) +
  scale_color_manual(values = color_arr)+  #设置填充的颜色
  labs(x=element_blank(), y=element_blank()) +
  theme_classic() + 
  theme(legend.position = "none", 
        panel.grid = element_blank(),
        axis.line = element_line(color = "black", size = 0.1),
        axis.text = element_text(size=17, family = "Arial", color = c("black")))
p7 <- scatterplot + geom_hdr_lines(linewidth=0.3, probs = c(0.4), 
                                   method = "kde", xlim=c(-4, 3), ylim=c(-3, 4))
file_path_name = sprintf("../results_bins/DRIVER/%s/%s/%s/scatter_partialR2_secondDerivative.emf", dataset, cs_name, gam_fig_filename)
emf(file_path_name, width=4, height=4, emfPlus=FALSE)
print(p7)
dev.off()


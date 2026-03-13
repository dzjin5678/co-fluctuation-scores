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

# ratioOfMeans, meanOfRatios, norm_zscore, norm_regress_beta, norm_regress_bAndBeta, leaveOneRegionOut
cs_name <- "1_region_level_ratioOfMeans"
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
gam_data_filename=sprintf("../results_bins/DRIVER/%s/%s/bins_gam_data_zscored_bins_rss_high_fd-%s_num_of_bins-%s_atlas-%s_regress_cofounds-%s_gsr-%s.csv", 
                          dataset, cs_name, high_fd, num_of_bins, atlas, regress_cofounds, pipeline);
gam_amp_filename=sprintf("../results_bins/mixedGAMs/%s/%s/bins_gam_amp_zscored_bins_rss_high_fd-%s_num_of_bins-%s_atlas-%s_regress_cofounds-%s_gsr-%s.csv", 
                         dataset, cs_name, high_fd, num_of_bins, atlas, regress_cofounds, pipeline);
gam_smooths_filename=sprintf("../results_bins/mixedGAMs/%s/%s/bins_gam_smooths_zscored_bins_rss_high_fd-%s_num_of_bins-%s_atlas-%s_regress_cofounds-%s_gsr-%s.csv", 
                             dataset, cs_name, high_fd, num_of_bins, atlas, regress_cofounds, pipeline);
gam_estimated_smooths_filename=sprintf("../results_bins/mixedGAMs/%s/%s/bins_gam_estimated_smooths_zscored_bins_rss_high_fd-%s_num_of_bins-%s_atlas-%s_regress_cofounds-%s_gsr-%s.csv", 
                                       dataset, cs_name, high_fd, num_of_bins, atlas, regress_cofounds, pipeline);
gam_fig_filename=sprintf("bins_gam_data_zscored_bins_rss_high_fd-%s_num_of_bins-%s_atlas-%s_regress_cofounds-%s_gsr-%s", 
                          high_fd, num_of_bins, atlas, regress_cofounds, pipeline);
save_dir <- sprintf("../results_bins/mixedGAMs/%s/%s/%s", dataset, cs_name, gam_fig_filename)
dir.create(save_dir)

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


# color_arr <- array(data = NA, dim = 8, dimnames = NULL)
# color_arr[1] <- rgb(128/255, 0/255, 128/255) # VIS
# color_arr[2] <- rgb(120/255, 154/255, 192/255) # SMN
# color_arr[3] <- rgb(64/255, 152/255, 50/255) # DAN
# color_arr[4] <- rgb(224/255, 102/255, 254/255) # SVAN
# color_arr[5] <- rgb(169/255, 169/255, 169/255) # LIM
# color_arr[6] <- rgb(238/255, 185/255, 67/255) # CONT
# color_arr[7] <- rgb(255/255, 0/255, 0/255) # DMN
# color_arr[8] <- rgb(0/255, 0/255, 128/255) # TP



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



################################ raw trajectories ##############################
# 89 "#80cdc1"
# 45 "#FBDC9D"
# plot single region.

# for(row in c(45:45)){
#   region <- sprintf("V%s", row)
#   V1.smooth <- gam.smooths.schaefer200x17 %>% filter(label == region)
#   file_path_name = sprintf("../../../results_bins/mixedGAMs/%s/region_level/ALL_SESSION_%s/%s_%s_%s_%s.jpeg", dataset, pipeline, session, atlas, pipeline, region)
#   jpeg(file = file_path_name,
#        width = 3000, height = 3000, units = "px", res = 2000)
#   gam_data <- gam.data %>% select(region, raw_rss)
#   colnames(gam_data)[1] <- "region_roi"
#   scatterplot <- ggplot(data=gam_data, aes(x = raw_rss, y = region_roi)) +
#     geom_point(color="#FBDC9D", size = 0.02) +
#     geom_line(data = V1.smooth, aes(x = raw_rss, y = fitted), color="black", size = 0.3) +
#     labs(x=element_blank(), y=element_blank()) +
#     theme_classic() +
#     theme(
#       axis.ticks=element_blank(),
#       axis.text = element_text(size=5, family = "Arial", color = c("black")),
#       axis.title.x=element_text(size=2, family ="Arial", color = c("darkgray")),
#       axis.title.y=element_text(size=2, family ="Arial", color = c("darkgray")),
#       axis.line = element_line(size=.12))
#   print(scatterplot)
#   dev.off()
# }

# for(row in c(45:45)){
#   region <- sprintf("V%s", row)
#   V1.smooth <- gam.smooths.schaefer200x17 %>% filter(label == region)
#   file_path_name = sprintf("../../../results_bins/mixedGAMs/%s/region_level/ALL_SESSION_%s/%s_%s_%s_%s.emf", dataset, pipeline, session, atlas, pipeline, region)
#   gam_data <- gam.data %>% select(region, raw_rss)
#   colnames(gam_data)[1] <- "region_roi"
#   scatterplot <- ggplot(data=gam_data, aes(x = raw_rss, y = region_roi)) +
#     geom_point(color="#FBDC9D", size = 1.5) +
#     geom_line(data = V1.smooth, aes(x = raw_rss, y = fitted), color="black", size = 1) +
#     labs(x=element_blank(), y=element_blank()) +
#     theme_classic() +
#     theme(
#       axis.text = element_text(size=17, family = "Arial", color = c("black")),
#       axis.line = element_line(size=0.22))
#   emf(file_path_name, width=4, height=4, emfPlus=FALSE) # 打开 emf 设备
#   print(scatterplot)                                           # 画图
#   dev.off()                                          # 关闭设备
# }



############################# fitted trajectories ##############################
# display trajectory by yeo8.
min_limit <- min(gam.smoothestimates.schaefer200x17$est);
max_limit <- max(gam.smoothestimates.schaefer200x17$est);

for(subnet_id in c(1:8)){
  subnet_data <- gam.smoothestimates.schaefer200x17[gam.smoothestimates.schaefer200x17$yeo8_label == subnet_id, ]
  file_path_name = sprintf("../results_bins/mixedGAMs/%s/%s/%s/%s.jpeg", dataset, cs_name, gam_fig_filename, as.character(subnet_id))
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
  file_path_name = sprintf("../results_bins/mixedGAMs/%s/%s/%s/%s.emf", dataset, cs_name, gam_fig_filename, as.character(subnet_id))
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


# display all trajectories colored by partial R2.
q_r2 <- quantile(gam.smoothestimates.schaefer200x17$GAM.age.partialR2, probs = c(0.25, 0.75), na.rm = FALSE)
df <- gam.smoothestimates.schaefer200x17
df_filtered <- df %>% filter(GAM.age.partialR2 <= q_r2[1] | GAM.age.partialR2 >= q_r2[2])

color_full <- gam.smoothestimates.schaefer200x17$GAM.age.partialR2
color <- df_filtered$GAM.age.partialR2

file_path_name = sprintf("../results_bins/mixedGAMs/%s/%s/%s/all_primary_effect.jpeg", dataset, cs_name, gam_fig_filename)

jpeg(file = file_path_name, width=1600, height=1500, units = "px",res =105)
ggplot(df_filtered, aes(zscored_bins_rss, est, group=index, color=color)) +
  geom_line(size=2.0, alpha = 0.8) +
  paletteer::scale_color_paletteer_c("ggthemes::Orange-Blue Diverging", direction = -1, limits = c(min(color_full), max(color_full))) +
  theme_classic() +
  ylab(NULL) +
  xlab(NULL) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size=70, family = "Arial", color = c("black")),
        axis.title = element_text(size=65, family = "Arial", color = c("black")),
        axis.line = element_line(size=.22), axis.ticks = element_line(size=.22)) +
  scale_y_continuous(limits = c(min_limit, max_limit), breaks = c(-0.02, 0, 0.02))
dev.off()


q_r2 <- quantile(gam.smoothestimates.schaefer200x17$GAM.age.partialR2, probs = c(0.25, 0.75), na.rm = FALSE)
df <- gam.smoothestimates.schaefer200x17
df_filtered <- df %>% filter(GAM.age.partialR2 <= q_r2[1] | GAM.age.partialR2 >= q_r2[2])
color_full <- gam.smoothestimates.schaefer200x17$GAM.age.partialR2
color <- df_filtered$GAM.age.partialR2
file_path_name = sprintf("../results_bins/mixedGAMs/%s/%s/%s/all_primary_effect.emf", dataset, cs_name, gam_fig_filename)
p <- ggplot(df_filtered, aes(zscored_bins_rss, est, group=index, color=color)) +
  geom_line(size=0.5, alpha = 0.8) +
  paletteer::scale_color_paletteer_c("ggthemes::Orange-Blue Diverging", direction = -1, limits = c(min(color_full), max(color_full))) +
  theme_classic() +
  ylab(NULL) +
  xlab(NULL) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size=17, family = "Arial", color = c("black")),
        axis.line = element_line(size=.22), axis.ticks = element_line(size=.22)) +
  scale_y_continuous(limits = c(min_limit, max_limit), breaks = c(-0.02, 0, 0.02))
emf(file_path_name, width=3, height=3, emfPlus=FALSE)
print(p)
dev.off()



# display all trajectories colored by second derivations

q_r2 <- quantile(gam.smoothestimates.schaefer200x17$mean.derivative_2l, probs = c(0.25, 0.75), na.rm = FALSE)
df <- gam.smoothestimates.schaefer200x17
df_filtered <- df %>% filter(mean.derivative_2l <= q_r2[1] | mean.derivative_2l >= q_r2[2])

color_full <- gam.smoothestimates.schaefer200x17$mean.derivative_2l
color <- df_filtered$mean.derivative_2l
file_path_name = sprintf("../results_bins/mixedGAMs/%s/%s/%s/all_second_effect.jpeg", dataset, cs_name, gam_fig_filename)
jpeg(file = file_path_name, width=1600, height=1500, units = "px",res =105)
ggplot(df_filtered, aes(zscored_bins_rss, est, group=index, color=color)) +
  geom_line(size=2, alpha = 0.8) +
  # paletteer::scale_color_paletteer_c("ggthemes::Orange-Blue-White Diverging", direction = 1, limits = c(min(color_full), max(color_full))) +
  paletteer::scale_color_paletteer_c("grDevices::Geyser", direction = 1, limits = c(min(color), max(color))) +
  theme_classic() +
  ylab(NULL) +
  xlab(NULL) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size=70, family = "Arial", color = c("black")),
        axis.title = element_text(size=65, family = "Arial", color = c("black")),
        axis.line = element_line(size=.22), axis.ticks = element_line(size=.22)) +
  scale_y_continuous(limits = c(min_limit, max_limit), breaks = c(-0.02, 0, 0.02))
dev.off()


q_r2 <- quantile(gam.smoothestimates.schaefer200x17$mean.derivative_2l, probs = c(0.25, 0.75), na.rm = FALSE)
df <- gam.smoothestimates.schaefer200x17
df_filtered <- df %>% filter(mean.derivative_2l <= q_r2[1] | mean.derivative_2l >= q_r2[2])
color_full <- gam.smoothestimates.schaefer200x17$mean.derivative_2l
color <- df_filtered$mean.derivative_2l
file_path_name = sprintf("../results_bins/mixedGAMs/%s/%s/%s/all_second_effect.emf", dataset, cs_name, gam_fig_filename)
p <- ggplot(df_filtered, aes(zscored_bins_rss, est, group=index, color=color)) +
  geom_line(size=0.5, alpha = 0.8) +
  paletteer::scale_color_paletteer_c("grDevices::Geyser", direction = 1, limits = c(min(color), max(color))) +
  theme_classic() +
  ylab(NULL) +
  xlab(NULL) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size=17, family = "Arial", color = c("black")),
        axis.line = element_line(size=.22), axis.ticks = element_line(size=.22)) +
  scale_y_continuous(limits = c(min_limit, max_limit), breaks = c(-0.02, 0, 0.02))
emf(file_path_name, width=3, height=3, emfPlus=FALSE)
print(p)
dev.off()


################################ amplitude effect ##############################
# effect_data <- read.csv(sprintf("../../../results_bins/mixedGAMs/%s/region_level/ALL_SESSION_%s/gam_amp_%s_%s_raw_rss.csv", 
#                                 dataset, pipeline, session, atlas))
effect_data <- gam.amp.schaefer200x17
some_data = tibble(
  region = schaefer200x17_roi_names,
  p = as.vector(effect_data$GAM.age.partialR2)
)
max_lim <- max(some_data$p)
min_lim <- min(some_data$p)

# four map.
for (lr in c("left", "right")){
  for (lm in c("lateral", "medial")){
    file_path_name = sprintf("../results_bins/mixedGAMs/%s/%s/%s/all_primary_effect_%s_%s.jpeg", dataset, cs_name, gam_fig_filename, lr, lm)
    jpeg(file = file_path_name, width=1500, height=1000, units = "px",res =800)
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

# single maps.
file_path_name = sprintf("../results_bins/mixedGAMs/%s/%s/%s/all_primary_effect_brain_map.jpeg", dataset, cs_name, gam_fig_filename)
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


# # plot one region (fig. 2a).
# one_data = tibble(
#   region = schaefer200x17_roi_names[c(45)],
#   p = as.vector(effect_data$GAM.age.partialR2[45])
# )
# max_lim <- max(some_data$p)
# min_lim <- min(some_data$p)
# jpeg(file = sprintf("../../../results_bins/mixedGAMs/%s/region_level/ALL_SESSION_%s/%s_amplitude_effect_%s_%s_zscored_bins_rss_V45.jpeg",
#                     dataset, pipeline, session, atlas, pipeline), width = 1500, height = 1500, units = "px",res = 800)
# ggseg(.data = one_data, atlas = schaefer17_200, mapping=aes(fill=p, colour=I("gray")), position = c("stacked"), hemisphere="left", view="medial") +
#   theme_void() +
#   labs(fill="") +
#   paletteer::scale_fill_paletteer_c("ggthemes::Orange-Blue Diverging", na.value="transparent", direction = -1,
#                                     limits = c(min_lim, max_lim), oob = squish) +
#   theme(legend.text = element_text(family = "Arial", size = 9, color = c("black")), legend.position = "no")
# dev.off()



############################# second derivations ###############################
second_derivatives_data <- gam.amp.schaefer200x17
some_data = tibble(
  region = schaefer200x17_roi_names,
  p = as.vector(second_derivatives_data$mean.derivative_2l)
)
max_lim <- max(some_data$p)
min_lim <- min(some_data$p)
# single map.
for (lr in c("left", "right")){
  for (lm in c("lateral", "medial")){
    file_path_name = sprintf("../results_bins/mixedGAMs/%s/%s/%s/all_second_effect_%s_%s.jpeg", dataset, cs_name, gam_fig_filename, lr, lm)
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
file_path_name = sprintf("../results_bins/mixedGAMs/%s/%s/%s/all_second_effect_brain_map.jpeg", dataset, cs_name, gam_fig_filename)
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


############## partial R2 and second derivatives in box plot ###################
# effect <- effect_data$GAM.age.partialR2
# Data1 = data.frame(Group = yeo8,Value = effect)
# Data1$Group = factor(Data1$Group,levels = c(1,2,3,4,5,6,7,8))
# P1 <- ggplot(Data1,aes(x=Group,y=Value,fill=Group))+ #”fill=“设置填充颜色
#   stat_boxplot(geom = "errorbar",width=0.15,aes(color="black"))+ #由于自带的箱形图没有胡须末端没有短横线，使用误差条的方式补上
#   geom_boxplot(size=0.2,fill="white",outlier.fill="white",outlier.color="white")+ #size设置箱线图的边框线和胡须的线宽度，fill设置填充颜色，outlier.fill和outlier.color设置异常点的属性
#   geom_jitter(aes(fill=Group), width=0.15, shape = 21, size=1.5, alpha=0.8, stroke=NA)+ #设置为向水平方向抖动的散点图，width指定了向水平方向抖动，不改变纵轴的值
#   scale_fill_manual(values = color_arr)+  #设置填充的颜色
#   scale_color_manual(values=c("black","black","black","black","black","black","black","black"))+ #设置散点图的圆圈的颜色为黑色
#   theme_bw()+ #背景变为白色
#   theme(legend.position="none", #不需要图例
#         axis.text.x=element_text(colour="black",family="Arial",size=11,face="plain", angle=90), #设置x轴刻度标签的字体属性
#         axis.text.y=element_text(colour="black",family="Arial",size=11,face="plain"), #设置x轴刻度标签的字体属性
#         axis.title.y=element_text(colour="black",family="Arial",size = 11,face="plain"), #设置y轴的标题的字体属性
#         axis.title.x=element_text(colour="black",family="Arial",size = 11,face="plain"), #设置x轴的标题的字体属性
#         # plot.title = element_text(family="Arial",size=11,face="bold",hjust = 0.5), #设置总标题的字体属性
#         panel.grid.major = element_blank(), #不显示网格线
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_line(colour="gray", linewidth=0.2),
#         axis.ticks.x = element_line(colour="gray", linewidth=0.1),
#         axis.ticks.y = element_line(colour="gray", linewidth=0.1))+
#   ylab(NULL) +
#   xlab(NULL) +
#   scale_x_discrete(labels = c("VIS", "SMN", "DAN", "SVAN", "LIMB", "CONT", "DMN", "TP"))
# P1
# jpeg(file = sprintf("../../../results_bins/mixedGAMs/%s/region_level/ALL_SESSION_%s/stat_amplitude_effect_%s_%s.jpg", dataset, pipeline, atlas, pipeline),
#      width =1200,height = 1200,units = "px",res =450)
# print(P1)
# dev.off()
# 
# scheafer200x17_mat <- readMat("../../../../atlas/Schaefer200x17.mat");
# effect <- effect_data$mean.derivative_2l
# yeo8 <- scheafer200x17_mat$lab17to8
# Data1 = data.frame(Group = yeo8,Value = effect)
# Data1$Group = factor(Data1$Group,levels = c(1,2,3,4,5,6,7,8))
# P1 <- ggplot(Data1,aes(x=Group,y=Value,fill=Group))+ #”fill=“设置填充颜色
#   stat_boxplot(geom = "errorbar",width=0.15,aes(color="black"))+ #由于自带的箱形图没有胡须末端没有短横线，使用误差条的方式补上
#   geom_boxplot(size=0.2,fill="white",outlier.fill="white",outlier.color="white")+ #size设置箱线图的边框线和胡须的线宽度，fill设置填充颜色，outlier.fill和outlier.color设置异常点的属性
#   geom_jitter(aes(fill=Group), width=0.15, shape = 21, size=1.5, alpha=0.8, stroke=NA)+ #设置为向水平方向抖动的散点图，width指定了向水平方向抖动，不改变纵轴的值
#   scale_fill_manual(values = color_arr)+  #设置填充的颜色
#   scale_color_manual(values=c("black","black","black","black","black","black","black","black"))+ #设置散点图的圆圈的颜色为黑色
#   theme_bw()+ #背景变为白色
#   theme(legend.position="none", #不需要图例
#         axis.text.x=element_text(colour="black",family="Arial",size=11,face="plain", angle=90), #设置x轴刻度标签的字体属性
#         axis.text.y=element_text(colour="black",family="Arial",size=11,face="plain"), #设置x轴刻度标签的字体属性
#         axis.title.y=element_text(colour="black",family="Arial",size = 11,face="plain"), #设置y轴的标题的字体属性
#         axis.title.x=element_text(colour="black",family="Arial",size = 11,face="plain"), #设置x轴的标题的字体属性
#         # plot.title = element_text(family="Arial",size=11,face="bold",hjust = 0.5), #设置总标题的字体属性
#         panel.grid.major = element_blank(), #不显示网格线
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_line(colour="gray", linewidth=0.2),
#         axis.ticks.x = element_line(colour="gray", linewidth=0.1),
#         axis.ticks.y = element_line(colour="gray", linewidth=0.1))+
#   ylab(NULL) +
#   xlab(NULL) +
#   scale_x_discrete(labels = c("VIS", "SMN", "DAN", "SVAN", "LIMB", "CONT", "DMN", "TP"))
# P1
# jpeg(file = sprintf("../../../results_bins/mixedGAMs/%s/region_level/ALL_SESSION_%s/stat_second_derivatives_%s_%s.jpg", dataset, pipeline, atlas, pipeline),
#      width =1200,height = 1200,units = "px",res =450)
# print(P1)
# dev.off()



########################## partial R2 and second derivative ####################
plot_data <- gam.amp.schaefer200x17 %>% select(GAM.age.partialR2, mean.derivative_2l)
plot_data <- scale(plot_data)
plot_data <- as.data.frame(plot_data)
plot_data[3] <- yeo8
colnames(plot_data)[1] <- "partialR2"
colnames(plot_data)[2] <- "secondDerivative"
colnames(plot_data)[3] <- "Group"
plot_data$Group = factor(plot_data$Group,levels = c(1,2,3,4,5,6,7,8))

scatterplot <- ggplot(data=plot_data, aes(x = partialR2, y = secondDerivative , color = Group)) +
  geom_point(size = 0.2) +
  # geom_hline(yintercept = 0, color = "gray", linetype = "solid", size = 0.3, alpha = 0.5) +
  # geom_vline(xintercept = 0, color = "gray", linetype = "solid", size = 0.3, alpha = 0.5) +
  scale_color_manual(values = color_arr)+  #设置填充的颜色
  labs(x=element_blank(), y=element_blank()) +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        # panel.border = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.1),
        axis.text = element_text(size=8, family = "Arial", color = c("black")))
p7 <- scatterplot + geom_hdr_lines(linewidth=0.3, probs = c(0.41, 0.4),
                                   method = "kde", xlim=c(-4, 3), ylim=c(-3, 4))
file_path_name = sprintf("../results_bins/mixedGAMs/%s/%s/%s/scatter_partialR2_secondDerivative.jpeg", dataset, cs_name, gam_fig_filename)
jpeg(file = file_path_name, width = 3000, height = 3000, units = "px", res = 1500)
print(p7)
dev.off()


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
file_path_name = sprintf("../results_bins/mixedGAMs/%s/%s/%s/scatter_partialR2_secondDerivative.emf", dataset, cs_name, gam_fig_filename)
emf(file_path_name, width=4, height=4, emfPlus=FALSE)
print(p7)
dev.off()


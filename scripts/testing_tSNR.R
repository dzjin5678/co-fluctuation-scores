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

schaefer200x17_atlas <- read.csv(sprintf("../../atlas/schaefer200x17_atlas.CSV"))
schaefer200x17_label <- read.csv("G:/datasets/atlas/Schaefer/Schaefer2018_200Parcels_17Networks_order_FSLMNI152_2mm.Centroid_RAS.csv")
schaefer200x17_roi_names <- schaefer200x17_label$ROI.Name

# load tSNR.
hcp3t_unrelated100_tSNR_mean <- readMat("../outputs_bins/HCP/tSNR/hcp3t_unrelated100_tSNR_mean.mat")
hcp3t_unrelated100_tSNR_mean$low.tSNR.roi.idx

# load primary amplitude effects.
gam.amp.schaefer200x17 <- read.csv("../results_bins/DRIVER/HCP/1_region_level_ratioOfMeans/bins_gam_amp_zscored_bins_rss_high_fd-0_num_of_bins-20_atlas-schaefer200_regress_cofounds-0_gsr-proc_regress_fix_ts_mat.csv")
# plot primary amplitude effects after remove 10 regions.
effect_data <- gam.amp.schaefer200x17
effect_data$GAM.age.partialR2[as.vector(hcp3t_unrelated100_tSNR_mean$low.tSNR.roi.idx)] <- NA
some_data = tibble(
  region = schaefer200x17_roi_names,
  p = as.vector(effect_data$GAM.age.partialR2)
)
max_lim <- max(some_data$p)
min_lim <- min(some_data$p)
# four map.
for (lr in c("left", "right")){
  for (lm in c("lateral", "medial")){
    file_path_name = sprintf("../results_bins/tSNR/all_primary_effect_%s_%s.jpeg", lr, lm)
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


# load myelination maps.
myelin_data <- readMat("../results_bins/mechanism/myelin/myelin_data.mat")
myelin_data <- myelin_data$myelin.map.schaefer200x17
myelin_data[as.vector(hcp3t_unrelated100_tSNR_mean$low.tSNR.roi.idx)] <- NA
atlas <- "schaefer200x17"
some_data = tibble(
  region = schaefer200x17_roi_names, 
  p = as.vector(myelin_data)
)
max_lim <- max(some_data$p)
min_lim <- min(some_data$p)
d <- some_data$p
d <- d[complete.cases(d)]
m <- median(d)
file_path_name = sprintf("../results_bins/tSNR/all_myelin.jpeg")
jpeg(file = file_path_name, width=100*15, height=70*15, units = "px",res =20*15)
ggseg(.data = some_data, atlas = schaefer17_200, mapping=aes(fill=p, colour=I("lightgray")), position = c("stacked")) +
  theme_void() +
  labs(fill="") +
  scale_fill_gradient2(low = "#4995c6", mid = "white", high = "#d51d08", guide = "colourbar", aesthetics = "fill", name = NULL, midpoint = 1.3) +
  theme(legend.text = element_text(family = "Arial", size = 15, color = c("black")), legend.position = "bottom")
dev.off()

# four map.
for (lr in c("left", "right")){
  for (lm in c("lateral", "medial")){
    file_path_name = sprintf("../results_bins/tSNR/all_myelin_%s_%s.jpeg", lr, lm)
    jpeg(file = file_path_name, width=1500, height=1000, units = "px",res =800)
    p <- ggseg(.data = some_data, atlas = schaefer17_200, mapping=aes(fill=p, colour=I("lightgray")), position = c("stacked"), hemisphere=lr, view=lm) +
      theme_void() +
      labs(fill="") +
      scale_fill_gradient2(low = "#4995c6", mid = "white", high = "#d51d08", guide = "colourbar", aesthetics = "fill", name = NULL, midpoint = 1.3) +
      theme(legend.text = element_text(family = "Arial", size = 20, color = c("black")), legend.position = "none")
    p
    print(p)
    dev.off()
  }
}


# load pvalb-sst maps.
pvalb_sst_data <- readMat("../results_bins/mechanism/pvalb_sst/pvalb_sst_data.mat")
pvalb_sst_data <- pvalb_sst_data$abha.expression.data.pvalb.sst.schaefer200x17
pvalb_sst_data[as.vector(hcp3t_unrelated100_tSNR_mean$low.tSNR.roi.idx)] <- NA
atlas <- "schaefer200x17"
some_data = tibble(
  region = schaefer200x17_roi_names, 
  p = as.vector(pvalb_sst_data)
)
max_lim <- max(some_data$p)
min_lim <- min(some_data$p)
d <- some_data$p
d <- d[complete.cases(d)]
m <- median(d)
file_path_name = sprintf("../results_bins/tSNR/all_pvalb_sst.jpeg")
jpeg(file = file_path_name, width=100*15, height=70*15, units = "px",res =20*15)
ggseg(.data = some_data, atlas = schaefer17_200, mapping=aes(fill=p, colour=I("lightgray")), position = c("stacked")) +
  theme_void() +
  labs(fill="") +
  scale_fill_gradient2(low = "#4995c6", mid = "white", high = "#d51d08", guide = "colourbar", aesthetics = "fill", name = NULL, midpoint = 0) +
  theme(legend.text = element_text(family = "Arial", size = 15, color = c("black")), legend.position = "bottom")
dev.off()

# four map.
for (lr in c("left", "right")){
  for (lm in c("lateral", "medial")){
    file_path_name = sprintf("../results_bins/tSNR/all_pvalb_sst_%s_%s.jpeg", lr, lm)
    jpeg(file = file_path_name, width=1500, height=1000, units = "px",res =800)
    p <- ggseg(.data = some_data, atlas = schaefer17_200, mapping=aes(fill=p, colour=I("lightgray")), position = c("stacked"), hemisphere=lr, view=lm) +
      theme_void() +
      labs(fill="") +
      scale_fill_gradient2(low = "#4995c6", mid = "white", high = "#d51d08", guide = "colourbar", aesthetics = "fill", name = NULL, midpoint = 0) +
      theme(legend.text = element_text(family = "Arial", size = 20, color = c("black")), legend.position = "none")
    p
    print(p)
    dev.off()
  }
}


# load pvalb-sst maps.
pvalb_sst_data <- readMat("../results_bins/mechanism/pvalb_sst/pvalb_sst_data.mat")
pvalb_sst_data <- pvalb_sst_data$abha.expression.data.pvalb.sst.schaefer200x17
atlas <- "schaefer200x17"
some_data = tibble(
  region = schaefer200x17_roi_names, 
  p = as.vector(pvalb_sst_data)
)
max_lim <- max(some_data$p)
min_lim <- min(some_data$p)
d <- some_data$p
d <- d[complete.cases(d)]
m <- median(d)
file_path_name = sprintf("../results_bins/mechanism/all_pvalb_sst.jpeg")
jpeg(file = file_path_name, width=100*15, height=70*15, units = "px",res =20*15)
ggseg(.data = some_data, atlas = schaefer17_200, mapping=aes(fill=p, colour=I("lightgray")), position = c("stacked")) +
  theme_void() +
  labs(fill="") +
  scale_fill_gradient2(low = "#4995c6", mid = "white", high = "#d51d08", guide = "colourbar", aesthetics = "fill", name = NULL, midpoint = 0) +
  theme(legend.text = element_text(family = "Arial", size = 15, color = c("black")), legend.position = "bottom")
dev.off()

# four map.
for (lr in c("left", "right")){
  for (lm in c("lateral", "medial")){
    file_path_name = sprintf("../results_bins/mechanism/all_pvalb_sst_%s_%s.jpeg", lr, lm)
    jpeg(file = file_path_name, width=1500, height=1000, units = "px",res =800)
    p <- ggseg(.data = some_data, atlas = schaefer17_200, mapping=aes(fill=p, colour=I("lightgray")), position = c("stacked"), hemisphere=lr, view=lm) +
      theme_void() +
      labs(fill="") +
      scale_fill_gradient2(low = "#4995c6", mid = "white", high = "#d51d08", guide = "colourbar", aesthetics = "fill", name = NULL, midpoint = 0) +
      theme(legend.text = element_text(family = "Arial", size = 20, color = c("black")), legend.position = "none")
    p
    print(p)
    dev.off()
  }
}


library(ggsegGlasser)
library(ggsegSchaefer)
library(ggseg)
library(dplyr)
library(ggseg3d)
library(R.matlab)
library(ggplot2)
library(cifti)
library(stringr)
library(factoextra)
library(matrixStats)
library(scales)
library(Hmisc)
library(tidyr)
library(cocor)


# atlas data.
schaefer200x17_label <- read.csv("G:/datasets/atlas/Schaefer/Schaefer2018_200Parcels_17Networks_order_FSLMNI152_2mm.Centroid_RAS.csv")
schaefer200x17_roi_names <- schaefer200x17_label$ROI.Name

# glasser180_label <- read.csv("G:/datasets/atlas/HCPMMP/glasser_regions.csv")
# glasser360_label <- rbind(glasser180_label, glasser180_label)
# glasser360_roi_names <- glasser360_label$region


# # timescale.
# timescale_data <- readMat("D:/matlab_proj/MyMatlab/ets/results_bins/mechanism/timescale/neuromaps_meg_timescale_glasser.mat")
# timescale_datas <- timescale_data[["meg.timescale.glasser"]]
# timescale_datas_vec <- as.vector(timescale_datas)
# 
# some_data = tibble(region = glasser360_roi_names, p = timescale_datas_vec)
# save_path_name <- sprintf("D:/matlab_proj/MyMatlab/ets/results_bins/mechanism/timescale/timescale_glasser.jpeg")
# plot_map_glasser(some_data, save_path_name)




# 利用功能解释功能，有循环验证的嫌疑
# timescale_data <- readMat("D:/matlab_proj/MyMatlab/ets/results_bins/mechanism/timescale/neuromaps_meg_timescale_schaefer200x17.mat")
# timescale_datas <- timescale_data
# timescale_datas_vec <- as.vector(timescale_datas$meg.timescale.schaefer200x17)# myelin.
# atlas <- "schaefer200x17"
# some_data = tibble(
#   region = schaefer200x17_roi_names, 
#   p = timescale_datas_vec
# )
# max_lim <- max(some_data$p)
# min_lim <- min(some_data$p)
# jpeg(file = sprintf("D:/matlab_proj/MyMatlab/ets/results_bins/mechanism/timescale/timescale_%s.jpeg", atlas), width=100*15, height=70*15, units = "px",res =20*15)
# ggseg(.data = some_data, atlas = schaefer17_200, mapping=aes(fill=p, colour=I("gray")), position = c("stacked")) +
#   theme_void() +
#   labs(fill="") +
#   paletteer::scale_fill_paletteer_c("pals::ocean.matter", direction = -1, limits = c(min_lim, max_lim), oob = squish) +
#   theme(legend.text = element_text(family = "Arial", size = 15, color = c("black")))
# dev.off()




# timescale_data <- readMat("D:/matlab_proj/MyMatlab/ets/results_bins/mechanism/timescale/neuromaps_meg_timescale_glasser.mat")
# timescale_datas <- timescale_data[["meg.timescale.glasser"]]
# timescale_datas_vec <- as.vector(timescale_datas)# myelin.
# data <- read.csv("D:/matlab_proj/MyMatlab/ets/testing/contribution/all_map_HCP_glasser.csv", header = TRUE)
# atlas <- "glasser"
# some_data = tibble(
#   region = glasser360_roi_names, 
#   p = timescale_datas_vec
# )
# max_lim <- max(some_data$p)
# min_lim <- min(some_data$p)
# jpeg(file = sprintf("D:/matlab_proj/MyMatlab/ets/results_bins/mechanism/timescale/timescale_%s.jpeg", atlas), width=100*15, height=70*15, units = "px",res =20*15)
# ggseg(.data = some_data, atlas = glasser, mapping=aes(fill=p, colour=I("gray")), position = c("stacked")) +
#   theme_void() +
#   labs(fill="") +
#   paletteer::scale_fill_paletteer_c("pals::ocean.matter", direction = -1, limits = c(min_lim, max_lim), oob = squish) +
#   theme(legend.text = element_text(family = "Arial", size = 15, color = c("black")))
# dev.off()


lr <- "right" 
# "right", "left"
ml <- "lateral" 
# "medial", "lateral"
data <- read.csv("D:/matlab_proj/MyMatlab/ets/testing/mechanism/s_a_axis_shaefer200x17.csv", header = FALSE)
# S-A axis
atlas <- "schaefer200x17"
some_data = tibble(
  region = schaefer200x17_roi_names, 
  p = data$V1
)
max_lim <- max(some_data$p)
min_lim <- min(some_data$p)
jpeg(file = sprintf("D:/matlab_proj/MyMatlab/ets/results_bins/mechanism/sa_axis_%s_%s.jpeg", lr, ml), width=100*15, height=70*15, units = "px",res =20*15)
ggseg(.data = some_data, atlas = schaefer17_200, mapping=aes(fill=p, colour=I("lightgray")), 
      position = c("stacked"), hemisphere = c(lr) , view = c(ml)) +
  theme_void() +
  labs(fill="") +
  scale_fill_gradient2(low = "goldenrod1", mid = "white", high = "#6f1282", guide = "colourbar", aesthetics = "fill", name = NULL, midpoint = 0) +
  theme(legend.text = element_text(family = "Arial", size = 15, color = c("black")), legend.position = "bottom")
dev.off()


jpeg(file = sprintf("D:/matlab_proj/MyMatlab/ets/results_bins/mechanism/sa_axis_new.jpeg"), width=100*15, height=70*15, units = "px",res =20*15)
ggseg(.data = some_data, atlas = schaefer17_200, mapping=aes(fill=p, colour=I("lightgray")), 
      position = c("stacked")) +
  theme_void() +
  labs(fill="") +
  scale_fill_gradient2(low = "goldenrod1", mid = "white", high = "#6f1282", guide = "colourbar", aesthetics = "fill", name = NULL, midpoint = 0) +
  theme(legend.text = element_text(family = "Arial", size = 15, color = c("black")), legend.position = "bottom")
dev.off()


data <- read.csv("D:/matlab_proj/MyMatlab/ets/testing/contribution/all_map_HCP_schaefer200x17.csv", header = TRUE)
data <- readMat("D:/matlab_proj/MyMatlab/ets/results_bins/mechanism/timescale/neuromaps_myeline_schaefer200x17.mat")
data <- data[['myelin.schaefer200x17']]
# myelin.
# p = as.vector(data$myelin_map)
atlas <- "schaefer200x17"
some_data = tibble(
  region = schaefer200x17_roi_names, 
  p = as.vector(data)
)
max_lim <- max(some_data$p)
min_lim <- min(some_data$p)
d <- some_data$p
d <- d[complete.cases(d)]
m <- median(d)
jpeg(file = sprintf("D:/matlab_proj/MyMatlab/ets/results_bins/mechanism/myelin/myelin_%s.jpeg", atlas), width=100*15, height=70*15, units = "px",res =20*15)
ggseg(.data = some_data, atlas = schaefer17_200, mapping=aes(fill=p, colour=I("lightgray")), position = c("stacked")) +
  theme_void() +
  labs(fill="") +
  scale_fill_gradient2(low = "#4995c6", mid = "white", high = "#d51d08", guide = "colourbar", aesthetics = "fill", name = NULL, midpoint = 1.2) +
  theme(legend.text = element_text(family = "Arial", size = 15, color = c("black")), legend.position = "bottom")
dev.off()


data <- read.csv("D:/matlab_proj/MyMatlab/ets/testing/contribution/all_map_HCP_glasser.csv", header = TRUE)
# myelin.
atlas <- "glasser"
some_data = tibble(
  region = glasser360_roi_names, 
  p = as.vector(data$myelin_map)
)
max_lim <- max(some_data$p)
min_lim <- min(some_data$p)
jpeg(file = sprintf("D:/matlab_proj/MyMatlab/ets/results_bins/mechanism/myelin/myelin_%s.jpeg", atlas), width=100*15, height=70*15, units = "px",res =20*15)
ggseg(.data = some_data, atlas = glasser, mapping=aes(fill=p), position = c("stacked")) +
  theme_void() +
  labs(fill="") +
  paletteer::scale_fill_paletteer_c("pals::ocean.matter", direction = -1, limits = c(min_lim, max_lim), oob = squish) +
  theme(legend.text = element_text(family = "Arial", size = 15, color = c("black")), legend.position = "bottom")
dev.off()




data <- read.csv("D:/matlab_proj/MyMatlab/ets/testing/contribution/all_map_HCP_schaefer200x17.csv", header = TRUE)
# gene pc1.
atlas <- "schaefer200x17"
some_data = tibble(
  region = schaefer200x17_roi_names, 
  p = as.vector(data$pc1)
)
max_lim <- max(some_data$p)
min_lim <- min(some_data$p)
d <- some_data$p
d <- d[complete.cases(d)]
m <- mean(d)
jpeg(file = sprintf("D:/matlab_proj/MyMatlab/ets/results_bins/mechanism/gene/gene_pc1_%s.jpeg", atlas), width=100*15, height=70*15, units = "px",res =20*15)
ggseg(.data = some_data, atlas = schaefer17_200, mapping=aes(fill=p, colour=I("white")), position = c("stacked")) +
  theme_void() +
  labs(fill="") +
  scale_fill_gradient2(low = "#4995c6", mid = "white", high = "#ff1d08", guide = "colourbar", aesthetics = "fill", name = NULL, midpoint = m) +
  theme(legend.text = element_text(family = "Arial", size = 15, color = c("black")), legend.position = "bottom")
dev.off()




data <- read.csv("D:/matlab_proj/MyMatlab/ets/testing/contribution/all_map_HCP_schaefer200x17.csv", header = TRUE)
# gene pvalb-sst.
atlas <- "schaefer200x17"
some_data = tibble(
  region = schaefer200x17_roi_names, 
  p = as.vector(data$abha_expression_data_pvalb_sst)
)
d <- some_data$p
d <- d[complete.cases(d)]
m <- mean(d)
jpeg(file = sprintf("D:/matlab_proj/MyMatlab/ets/results_bins/mechanism/gene/pvalb_sst_difference_%s.jpeg", atlas), width=100*15, height=70*15, units = "px",res =20*15)
ggseg(.data = some_data, atlas = schaefer17_200, mapping=aes(fill=p, colour=I("white")), position = c("stacked")) +
  theme_void() +
  labs(fill="") +
  scale_fill_gradient2(low = "#4995c6", mid = "white", high = "#ff1d08", guide = "colourbar", aesthetics = "fill", name = NULL, midpoint = 0) +
  theme(legend.text = element_text(family = "Arial", size = 15, color = c("black")), legend.position = "bottom")
dev.off()




data <- read.csv("D:/matlab_proj/MyMatlab/ets/testing/contribution/all_map_HCP_glasser.csv", header = TRUE)
# gene pc1.
atlas <- "glaser"
some_data = tibble(
  region = glasser360_roi_names, 
  p = as.vector(data$pc1)
)
max_lim <- max(some_data$p)
min_lim <- min(some_data$p)
jpeg(file = sprintf("D:/matlab_proj/MyMatlab/ets/results_bins/mechanism/gene/gene_pc1_%s.jpeg", atlas), width=100*15, height=70*15, units = "px",res =20*15)
ggseg(.data = some_data, atlas = glasser, mapping=aes(fill=p), position = c("stacked")) +
  theme_void() +
  labs(fill="") +
  paletteer::scale_fill_paletteer_c("pals::ocean.matter", direction = -1, limits = c(min_lim, max_lim), oob = squish) +
  theme(legend.text = element_text(family = "Arial", size = 15, color = c("black")), legend.position = "bottom")
dev.off()




data <- read.csv("D:/matlab_proj/MyMatlab/ets/testing/contribution/all_map_HCP_glasser.csv", header = TRUE)
# gene pvalb-sst.
atlas <- "glasser"
some_data = tibble(
  region = glasser360_roi_names, 
  p = as.vector(data$abha_expression_data_pvalb_sst)
)
max_lim <- max(some_data$p)
min_lim <- min(some_data$p)
jpeg(file = sprintf("D:/matlab_proj/MyMatlab/ets/results_bins/mechanism/gene/pvalb_sst_difference_%s.jpeg", atlas), width=100*15, height=70*15, units = "px",res =20*15)
ggseg(.data = some_data, atlas = glasser, mapping=aes(fill=p), position = c("stacked")) +
  theme_void() +
  labs(fill="") +
  paletteer::scale_fill_paletteer_c("pals::ocean.matter", direction = -1, limits = c(min_lim, max_lim), oob = squish) +
  theme(legend.text = element_text(family = "Arial", size = 15, color = c("black")), legend.position = "bottom")
dev.off()


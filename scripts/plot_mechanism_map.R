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





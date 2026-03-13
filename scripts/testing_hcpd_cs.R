library(dplyr)
library(cifti)
library(ggplot2)
library(ggExtra)
library(R.matlab)
library(corrplot);
library(ComplexHeatmap);
library(colorRamp2);
library(ggsegGlasser)
library(ggsegSchaefer)
library(ggseg)
library(dplyr)
library(ggseg3d)
library(stringr)
library(factoextra)
library(matrixStats)
library(scales)
library(Hmisc)
library(tidyr)
library(cocor)
library(oompaBase)
library(ggthemes)
library(devEMF)


source("./dev/GAM_functions.R")


## work 1, figure development.

root_path<-"D:/matlab_proj/MyMatlab/ets/results_bins/DRIVER/HCPD/region_level/bins_20_age_effect_lme";

bins_top_middle_bottom_sa_sim_individual = readMat("../../outputs_bins/HCPD/bins_20_top_middle_bottom_sa_sim_individual_spearman.mat");
ages <- array(bins_top_middle_bottom_sa_sim_individual$age.6.18);
sexs <- array(bins_top_middle_bottom_sa_sim_individual$sex.6.18);
bins_1_2 <- array(bins_top_middle_bottom_sa_sim_individual$bin.1.2.sims.all);
bins_9_12 <- array(bins_top_middle_bottom_sa_sim_individual$bin.9.12.sims.all);
bins_19_20 <- array(bins_top_middle_bottom_sa_sim_individual$bin.19.20.sims.all);

gam_data_f_bins_1_2 <- as.data.frame(cbind(ages, bins_1_2, sexs));
gam_data_f_bins_9_12 <- as.data.frame(cbind(ages, bins_9_12, sexs));
gam_data_f_bins_19_20 <- as.data.frame(cbind(ages, bins_19_20, sexs));


# GAM.RESULTS <- gam.fit.smooth(gam_data_f, measure = "bins_1_2", smooth_var = "ages", covariates = "sexs",
#                               knots = 3, set_fx = TRUE, stats_only = FALSE) #run the gam.fit.smooth function
# gam.estimated.smooths.schaefer200x7 <- matrix(data=NA, ncol=4) 
# colnames(gam.estimated.smooths.schaefer200x7) <- c("ages","est","index","label")


## bin 1, trajectory.
GAM.SMOOTH <- gam.smooth.predict(gam_data_f_bins_1_2, measure = "bins_1_2", smooth_var = "ages", covariates = "sexs",
                                 knots = 3, set_fx = TRUE, increments = 1000) #run the gam.smooth.predict function
preddata <- as.data.frame(GAM.SMOOTH[2]) #get predicted.smooth df from function output

# jpeg(file = sprintf("%s/bins_1_2_spearman_new.jpeg", root_path), width = 3700, height = 3000, units = "px", res = 1500)
p <- ggplot(data=gam_data_f_bins_1_2, aes(x = ages, y = bins_1_2)) +
  geom_point(color="#bdbdbd", size = 1) +
  geom_line(data = preddata, aes(x = ages, y = fitted), color="darkred", size = 0.7) +
  labs(x="", y="") +
  theme_classic() +
  theme(
    axis.ticks=element_blank(),
    axis.text = element_text(colour="black",family="Arial",size=20, color = c("black")),
    axis.title.x=element_text(colour="black",family="Arial",size=20, color = c("black")),
    axis.title.y=element_text(colour="black",family="Arial",size=20, color = c("black")),
    axis.line = element_line(size=0.3)) +
  ylim(-0.8, 0.8) + 
  scale_x_continuous(breaks=c(6, 10, 15, 18), labels = c(6, 10, 15, 18), expand = c(0,.15))
  # scale_y_continuous(breaks=c(-0.8, 0, 0.8), labels = c(-0.8, 0, 0.8))
# print(p)
# dev.off()
file_path_name = sprintf("%s/bins_1_2_spearman_new.emf", root_path)
emf(file_path_name, width=3, height=3, emfPlus=FALSE)
print(p)
dev.off()



## bin 2, trajectory.
GAM.SMOOTH <- gam.smooth.predict(gam_data_f_bins_9_12, measure = "bins_9_12", smooth_var = "ages", covariates = "sexs",
                                 knots = 3, set_fx = TRUE, increments = 200) #run the gam.smooth.predict function
preddata <- as.data.frame(GAM.SMOOTH[2]) #get predicted.smooth df from function output
# jpeg(file = sprintf("%s/bins_9_12_spearman_new.jpeg", root_path), width = 3700, height = 3000, units = "px", res = 1500)
p <- ggplot(data=gam_data_f_bins_9_12, aes(x = ages, y = bins_9_12)) +
  geom_point(color="#bdbdbd", size = 1) +
  geom_line(data = preddata, aes(x = ages, y = fitted), color="darkred", size = 0.7) +
  labs(x="", y="") +
  theme_classic() +
  theme(
    axis.ticks=element_blank(),
    axis.text = element_text(colour="black",size=20, family = "Arial", color = c("black")),
    axis.title.x=element_text(colour="black",size=20, family ="Arial", color = c("black")),
    axis.title.y=element_text(colour="black",size=20, family ="Arial", color = c("black")),
    axis.line = element_line(size=0.3)) +
  ylim(-0.8, 0.8) + 
  scale_x_continuous(breaks=c(6, 10, 15, 18), labels = c(6, 10, 15, 18), expand = c(0,.15))
# print(p)
# dev.off()
file_path_name = sprintf("%s/bins_9_12_spearman_new.emf", root_path)
emf(file_path_name, width=3, height=3, emfPlus=FALSE)
print(p)
dev.off()


# ## bin 3, trajectory.
# GAM.SMOOTH <- gam.smooth.predict(gam_data_f_bins_19_20, measure = "bins_19_20", smooth_var = "ages", covariates = "sexs",
#                                  knots = 3, set_fx = TRUE, increments = 200) #run the gam.smooth.predict function
# preddata <- as.data.frame(GAM.SMOOTH[2]) #get predicted.smooth df from function output
# jpeg(file = sprintf("%s/bins_19_20_spearman_new.jpeg", root_path), width = 3700, height = 3000, units = "px", res = 1500)
# p <- ggplot(data=gam_data_f_bins_19_20, aes(x = ages, y = bins_19_20)) +
#   geom_point(color="#bdbdbd", size = 0.1) +
#   geom_line(data = preddata, aes(x = ages, y = fitted), color="darkred", size = 0.7) +
#   labs(x="Age", y="Similarity") +
#   theme_classic() +
#   theme(
#     axis.ticks=element_blank(),
#     axis.text = element_text(colour="black",size=11, family = "Arial", color = c("black")),
#     axis.title.x=element_text(colour="black",size=11, family ="Arial", color = c("black")),
#     axis.title.y=element_text(colour="black",size=11, family ="Arial", color = c("black")),
#     axis.line = element_line(size=0.3)) +
#   ylim(-0.8, 0.8) + 
#   scale_x_continuous(breaks=c(6, 10, 15, 18), labels = c(6, 10, 15, 18), expand = c(0,.15))
# print(p)
# dev.off()






schaefer200x17_label <- read.csv("G:/datasets/atlas/Schaefer/Schaefer2018_200Parcels_17Networks_order_FSLMNI152_2mm.Centroid_RAS.csv")
schaefer200x17_roi_names <- schaefer200x17_label$ROI.Name

## plot age effects.
age_effects <- readMat(sprintf("%s/combins_cofluctuation_score_age_effect.mat", root_path))

bin_index <- 1;
t_values <- age_effects$corr.t.combin.bin.roi;
t_values_top <- t_values[bin_index, 1:200];
p_values <- age_effects$corr.lme.p.combin.bin.roi;
p_values_top <- p_values[bin_index, 1:200];
t_values_top[p_values_top>0.001/200] <- NaN
some_data = tibble(
  region = schaefer200x17_roi_names,
  p = as.vector(t_values_top)
)
# max_lim <- max(some_data$p)
# min_lim <- min(some_data$p)
min_lim <- -15
max_lim <- 15

for (lr in c("left", "right")){
  for (lm in c("lateral", "medial")){
    jpeg(file = sprintf("%s/combins_cofluctuation_score_age_effect_bin_%i_%s_%s.png", root_path, bin_index, lr, lm),
         width=100*15, height=70*15, units = "px",res=600)
    p <- ggseg(.data = some_data, atlas = schaefer17_200, mapping=aes(fill=p, colour=I("lightgray")), position = c("stacked"), hemisphere=lr, view=lm) +
      theme_void() +
      labs(fill="") +
      paletteer::scale_fill_paletteer_c("pals::coolwarm", direction = 1, na.value="transparent", 
                                        limits = c(min_lim, max_lim), oob = squish, breaks=c(-15, 0, 15), labels=c("-15", "0", "15")) +
      theme(legend.text = element_text(family = "Arial", size = 20, color = c("black")), legend.position = "none")
    p
    print(p)
    dev.off()
  }
}
jpeg(file = sprintf("%s/combins_cofluctuation_score_age_effect_bin_%i.png", root_path, bin_index),
  width=100*20, height=70*20, units = "px",res =20*20)
p <- ggseg(.data = some_data, atlas = schaefer17_200, mapping=aes(fill=p, colour=I("lightgray")), position = c("stacked")) +
  theme_void() +
  labs(fill="") +
  paletteer::scale_fill_paletteer_c("pals::coolwarm", direction = 1, na.value="transparent", 
                                    limits = c(min_lim, max_lim), oob = squish, breaks=c(-15, 0, 15), labels=c("-15", "0", "15")) +
  theme(legend.text = element_text(family = "Arial", size = 20, color = c("black")), legend.position = "bottom")
p
print(p)
dev.off()




## association between age effect and SA rank.
age_effects <- readMat(sprintf("%s/combins_cofluctuation_score_age_effect.mat", root_path))
t_values <- age_effects$corr.t.combin.bin.roi;
bin_index <- 1;
t_values_top <- t_values[bin_index, 1:200];
age_effect <- t_values_top;
# p_values <- age_effects$corr.lme.p.combin.bin.roi;
# p_values_top <- p_values[1, 1:200];
# t_values_top[p_values_top>0.001/200] <- NaN
SAranks <- read.csv("D:/matlab_proj/MyMatlab/ets/testing/mechanism/s_a_axis_shaefer200x17.csv", header = F)

some_data = tibble(
  age_effect = age_effect,
  SAranks = SAranks$V1
)

# jpeg(file = sprintf("%s/combins_cofluctuation_score_age_effect_SAranks_bin_%i.png", root_path, bin_index),
#   width=100*20, height=100*20, units = "px",res =20*20)
p <- ggplot(some_data, aes(x = SAranks, y = age_effect, fill = SAranks)) +
  geom_point(aes(color = SAranks), shape = 21, size = 2.2) +
  scale_fill_gradient2(low = "goldenrod1", mid = "white", high = "#6f1282", guide = "colourbar", aesthetics = "fill", name = NULL, midpoint = 0) +
  scale_fill_gradient2(low = "goldenrod1", mid = "white", high = "#6f1282", guide = "colourbar", aesthetics = "color", name = NULL, midpoint = 0) +
  geom_smooth(method = 'lm', se = TRUE, fill = alpha(c("lightgray"), 1), col = "black", size = .25) +
  theme_classic() +
  ylab(NULL) +
  xlab(NULL) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size = 19, family = "Arial", color = c("black")), 
        axis.title = element_text(size = 19, family = "Arial", color = c("black")), 
        axis.line = element_line(size = .22), axis.ticks = element_line(size = .22)) + 
  scale_y_continuous(limits = c(-15, 10), breaks = c(-15, -10, -5, 0, 5, 10))
  
# p
# print(p)
# dev.off()
file_path_name=sprintf("%s/combins_cofluctuation_score_age_effect_SAranks_bin_%i.emf", root_path, bin_index)
emf(file_path_name, width=3.5, height=3.5, emfPlus=FALSE)
print(p)
dev.off()


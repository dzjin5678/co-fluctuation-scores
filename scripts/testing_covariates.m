new_gam_model = readtable("../results_bins/DRIVER/HCP/1_region_level_ratioOfMeans/" + ...
    "bins_gam_amp_zscored_bins_rss_high_fd-0_num_of_bins-20_atlas-schaefer200_regress_cofounds-0_gsr-proc_regress_fix_ts_mat.csv");
old_rest1lr_gam_model = readtable("../results_bins/DRIVER/HCP/region_level/gam_amp_REST1_LR_schaefer200x17_icafix_ts_raw_rss.csv");
old_rest1rl_gam_model = readtable("../results_bins/DRIVER/HCP/region_level/gam_amp_REST1_RL_schaefer200x17_icafix_ts_raw_rss.csv");
old_rest2lr_gam_model = readtable("../results_bins/DRIVER/HCP/region_level/gam_amp_REST2_LR_schaefer200x17_icafix_ts_raw_rss.csv");
old_rest2rl_gam_model = readtable("../results_bins/DRIVER/HCP/region_level/gam_amp_REST2_RL_schaefer200x17_icafix_ts_raw_rss.csv");


%% partial R2.
new_partial_r2 = new_gam_model.GAM_age_partialR2;
old_partial_r2 = (old_rest1lr_gam_model.GAM_age_partialR2 + ...
    old_rest1rl_gam_model.GAM_age_partialR2 + ...
    old_rest2lr_gam_model.GAM_age_partialR2 + ...
    old_rest2rl_gam_model.GAM_age_partialR2)/4;
[coef, pval] = corr(new_partial_r2, old_partial_r2, "type","Spearman");

x = zscore(new_partial_r2);
y = zscore(old_partial_r2);
[coef, pval] = corr(x, y, "type","Spearman");
figure("Position", [100 100 260 260]);
scatter(x, y, 40, 'filled', 'MarkerFaceAlpha',0.4);
hold on
pfit = polyfit(x,y,1);
x_fit = linspace(min(x),max(x),100);
y_fit = polyval(pfit,x_fit);
plot(x_fit,y_fit,'r','LineWidth',2);
ax=gca;hold on;
ax.FontSize=20;
set(gca,'Box','off');
save_filename=strcat("hcp3t_corr_new_old_partial_R2");
save_path=strcat("../results_bins/covariates");
print(gcf, '-dpng', '-r1000', strcat(save_path, "/", save_filename));
print(gcf, '-dmeta', strcat(save_path, "/", save_filename));
close;


%% mean second derivative.
new_mean_2th_derivative = new_gam_model.mean_derivative_2l;
old_mean_2th_derivative = (old_rest1lr_gam_model.mean_derivative_2l + ...
    old_rest1rl_gam_model.mean_derivative_2l + ...
    old_rest2lr_gam_model.mean_derivative_2l + ...
    old_rest2rl_gam_model.mean_derivative_2l)/4;
[coef, pval] = corr(new_mean_2th_derivative, old_mean_2th_derivative, "type","Spearman");


x = zscore(new_mean_2th_derivative);
y = zscore(old_mean_2th_derivative);
[coef, pval] = corr(x, y, "type","Spearman");
figure("Position", [100 100 260 260]);
scatter(x, y, 40, 'filled', 'MarkerFaceAlpha',0.4);
hold on
pfit = polyfit(x,y,1);
x_fit = linspace(min(x),max(x),100);
y_fit = polyval(pfit,x_fit);
plot(x_fit,y_fit,'r','LineWidth',2);
ax=gca;hold on;
ax.FontSize=20;
set(gca,'Box','off');
save_filename=strcat("hcp3t_corr_new_old_second_effect");
save_path=strcat("../results_bins/covariates");
print(gcf, '-dpng', '-r1000', strcat(save_path, "/", save_filename));
print(gcf, '-dmeta', strcat(save_path, "/", save_filename));
close;
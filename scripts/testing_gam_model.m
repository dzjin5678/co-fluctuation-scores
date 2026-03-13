%% Remove the lowest 5% tSNR regions
atlas="schaefer200";
save_dir=strcat("../results_bins/tSNR");
load(strcat(save_dir, strcat("/hcp3t_unrelated100_tSNR_", atlas)), "tSNR_all", "tSNR_mean", "low_tSNR_roi_idx");
roi_idx = 1:200;
roi_idx(low_tSNR_roi_idx)=[];


%% neurobiological.
new_gam_model = readtable("../results_bins/DRIVER/HCP/1_region_level_ratioOfMeans/" + ...
    "bins_gam_amp_zscored_bins_rss_high_fd-0_num_of_bins-20_atlas-schaefer200_regress_cofounds-0_gsr-proc_regress_fix_ts_mat.csv");
myelin_data = load("../results_bins/mechanism/myelin/myelin_data.mat");
pvalb_sst_data = load("../results_bins/mechanism/pvalb_sst/pvalb_sst_data.mat");

mean_derivative_2l = new_gam_model.mean_derivative_2l(roi_idx);
mean_curvature = new_gam_model.mean_curvature(roi_idx);
myelin = myelin_data.myelin_map_schaefer200x17(roi_idx);
pvalb_sst = pvalb_sst_data.abha_expression_data_pvalb_sst_schaefer200x17(roi_idx);

mean_derivative_2l_myelin = corr(mean_derivative_2l, myelin, ...
    "type","Spearman", 'Rows','pairwise');
mean_derivative_2l_pvalb_sst = corr(mean_derivative_2l, pvalb_sst, ...
    "type","Spearman", 'Rows','pairwise');

mean_curvature_myelin = corr(mean_curvature, myelin, ...
    "type","Spearman", 'Rows','pairwise');
mean_curvature_pvalb_sst = corr(mean_curvature, pvalb_sst, ...
    "type","Spearman", 'Rows','pairwise');


y = zscore(mean_derivative_2l(~isnan(myelin)));
x = zscore(myelin(~isnan(myelin)));
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
save_filename=strcat("hcp3t_corr_mean_derivative_2l_myelin");
save_path=strcat("../results_bins/tSNR");
print(gcf, '-dpng', '-r1000', strcat(save_path, "/", save_filename));
print(gcf, '-dmeta', strcat(save_path, "/", save_filename));
close;


y = zscore(mean_curvature(~isnan(myelin)));
x = zscore(myelin(~isnan(myelin)));
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
save_filename=strcat("hcp3t_corr_mean_curvature_myelin");
save_path=strcat("../results_bins/tSNR");
print(gcf, '-dpng', '-r1000', strcat(save_path, "/", save_filename));
print(gcf, '-dmeta', strcat(save_path, "/", save_filename));
close;

y = zscore(mean_curvature(~isnan(pvalb_sst)));
x = zscore(pvalb_sst(~isnan(pvalb_sst)));
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
save_filename=strcat("hcp3t_corr_mean_curvature_pvalb_sst");
save_path=strcat("../results_bins/tSNR");
print(gcf, '-dpng', '-r1000', strcat(save_path, "/", save_filename));
print(gcf, '-dmeta', strcat(save_path, "/", save_filename));
close;

y = zscore(mean_derivative_2l(~isnan(pvalb_sst)));
x = zscore(pvalb_sst(~isnan(pvalb_sst)));
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
save_filename=strcat("hcp3t_corr_mean_derivative_2l_pvalb_sst");
save_path=strcat("../results_bins/tSNR");
print(gcf, '-dpng', '-r1000', strcat(save_path, "/", save_filename));
print(gcf, '-dmeta', strcat(save_path, "/", save_filename));
close;
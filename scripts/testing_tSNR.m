% testing tSNR.
addpath(genpath("D:\Program Files\MATLAB\R2023a\toolbox\dzjin/fieldtrip-master/external/freesurfer"));
addpath(genpath("D:\Program Files\MATLAB\R2023a\toolbox\dzjin\spm12"));
addpath(genpath("D:\Program Files\MATLAB\R2023a\toolbox\dzjin/brainstat_matlab/io"));
addpath(genpath("D:\Program Files\MATLAB\R2023a\toolbox\dzjin/BrainSpace-0.1.10"));
addpath(genpath("D:\Program Files\MATLAB\R2023a\toolbox\dzjin\gifti-main"));
addpath(genpath("../../fcn"));

X = ft_read_cifti(char(strcat("../../atlas/hcp_fslr32k_cifti/Schaefer2018_200Parcels_17Networks_order.dlabel.nii")), 'mpname', 'array');
vertex_to_roi = X.parcels;
vertex_to_roi(vertex_to_roi==0)=NaN;


%% group_level tSNR.
atlas="schaefer400x17";
data_dir=strcat("G:\datasets\HCP\hcp_3T_unrelated_100_tSNR\", atlas);
load("../../data/hcp/subj_ids.mat");
tSNR_all = nan(100, 4, 400);
for subj_i=1:100
    disp(subj_i);
    X = ft_read_cifti(char(strcat(data_dir, "/", string(subj_ids(subj_i, :)), "_rfMRI_REST1_LR_tSNR_", atlas, ".pscalar.nii")), 'mpname', 'array');
    tSNR_all(subj_i, 1, :) = X.pscalar;
    X = ft_read_cifti(char(strcat(data_dir, "/", string(subj_ids(subj_i, :)), "_rfMRI_REST1_RL_tSNR_", atlas, ".pscalar.nii")), 'mpname', 'array');
    tSNR_all(subj_i, 2, :) = X.pscalar;
    X = ft_read_cifti(char(strcat(data_dir, "/", string(subj_ids(subj_i, :)), "_rfMRI_REST2_LR_tSNR_", atlas, ".pscalar.nii")), 'mpname', 'array');
    tSNR_all(subj_i, 3, :) = X.pscalar;
    X = ft_read_cifti(char(strcat(data_dir, "/", string(subj_ids(subj_i, :)), "_rfMRI_REST2_RL_tSNR_", atlas, ".pscalar.nii")), 'mpname', 'array');
    tSNR_all(subj_i, 4, :) = X.pscalar;
end
tSNR_mean = squeeze(mean(squeeze(mean(tSNR_all))));

% % plot tSNR.
% plot_data = vertex_to_roi;
% for roi_i=1:200
%     plot_data(plot_data==roi_i)=tSNR_mean(roi_i);
% end
% save_dir=strcat("../results_bins/DRIVER/HCP/tSNR");
% save_filename = "hcp37_tSNR";
% plot_my_fslr32k_v_parcel(plot_data, turbo, strcat(save_dir, "/", save_filename));

[~,idxsort] = sort(tSNR_mean, 'ascend');
low_tSNR_roi_idx=idxsort(1:10);
save_dir=strcat("../results_bins/tSNR");
save(strcat(save_dir, strcat("/hcp3t_unrelated100_tSNR_", atlas)), "tSNR_all", "tSNR_mean", "low_tSNR_roi_idx");


%% group_level tSNR.
atlas="glasser360";
data_dir=strcat("G:\datasets\HCP\hcp_3T_unrelated_100_tSNR\", atlas);
load("../../data/hcp/subj_ids.mat");
tSNR_all = nan(100, 4, 360);
for subj_i=1:100
    disp(subj_i);
    X = ft_read_cifti(char(strcat(data_dir, "/", string(subj_ids(subj_i, :)), "_rfMRI_REST1_LR_tSNR_", atlas, ".pscalar.nii")), 'mpname', 'array');
    tSNR_all(subj_i, 1, :) = X.pscalar;
    X = ft_read_cifti(char(strcat(data_dir, "/", string(subj_ids(subj_i, :)), "_rfMRI_REST1_RL_tSNR_", atlas, ".pscalar.nii")), 'mpname', 'array');
    tSNR_all(subj_i, 2, :) = X.pscalar;
    X = ft_read_cifti(char(strcat(data_dir, "/", string(subj_ids(subj_i, :)), "_rfMRI_REST2_LR_tSNR_", atlas, ".pscalar.nii")), 'mpname', 'array');
    tSNR_all(subj_i, 3, :) = X.pscalar;
    X = ft_read_cifti(char(strcat(data_dir, "/", string(subj_ids(subj_i, :)), "_rfMRI_REST2_RL_tSNR_", atlas, ".pscalar.nii")), 'mpname', 'array');
    tSNR_all(subj_i, 4, :) = X.pscalar;
end
tSNR_mean = squeeze(mean(squeeze(mean(tSNR_all))));

[~,idxsort] = sort(tSNR_mean, 'ascend');
low_tSNR_roi_idx=idxsort(1:10);
save_dir=strcat("../results_bins/tSNR");
save(strcat(save_dir, strcat("/hcp3t_unrelated100_tSNR_", atlas)), "tSNR_all", "tSNR_mean", "low_tSNR_roi_idx");


%% Remove the lowest 5% tSNR regions
atlas="schaefer200";
save_dir=strcat("../results_bins/tSNR");
load(strcat(save_dir, strcat("/hcp3t_unrelated100_tSNR_", atlas)), "tSNR_all", "tSNR_mean", "low_tSNR_roi_idx");
roi_idx = 1:200;
roi_idx(low_tSNR_roi_idx)=[];
tSNR_mean(roi_idx)=NaN;

% plot tSNR.
plot_data = vertex_to_roi;
for roi_i=1:200
    plot_data(plot_data==roi_i)=tSNR_mean(roi_i);
end
save_dir=strcat("../results_bins/tSNR");
save_filename = "hcp37_tSNR_remove";
plot_my_fslr32k_v_parcel(plot_data, turbo, strcat(save_dir, "/", save_filename));


%% neurobiological.
new_gam_model = readtable("../results_bins/DRIVER/HCP/1_region_level_ratioOfMeans/" + ...
    "bins_gam_amp_zscored_bins_rss_high_fd-0_num_of_bins-20_atlas-schaefer200_regress_cofounds-0_gsr-proc_regress_fix_ts_mat.csv");
myelin_data = load("../results_bins/mechanism/myelin/myelin_data.mat");
pvalb_sst_data = load("../results_bins/mechanism/pvalb_sst/pvalb_sst_data.mat");

partial_r2 = new_gam_model.GAM_age_partialR2(roi_idx);
myelin = myelin_data.myelin_map_schaefer200x17(roi_idx);
pvalb_sst = pvalb_sst_data.abha_expression_data_pvalb_sst_schaefer200x17(roi_idx);

partial_r2_myelin = corr(partial_r2, myelin, ...
    "type","Spearman", 'Rows','pairwise');
partial_r2_pvalb_sst = corr(partial_r2, pvalb_sst, ...
    "type","Spearman", 'Rows','pairwise');

y = zscore(partial_r2(~isnan(myelin)));
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
save_filename=strcat("hcp3t_corr_partial_r2_myelin");
save_path=strcat("../results_bins/tSNR");
print(gcf, '-dpng', '-r1000', strcat(save_path, "/", save_filename));
print(gcf, '-dmeta', strcat(save_path, "/", save_filename));
close;


y = zscore(partial_r2(~isnan(pvalb_sst)));
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
save_filename=strcat("hcp3t_corr_partial_r2_pvalb_sst");
save_path=strcat("../results_bins/tSNR");
print(gcf, '-dpng', '-r1000', strcat(save_path, "/", save_filename));
print(gcf, '-dmeta', strcat(save_path, "/", save_filename));
close;



%% bins 20 corr colormap.
maxcolor    = [247/255,126/255,105/255];
mediancolor = [255/255 255/255 255/255];
mincolor    = [52/255,96/255,141/255];
ColorMapSize = 50;
int1 = zeros(ColorMapSize,3);
int2 = zeros(ColorMapSize,3);
for k=1:3
    int1(:,k) = linspace(mincolor(k), mediancolor(k), ColorMapSize);
    int2(:,k) = linspace(mediancolor(k), maxcolor(k), ColorMapSize);
end
colormap_corrs = [int1(1:end-1,:); int2];

session_names = ['REST1_LR'; 'REST1_RL'; 'REST2_LR'; 'REST2_RL'];
num_of_bins=20;
pipeline='proc_regress_fix_ts_mat';
atlas='schaefer200';
high_fd=0;
regress_cofounds=0;
cs_name = "ets_region_level_contribution_ratioOfMeans_all";
rest1_lr = load(['../outputs_bins/HCP/bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(1, :) '_' atlas '_hfd_' num2str(high_fd) '_regressCofounds_' num2str(regress_cofounds) '.mat'], cs_name);
rest1_rl = load(['../outputs_bins/HCP/bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(2, :) '_' atlas '_hfd_' num2str(high_fd) '_regressCofounds_' num2str(regress_cofounds) '.mat'], cs_name);
rest2_lr = load(['../outputs_bins/HCP/bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(3, :) '_' atlas '_hfd_' num2str(high_fd) '_regressCofounds_' num2str(regress_cofounds) '.mat'], cs_name);
rest2_rl = load(['../outputs_bins/HCP/bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(4, :) '_' atlas '_hfd_' num2str(high_fd) '_regressCofounds_' num2str(regress_cofounds) '.mat'], cs_name);
load("../outputs_bins/HCP/physio_phase/BOLD_PLVs_proc_regress_fix_ts_mat_schaefer200.mat", "physios_avg1", "physios_ds");
HR_all = permute(physios_avg1, [2 3 1]);
BR_all = permute(physios_ds, [2 3 1]);
is_physio_exists = zeros(100, 4);
for subj_i=1:100
    for ses_i=1:4
        if sum(isnan(BR_all(subj_i, ses_i, :)))==0 && sum(isnan(HR_all(subj_i, ses_i, :)))==0
            is_physio_exists(subj_i, ses_i)=1;
        end
    end
end
bins_corr_rest1_lr = corr(squeeze(mean(rest1_lr.(cs_name)(is_physio_exists(:, 1)==1, roi_idx, :), 1)));
bins_corr_rest1_rl = corr(squeeze(mean(rest1_rl.(cs_name)(is_physio_exists(:, 2)==1, roi_idx, :), 1)));
bins_corr_rest2_lr = corr(squeeze(mean(rest2_lr.(cs_name)(is_physio_exists(:, 3)==1, roi_idx, :), 1)));
bins_corr_rest2_rl = corr(squeeze(mean(rest2_rl.(cs_name)(is_physio_exists(:, 4)==1, roi_idx, :), 1)));
bins_corr = (bins_corr_rest1_lr + bins_corr_rest1_rl + bins_corr_rest2_lr + bins_corr_rest2_rl)./4;
% plot.
figure("Position", [100, 100, 250, 240]);
imagesc(bins_corr);
colormap(colormap_corrs);
clim([-1 1]);
set(gca, 'XTick', [1, 5, 10, 15, 20], 'TickDir', 'None', ...
    'XTickLabel', ["95-100%", "75-80%", "50-55%", "25-30%", "0-5%"], "XTickLabelRotation", 45);
set(gca, 'YTick', [1, 5, 10, 15, 20], 'TickDir', 'None', ...
    'YTickLabel', ["95-100%", "75-80%", "50-55%", "25-30%", "0-5%"], "YTickLabelRotation", 0);
% set(gca, 'XTick', [1, 5, 10, 15, 20], 'TickDir', 'None', ...
%     'XTickLabel', ["", "", "", "", ""]);
% set(gca, 'YTick', [1, 5, 10, 15, 20], 'TickDir', 'None', ...
%     'YTickLabel', ["", "", "", "", ""]);
ax=gca;hold on;
ax.FontSize=12;
save_filename=strcat("bins_corr_high_fd-",num2str(high_fd), ...
    "_num_of_bins-",num2str(num_of_bins), ...
    "_atlas-",atlas, ...
    "_regress_cofounds-",num2str(regress_cofounds), ...
    "_gsr-", pipeline);
save_path=strcat("../results_bins/tSNR");
print(gcf, '-dpng', '-r1000', strcat(save_path, "/", save_filename));
print(gcf, '-dmeta', strcat(save_path, "/", save_filename));
close;

%% SA axis aligment.
% cal group-level cs.
roi_contri_group_rest1_lr = squeeze(mean(rest1_lr.(cs_name)(is_physio_exists(:, 1)==1, :, :), 1));
roi_contri_group_rest1_rl = squeeze(mean(rest1_rl.(cs_name)(is_physio_exists(:, 2)==1, :, :), 1));
roi_contri_group_rest2_lr = squeeze(mean(rest2_lr.(cs_name)(is_physio_exists(:, 3)==1, :, :), 1));
roi_contri_group_rest2_rl = squeeze(mean(rest2_rl.(cs_name)(is_physio_exists(:, 4)==1, :, :), 1));
roi_contri_group = roi_contri_group_rest1_lr + roi_contri_group_rest1_rl + ...
roi_contri_group_rest2_lr + roi_contri_group_rest2_rl;
roi_contri_group = roi_contri_group./4;
% read sa axis.
brain_map_sa_axis = readmatrix("../results_bins/mechanism/SA_Axis/s_a_axis_shaefer200x17.csv");
brain_map_sa_axis = brain_map_sa_axis(roi_idx);
% cal corr between cs and sa.
coefs = nan(num_of_bins, 1);
for bin_i=1:num_of_bins
    [coef, pval] = corr(brain_map_sa_axis, roi_contri_group(roi_idx, bin_i), "type", "Spearman");
    disp(coef);
    coefs(bin_i) = coef;
end
% file_name = strcat(cs_name, "_sa_coefs_bins_", num2str(num_of_bins), "_bin_hcp_", pipeline, "_scan_", atlas, "_hfd_0.mat");
% save(strcat("../results_bins/DRIVER/HCP/",save_dir_name, "/", file_name), "roi_contri_group", "coefs");
% plot results.
coefs_abs = abs(coefs);
min_idx = find(coefs_abs==min(coefs_abs(1:num_of_bins/2)));
x_1 = linspace(1, min_idx, min_idx);
x_2 = linspace(min_idx, num_of_bins, num_of_bins-min_idx+1);
figure("Position", [100 100 500 360]);
plot([1, 20], [0, 0], '-o', 'Color', [235/255 235/255 235/255], 'LineWidth', 0.8, 'MarkerFaceColor',[218,165,32]./255, 'MarkerSize', 0.001); hold on;
plot(x_1, coefs(1:min_idx), '-o', 'Color', [189,189,189]./255, 'LineWidth', 1.5, ...
'MarkerFaceColor',[218,165,32]./255, 'MarkerSize',6, 'MarkerEdgeColor','none'); hold on;
plot(x_2, coefs(min_idx:num_of_bins), '-o', 'Color', [189,189,189]./255, 'LineWidth', 1.5, ...
'MarkerFaceColor',[111 18 130]./255, 'MarkerSize',6, 'MarkerEdgeColor','none'); hold on;
ylim([-1 1]);
xlim([0, num_of_bins+1]);
xticks([1,5,10,15,20]);
% xticklabels([" ", " ", " ", " ", " "]);
xticklabels(["95-100%", "75-80%", "50-55%", "25-30%", "0-5%"]);

yticks(-1:0.2:1);
ax=gca;hold on;
ax.FontSize=16;
ax.XTickLabelRotation=45;
set(gca,'Box','off');
save_filename=strcat("bins_cs_sa_high_fd-",num2str(high_fd), ...
"_num_of_bins-",num2str(num_of_bins), ...
"_atlas-",atlas, ...
"_regress_cofounds-",num2str(regress_cofounds), ...
"_gsr-", pipeline);
save_path=strcat("../results_bins/tSNR");
print(gcf, '-dpng', '-r1000', strcat(save_path, "/", save_filename));
print(gcf, '-dmeta', strcat(save_path, "/", save_filename));
close;











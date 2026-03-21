% testing tSNR.
addpath(genpath("D:\Program Files\MATLAB\R2023a\toolbox\dzjin/fieldtrip-master/external/freesurfer"));
addpath(genpath("D:\Program Files\MATLAB\R2023a\toolbox\dzjin\spm12"));
addpath(genpath("D:\Program Files\MATLAB\R2023a\toolbox\dzjin/brainstat_matlab/io"));
addpath(genpath("D:\Program Files\MATLAB\R2023a\toolbox\dzjin/BrainSpace-0.1.10"));
addpath(genpath("D:\Program Files\MATLAB\R2023a\toolbox\dzjin\gifti-main"));
addpath(genpath("../../fcn"));

% subject ids.
load('../../data/hcp/subj_ids.mat');
num_of_bins = 20;
pipeline = 'proc_regress_fix_ts_mat'; % proc_regress_fixWglob_ts_mat, proc_regress_fix_ts_mat
atlas = 'schaefer200'; % schaefer200, schaefer400, hcpmmp
num_of_rois = 200;
regressCofounds = 1; % regress global signal, fd, respiration, heart rate.
cs_name="ets_region_level_contribution_ratioOfMeans_all";
% read sa axis.
brain_map_sa_axis = readmatrix("../results_bins/mechanism/SA_Axis/s_a_axis_shaefer200x17.csv");
load("../../atlas/Schaefer200x17.mat");
[~,idxsort] = sort(lab17to8);


%% load physio data.
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


%{
  individual-level cs.
%}
%% load cs data.
session_names = ['REST1_LR'; 'REST1_RL'; 'REST2_LR'; 'REST2_RL'];
rest1_lr = load(['../outputs_bins/HCP/bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(1, :) '_' atlas '_hfd_0_regressCofounds_' num2str(regressCofounds) '.mat'], cs_name);
rest1_rl = load(['../outputs_bins/HCP/bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(2, :) '_' atlas '_hfd_0_regressCofounds_' num2str(regressCofounds) '.mat'], cs_name);
rest2_lr = load(['../outputs_bins/HCP/bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(3, :) '_' atlas '_hfd_0_regressCofounds_' num2str(regressCofounds) '.mat'], cs_name);
rest2_rl = load(['../outputs_bins/HCP/bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(4, :) '_' atlas '_hfd_0_regressCofounds_' num2str(regressCofounds) '.mat'], cs_name);

% organize cs data.
cs_all = nan(4, 100, 200, 20);
cs_all(1, :, :, :) = rest1_lr.(cs_name);
cs_all(2, :, :, :) = rest1_rl.(cs_name);
cs_all(3, :, :, :) = rest2_lr.(cs_name);
cs_all(4, :, :, :) = rest2_rl.(cs_name);
cs_all(is_physio_exists.'==0, :)=NaN;
cs_avg4Sess = squeeze(mean(cs_all, 1, 'omitmissing'));


%% ind-level cs and SA Axis.
figure("Position", [100 100 500 360]);
plot([1, 20], [0, 0], '-o', 'Color', [0/255 0/255 0/255], 'LineWidth', 0.5, 'MarkerFaceColor',[218,165,32]./255, 'MarkerSize', 0.001); hold on;
coefs = nan(100, num_of_bins);
for subj_i=1:100
    % cal corr between cs and sa.
    for bin_i=1:num_of_bins
        [coef, pval] = corr(brain_map_sa_axis, squeeze(cs_avg4Sess(subj_i, :, bin_i)).', "type", "Spearman");
        coefs(subj_i, bin_i) = coef;
    end
    x_1 = linspace(1, num_of_bins, num_of_bins);
    plot(x_1, coefs(subj_i, :), '-o', 'Color', [189,189,189]./255, 'LineWidth', 0.5, ...
        'MarkerFaceColor',[99,99,99]./255, 'MarkerSize',2, 'MarkerEdgeColor','none'); hold on;
end
ylim([-1 1]);
xlim([0, num_of_bins+1]);
xticks([1,5,10,15,20]);
xticklabels([" ", " ", " ", " ", " "]);
xticklabels(["95-100%", "75-80%", "50-55%", "25-30%", "0-5%"]);
yticks(-1:0.2:1);
% ylabel("Similarity between " + newline + "Cofluctuation Score and SA Rank (r)");
ax=gca;hold on;
ax.FontSize=13;
ax.XTickLabelRotation=45;
set(gca,'Box','off');
% set(gca,'Box','off', 'TickDir', 'none');
file_name=strcat("hcp_all_sess_spearman_cs_sa_axis_raw_bins_", num2str(num_of_bins), "_bin_hcp_", pipeline, "_scan_", atlas, "_hfd_0_ind");
save_path_filename = strcat("../results_bins/DRIVER/HCP/1_region_level_ratioOfMeans/", file_name);
print(gcf, '-dpng', '-r600', save_path_filename);
print(gcf, '-dmeta', save_path_filename);
close;
save(save_path_filename, "coefs");


%% ind-level cs and SMN-VIS Axis.
hcp_group_grads = load("../../func_gradient/hcp_results/parcel_nfc_grad/proc_regress_fix_ts_mat_schaefer200_hfd_0_nfc_grad.mat");
hcp_group_grad_2 = hcp_group_grads.func_gradient(:, 2);
hcp_group_grad_2 = zscore(hcp_group_grad_2);
figure("Position", [100 100 500 360]);
plot([1, 20], [0, 0], '-o', 'Color', [0/255 0/255 0/255], 'LineWidth', 0.5, 'MarkerFaceColor',[218,165,32]./255, 'MarkerSize', 0.001); hold on;
coefs = nan(100, num_of_bins);
for subj_i=1:100
    % cal corr between cs and sa.
    for bin_i=1:num_of_bins
        [coef, pval] = corr(hcp_group_grad_2, squeeze(cs_avg4Sess(subj_i, :, bin_i)).', "type", "Spearman");
        coefs(subj_i, bin_i) = coef;
    end
    x_1 = linspace(1, num_of_bins, num_of_bins);
    plot(x_1, coefs(subj_i, :), '-o', 'Color', [189,189,189]./255, 'LineWidth', 0.5, ...
        'MarkerFaceColor',[99,99,99]./255, 'MarkerSize',2, 'MarkerEdgeColor','none'); hold on;
end
ylim([-1 1]);
xlim([0, num_of_bins+1]);
xticks([1,5,10,15,20]);
xticklabels([" ", " ", " ", " ", " "]);
xticklabels(["95-100%", "75-80%", "50-55%", "25-30%", "0-5%"]);
yticks(-1:0.2:1);
% ylabel("Similarity between " + newline + "Cofluctuation Score and SA Rank (r)");
ax=gca;hold on;
ax.FontSize=13;
ax.XTickLabelRotation=45;
set(gca,'Box','off');
% set(gca,'Box','off', 'TickDir', 'none');
file_name=strcat("hcp_all_sess_spearman_cs_smn_vis_axis_raw_bins_", num2str(num_of_bins), "_bin_hcp_", pipeline, "_scan_", atlas, "_hfd_0_ind");
save_path_filename = strcat("../results_bins/DRIVER/HCP/1_region_level_ratioOfMeans/", file_name);
print(gcf, '-dpng', '-r600', save_path_filename);
print(gcf, '-dmeta', save_path_filename);
close;
save(save_path_filename, "coefs");


%% ind-level cs and ind-level primary and secondary functional gradient.
grad_i=1;
figure("Position", [100 100 500 360]);
plot([1, 20], [0, 0], '-o', 'Color', [0/255 0/255 0/255], 'LineWidth', 0.5, 'MarkerFaceColor',[218,165,32]./255, 'MarkerSize', 0.001); hold on;
coefs = nan(100, num_of_bins);
for subj_i=1:100
    % load ind-level primary functional gradient.
    load(strcat("../../func_gradient/hcp_results/parcel_nfc_grad_ind/",string(subj_ids(subj_i, :)),"_proc_regress_fixWglob_ts_mat_schaefer200_hfd_0_nfc_grad.mat"), "func_gradient");
    grad = func_gradient(:, grad_i);
    if corr(grad, brain_map_sa_axis) < 0
        grad = -grad;
    end
    grad = zscore(grad);
    % cal corr between cs and ind-level primary gradient.
    for bin_i=1:num_of_bins
        [coef, pval] = corr(grad, squeeze(cs_avg4Sess(subj_i, :, bin_i)).', "type", "Spearman");
        coefs(subj_i, bin_i) = coef;
    end
    x_1 = linspace(1, num_of_bins, num_of_bins);
    plot(x_1, coefs(subj_i, :), '-o', 'Color', [189,189,189]./255, 'LineWidth', 0.5, ...
        'MarkerFaceColor',[99,99,99]./255, 'MarkerSize',2, 'MarkerEdgeColor','none'); hold on;
end
ylim([-1 1]);
xlim([0, num_of_bins+1]);
xticks([1,5,10,15,20]);
xticklabels([" ", " ", " ", " ", " "]);
xticklabels(["95-100%", "75-80%", "50-55%", "25-30%", "0-5%"]);
yticks(-1:0.2:1);
% ylabel("Similarity between " + newline + "Cofluctuation Score and SA Rank (r)");
ax=gca;hold on;
ax.FontSize=13;
ax.XTickLabelRotation=45;
set(gca,'Box','off');
% set(gca,'Box','off', 'TickDir', 'none');
file_name=strcat("hcp_all_sess_spearman_cs_grad",num2str(grad_i),"_axis_raw_bins_", num2str(num_of_bins), "_bin_hcp_", pipeline, "_scan_", atlas, "_hfd_0_ind");
save_path_filename = strcat("../results_bins/DRIVER/HCP/1_region_level_ratioOfMeans/", file_name);
print(gcf, '-dpng', '-r600', save_path_filename);
print(gcf, '-dmeta', save_path_filename);
close;
save(save_path_filename, "coefs");


%{
  weighted reconstruction of nFC.
%}
%% load bin-level ets data.
variable_name="ets_mean_bins_all";
session_names = ['REST1_LR'; 'REST1_RL'; 'REST2_LR'; 'REST2_RL'];
rest1_lr = load(['../outputs_bins/HCP/bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(1, :) '_' atlas '_hfd_0_regressCofounds_' num2str(regressCofounds) '.mat'], variable_name);
rest1_rl = load(['../outputs_bins/HCP/bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(2, :) '_' atlas '_hfd_0_regressCofounds_' num2str(regressCofounds) '.mat'], variable_name);
rest2_lr = load(['../outputs_bins/HCP/bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(3, :) '_' atlas '_hfd_0_regressCofounds_' num2str(regressCofounds) '.mat'], variable_name);
rest2_rl = load(['../outputs_bins/HCP/bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(4, :) '_' atlas '_hfd_0_regressCofounds_' num2str(regressCofounds) '.mat'], variable_name);
% organize ets_mean_bin data.
ets_mean_bin_all = nan(4, 100, 20, 19900);
ets_mean_bin_all(1, :, :, :) = rest1_lr.(variable_name);
ets_mean_bin_all(2, :, :, :) = rest1_rl.(variable_name);
ets_mean_bin_all(3, :, :, :) = rest2_lr.(variable_name);
ets_mean_bin_all(4, :, :, :) = rest2_rl.(variable_name);
for subj_i=1:100
    for ses_i=1:4
        if is_physio_exists(subj_i, ses_i)==0
            ets_mean_bin_all(ses_i, subj_i, :, :)=NaN;
        end
    end
end
ets_mean_bin_avg4Sess = squeeze(mean(ets_mean_bin_all, 1, 'omitmissing'));

file_name=strcat("hcp_all_sess_spearman_cs_grad1_axis_raw_bins_", num2str(num_of_bins), "_bin_hcp_", pipeline, "_scan_", atlas, "_hfd_0_ind");
save_path_filename = strcat("../results_bins/DRIVER/HCP/1_region_level_ratioOfMeans/", file_name);
load(save_path_filename, "coefs");

for subj_i=1:100
    for bin_i=1:num_of_bins
        ets_mean_bin_avg4Sess(subj_i, bin_i, :) = abs(coefs(subj_i, bin_i)).*ets_mean_bin_avg4Sess(subj_i, bin_i, :);
    end
end
ets_mean_bin_avg4Sess = squeeze(mean(ets_mean_bin_avg4Sess, 1, 'omitmissing'));
ets_mean_bin_avg4Sess = squeeze(mean(ets_mean_bin_avg4Sess, 1, 'omitmissing'));

% Index of upper triabgle in FC matrix.
[u,v] = find(triu(ones(num_of_rois),1));
idx = (v - 1)*num_of_rois + u;
weighted_reconstruct_nfc = zeros(200, 200);
weighted_reconstruct_nfc(idx) = ets_mean_bin_avg4Sess;
weighted_reconstruct_nfc = weighted_reconstruct_nfc + weighted_reconstruct_nfc.';

plot_nfc_schaefer200(weighted_reconstruct_nfc, lab17to8, idxsort, net17to8, 200, 1, -0.4, 0.4);
print(gcf, '-dmeta', strcat("../results_bins/DRIVER/HCP/1_region_level_ratioOfMeans/hcp3t_weighted_reconstruct_nfc"));
print(gcf, '-dpng', strcat("../results_bins/DRIVER/HCP/1_region_level_ratioOfMeans/hcp3t_weighted_reconstruct_nfc"));
close;


%% cal and plot grad.
grad = GradientMaps('kernel','cs','approach', 'dm');
grad = grad.fit(weighted_reconstruct_nfc, 'Sparsity', 90);
func_gradient = grad.gradients{1}(:, 1:10);
variance_explained = grad.lambda{1};
grad_1=func_gradient(:, 2);
grad_1=zscore(grad_1);

X = ft_read_cifti(char(strcat("../../atlas/hcp_fslr32k_cifti/Schaefer2018_200Parcels_17Networks_order.dlabel.nii")), 'mpname', 'array');
vertex_to_roi = X.parcels;
vertex_to_roi(vertex_to_roi==0)=NaN;

% plot grad.
plot_data = vertex_to_roi;
for roi_i=1:200
    plot_data(plot_data==roi_i)=grad_1(roi_i);
end
save_dir=strcat("../results_bins/DRIVER/HCP/1_region_level_ratioOfMeans/");
save_filename = "hcp3t_weighted_reconstruct_nfc_grad_2";
plot_my_fslr32k_v_parcel(plot_data, turbo, strcat(save_dir, "/", save_filename));
save(strcat(save_dir, "/hcp3t_weighted_reconstruct_nfc_grads"), "variance_explained", "func_gradient", "weighted_reconstruct_nfc");


%% plot hcp3t group nfc.
hcp_group_nfc = load("../../func_gradient/hcp_results/parcel_nfc_grad/proc_regress_fix_ts_mat_schaefer200_hfd_0_nfc_grad.mat");
hcp_group_nfc = hcp_group_nfc.nfc_group;
hcp_group_nfc = tanh(hcp_group_nfc);
plot_nfc_schaefer200(hcp_group_nfc, lab17to8, idxsort, net17to8, 200, 1, -1, 1);
print(gcf, '-dmeta', strcat("../results_bins/DRIVER/HCP/1_region_level_ratioOfMeans/hcp3t_c_nfc"));
print(gcf, '-dpng', strcat("../results_bins/DRIVER/HCP/1_region_level_ratioOfMeans/hcp3t_hcp_group_nfc_nfc"));
close;


%% corr nfc.
[u,v] = find(triu(ones(num_of_rois),1));
idx = (v - 1)*num_of_rois + u;
[coef, pval] = corr(weighted_reconstruct_nfc(idx), hcp_group_nfc(idx), "type","Spearman");
x = weighted_reconstruct_nfc(idx);
y = hcp_group_nfc(idx);
figure("Position", [100 100 360 360]);
scatter(x, y, 40, 'filled', 'MarkerFaceAlpha',0.4);
hold on
pfit = polyfit(x,y,1);
x_fit = linspace(min(x),max(x),100);
y_fit = polyval(pfit,x_fit);
plot(x_fit,y_fit,'r','LineWidth',2);
ax=gca;hold on;
ax.FontSize=15;
set(gca,'Box','off');
save_filename=strcat("hcp3t_nfc_reconstructed_nfc");
save_path=strcat("../results_bins/DRIVER/HCP/1_region_level_ratioOfMeans");
print(gcf, '-dpng', '-r1000', strcat(save_path, "/", save_filename));
print(gcf, '-dmeta', strcat(save_path, "/", save_filename));
close;


%% corr grad.
hcp_group_grads = load("../../func_gradient/hcp_results/parcel_nfc_grad/proc_regress_fix_ts_mat_schaefer200_hfd_0_nfc_grad.mat");
grad_i=1;
x = zscore(func_gradient(:, grad_i));
y = zscore(hcp_group_grads.func_gradient(:, grad_i));
[coef, pval] = corr(x, y, "type","Spearman");
figure("Position", [100 100 360 360]);
scatter(x, y, 40, 'filled', 'MarkerFaceAlpha',0.4);
hold on
pfit = polyfit(x,y,1);
x_fit = linspace(min(x),max(x),100);
y_fit = polyval(pfit,x_fit);
plot(x_fit,y_fit,'r','LineWidth',2);
ax=gca;hold on;
ax.FontSize=15;
set(gca,'Box','off');
save_filename=strcat("hcp3t_corr_grad_", num2str(grad_i));
save_path=strcat("../results_bins/DRIVER/HCP/1_region_level_ratioOfMeans");
print(gcf, '-dpng', '-r1000', strcat(save_path, "/", save_filename));
print(gcf, '-dmeta', strcat(save_path, "/", save_filename));
close;


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
session_names = ['REST1_LR'; 'REST1_RL'; 'REST2_LR'; 'REST2_RL'];


%% load local_global_lead_lag_all.
var_name="rss_local_global_lead_lag_all";
rest1_lr = load(['../outputs_bins\HCP\bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(1, :) '_' atlas '_hfd_0_regressCofounds_' num2str(regressCofounds) '.mat'], var_name);
rest1_rl = load(['../outputs_bins\HCP\bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(2, :) '_' atlas '_hfd_0_regressCofounds_' num2str(regressCofounds) '.mat'], var_name);
rest2_lr = load(['../outputs_bins\HCP\bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(3, :) '_' atlas '_hfd_0_regressCofounds_' num2str(regressCofounds) '.mat'], var_name);
rest2_rl = load(['../outputs_bins\HCP\bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(4, :) '_' atlas '_hfd_0_regressCofounds_' num2str(regressCofounds) '.mat'], var_name);
local_global_lead_lag_all = nan(100, 4, num_of_rois, 2);
local_global_lead_lag_all(:, 1, :, :) = rest1_lr.(var_name);
local_global_lead_lag_all(:, 2, :, :) = rest1_rl.(var_name);
local_global_lead_lag_all(:, 3, :, :) = rest2_lr.(var_name);
local_global_lead_lag_all(:, 4, :, :) = rest2_rl.(var_name);

local_global_lead_lag_all = mean(squeeze(mean(squeeze(local_global_lead_lag_all(:, :, :, 1)), 'omitmissing')), 'omitmissing');

X = ft_read_cifti(char(strcat("../../atlas/hcp_fslr32k_cifti/Schaefer2018_200Parcels_17Networks_order.dlabel.nii")), 'mpname', 'array');
vertex_to_roi = X.parcels;
vertex_to_roi(vertex_to_roi==0)=NaN;

% plot phase.
plot_data = vertex_to_roi;
for roi_i=1:200
    plot_data(plot_data==roi_i)=local_global_lead_lag_all(roi_i);
end
save_dir=strcat("../results_bins/local_global_lead_lag/");
save_filename = "hcp3t_rss_local_global_lead_lag";
plot_my_fslr32k_v_parcel(plot_data, turbo, strcat(save_dir, "/", save_filename));

save_dir=strcat("../results_bins/local_global_lead_lag/");
save_filename = "hcp3t_rss_local_global_lead_lag_hist";
figure("Position", [100 100 300 230]);
histogram(local_global_lead_lag_all);
xlim([-2.5 2.5]);
ax=gca;hold on;
ax.FontSize=13;
% ax.XTickLabelRotation=45;
set(gca,'Box','off');
print('-dpng', '-r600', strcat(save_dir, "/", save_filename));
print('-dmeta', strcat(save_dir, "/", save_filename));
close;



%% load local_global_lead_lag_all.
var_name="cs_local_global_lead_lag_all";
rest1_lr = load(['../outputs_bins\HCP\bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(1, :) '_' atlas '_hfd_0_regressCofounds_' num2str(regressCofounds) '.mat'], var_name);
rest1_rl = load(['../outputs_bins\HCP\bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(2, :) '_' atlas '_hfd_0_regressCofounds_' num2str(regressCofounds) '.mat'], var_name);
rest2_lr = load(['../outputs_bins\HCP\bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(3, :) '_' atlas '_hfd_0_regressCofounds_' num2str(regressCofounds) '.mat'], var_name);
rest2_rl = load(['../outputs_bins\HCP\bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(4, :) '_' atlas '_hfd_0_regressCofounds_' num2str(regressCofounds) '.mat'], var_name);
local_global_lead_lag_all = nan(100, 4, num_of_rois, 2);
local_global_lead_lag_all(:, 1, :, :) = rest1_lr.(var_name);
local_global_lead_lag_all(:, 2, :, :) = rest1_rl.(var_name);
local_global_lead_lag_all(:, 3, :, :) = rest2_lr.(var_name);
local_global_lead_lag_all(:, 4, :, :) = rest2_rl.(var_name);

local_global_lead_lag_all = mean(squeeze(mean(squeeze(local_global_lead_lag_all(:, :, :, 1)), 'omitmissing')), 'omitmissing');

X = ft_read_cifti(char(strcat("../../atlas/hcp_fslr32k_cifti/Schaefer2018_200Parcels_17Networks_order.dlabel.nii")), 'mpname', 'array');
vertex_to_roi = X.parcels;
vertex_to_roi(vertex_to_roi==0)=NaN;

% plot phase.
plot_data = vertex_to_roi;
for roi_i=1:200
    plot_data(plot_data==roi_i)=local_global_lead_lag_all(roi_i);
end
save_dir=strcat("../results_bins/local_global_lead_lag/");
save_filename = "hcp3t_cs_local_global_lead_lag";
plot_my_fslr32k_v_parcel(plot_data, turbo, strcat(save_dir, "/", save_filename));

save_dir=strcat("../results_bins/local_global_lead_lag/");
save_filename = "hcp3t_cs_local_global_lead_lag_hist";
figure("Position", [100 100 300 230]);
histogram(local_global_lead_lag_all);
xlim([-2.5 2.5]);
ax=gca;hold on;
ax.FontSize=13;
% ax.XTickLabelRotation=45;
set(gca,'Box','off');
print('-dpng', '-r600', strcat(save_dir, "/", save_filename));
print('-dmeta', strcat(save_dir, "/", save_filename));
close;
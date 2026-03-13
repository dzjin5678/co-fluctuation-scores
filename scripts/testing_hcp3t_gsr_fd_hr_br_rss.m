% subject ids.
load('../../data/hcp/subj_ids.mat');
num_of_bins = 20;
pipeline = 'proc_regress_fix_ts_mat'; % proc_regress_fixWglob_ts_mat, proc_regress_fix_ts_mat
atlas = 'schaefer200'; % schaefer200, schaefer400, hcpmmp
num_of_rois = 200;
regressCofounds = 0; % regress global signal, fd, respiration, heart rate.
session_names = ['REST1_LR'; 'REST1_RL'; 'REST2_LR'; 'REST2_RL'];


%% load phsio data.
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


%% load fd and gs.
fd_all = nan(100, 4, 1191);
gs_all = nan(100, 4, 1191);
for subj_i=1:100
    for ses_i=1:4
        confounds = readmatrix(['G:\datasets\HCP\regressors_cp/' subj_ids(subj_i, :) '/rfMRI_' session_names(1, :) '/' subj_ids(subj_i, :) '_rfMRI_' session_names(1, :) '_confounds.csv']);
        fd = confounds(:, 7);fd = fd(5:end-5);
        gs = confounds(:, 10);gs = gs(5:end-5);
        fd_all(subj_i, ses_i, :) = fd;
        gs_all(subj_i, ses_i, :) = gs;
    end
end


%% load global rss.
var_name="ets_global_rss_all";
rest1_lr = load(['../outputs_bins/HCP/bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(1, :) '_' atlas '_hfd_1_regressCofounds_' num2str(regressCofounds) '.mat'], var_name);
rest1_rl = load(['../outputs_bins/HCP/bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(2, :) '_' atlas '_hfd_1_regressCofounds_' num2str(regressCofounds) '.mat'], var_name);
rest2_lr = load(['../outputs_bins/HCP/bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(3, :) '_' atlas '_hfd_1_regressCofounds_' num2str(regressCofounds) '.mat'], var_name);
rest2_rl = load(['../outputs_bins/HCP/bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(4, :) '_' atlas '_hfd_1_regressCofounds_' num2str(regressCofounds) '.mat'], var_name);
global_rss_all = nan(100, 4, 1191);
global_rss_all(:, 1, :) = rest1_lr.(var_name);
global_rss_all(:, 2, :) = rest1_rl.(var_name);
global_rss_all(:, 3, :) = rest2_lr.(var_name);
global_rss_all(:, 4, :) = rest2_rl.(var_name);


%% corr.
global_rss_fd_corr = nan(100, 4);
global_rss_gs_corr = nan(100, 4);
global_rss_hr_corr = nan(100, 4);
global_rss_br_corr = nan(100, 4);
for subj_i=1:100
    for ses_i=1:4
        global_rss = squeeze(global_rss_all(subj_i, ses_i, :));
        fd = squeeze(fd_all(subj_i, ses_i, :));
        gs = squeeze(gs_all(subj_i, ses_i, :));
        hr = squeeze(HR_all(subj_i, ses_i, :));
        br = squeeze(BR_all(subj_i, ses_i, :));
        
        if is_physio_exists(subj_i, ses_i)==1
            global_rss_fd_corr(subj_i, ses_i) = corr(global_rss, fd, "type", "Spearman");
            global_rss_gs_corr(subj_i, ses_i) = corr(global_rss, gs, "type", "Spearman");
            global_rss_hr_corr(subj_i, ses_i) = corr(global_rss, hr, "type", "Spearman");
            global_rss_br_corr(subj_i, ses_i) = corr(global_rss, br, "type", "Spearman");
        end
    end
end
global_rss_br_corr=mean(global_rss_br_corr, "omitmissing");
global_rss_gs_corr=mean(global_rss_gs_corr, "omitmissing");
global_rss_hr_corr=mean(global_rss_hr_corr, "omitmissing");
global_rss_br_corr=mean(global_rss_br_corr, "omitmissing");

save("../results_bins/noise/hcp3t_gsr_fd_hr_br_rss", ...
    "global_rss_br_corr", ...
    "global_rss_gs_corr", ...
    "global_rss_hr_corr", ...
    "global_rss_br_corr");



%%
%% hr and br.
rest1_lr = load(['../outputs_bins\HCP\bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(1, :) '_' atlas '_hfd_0_regressCofounds_' num2str(regressCofounds) '.mat'], ...
 "ets_global_rss_all", "high_fd_label_all");
rest1_rl = load(['../outputs_bins\HCP\bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(2, :) '_' atlas '_hfd_0_regressCofounds_' num2str(regressCofounds) '.mat'], ...
 "ets_global_rss_all", "high_fd_label_all");
rest2_lr = load(['../outputs_bins\HCP\bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(3, :) '_' atlas '_hfd_0_regressCofounds_' num2str(regressCofounds) '.mat'], ...
 "ets_global_rss_all", "high_fd_label_all");
rest2_rl = load(['../outputs_bins\HCP\bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(4, :) '_' atlas '_hfd_0_regressCofounds_' num2str(regressCofounds) '.mat'], ...
 "ets_global_rss_all", "high_fd_label_all");
load("../outputs_bins/HCP/physio_phase/BOLD_PLVs_proc_regress_fix_ts_mat_schaefer200.mat", "physios_avg1", "physios_ds");
HR_all = permute(physios_avg1, [2 3 1]);
BR_all = permute(physios_ds, [2 3 1]);

ets_rss_all = nan(100, 4, 1191);
ets_rss_all(:, 1, :) = rest1_lr.ets_global_rss_all;
ets_rss_all(:, 2, :) = rest1_rl.ets_global_rss_all;
ets_rss_all(:, 3, :) = rest2_lr.ets_global_rss_all;
ets_rss_all(:, 4, :) = rest2_rl.ets_global_rss_all;

high_fd_label_all = nan(100, 4, 1191);
high_fd_label_all(:, 1, :) = rest1_lr.high_fd_label_all;
high_fd_label_all(:, 2, :) = rest1_rl.high_fd_label_all;
high_fd_label_all(:, 3, :) = rest2_lr.high_fd_label_all;
high_fd_label_all(:, 4, :) = rest2_rl.high_fd_label_all;

br_corr_all = [];
hr_corr_all = [];
for subj_i=1:100
    br_corr_subj_i=[];
    hr_corr_subj_i=[];
    for ses_i=1:4
        % if is_physio_exists(subj_i, ses_i)==1
        if sum(isnan(BR_all(subj_i, ses_i, :)))==0 && sum(isnan(HR_all(subj_i, ses_i, :)))==0
            frame_label = squeeze(high_fd_label_all(subj_i, ses_i, :)==0);
            a = squeeze(ets_rss_all(subj_i, ses_i, 1:sum(frame_label)));
            b = squeeze(BR_all(subj_i, ses_i, frame_label));
            br_ets_rss_corr = corr(squeeze(ets_rss_all(subj_i, ses_i, 1:sum(frame_label))), squeeze(BR_all(subj_i, ses_i, frame_label)));
            br_corr_subj_i = [br_corr_subj_i br_ets_rss_corr];
            hr_ets_rss_corr = corr(squeeze(ets_rss_all(subj_i, ses_i, 1:sum(frame_label))), squeeze(HR_all(subj_i, ses_i, frame_label)));
            hr_corr_subj_i = [hr_corr_subj_i hr_ets_rss_corr];
        end
    end
    br_corr_all = [br_corr_all mean(br_corr_subj_i)];
    hr_corr_all = [hr_corr_all mean(hr_corr_subj_i)];
end

figure('Position', [100 100 500 500]);
boxplot(br_corr_all);
ylim([-1 1]);
ax=gca;hold on;
ax.FontSize=20;
print(gcf, '-dpng', '-r600', "../results_bins/noise/HCP3T_br");
close;

figure('Position', [100 100 500 500]);
boxplot(hr_corr_all);
ylim([-1 1]);
ax=gca;hold on;
ax.FontSize=20;
print(gcf, '-dpng', '-r600', "../results_bins/noise/HCP3T_hr");
close;

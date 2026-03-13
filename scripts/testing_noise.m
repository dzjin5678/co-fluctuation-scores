global_rss_rest1_lr = load("../outputs_bins/HCP/bins_20_proc_regress_fix_ts_mat_scan_REST1_LR_schaefer200_hfd_1_regressCofounds_0.mat", "ets_global_rss_all");
global_rss_rest1_rl = load("../outputs_bins/HCP/bins_20_proc_regress_fix_ts_mat_scan_REST1_RL_schaefer200_hfd_1_regressCofounds_0.mat", "ets_global_rss_all");
global_rss_rest2_lr = load("../outputs_bins/HCP/bins_20_proc_regress_fix_ts_mat_scan_REST2_LR_schaefer200_hfd_1_regressCofounds_0.mat", "ets_global_rss_all");
global_rss_rest2_rl = load("../outputs_bins/HCP/bins_20_proc_regress_fix_ts_mat_scan_REST2_RL_schaefer200_hfd_1_regressCofounds_0.mat", "ets_global_rss_all");
global_rss_all = nan(4, 100, 1191);
global_rss_all(1, :, :) = global_rss_rest1_lr.ets_global_rss_all;
global_rss_all(2, :, :) = global_rss_rest1_rl.ets_global_rss_all;
global_rss_all(3, :, :) = global_rss_rest2_lr.ets_global_rss_all;
global_rss_all(4, :, :) = global_rss_rest2_rl.ets_global_rss_all;

physio_data = load("../outputs_bins/HCP/physio_phase/ETS_REGION_RSS_RV_PLVs_proc_regress_fix_ts_mat.mat");

fd_corr_all = nan(100, 4);
gs_corr_all = nan(100, 4);
ds_corr_all = nan(100, 4);
avg1_corr_all = nan(100, 4);

load("../../data/hcp/subj_ids.mat");
scan_names = ["REST1_LR", "REST1_RL", "REST2_LR", "REST2_RL"];
for subj_i=1:100
    for scan_i=1:4
        scan = char(scan_names(scan_i));
        subj_i_scan_i_global_rss = squeeze(global_rss_all(scan_i, subj_i, :));
        % load fd and gs.
        confounds = readmatrix(['G:/datasets/HCP/regressors_cp/' subj_ids(subj_i, :) '/rfMRI_' scan '/' subj_ids(subj_i, :) '_rfMRI_' scan '_confounds.csv']);
        fd = confounds(:, 7);fd = fd(5:end-5);
        gs = confounds(:, 10);gs = gs(5:end-5);
        ds = physio_data.physios_ds(:, subj_i, scan_i);
        avg1 = physio_data.physios_avg1(:, subj_i, scan_i);

        if sum(isnan(ds))>1 || sum(isnan(avg1))>1
            continue;
        else
            fd_corr_all(subj_i, scan_i) = corr(subj_i_scan_i_global_rss, fd, "type","Spearman");
            gs_corr_all(subj_i, scan_i) = corr(subj_i_scan_i_global_rss, gs, "type","Spearman");
            ds_corr_all(subj_i, scan_i) = corr(subj_i_scan_i_global_rss, ds, "type","Spearman");
            avg1_corr_all(subj_i, scan_i) = corr(subj_i_scan_i_global_rss, avg1, "type","Spearman");
        end
    end
end


fd_corr_mean_all = mean(fd_corr_all, 2, "omitmissing"); fd_corr_mean_all = fd_corr_mean_all(~isnan(fd_corr_mean_all));
gs_corr_mean_all = mean(gs_corr_all, 2, "omitmissing"); gs_corr_mean_all = gs_corr_mean_all(~isnan(gs_corr_mean_all));
ds_corr_mean_all = mean(ds_corr_all, 2, "omitmissing"); ds_corr_mean_all = ds_corr_mean_all(~isnan(ds_corr_mean_all));
avg1_corr_mean_all = mean(avg1_corr_all, 2, "omitmissing"); avg1_corr_mean_all = avg1_corr_mean_all(~isnan(avg1_corr_mean_all));
save("../results_bins/noise/fd_gs_ds_avg1_global_rss_corr.mat", ...
    "fd_corr_mean_all", ...
    "gs_corr_mean_all", ...
    "ds_corr_mean_all", ...
    "avg1_corr_mean_all");

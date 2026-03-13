clear
close all
clc


%% default settings.
% high_fd = 0; % remove high head motion frames.
% num_of_bins = 20; % 20 bins.
% atlas_i=1; % schaefer200x17.
% regress_cofounds=0; % not regress fd, gs, hr and br.
% gsr=0; % no global signal regression.


%% loop settings.

for high_fd = [1, 0] % 是否移除高头动帧，1保留，0移除
    for num_of_bins=[10, 20, 40]
        for atlas_i=1:3
            for regress_cofounds=[1, 0] % 1: regress, 0: skip.
                for gsr=[1, 0] % 1: with GSR, 0: no GSR.
                    loop_settings(high_fd, num_of_bins, atlas_i, regress_cofounds, gsr);
                end
            end
        end
    end
end


function loop_settings(high_fd, num_of_bins, atlas_i, regress_cofounds, gsr)
    % extract atlas name and number of regions.
    atlas_names=["schaefer200", "schaefer400", "hcpmmp"];
    atlas_rois=[200, 400, 360];
    atlas = char(atlas_names(atlas_i));
    num_of_rois = atlas_rois(atlas_i);
    % Index of upper triabgle in FC matrix.
    [u,v] = find(triu(ones(num_of_rois),1));
    idx = (v - 1)*num_of_rois + u;
    
    % extract pipeline name,
    if gsr==1
        pipeline = 'proc_regress_fixWglob_ts_mat';
    elseif gsr==0
        pipeline = 'proc_regress_fix_ts_mat';
    end
    
    % load subject ids and phenotypes.
    load('../../data/hcp/subj_ids.mat');
    phenotype_data = readtable("../../data/hcp/unrestricted_dzjin_7_30_2023_20_45_23.csv");
    
    % four sessions.
    session_names = ['REST1_LR'; 'REST1_RL'; 'REST2_LR'; 'REST2_RL'];
    
    
    % auto corrs of co-fluctuation score (and its variants). averaged across four sessions.
    cal_and_plot_bin_corr(high_fd, num_of_bins, pipeline, atlas, regress_cofounds, "ets_region_level_contribution_ratioOfMeans_all", "1_region_level_ratioOfMeans");
    cal_and_plot_bin_corr(high_fd, num_of_bins, pipeline, atlas, regress_cofounds, "ets_region_level_contribution_meanOfRatios_all", "2_region_level_meanOfRatios");
    cal_and_plot_bin_corr(high_fd, num_of_bins, pipeline, atlas, regress_cofounds, "ets_region_level_contribution_leaveOneRegionOut_all", "3_region_level_leaveOneRegionOut");
    cal_and_plot_bin_corr(high_fd, num_of_bins, pipeline, atlas, regress_cofounds, "ets_region_level_contribution_regressOneRegionOut_all", "4_region_level_regressOneRegionOut");
    cal_and_plot_bin_corr(high_fd, num_of_bins, pipeline, atlas, regress_cofounds, "ets_region_level_contribution_norm_regress_bAndBeta_all", "5_region_level_norm_regress_bAndBeta");
    cal_and_plot_bin_corr(high_fd, num_of_bins, pipeline, atlas, regress_cofounds, "ets_region_level_contribution_norm_regress_beta_all", "5_region_level_norm_regress_beta");
    cal_and_plot_bin_corr(high_fd, num_of_bins, pipeline, atlas, regress_cofounds, "ets_region_level_contribution_norm_zscore_all", "6_region_level_norm_zscore");
    cal_and_plot_bin_corr(high_fd, num_of_bins, pipeline, atlas, regress_cofounds, "ets_region_rss_bins_all", "7_region_level_means");
    
    
    % association between cs (and its variants) and sa. averaged across four sessions.
    cal_and_plot_cs_sa_association(high_fd, num_of_bins, pipeline, atlas, regress_cofounds, "ets_region_level_contribution_ratioOfMeans_all", "1_region_level_ratioOfMeans");
    cal_and_plot_cs_sa_association(high_fd, num_of_bins, pipeline, atlas, regress_cofounds, "ets_region_level_contribution_meanOfRatios_all", "2_region_level_meanOfRatios");
    cal_and_plot_cs_sa_association(high_fd, num_of_bins, pipeline, atlas, regress_cofounds, "ets_region_level_contribution_leaveOneRegionOut_all", "3_region_level_leaveOneRegionOut");
    cal_and_plot_cs_sa_association(high_fd, num_of_bins, pipeline, atlas, regress_cofounds, "ets_region_level_contribution_regressOneRegionOut_all", "4_region_level_regressOneRegionOut");
    cal_and_plot_cs_sa_association(high_fd, num_of_bins, pipeline, atlas, regress_cofounds, "ets_region_level_contribution_norm_regress_bAndBeta_all", "5_region_level_norm_regress_bAndBeta");
    cal_and_plot_cs_sa_association(high_fd, num_of_bins, pipeline, atlas, regress_cofounds, "ets_region_level_contribution_norm_regress_beta_all", "5_region_level_norm_regress_beta");
    cal_and_plot_cs_sa_association(high_fd, num_of_bins, pipeline, atlas, regress_cofounds, "ets_region_level_contribution_norm_zscore_all", "6_region_level_norm_zscore");
    cal_and_plot_cs_sa_association(high_fd, num_of_bins, pipeline, atlas, regress_cofounds, "ets_region_rss_bins_all", "7_region_level_means");
    
    
    % prepare data for gam data. averaged across four sessions.
    prepare_all_ses_data_for_gam(high_fd, num_of_bins, pipeline, atlas, regress_cofounds, "ets_region_level_contribution_ratioOfMeans_all", "1_region_level_ratioOfMeans", phenotype_data, subj_ids, num_of_rois);
    prepare_all_ses_data_for_gam(high_fd, num_of_bins, pipeline, atlas, regress_cofounds, "ets_region_level_contribution_meanOfRatios_all", "2_region_level_meanOfRatios", phenotype_data, subj_ids, num_of_rois);
    prepare_all_ses_data_for_gam(high_fd, num_of_bins, pipeline, atlas, regress_cofounds, "ets_region_level_contribution_leaveOneRegionOut_all", "3_region_level_leaveOneRegionOut", phenotype_data, subj_ids, num_of_rois);
    prepare_all_ses_data_for_gam(high_fd, num_of_bins, pipeline, atlas, regress_cofounds, "ets_region_level_contribution_regressOneRegionOut_all", "4_region_level_regressOneRegionOut", phenotype_data, subj_ids, num_of_rois);
    prepare_all_ses_data_for_gam(high_fd, num_of_bins, pipeline, atlas, regress_cofounds, "ets_region_level_contribution_norm_regress_bAndBeta_all", "5_region_level_norm_regress_bAndBeta", phenotype_data, subj_ids, num_of_rois);
    prepare_all_ses_data_for_gam(high_fd, num_of_bins, pipeline, atlas, regress_cofounds, "ets_region_level_contribution_norm_regress_beta_all", "5_region_level_norm_regress_beta", phenotype_data, subj_ids, num_of_rois);
    prepare_all_ses_data_for_gam(high_fd, num_of_bins, pipeline, atlas, regress_cofounds, "ets_region_level_contribution_norm_zscore_all", "6_region_level_norm_zscore", phenotype_data, subj_ids, num_of_rois);
    prepare_all_ses_data_for_gam(high_fd, num_of_bins, pipeline, atlas, regress_cofounds, "ets_region_rss_bins_all", "7_region_level_means", phenotype_data, subj_ids, num_of_rois);
end


%% functions.
function cal_and_plot_bin_corr(high_fd, num_of_bins, pipeline, atlas, regress_cofounds, cs_name, save_dir_name)
    % bins 20 corr colormap.
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

    mkdir(strcat("../results_bins/DRIVER/HCP/", save_dir_name));
    session_names = ['REST1_LR'; 'REST1_RL'; 'REST2_LR'; 'REST2_RL'];
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
    bins_corr_rest1_lr = corr(squeeze(mean(rest1_lr.(cs_name)(is_physio_exists(:, 1)==1, :, :), 1)));
    bins_corr_rest1_rl = corr(squeeze(mean(rest1_rl.(cs_name)(is_physio_exists(:, 2)==1, :, :), 1)));
    bins_corr_rest2_lr = corr(squeeze(mean(rest2_lr.(cs_name)(is_physio_exists(:, 3)==1, :, :), 1)));
    bins_corr_rest2_rl = corr(squeeze(mean(rest2_rl.(cs_name)(is_physio_exists(:, 4)==1, :, :), 1)));
    bins_corr = (bins_corr_rest1_lr + bins_corr_rest1_rl + bins_corr_rest2_lr + bins_corr_rest2_rl)./4;
    % plot.
    figure("Position", [100, 100, 250, 240]);
    imagesc(bins_corr);
    colormap(colormap_corrs);
    clim([-1 1]);
    set(gca, 'XTick', [1, 5, 10, 15, 20], 'TickDir', 'None', ...
        'XTickLabel', ["", "", "", "", ""]);
    set(gca, 'YTick', [1, 5, 10, 15, 20], 'TickDir', 'None', ...
        'YTickLabel', ["", "", "", "", ""]);
    ax=gca;hold on;
    ax.FontSize=13;
    save_filename=strcat("bins_corr_high_fd-",num2str(high_fd), ...
        "_num_of_bins-",num2str(num_of_bins), ...
        "_atlas-",atlas, ...
        "_regress_cofounds-",num2str(regress_cofounds), ...
        "_gsr-", pipeline);
    save_path=strcat("../results_bins/DRIVER/HCP/",save_dir_name);
    print(gcf, '-dpng', '-r1000', strcat(save_path, "/", save_filename));
    print(gcf, '-dmeta', strcat(save_path, "/", save_filename));
    close;
end


function cal_and_plot_cs_sa_association(high_fd, num_of_bins, pipeline, atlas, regress_cofounds, cs_name, save_dir_name)
    % cal group-level cs.
    session_names = ['REST1_LR'; 'REST1_RL'; 'REST2_LR'; 'REST2_RL'];
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

    roi_contri_group_rest1_lr = squeeze(mean(rest1_lr.(cs_name)(is_physio_exists(:, 1)==1, :, :), 1));
    roi_contri_group_rest1_rl = squeeze(mean(rest1_rl.(cs_name)(is_physio_exists(:, 2)==1, :, :), 1));
    roi_contri_group_rest2_lr = squeeze(mean(rest2_lr.(cs_name)(is_physio_exists(:, 3)==1, :, :), 1));
    roi_contri_group_rest2_rl = squeeze(mean(rest2_rl.(cs_name)(is_physio_exists(:, 4)==1, :, :), 1));
    roi_contri_group = roi_contri_group_rest1_lr + roi_contri_group_rest1_rl + ...
        roi_contri_group_rest2_lr + roi_contri_group_rest2_rl;
    roi_contri_group = roi_contri_group./4;
    % read sa axis.
    brain_map_sa_axis = readmatrix("./mechanism/s_a_axis_shaefer200x17.csv");
    % cal corr between cs and sa.
    coefs = nan(num_of_bins, 1);
    for bin_i=1:num_of_bins
        [coef, pval] = corr(brain_map_sa_axis, roi_contri_group(:, bin_i), "type", "Spearman");
        disp(coef);
        coefs(bin_i) = coef;
    end
    file_name = strcat(cs_name, "_sa_coefs_bins_", num2str(num_of_bins), "_bin_hcp_", pipeline, "_scan_", atlas, "_hfd_0.mat");
    save(strcat("../results_bins/DRIVER/HCP/",save_dir_name, "/", file_name), "roi_contri_group", "coefs");
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
    xticklabels([" ", " ", " ", " ", " "]);
    yticks(-1:0.2:1);
    ax=gca;hold on;
    ax.FontSize=13;
    ax.XTickLabelRotation=45;
    set(gca,'Box','off');
    save_filename=strcat("bins_cs_sa_high_fd-",num2str(high_fd), ...
        "_num_of_bins-",num2str(num_of_bins), ...
        "_atlas-",atlas, ...
        "_regress_cofounds-",num2str(regress_cofounds), ...
        "_gsr-", pipeline);
    save_path=strcat("../results_bins/DRIVER/HCP/",save_dir_name);
    print(gcf, '-dpng', '-r1000', strcat(save_path, "/", save_filename));
    print(gcf, '-dmeta', strcat(save_path, "/", save_filename));
    close;
end


function prepare_all_ses_data_for_gam(high_fd, num_of_bins, pipeline, atlas, regress_cofounds, cs_name, save_dir_name, phenotype_data, subj_ids, num_of_rois)
    % load cs data.
    % cs averaged scross four sessions for each subject.
    session_names = ['REST1_LR'; 'REST1_RL'; 'REST2_LR'; 'REST2_RL'];
    rest1_lr = load(['../outputs_bins/HCP/bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(1, :) '_' atlas '_hfd_' num2str(high_fd) '_regressCofounds_' num2str(regress_cofounds) '.mat'], ...
        cs_name, 'ets_global_rss_bins_all', 'fd_bins_all', 'gs_bins_all', 'br_bins_all', 'hr_bins_all');
    rest1_rl = load(['../outputs_bins/HCP/bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(2, :) '_' atlas '_hfd_' num2str(high_fd) '_regressCofounds_' num2str(regress_cofounds) '.mat'], ...
        cs_name, 'ets_global_rss_bins_all', 'fd_bins_all', 'gs_bins_all', 'br_bins_all', 'hr_bins_all');
    rest2_lr = load(['../outputs_bins/HCP/bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(3, :) '_' atlas '_hfd_' num2str(high_fd) '_regressCofounds_' num2str(regress_cofounds) '.mat'], ...
        cs_name, 'ets_global_rss_bins_all', 'fd_bins_all', 'gs_bins_all', 'br_bins_all', 'hr_bins_all');
    rest2_rl = load(['../outputs_bins/HCP/bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(4, :) '_' atlas '_hfd_' num2str(high_fd) '_regressCofounds_' num2str(regress_cofounds) '.mat'], ...
        cs_name, 'ets_global_rss_bins_all', 'fd_bins_all', 'gs_bins_all', 'br_bins_all', 'hr_bins_all');
    % load phsio data.
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
    % load tSNR.
    load(strcat("../results_bins/tSNR/hcp3t_unrelated100_tSNR_", atlas, ".mat"));
    tSNR_avg4Sess = squeeze(mean(tSNR_all, 2, 'omitmissing'));

    cs_all = nan(4, 100, 200, 20);
    cs_all(1, :, :, :) = rest1_lr.(cs_name);
    cs_all(2, :, :, :) = rest1_rl.(cs_name);
    cs_all(3, :, :, :) = rest2_lr.(cs_name);
    cs_all(4, :, :, :) = rest2_rl.(cs_name);
    cs_avg4Sess = squeeze(mean(cs_all, 1, 'omitmissing'));
    
    bins_rss_ALL = nan(4, 100, 20);
    bins_rss_ALL(1, :, :) = transpose(zscore(rest1_lr.ets_global_rss_bins_all.'));
    bins_rss_ALL(2, :, :) = transpose(zscore(rest1_rl.ets_global_rss_bins_all.'));
    bins_rss_ALL(3, :, :) = transpose(zscore(rest2_lr.ets_global_rss_bins_all.'));
    bins_rss_ALL(4, :, :) = transpose(zscore(rest2_rl.ets_global_rss_bins_all.'));
    bins_rss_avg4Sess = squeeze(mean(bins_rss_ALL, 1, 'omitmissing'));
    
    mean_fd_ALL = nan(4, 100, 20);
    mean_fd_ALL(1, :, :) = transpose(zscore(rest1_lr.fd_bins_all.'));
    mean_fd_ALL(2, :, :) = transpose(zscore(rest1_rl.fd_bins_all.'));
    mean_fd_ALL(3, :, :) = transpose(zscore(rest2_lr.fd_bins_all.'));
    mean_fd_ALL(4, :, :) = transpose(zscore(rest2_rl.fd_bins_all.'));
    mean_fd_avg4Sess = squeeze(mean(mean_fd_ALL, 1, 'omitmissing'));

    mean_gs_ALL = nan(4, 100, 20);
    mean_gs_ALL(1, :, :) = transpose(zscore(rest1_lr.gs_bins_all.'));
    mean_gs_ALL(2, :, :) = transpose(zscore(rest1_rl.gs_bins_all.'));
    mean_gs_ALL(3, :, :) = transpose(zscore(rest2_lr.gs_bins_all.'));
    mean_gs_ALL(4, :, :) = transpose(zscore(rest2_rl.gs_bins_all.'));
    mean_gs_avg4Sess = squeeze(mean(mean_gs_ALL, 1, 'omitmissing'));

    mean_br_ALL = nan(4, 100, 20);
    mean_br_ALL(1, :, :) = transpose(zscore(rest1_lr.br_bins_all.'));
    mean_br_ALL(2, :, :) = transpose(zscore(rest1_rl.br_bins_all.'));
    mean_br_ALL(3, :, :) = transpose(zscore(rest2_lr.br_bins_all.'));
    mean_br_ALL(4, :, :) = transpose(zscore(rest2_rl.br_bins_all.'));
    mean_br_avg4Sess = squeeze(mean(mean_br_ALL, 1, 'omitmissing'));

    mean_hr_ALL = nan(4, 100, 20);
    mean_hr_ALL(1, :, :) = transpose(zscore(rest1_lr.hr_bins_all.'));
    mean_hr_ALL(2, :, :) = transpose(zscore(rest1_rl.hr_bins_all.'));
    mean_hr_ALL(3, :, :) = transpose(zscore(rest2_lr.hr_bins_all.'));
    mean_hr_ALL(4, :, :) = transpose(zscore(rest2_rl.hr_bins_all.'));
    mean_hr_avg4Sess = squeeze(mean(mean_hr_ALL, 1, 'omitmissing'));
    
    % Contribtions of each ROI, bin_label, roi_label, sex
    gam_data = nan(100*num_of_bins, num_of_rois+7);
    idx = 1;
    for subj_i=1:100
        if sum(is_physio_exists(subj_i, :))==0
            continue;
        end
        bins_rss_subj_i = bins_rss_avg4Sess(subj_i, :);
        for bin_i=1:num_of_bins
            for roi_i=1:num_of_rois
                % Contribution
                gam_data(idx, roi_i) = cs_avg4Sess(subj_i, roi_i, bin_i);            
            end
            % averaged_and_zscored_raw_rss
            gam_data(idx, num_of_rois+1) = bins_rss_subj_i(bin_i);
            % sex_label
            subj_i_row_id = phenotype_data.('Subject') == str2num(subj_ids(subj_i, :));
            % if phenotype_data.("Gender")(subj_i_row_id)=='M'
            if phenotype_data{subj_i_row_id, "Gender"}{1,1}=='M'
                gam_data(idx, num_of_rois+2) = 0;
            else
                gam_data(idx, num_of_rois+2) = 1;
            end
            % subj_label
            gam_data(idx, num_of_rois+3) = subj_i;
            % mean fd, gs, br and hr with each bin.
            gam_data(idx, num_of_rois+4) = mean_fd_avg4Sess(subj_i, bin_i);
            gam_data(idx, num_of_rois+5) = mean_gs_avg4Sess(subj_i, bin_i);
            gam_data(idx, num_of_rois+6) = mean_br_avg4Sess(subj_i, bin_i);
            gam_data(idx, num_of_rois+7) = mean_hr_avg4Sess(subj_i, bin_i);
            % age
            age_char_arr = phenotype_data{subj_i_row_id, "Age"}{1,1};
            gam_data(idx, num_of_rois+8) = str2double(age_char_arr(1:2));
            gam_data(idx, num_of_rois+9) = tSNR_avg4Sess(subj_i, roi_i);
            idx = idx + 1;
        end
    end
    gam_data(idx:100*num_of_bins, :)=[];
    % save data in results_bins.
    save_filename=strcat("bins_gam_data_zscored_bins_rss_high_fd-",num2str(high_fd), ...
        "_num_of_bins-",num2str(num_of_bins), ...
        "_atlas-",atlas, ...
        "_regress_cofounds-",num2str(regress_cofounds), ...
        "_gsr-", pipeline, ".csv");
    save_path=strcat("../results_bins/DRIVER/HCP/",save_dir_name);
    writematrix(gam_data, strcat(save_path, "/", save_filename));

    % mkdir(strcat("../results_bins/DRIVER/HCP/",save_dir_name, "/", pipeline));
    % schaefer200x17_names = readtable("G:/datasets/atlas/Schaefer/Schaefer2018_200Parcels_17Networks_order_FSLMNI152_1mm.Centroid_RAS.csv");
    % schaefer200x17_names = schaefer200x17_names.ROIName;
    % schaefer200x17_names = strrep(schaefer200x17_names, "17Networks_", "");
    % schaefer200x17_names = strrep(schaefer200x17_names, "_", "-");
    % for roi_i=1:100
    %     figure('Position', [100 100 500 500]);
    %     scatter(gam_data(:, 201), gam_data(:, roi_i));
    %     title(schaefer200x17_names(roi_i));
    %     xlabel("GLOBAL RSS (zscore)");
    %     ylabel("REGIONAL");
    %     a = round(min(min(gam_data(:, 1:200))), 1);
    %     b = round(max(max(gam_data(:, 1:200))), 1);
    %     ylim([round(min(min(gam_data(:, 1:200))), 2), round(max(max(gam_data(:, 1:200))), 2)]);
    %     ax=gca;hold on;
    %     ax.FontSize=20;
    %     file_name=strcat("hcp_all_sess_regional_global_rss_raw_bins_", num2str(num_of_bins), "_bin_hcp_", pipeline, "_scan_", atlas, "_hfd_0_region_", num2str(roi_i));
    %     save_path_filename = strcat("../results_bins/DRIVER/HCP/",save_dir_name, "/", pipeline, "/", file_name);
    %     print(gcf, '-dpng', '-r600', save_path_filename);
    %     % print(gcf, '-dmeta', save_path_filename);
    %     close;
    % end
end


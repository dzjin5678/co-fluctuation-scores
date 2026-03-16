clear
close all
clc
addpath(genpath('fcn'));

high_fd = 0; % 是否移除高头动帧，1保留，0移除
load("../data/hcp7t/hcp7t_subj_sessions_all.mat", "subj_sessions_all");
session_names=["rfMRI_REST1_7T_PA","rfMRI_REST2_7T_AP","rfMRI_REST3_7T_PA","rfMRI_REST4_7T_AP", ...
    "tfMRI_MOVIE1_7T_AP","tfMRI_MOVIE2_7T_PA","tfMRI_MOVIE3_7T_PA","tfMRI_MOVIE4_7T_AP"];
session_timepoints=[900,900,900,900,921,918,915,901];
root_path = "/mnt/system_v2022/my_data/dzjin_data/hcp7t";
atlas_names=["schaefer200x17", "schaefer400x17", "hcpmmp"];
atlas_rois=[200, 400, 360];

for session_idx=1:8
    subj_sessions_idx = subj_sessions_all(:, session_idx)==1;
    subj_list = subj_sessions_all(subj_sessions_idx, 9);
    num_of_subjects = size(subj_list, 1);

    session_name = session_names(session_idx);
    num_of_timepoints = session_timepoints(session_idx);
    for num_of_bins=[10, 20, 40]
        for atlas_i=1:3
            atlas=atlas_names(atlas_i);
            num_of_rois=atlas_rois(atlas_i);
            % Index of upper triabgle in FC matrix.
            [u,v] = find(triu(ones(num_of_rois),1));
            idx = (v - 1)*num_of_rois + u;

            %% define variable.
            % zscored bold ts.
            bold_zscored_all = nan(num_of_subjects, num_of_rois, num_of_timepoints);
            % bin-level ets.
            ets_mean_bins_all = nan(num_of_subjects, num_of_bins, size(idx, 1));

            % the rss timeseries of whole-connectome.
            ets_global_rss_all = nan(num_of_subjects, num_of_timepoints);
            ets_global_rss_bins_all = nan(num_of_subjects, num_of_bins);
            % the rss timeseries of each region (co-fluc with other regions).
            ets_region_rss_all = nan(num_of_subjects, num_of_rois, num_of_timepoints);
            ets_region_rss_bins_all = nan(num_of_subjects, num_of_rois, num_of_bins);
            high_fd_label_all = nan(num_of_subjects, num_of_timepoints);
            % the contribution of each region (co-fluc with other regions) to whole-brain co-fluc amp.
            ets_region_level_contribution_ratioOfMeans_all = nan(num_of_subjects, num_of_rois, num_of_bins);
            ets_region_level_contribution_meanOfRatios_all = nan(num_of_subjects, num_of_rois, num_of_bins);        
            % the normalization of regional level co-fluctuation.
            ets_region_level_contribution_norm_zscore_all = nan(num_of_subjects, num_of_rois, num_of_bins);
            ets_region_level_contribution_norm_regress_bAndBeta_all = nan(num_of_subjects, num_of_rois, num_of_bins);
            ets_region_level_contribution_norm_regress_beta_all = nan(num_of_subjects, num_of_rois, num_of_bins);
            % leave one region out.
            ets_region_level_contribution_leaveOneRegionOut_all = nan(num_of_subjects, num_of_rois, num_of_bins);
            % regress one region out.
            ets_region_level_contribution_regressOneRegionOut_all = nan(num_of_subjects, num_of_rois, num_of_bins);

            % local-global lead-lag analysis.
            local_global_lead_lag_all = nan(num_of_subjects, num_of_rois, 2);
            % eta analysis.
            template = struct( ...
                'trigIdx', [], ...
                'nEvents', [], ...
                'threshold', [], ...
                'tau', [], ...
                'X', [], ...
                'etaMean', [], ...
                'etaSEM', [] ...
            );
            eta_out_all = repmat(template, num_of_subjects, 1);

            %% main function.
            for subj_i=1:size(subj_list, 1)
                subj_id=num2str(subj_list(subj_i));
                disp(subj_id);
                load(strcat(root_path, "/hcp7t_roi_ts_", atlas, "/", subj_id, "_", session_name, "_roi_ts_", atlas, ".mat"), "roi_ts");
                timeseries=roi_ts;
                if high_fd == 0
                    fds = readmatrix(strcat(root_path, "/hcp7t_fd/", subj_id, "_", session_name, "_fwd_abssum.txt"));
                    high_fd_label = (fds > 0.2);
                    high_fd_label_all(subj_i, 1:size(timeseries, 2)) = high_fd_label;
                    timeseries_no_high_fd = timeseries;
                    timeseries_no_high_fd(:, high_fd_label) = [];
                    timeseries = timeseries_no_high_fd;
                    timeseries = timeseries.';
                end
                bin_size = floor(size(timeseries, 1)/num_of_bins);
                timeseries = zscore(timeseries);
                bold_zscored_all(subj_i, :, 1:size(timeseries, 1))=timeseries.';
                ets = timeseries(:, u).*timeseries(:, v);
                ets_rss = sum(ets.^2, 2).^0.5;
                ets_rss_orig = ets_rss;
                % for leave one region out.
                ets_rss_r2 = sum(ets.^2, 2);
                
                ets_global_rss_all(subj_i, 1:size(timeseries, 1)) = ets_rss;
                [~,idxsort] = sort(ets_rss, 'descend');

                % ets mean bins.
                for bin_i=1:num_of_bins
                    bin_idx = idxsort((bin_i-1)*bin_size+1:bin_i*bin_size);
                    ets_mean_bin_i = squeeze(mean(ets(bin_idx, :)));
                    ets_mean_bins_all(subj_i, bin_i, :) = ets_mean_bin_i;
                end

                ets_region_rss_subj_i = nan(num_of_rois, size(timeseries, 1));
                % CS. (ratio of means)
                for roi_i=1:num_of_rois
                    edges_related_roi_i = (u==roi_i) | (v==roi_i);
                    edges_idx_related_roi_i = find(edges_related_roi_i==1);
                    ets_related_roi_i = ets(:, edges_idx_related_roi_i);
                    % region level RSS.
                    ets_related_roi_i_rss = sum(ets_related_roi_i.^2, 2).^0.5;
                    ets_region_rss_all(subj_i, roi_i, 1:size(timeseries, 1)) = ets_related_roi_i_rss;
                    ets_region_rss_subj_i(roi_i, :) = ets_related_roi_i_rss;
                    for bin_i=1:num_of_bins
                        bin_idx = idxsort((bin_i-1)*bin_size+1:bin_i*bin_size);
                        rss_d = mean(ets_rss(bin_idx));
                        rss_n = mean(ets_related_roi_i_rss(bin_idx));
                        ets_region_level_contribution_ratioOfMeans_all(subj_i, roi_i, bin_i) = rss_n/rss_d;
                        ets_global_rss_bins_all(subj_i, bin_i)=rss_d;
                        ets_region_rss_bins_all(subj_i, roi_i, bin_i)=rss_n;
                    end

                    % local-global lead-lag analysis.
                    [peak_lag, peak_corr] = lead_lag_analysis(ets_related_roi_i_rss, ets_rss);
                    local_global_lead_lag_all(subj_i, roi_i, 1)=peak_lag;
                    local_global_lead_lag_all(subj_i, roi_i, 2)=peak_corr;
                end
                ets_region_rss_subj_i = ets_region_rss_subj_i.';
                % Inputs:
                % rssG: [T x 1]
                % rssN: [T x K] (e.g., K=7 or 17 networks)
                preTR = 20;
                postTR = 20;
                topPct = 95;      % top 5%
                minSepTR = 5;     % minimum separation
                doZscore = true;
                out = eta_network_rss(ets_rss, ets_region_rss_subj_i, preTR, postTR, topPct, minSepTR, doZscore);
                eta_out_all(subj_i) = out;

                % CS. (mean of ratios).
                for roi_i=1:num_of_rois
                    % region level RSS.
                    ets_related_roi_i_rss = squeeze(ets_region_rss_all(subj_i, roi_i, 1:size(timeseries, 1)));
                    for bin_i=1:num_of_bins
                        bin_idx = idxsort((bin_i-1)*bin_size+1:bin_i*bin_size);
                        local_global_ratios = ets_related_roi_i_rss(bin_idx) ./ ets_rss(bin_idx);
                        ets_region_level_contribution_meanOfRatios_all(subj_i, roi_i, bin_i) = mean(local_global_ratios);
                    end
                end

                % CS. (zscore).
                for bin_i=1:num_of_bins
                    bin_idx = idxsort((bin_i-1)*bin_size+1:bin_i*bin_size);
                    % region level RSS within Bin i.
                    ets_related_bin_i_roi_rss = squeeze(ets_region_rss_all(subj_i, :, bin_idx));
                    ets_related_bin_i_roi_rss = mean(ets_related_bin_i_roi_rss, 2);
                    ets_region_level_contribution_norm_zscore_all(subj_i, :, bin_i) = zscore(ets_related_bin_i_roi_rss);
                end

                % CS. (regress, both b and beta).
                ets_related_roi_rss = squeeze(ets_region_rss_all(subj_i, :, 1:size(timeseries, 1)));
                ets_related_roi_rss = ets_related_roi_rss.';
                Xcov = [ones(size(ets_rss, 1), 1), ets_rss];
                beta = Xcov \ ets_related_roi_rss;
                ets_related_roi_rss_regress = ets_related_roi_rss - Xcov * beta;
                ets_related_roi_rss_regress = ets_related_roi_rss_regress.';
                for roi_i=1:num_of_rois
                    ets_related_roi_i_rss = squeeze(ets_related_roi_rss_regress(roi_i, 1:size(timeseries, 1)));
                    for bin_i=1:num_of_bins
                        bin_idx = idxsort((bin_i-1)*bin_size+1:bin_i*bin_size);
                        rss_n = mean(ets_related_roi_i_rss(bin_idx));
                        ets_region_level_contribution_norm_regress_bAndBeta_all(subj_i, roi_i, bin_i) = rss_n;
                    end
                end

                % CS. (regress, only beta).
                ets_related_roi_rss = squeeze(ets_region_rss_all(subj_i, :, 1:size(timeseries, 1)));
                ets_related_roi_rss = ets_related_roi_rss.';
                Xcov = [ones(size(ets_rss, 1), 1), ets_rss];
                beta = Xcov \ ets_related_roi_rss;
                ets_related_roi_rss_regress = ets_related_roi_rss - Xcov(:,2) * beta(2,:);
                ets_related_roi_rss_regress = ets_related_roi_rss_regress.';
                for roi_i=1:num_of_rois
                    ets_related_roi_i_rss = squeeze(ets_related_roi_rss_regress(roi_i, 1:size(timeseries, 1)));
                    for bin_i=1:num_of_bins
                        bin_idx = idxsort((bin_i-1)*bin_size+1:bin_i*bin_size);
                        rss_n = mean(ets_related_roi_i_rss(bin_idx));
                        ets_region_level_contribution_norm_regress_beta_all(subj_i, roi_i, bin_i) = rss_n;
                    end
                end

                % CS. (leave one region out).
                for roi_i=1:num_of_rois
                    edges_related_roi_i = (u==roi_i) | (v==roi_i);
                    edges_idx_related_roi_i = find(edges_related_roi_i==1);
                    ets_related_roi_i = ets(:, edges_idx_related_roi_i);
                    % region level RSS.
                    ets_related_roi_i_rss_r2 = sum(ets_related_roi_i.^2, 2);
                    for bin_i=1:num_of_bins
                        bin_idx = idxsort((bin_i-1)*bin_size+1:bin_i*bin_size);
                        rss_d = mean(ets_rss_r2(bin_idx));
                        rss_n = mean(ets_related_roi_i_rss_r2(bin_idx));
                        ets_region_level_contribution_leaveOneRegionOut_all(subj_i, roi_i, bin_i) = log(rss_n/(rss_d-rss_n));
                        ets_global_rss_r2_bins_all(subj_i, bin_i)=rss_d;
                        ets_region_rss_r2_bins_all(subj_i, roi_i, bin_i)=rss_n;
                    end
                end

                % CS. (regress one region out).
                for roi_i=1:num_of_rois
                    ets_related_roi_i_rss = squeeze(ets_region_rss_all(subj_i, roi_i, 1:size(timeseries, 1)));
                    % regress local on global.
                    Xcov = [ones(size(ets_related_roi_i_rss, 1), 1), ets_related_roi_i_rss];
                    beta = Xcov \ ets_rss;
                    ets_rss_regress_local = ets_rss - Xcov(:,2) * beta(2,:);
                    ets_rss_regress_local = ets_rss_regress_local.';
                    for bin_i=1:num_of_bins
                        bin_idx = idxsort((bin_i-1)*bin_size+1:bin_i*bin_size);
                        rss_d = mean(ets_rss_regress_local(bin_idx));
                        rss_n = mean(ets_related_roi_i_rss(bin_idx));
                        ets_region_level_contribution_regressOneRegionOut_all(subj_i, roi_i, bin_i) = rss_n/rss_d;
                    end
                end

            end
            %% save results.
            file_name=strcat("./outputs_bins/HCP7T/bins_", num2str(num_of_bins), "_session_", session_name, "_", atlas, "_hfd_", num2str(high_fd), ".mat");
            save(file_name, ...
                'subj_list', ...
                'bold_zscored_all', ...
                'ets_mean_bins_all', ...
                'ets_global_rss_all', ...
                'ets_global_rss_bins_all', ...
                'ets_region_rss_all', ...
                'ets_region_rss_bins_all', ...
                'high_fd_label_all', ...
                'ets_region_level_contribution_ratioOfMeans_all', ...
                'ets_region_level_contribution_meanOfRatios_all', ...
                'ets_region_level_contribution_norm_zscore_all', ...
                'ets_region_level_contribution_norm_regress_beta_all', ...
                'ets_region_level_contribution_norm_regress_bAndBeta_all', ...
                'ets_region_level_contribution_leaveOneRegionOut_all', ...
                'ets_region_level_contribution_regressOneRegionOut_all', ...
                'local_global_lead_lag_all', 'eta_out_all', ...
                '-v7.3');
        end
    end
end

function [peak_lag, peak_corr] = lead_lag_analysis(rss_local, rss_global)
    max_lag = 20;
    local_z = zscore(rss_local);
    global_z = zscore(rss_global);
    [xc,lags] = xcorr(local_z, global_z, max_lag, 'coeff');
    [~, idx] = max(abs(xc));
    peak_lag = lags(idx);
    peak_corr = xc(idx);
    % disp(['Peak lag: ', num2str(peak_lag)]);
    % disp(['Peak corr: ', num2str(peak_corr)]);
end


function out = eta_network_rss(rssG, rssN, preTR, postTR, topPct, minSepTR, doZscore)
    % ETA of network-level RSS aligned to high-amplitude RSS_GLOBAL events.
    %
    % Inputs:
    %   rssG      : [T x 1] global RSS time series
    %   rssN      : [T x K] network-level RSS time series (K networks)
    %   preTR     : number of TRs before event (e.g., 20)
    %   postTR    : number of TRs after event (e.g., 20)
    %   topPct    : trigger percentile (e.g., 95 means top 5% events)
    %   minSepTR  : minimum separation between events in TRs (e.g., 5 or 10)
    %   doZscore  : true/false, whether to z-score rssG and rssN
    %
    % Output struct 'out' contains triggers, ETA, etc.
    
    rssG = rssG(:);
    [T, K] = size(rssN);

    if doZscore
        rssGz = zscore(rssG);
        rssNz = zscore(rssN);
    else
        rssGz = rssG;
        rssNz = rssN;
    end

    % 1) define threshold
    thr = prctile(rssGz, topPct);

    % 2) find candidate points above threshold
    above = rssGz >= thr;

    % 3) keep local maxima among above-threshold points
    % local maxima: rssGz(t) > rssGz(t-1) and rssGz(t) >= rssGz(t+1)
    % handle boundaries safely
    isPeak = false(T,1);
    for t = 2:T-1
        if above(t) && rssGz(t) > rssGz(t-1) && rssGz(t) >= rssGz(t+1)
            isPeak(t) = true;
        end
    end
    candIdx = find(isPeak);

    % 4) enforce minimum separation between events
    if isempty(candIdx)
        warning('No events found. Consider lowering threshold (topPct) or checking rssG.');
        out = struct();
        return;
    end

    % Greedy selection: keep peaks in descending amplitude order, reject neighbors within minSepTR
    [~, order] = sort(rssGz(candIdx), 'descend');
    selected = [];
    for ii = 1:numel(order)
        t0 = candIdx(order(ii));
        if isempty(selected) || all(abs(selected - t0) > minSepTR)
            selected(end+1,1) = t0; %#ok<AGROW>
        end
    end
    trigIdx = sort(selected);

    % 5) remove triggers too close to edges for window extraction
    good = trigIdx > preTR & trigIdx <= (T - postTR);
    trigIdx = trigIdx(good);

    nEvt = numel(trigIdx);
    if nEvt == 0
        warning('All events too close to edges. Reduce preTR/postTR or check triggers.');
        out = struct();
        return;
    end

    % 6) extract windows
    winLen = preTR + postTR + 1;
    X = nan(winLen, nEvt, K);
    for i = 1:nEvt
        t0 = trigIdx(i);
        seg = rssNz((t0-preTR):(t0+postTR), :); % [winLen x K]
        X(:, i, :) = seg;
    end

    % 7) ETA mean and SEM across events
    etaMean = squeeze(mean(X, 2, 'omitnan')); % [winLen x K]
    etaSEM  = squeeze(std(X, 0, 2, 'omitnan')) / sqrt(nEvt); % [winLen x K]

    % time axis in TR units (0 = event time)
    tau = (-preTR:postTR)';

    out = struct();
    out.trigIdx  = trigIdx;
    out.nEvents  = nEvt;
    out.threshold= thr;
    out.tau      = tau;
    out.X        = X;        % [winLen x nEvt x K]
    out.etaMean  = etaMean;  % [winLen x K]
    out.etaSEM   = etaSEM;   % [winLen x K]
end


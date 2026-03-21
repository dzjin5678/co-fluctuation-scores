clear
close all
clc
addpath(genpath('fcn'));

num_of_bins = 20;
root_path = '/mnt/system_v2022/my_data/dzjin_data/hcp';
high_fd = 1;    % 是否移除高头动帧，1保留，0移除

% load subject id and atlas data.
load('../data/hcp/subj_ids.mat');
% load respiration and heart rate data.
load("./outputs_bins/HCP/physio_phase/BOLD_PLVs_proc_regress_fix_ts_mat_schaefer200.mat", "physios_ds", "physios_avg1");


atlas_names=["schaefer200"];
atlas_rois=[200];
scan_names = ["REST1_LR", "REST1_RL", "REST2_LR", "REST2_RL"];
%% settings.
for atlas_i=1:1
    atlas = char(atlas_names(atlas_i));
    num_of_rois = atlas_rois(atlas_i);
    % Index of upper triabgle in FC matrix.
    [u,v] = find(triu(ones(num_of_rois),1));
    idx = (v - 1)*num_of_rois + u;

    for gsr=[1, 0]
        % gsr = 1; % 1: with GSR, 0: no GSR.
        if gsr==1
            dataset = 'proc_regress_fixWglob_ts_mat';
        elseif gsr==0
            dataset = 'proc_regress_fix_ts_mat';
        end

        for scan_i=1:4
            scan = char(scan_names(scan_i));
            
            %% define variable.
            % ts phase shift.
            ets_ps_global_rss_all = nan(100, 1200);
            ets_ps_region_rss_all = nan(100, num_of_rois, 1200);
            ets_ps_global_rss_bins_all = nan(100, num_of_bins);
            ets_ps_region_rss_bins_all = nan(100, num_of_rois, num_of_bins);
            ets_ps_region_level_contribution_ratioOfMeans_all = nan(100, num_of_rois, num_of_bins);

            %% main function.
            for subj_i=1:100
                disp(['Index of subject: ' num2str(subj_i) ' !']);
                % skip scan without physio data.
                if sum(isnan(physios_ds(:, subj_i, scan_i)))>1 | sum(isnan(physios_avg1(:, subj_i, scan_i)))>1
                    continue;
                end
                load([root_path '/' dataset '_' atlas '/hcp_' subj_ids(subj_i, :) '_rfMRI_' scan '_' atlas '.mat']);
                timeseries = ts_roi.';

                % phase shift.
                timeseries_phase_shift = zeros(size(timeseries));
                for r=1:num_of_rois
                    timeseries_phase_shift(:, r) = phase_randomize(timeseries(:, r));
                end

                bin_size = floor(size(timeseries, 1)/num_of_bins);
                timeseries_p_s = zscore(timeseries_phase_shift);
                % from bold to ets and rss.
                bin_size_ps = floor(size(timeseries_p_s, 1)/num_of_bins);
                timeseries_p_s = zscore(timeseries_p_s);
                ets_ps = timeseries_p_s(:, u).*timeseries_p_s(:, v);
                ets_rss_ps = sum(ets_ps.^2, 2).^0.5;
                ets_ps_global_rss_all(subj_i, 1:size(timeseries_p_s, 1)) = ets_rss_ps;
                [~,idxsort_ps] = sort(ets_rss_ps, 'descend');
                % CS. (ratio of means).
                for roi_i=1:num_of_rois
                    edges_related_roi_i = (u==roi_i) | (v==roi_i);
                    edges_idx_related_roi_i = find(edges_related_roi_i==1);
                    ets_ps_related_roi_i = ets_ps(:, edges_idx_related_roi_i);
                    % region level RSS.
                    ets_ps_related_roi_i_rss = sum(ets_ps_related_roi_i.^2, 2).^0.5;
                    ets_ps_region_rss_all(subj_i, roi_i, 1:size(timeseries_p_s, 1)) = ets_ps_related_roi_i_rss;
                    for bin_i=1:num_of_bins
                        bin_idx = idxsort_ps((bin_i-1)*bin_size+1:bin_i*bin_size);
                        rss_d = mean(ets_rss_ps(bin_idx));
                        rss_n = mean(ets_ps_related_roi_i_rss(bin_idx));
                        ets_ps_region_level_contribution_ratioOfMeans_all(subj_i, roi_i, bin_i) = rss_n/rss_d;
                        ets_ps_global_rss_bins_all(subj_i, bin_i)=rss_d;
                        ets_ps_region_rss_bins_all(subj_i, roi_i, bin_i)=rss_n;
                    end
                end
            end
            
            %% save results.            
            file_name = ['./outputs_bins/HCP/phase_random/phase_shift_1_bins_' num2str(num_of_bins) '_' dataset '_scan_', scan , '_', atlas, '_hfd_' num2str(high_fd) '.mat'];
            save(file_name, 'ets_ps_global_rss_all', 'ets_ps_region_rss_all', 'ets_ps_global_rss_bins_all', 'ets_ps_region_rss_bins_all', ...
                'ets_ps_region_level_contribution_ratioOfMeans_all');
        end
    end
end


function x_rand = phase_randomize(x)
    % PHASE_RANDOMIZE  Phase randomization of a real-valued time series
    % 保持功率谱，随机相位，返回实值 surrogate
    
    x = x(:);
    N = length(x);
    
    % 去均值（推荐）
    x = x - mean(x);
    
    % FFT
    X = fft(x);
    A = abs(X);
    
    % 初始化随机相位
    phi_rand = zeros(N,1);
    
    % 正频率索引
    idx = 2:floor(N/2);
    
    % 随机相位
    phi_rand(idx) = 2*pi*rand(length(idx),1);
    
    % 构造共轭对称
    phi_rand(N-idx+2) = -phi_rand(idx);
    
    % DC 分量
    phi_rand(1) = angle(X(1));
    
    % Nyquist（仅当 N 为偶数）
    if mod(N,2) == 0
        phi_rand(N/2+1) = angle(X(N/2+1));
    end
    
    % 重构频谱
    X_rand = A .* exp(1i * phi_rand);
    
    % 逆变换
    x_rand = real(ifft(X_rand));
    
end
    

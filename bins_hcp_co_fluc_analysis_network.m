clear
close all
clc
addpath(genpath('fcn'));

num_of_bins = 20;
root_path = '/mnt/system_v2022/my_data/dzjin_data/hcp';
high_fd = 1;    % 是否移除高头动帧，1保留，0移除

% load subject id and atlas data.
load('../data/hcp/subj_ids.mat');
load("../atlas/Schaefer200x17.mat");
atlas_names=["schaefer200"];
atlas_rois=[200];
scan_names = ["REST1_LR", "REST1_RL", "REST2_LR", "REST2_RL"];
for atlas_i=1:1
    atlas = char(atlas_names(atlas_i));
    num_of_rois = atlas_rois(atlas_i);
    % Index of upper triabgle in FC matrix.
    [u,v] = find(triu(ones(num_of_rois),1));
    idx = (v - 1)*num_of_rois + u;

    for gsr=[1, 0]
        % 1: with GSR, 0: no GSR.
        if gsr==1
            dataset = 'proc_regress_fixWglob_ts_mat_';
        elseif gsr==0
            dataset = 'proc_regress_fix_ts_mat_';
        end

        for scan_i=1:4
            scan = char(scan_names(scan_i));
            % define variable.
            % the rss timeseries of whole-connectome.
            ets_global_rss_all = nan(100, 1191);
            % the rss timeseries of each network (co-fluc with other networks).
            ets_network_rss_all = nan(100, 16, 1191);
            % main function.
            for subj_i=1:100
                disp(['Index of subject: ' num2str(subj_i) ' !']);
                load([root_path '/' dataset atlas '/hcp_' subj_ids(subj_i, :) '_rfMRI_' scan '_' atlas '.mat']);
                timeseries = ts_roi.';
                timeseries = timeseries(5:end-5, :);
            
                % scrap frames with high head motion.
                confounds = readmatrix([root_path '/regressors_cp/' subj_ids(subj_i, :) '/rfMRI_' scan '/' subj_ids(subj_i, :) '_rfMRI_' scan '_confounds.csv']);
                if high_fd == 0
                    fd = confounds(:, 7);fd = fd(5:end-5, :);
                    high_fd_label = (fd > 0.2);
                    timeseries(high_fd_label, :) = [];
                end
                % from bold to ets and rss.
                bin_size = floor(size(timeseries, 1)/num_of_bins);
                timeseries = zscore(timeseries);
                ets = timeseries(:, u).*timeseries(:, v);
                ets_rss = sum(ets.^2, 2).^0.5;
                ets_global_rss_all(subj_i, 1:size(timeseries, 1)) = ets_rss;
                
                % network level RSS.
                for network_i=1:16
                    roi_within_network_i=find(lab16==network_i);
                    [~, idx_u]=ismember(u, roi_within_network_i);
                    [~, idx_v]=ismember(v, roi_within_network_i);
                    edges_related_network_i = (idx_u~=0) | (idx_v~=0);
                    edges_related_network_i = find(edges_related_network_i==1);
                    ets_related_network_i = ets(:, edges_related_network_i);
                    % network level RSS.
                    ets_related_network_i_rss = sum(ets_related_network_i.^2, 2).^0.5;
                    ets_network_rss_all(subj_i, network_i, 1:size(timeseries, 1)) = ets_related_network_i_rss;
                end
            end
            % save results.
            file_name = ['./outputs_bins/HCP/N_G_RSS_' num2str(num_of_bins) '_' dataset '_scan_', scan , '_', atlas, '_hfd_' num2str(high_fd) '.mat'];
            save(file_name, 'ets_global_rss_all', 'ets_network_rss_all');
        end
    end
end


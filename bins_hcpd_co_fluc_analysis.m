addpath(genpath('fcn'));
addpath(genpath('/usr/local/R2021a/toolbox/dzjin/2019_03_03_BCT'));

root_path="/mnt/system_v2022/my_data/dzjin_data/hcpd";
atlas_names=["schaefer200x17", "schaefer400x17", "glasser360"];
atlas_rois=[200, 400, 360];

for num_of_bins=[10, 20]
    for fd_threshold=[2,3]
        for atlas_i=1:3
            atlas = atlas_names(atlas_i);
            num_of_rois = atlas_rois(atlas_i);
            % Index of upper triabgle in FC matrix.
            [u,v] = find(triu(ones(num_of_rois),1));
            idx = (v - 1)*num_of_rois + u;
            % search files.
            files = dir(strcat(root_path, "/minimal_preprocessing_pipeline_ts_mat_", atlas, "/hcpd_*_V1_MR_rfMRI_REST*_", atlas, ".mat"));
            disp(size(files));
            num_of_scans = size(files, 1);
            num_of_timepoints = 478;

            %% define variable.
            subj_scan_ids = [];
            % zscored bold ts.
            bold_zscored_all = nan(num_of_scans, num_of_rois, num_of_timepoints);
            % the rss timeseries of whole-connectome.
            ets_global_rss_all = nan(num_of_scans, num_of_timepoints);
            ets_global_rss_bins_all = nan(num_of_scans, num_of_bins);
            % the rss timeseries of each region (co-fluc with other regions).
            ets_region_rss_all = nan(num_of_scans, num_of_rois, num_of_timepoints);
            ets_region_rss_bins_all = nan(num_of_scans, num_of_rois, num_of_bins);
            high_fd_label_all = nan(num_of_scans, num_of_timepoints);
            % the contribution of each region (co-fluc with other regions) to whole-brain co-fluc amp.
            ets_region_level_contribution_ratioOfMeans_all = nan(num_of_scans, num_of_rois, num_of_bins);

            %% main function.
            for scan_i=1:num_of_scans
                timeseries = load(strcat(files(scan_i).folder, "/", files(scan_i).name));
                timeseries = timeseries.ts_roi;
                timeseries = timeseries.';

                % load fd data.
                folder_arr = strsplit(files(scan_i).name, '_');
                subj_id = folder_arr(2);
                subj_id=string(subj_id{1,1});
                scan_id = extractBetween(files(scan_i).name, "_V1_MR_", strcat("_", atlas, ".mat"));
                scan_id=string(scan_id{1,1});
                subj_scan_ids = [subj_scan_ids, strcat(subj_id, scan_id)];
                fd_subj_i_scan_i = load(strcat("/data/HCP-D/fmriresults01/", subj_id, "_V1_MR/MNINonLinear/Results/", scan_id, "/Movement_RelativeRMS.txt"));
                high_fd_label = (fd_subj_i_scan_i > (fd_threshold/10));
                % some scan's length is less than 478.
                high_fd_label_all(scan_i, 1:size(high_fd_label, 1)) = high_fd_label;

                % remove high head motion frames.
                timeseries(high_fd_label, :) = [];

                % bins.
                bin_size = floor(size(timeseries, 1)/num_of_bins);
                timeseries = zscore(timeseries);
                % save zescored bold ts.
                bold_zscored_all(scan_i, :, 1:size(timeseries, 1))=timeseries.';
                ets = timeseries(:, u).*timeseries(:, v);
                ets_rss = sum(ets.^2, 2).^0.5;
                ets_global_rss_all(scan_i, 1:size(timeseries, 1)) = ets_rss;
                [~,idxsort] = sort(ets_rss, 'descend');

                % CS.
                for roi_i=1:num_of_rois
                    edges_related_roi_i = (u==roi_i) | (v==roi_i);
                    edges_idx_related_roi_i = find(edges_related_roi_i==1);
                    ets_related_roi_i = ets(:, edges_idx_related_roi_i);
                    % region level RSS.
                    ets_related_roi_i_rss = sum(ets_related_roi_i.^2, 2).^0.5;
                    ets_region_rss_all(scan_i, roi_i, 1:size(timeseries, 1)) = ets_related_roi_i_rss;
                    for bin_i=1:num_of_bins
                        bin_idx = idxsort((bin_i-1)*bin_size+1:bin_i*bin_size);
                        rss_d = mean(ets_rss(bin_idx));
                        rss_n = mean(ets_related_roi_i_rss(bin_idx));
                        ets_region_level_contribution_ratioOfMeans_all(scan_i, roi_i, bin_i) = rss_n/rss_d;
                        ets_global_rss_bins_all(scan_i, bin_i)=rss_d;
                        ets_region_rss_bins_all(scan_i, roi_i, bin_i)=rss_n;
                    end
                end
            end

            %% save results.
            file_name = strcat("./outputs_bins/HCPD/bins_", num2str(num_of_bins), "_", atlas, "_default_fd", num2str(fd_threshold), "mm.mat");
            save(file_name, ...
                'subj_scan_ids', ...
                'bold_zscored_all', ...
                'ets_global_rss_all', ...
                'ets_global_rss_bins_all', ...
                'ets_region_rss_all', ...
                'ets_region_rss_bins_all', ...
                'high_fd_label_all', ...
                'ets_region_level_contribution_ratioOfMeans_all', ...
                '-v7.3');
        end
    end
end


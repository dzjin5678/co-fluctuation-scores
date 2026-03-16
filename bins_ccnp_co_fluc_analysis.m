addpath(genpath('fcn'));

% atlas settings.
atlas = 'schaefer200x17';
load('../atlas/Schaefer200x17.mat');
load('../atlas/scheafer200x17to8_subnet_edge_idx.mat');
num_of_rois = 200;
[u,v] = find(triu(ones(num_of_rois),1));
idx = (v - 1)*num_of_rois + u;
num_of_bins=10;
fd_perc=60;

for site_name=["PEK", "CKG"]
    data = load(strcat("../data/ccnp/ccnp_", lower(site_name), "_session_ids.mat"), "session_ids");
    session_ids = data.session_ids;
    num_of_sessions = size(session_ids, 1);

    for fd_thres=[20, 30];
        output_file_name=strcat(site_name, "_bins_", num2str(num_of_bins), "_thres", num2str(fd_thres), "_perc", num2str(fd_perc));

        %% define variable.
        num_of_timepoints = 360;
        subj_session_ids = [];
        % zscored bold ts.
        bold_zscored_all = nan(num_of_sessions, num_of_rois, num_of_timepoints);
        % the rss timeseries of whole-connectome.
        ets_global_rss_all = nan(num_of_sessions, num_of_timepoints);
        ets_global_rss_bins_all = nan(num_of_sessions, num_of_bins);
        % the rss timeseries of each region (co-fluc with other regions).
        ets_region_rss_all = nan(num_of_sessions, num_of_rois, num_of_timepoints);
        ets_region_rss_bins_all = nan(num_of_sessions, num_of_rois, num_of_bins);
        % high_fd_label_all = nan(num_of_sessions, num_of_timepoints);
        % the contribution of each region (co-fluc with other regions) to whole-brain co-fluc amp.
        ets_region_level_contribution_ratioOfMeans_all = nan(num_of_sessions, num_of_rois, num_of_bins);

        % main.
        for sess_i=1:num_of_sessions
            disp(session_ids(sess_i));
            % you need to search corresponding ts files, at least two session.
            [timeseries, mean_fd] = extract_and_control_ts(session_ids, sess_i, fd_thres, fd_perc);
            if sum(isnan(timeseries))>0
                continue;
            end
            subj_session_ids=[subj_session_ids, string(session_ids(sess_i))];
            disp(size(timeseries));
            bin_size = floor(size(timeseries, 1)/num_of_bins);
            timeseries = zscore(timeseries);
            % save zescored bold ts.
            bold_zscored_all(sess_i, :, 1:size(timeseries, 1))=timeseries.';
            ets = timeseries(:, u).*timeseries(:, v);
            ets_rss = sum(ets.^2, 2).^0.5;
            ets_global_rss_all(sess_i, 1:size(timeseries, 1)) = ets_rss;
            [~,idxsort] = sort(ets_rss, 'descend');

            % ROI RSS. and region-level dreiver.
            for roi_i=1:num_of_rois
                edges_related_roi_i = (u==roi_i) | (v==roi_i);
                edges_idx_related_roi_i = find(edges_related_roi_i==1);
                ets_related_roi_i = ets(:, edges_idx_related_roi_i);
                % region level RSS.
                ets_related_roi_i_rss = sum(ets_related_roi_i.^2, 2).^0.5;
                % CS.
                for roi_i=1:num_of_rois
                    edges_related_roi_i = (u==roi_i) | (v==roi_i);
                    edges_idx_related_roi_i = find(edges_related_roi_i==1);
                    ets_related_roi_i = ets(:, edges_idx_related_roi_i);
                    % region level RSS.
                    ets_related_roi_i_rss = sum(ets_related_roi_i.^2, 2).^0.5;
                    ets_region_rss_all(sess_i, roi_i, 1:size(timeseries, 1)) = ets_related_roi_i_rss;
                    for bin_i=1:num_of_bins
                        bin_idx = idxsort((bin_i-1)*bin_size+1:bin_i*bin_size);
                        rss_d = mean(ets_rss(bin_idx));
                        rss_n = mean(ets_related_roi_i_rss(bin_idx));
                        ets_region_level_contribution_ratioOfMeans_all(sess_i, roi_i, bin_i) = rss_n/rss_d;
                        ets_global_rss_bins_all(sess_i, bin_i)=rss_d;
                        ets_region_rss_bins_all(sess_i, roi_i, bin_i)=rss_n;
                    end
                end
            end
        end
        file_name = strcat("./outputs_bins/CCNP/", output_file_name, ".mat");
        save(file_name, ...
            'subj_session_ids', ...
            'bold_zscored_all', ...
            'ets_global_rss_all', ...
            'ets_global_rss_bins_all', ...
            'ets_region_rss_all', ...
            'ets_region_rss_bins_all', ...
            'ets_region_level_contribution_ratioOfMeans_all', ...
            '-v7.3');
    end
end


function [ts_all, mean_fd] = extract_and_control_ts(session_ids, sesssion_i, fd_thres, fd_perc)

    num_of_timepoints_orig = 0;
    num_of_timepoints_controlled = 0;
    mean_fd = 0;

    ts_all = [];
    session_char = char(session_ids(sesssion_i, :));
    site_name = lower(session_char(9:11));

    files=dir(strcat("/mnt/system_v2022/my_data/dzjin_data/ccnp/ccnp_",site_name,"_roi_ts/", session_ids(sesssion_i, :), "_run_rest*_roi_ts_schaefer200x17.mat"));
    if size(files, 1)<2
        ts_all = nan(1, 1);
        mean_fd=999;
        return;
    end

    for f_i=1:size(files, 1)
        % roi ts.
        ts_rest1 = load(strcat(files(f_i).folder, "/", files(f_i).name));
        ts_rest1 = ts_rest1.roi_ts;
        num_of_timepoints_orig = num_of_timepoints_orig + size(ts_rest1, 2);
        
        % remove high head motion frames.
        fd_file_path=strcat(files(f_i).folder, "/", files(f_i).name);
        fd_file_path = strrep(fd_file_path, strcat("ccnp_",site_name,"_roi_ts"), strcat("ccnp_",site_name,"_fd"));
        fd_file_path = strrep(fd_file_path, "roi_ts_schaefer200x17.mat", "rest_fwd_abssum.1D");
        fd_1 = readmatrix(fd_file_path, "FileType", "delimitedtext");
        fd_1 = fd_1(:, 2);
        mean_fd = (mean_fd+mean(fd_1(fd_1 <= (fd_thres/100))))/2;
        ts_rest1(:, fd_1 > (fd_thres/100))=[];
        num_of_timepoints_controlled = num_of_timepoints_controlled + size(ts_rest1, 2);

        ts_rest1 = ts_rest1.';
        ts_rest1 = zscore(ts_rest1);
        if f_i==1
            ts_all = ts_rest1;
        else
            ts_all = vertcat(ts_all, ts_rest1);
        end
    end

    if num_of_timepoints_controlled~=0
        if (num_of_timepoints_controlled/num_of_timepoints_orig) < (fd_perc/100)
            ts_all = nan(1, 1);
            mean_fd=999;
        end
    else
        ts_all = nan(1, 1);
        mean_fd=999;
    end
end


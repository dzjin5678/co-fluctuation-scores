num_of_rois=200;
% Index of upper triabgle in FC matrix.
[u,v] = find(triu(ones(num_of_rois),1));
idx = (v - 1)*num_of_rois + u;
num_of_bins=20;

% load BOLD timeseries.
% timeseries: timepoints * num_of_rois

% calcualte bin size.
bin_size = floor(size(timeseries, 1)/num_of_bins);
% zsore timeseries.
timeseries = zscore(timeseries);
% calculate edge timeseries.
ets = timeseries(:, u).*timeseries(:, v);
% calculate GLOBAL RSS timeseries.
ets_rss = sum(ets.^2, 2).^0.5;
% sort timepoints according to GLOABL RSS.
[~,idxsort] = sort(ets_rss, 'descend');

cofluc_score(roi_i, bin_i) = nan(num_of_rois, num_of_bins);
% calculate REGION RSS and cofluctuation score.
for roi_i=1:num_of_rois
    edges_related_roi_i = (u==roi_i) | (v==roi_i);
    edges_idx_related_roi_i = find(edges_related_roi_i==1);
    ets_related_roi_i = ets(:, edges_idx_related_roi_i);
    % region level RSS.
    ets_related_roi_i_rss = sum(ets_related_roi_i.^2, 2).^0.5;

    for bin_i=1:num_of_bins
        % whole brain idxsort !!!
        bin_idx = idxsort((bin_i-1)*bin_size+1:bin_i*bin_size);
        rss_d = mean(ets_rss(bin_idx));
        rss_n = mean(ets_related_roi_i_rss(bin_idx));
        cofluc_score(roi_i, bin_i) = rss_n/rss_d;
    end
end
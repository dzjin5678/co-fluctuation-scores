%% dark red and blue corlormap.
maxcolor    = [219/255 49/255 36/255]; % darkred
mediancolor = [1 1 1]; % white   
mincolor    = [69/255 123/255 157/255]; % darkblue      
ColorMapSize = 90;
int1 = zeros(ColorMapSize,3); 
int1_0 = zeros(ColorMapSize,3); 
int2 = zeros(ColorMapSize,3);
for k=1:3
    int1(:,k) = linspace(mincolor(k), mediancolor(k), ColorMapSize);
    int1_0(:,k) = linspace(mediancolor(k), mincolor(k), ColorMapSize);
    int2(:,k) = linspace(mediancolor(k), maxcolor(k), ColorMapSize);
end
meep = [int1(1:end-1,:); int2];
meep_blue = [int1(1:end-1,:)];
meep_blue_reverse = [int1_0(1:end-1,:)];
meep_red = [int2];
meep_orig = meep;
% meep(45:60, :) = [];


%% sktetch of cs matrix.
load("../outputs_bins/HCP/bins_20_bin_hcp_icafixgs_ts_scan_REST1_LR_schaefer200x17_hfd_0.mat", "ets_region_level_contribution_all");
ets_region_level_contribution_group = squeeze(mean(ets_region_level_contribution_all, 1));


%%
figure("Position", [100 100 200 420]);
imagesc(ets_region_level_contribution_group(1:100, :));
colormap(meep_orig);
hold on;
line_width=1;

line_color = [100 100 100]./255;

plot([0.5,20 + 0.5],(20 + 0.5)*ones(1,2),'k', 'color',line_color, LineWidth=line_width);hold on;
plot([0.5,20 + 0.5],(21 + 0.5)*ones(1,2),'k', 'color',line_color, LineWidth=line_width);hold on;
plot([0.5, 0.5], [20.5, 21.5],'k', 'color',line_color, LineWidth=line_width);hold on;
plot([20.5, 20.5], [20.5, 21.5],'k', 'color',line_color, LineWidth=line_width);hold on;

plot([0.5,20 + 0.5],(70 + 0.5)*ones(1,2),'k', 'color',line_color, LineWidth=line_width);hold on;
plot([0.5,20 + 0.5],(71 + 0.5)*ones(1,2),'k', 'color',line_color, LineWidth=line_width);hold on;
plot([0.5, 0.5], [70.5, 71.5],'k', 'color',line_color, LineWidth=line_width);hold on;
plot([70.5, 70.5], [70.5, 71.5],'k', 'color',line_color, LineWidth=line_width);hold on;

red = [222,45,38]./255;
light_red = [252,146,114]./255;
blue = [49,130,189]./255;

red = line_color;

plot([0.5, 0.5],[0.5, 100.5],'k', 'color',red, LineWidth=line_width);hold on;
plot([1.5, 1.5],[0.5, 100.5],'k', 'color',red, LineWidth=line_width);hold on;
plot([0.5, 1.5],[0.5, 0.5],'k', 'color',red, LineWidth=line_width);hold on;
plot([0.5, 1.5],[100.5, 100.5],'k', 'color',red, LineWidth=line_width);hold on;

plot([9.5, 9.5],[0.5, 100.5],'k', 'color',red, LineWidth=line_width);hold on;
plot([10.5, 10.5],[0.5, 100.5],'k', 'color',red, LineWidth=line_width);hold on;
plot([9.5, 10.5],[0.5, 0.5],'k', 'color',red, LineWidth=line_width);hold on;
plot([9.5, 10.5],[100.5, 100.5],'k', 'color',red, LineWidth=line_width);hold on;

plot([19.5, 19.5],[0.5, 100.5],'k', 'color',red, LineWidth=line_width);hold on;
plot([20.5, 20.5],[0.5, 100.5],'k', 'color',red, LineWidth=line_width);hold on;
plot([19.5, 20.5],[0.5, 0.5],'k', 'color',red, LineWidth=line_width);hold on;
plot([19.5, 20.5],[100.5, 100.5],'k', 'color',red, LineWidth=line_width);hold on;

ax=gca;hold on;
ax.XColor='none';
ax.YColor='none';
set(gca,'Box','off');
print(gcf, '-dpng', '-r600', './new_sketch/cs');
print(gcf, '-dmeta', './new_sketch/cs');
close;


%% settings.
% Index of upper triabgle in FC matrix.
num_of_rois = 200;
[u,v] = find(triu(ones(num_of_rois),1));
idx = (v - 1)*num_of_rois + u;

data_path = 'G:\datasets\HCP\hcp_icafixgs_ts_mat_schaefer200x17\100307_rfMRI_REST1_LR.mat';
load(data_path);

dan_roi_idx = horzcat(29:39, 131:142).';
dan_roi_1_idx = 29;

dan_related_edge_idx = zeros(size(idx, 1), 1);
dan_roi_1_related_edge_idx = zeros(size(idx, 1), 1);
for i=1:size(idx, 1)
    if (sum(dan_roi_idx==u(i)) > 0) && (sum(dan_roi_idx==v(i)) > 0)
        dan_related_edge_idx(i) = 1;
    end

    if dan_roi_1_idx==u(i)
        dan_roi_1_related_edge_idx(i) = 1;
    end

    if dan_roi_1_idx==v(i)
        dan_roi_1_related_edge_idx(i) = 1;
    end
    
end

% timeseries = timeseries(201:600, :);
timeseries = zscore(timeseries);

% fc = corr(timeseries);
% imagesc(fc);
% colormap(meep);

ets = timeseries(:, u).*timeseries(:, v);
ets_rss = sum(ets.^2, 2).^0.5;
[~,idxsort] = sort(ets_rss, 'descend');

ets_bin_1 = ets(idxsort(1:60), :);
ets_rss_bin_1 = ets_rss(idxsort(1:60));

ets_bin_2 = ets(idxsort(61:120), :);
ets_rss_bin_2 = ets_rss(idxsort(61:120));

ets_bin_20 = ets(idxsort(1141:1200), :);
ets_rss_bin_20 = ets_rss(idxsort(1141:1200));

% ets_dan_reltaed = ets(:, find(dan_related_edge_idx));
% ets_dan_reltaed_bin_1 = ets_dan_reltaed(idxsort(1:60), :);
% ets_dan_reltaed_rss = sum(ets_dan_reltaed.^2, 2).^0.5;
% ets_dan_reltaed_rss_bin_1 = ets_dan_reltaed_rss(idxsort(1:60));

ets_dan_roi_1_reltaed = ets(:, find(dan_roi_1_related_edge_idx));
% ets_dan_roi_1_reltaed_bin_1 = ets_dan_roi_1_reltaed(idxsort(1:60), :);
ets_dan_roi_1_reltaed_rss = sum(ets_dan_roi_1_reltaed.^2, 2).^0.5;
[~,idxsort_dan_roi_1] = sort(ets_dan_roi_1_reltaed_rss, 'descend');


%% ets.
figure("Position", [100 100 500 60]);
% plot(timeseries(:, 8), "Color", [255/255 80/255 80/255], "LineWidth", 1); hold on;
plot(timeseries(:, 8), "Color", [253,174,107]./255, "LineWidth", 1); hold on;
plot([0.5,1200 + 0.5], [0 0],'k', 'color',[0,0,0]./255, 'LineWidth', 0.3);
ax=gca;hold on;
ax.XColor='none';
ax.YColor='none';
set(gca,'Box','off');
print(gcf, '-dpng', '-r600', './new_sketch/bold_roi_8');
print(gcf, '-dmeta', './new_sketch/bold_roi_8');
close;

figure("Position", [100 100 500 60]);
% plot(timeseries(:, 9), "Color", [0/255 112/255 192/255], "LineWidth", 1); hold on;
plot(timeseries(:, 9), "Color", [128,205,193]./255, "LineWidth", 1); hold on;
plot([0.5,1200 + 0.5], [0 0],'k', 'color',[0,0,0]./255, 'LineWidth', 0.3);
plot([0.5,1200 + 0.5], [0 0],'k', 'color',[0,0,0]./255, 'LineWidth', 0.3);
ax=gca;hold on;
ax.XColor='none';
ax.YColor='none';
set(gca,'Box','off');
print(gcf, '-dpng', '-r600', './new_sketch/bold_roi_9');
print(gcf, '-dmeta', './new_sketch/bold_roi_9');
close;

figure("Position", [100 100 500 60]);
plot(ets(:, 36), "Color", [100 100 100]./255, "LineWidth", 1);  hold on;
plot([0.5,1200 + 0.5], [0 0],'k', 'color',[0,0,0]./255, 'LineWidth', 0.3);
ax=gca;hold on;
ax.XColor='none';
ax.YColor='none';
set(gca,'Box','off');
print(gcf, '-dpng', '-r600', './new_sketch/ets_edge_36');
print(gcf, '-dmeta', './new_sketch/ets_edge_36');
close;

figure("Position", [100 100 500 200]);
imagesc(ets.');
colormap(meep);hold on;
line_color = meep_red(70, :);
plot([1, 1200], [1, 1], "Color", line_color, 'LineWidth', 2);hold on;
plot([1, 1], [1, 19900], "Color", line_color, 'LineWidth', 3);hold on;
plot([1200, 1200], [1, 19900], "Color", line_color, 'LineWidth', 3);hold on;
plot([1, 1200], [19900, 19900], "Color", line_color, 'LineWidth', 2);hold on;
line_color = meep_blue(20, :);
plot([1, 1200], [8000, 8000], "Color", line_color, 'LineWidth', 2);hold on;
plot([1, 1], [8000, 12000], "Color", line_color, 'LineWidth', 3.5);hold on;
plot([1200, 1200], [8000, 12000], "Color", line_color, 'LineWidth', 3.5);hold on;
plot([1, 1200], [12000, 12000], "Color", line_color, 'LineWidth', 2);hold on;
ax=gca;hold on;
% xticks([1, 1200]);
% xticklabels([1, 1200]);
% yticks([1, 19900]);
% yticklabels([1, 19900]);
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
% ylim([0 500]);
ax.FontSize=21;
print(gcf, '-dpng', '-r600', './new_sketch/ets');
print(gcf, '-dmeta', './new_sketch/ets');
close;


%% rss timeserries.
figure("Position", [100 100 500 180]);
h1 = plot(ets_rss, "Color", meep_red(43, :), "LineWidth", 1); 
ax=gca;hold on;
set(gca,'Box','off', 'TickDir', 'none');
% xticks([0, 1200]);
% xticklabels([0, 1200]);
% yticks([0, 500]);
% yticklabels([0, 500]);
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
ylim([0 500]);
ax.FontSize=23;
print(gcf, '-dpng', '-r600', './new_sketch/ets_rss_global');
print(gcf, '-dmeta', './new_sketch/ets_rss_global');
close;

figure("Position", [100 100 500 180]);
h2 = plot(ets_dan_roi_1_reltaed_rss, "Color", meep_blue(1, :), "LineWidth", 1); 
ax=gca;hold on;
set(gca,'Box','off', 'TickDir', 'none');
% xticks([0, 1200]);
% xticklabels([0, 1200]);
% yticks([0, 100]);
% yticklabels([0, 100]);
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
ylim([0 100]);
ax.FontSize=23;
print(gcf, '-dpng', '-r600', './new_sketch/ets_rss_local');
print(gcf, '-dmeta', './new_sketch/ets_rss_local');
close;


%% rss distributuion.
figure("Position", [100 100 300 200]);
histogram(ets_rss, 30, "FaceColor", [189/255 189/255 189/255]);
ax=gca;hold on;
% xticks([0, 400, 800, 1200]);
% xticklabels([0, 400, 800, 1200]);
% yticks([0, 250, 500]);
% yticklabels([0, 250, 500]);
% ylim([0 500]);
ax.FontSize=17;
% ax.XColor='none';
ax.YColor='none';
set(gca,'Box','off')
print(gcf, '-dpng', '-r600', './new_sketch/ets_rss_global_distribution');
close;


%% rss imagesc.
figure("Position", [100 100 700 40]);
imagesc(ets_rss.');
colormap(meep_red);
hold on;
plot((1)*ones(1,2),[0.5,1 + 0.5],'k', 'color',[0,0,0]./255, 'LineWidth', 0.3);
plot((1200)*ones(1,2),[0.5,1 + 0.5],'k', 'color',[0,0,0]./255, 'LineWidth', 0.3);
plot([0.5,1200 + 0.5], [0.5, 0.5],'k', 'color',[0,0,0]./255, 'LineWidth', 0.3);
plot([0.5,1200 + 0.5], [1.5, 1.5],'k', 'color',[0,0,0]./255, 'LineWidth', 0.3);
ax=gca;hold on;
ax.XColor='none';
ax.YColor='none';
print(gcf, '-dpng', '-r600', './new_sketch/ets_rss_global_imagesc');
print(gcf, '-dmeta', './new_sketch/ets_rss_global_imagesc');
close;

% figure("Position", [100 100 700 40]);
% imagesc(ets_rss(idxsort).');
% colormap(meep_red);
% hold on;
% plot((1)*ones(1,2),[0.5,1 + 0.5],'k', 'color',[0,0,0]./255, 'LineWidth', 0.3);
% for bin_i=1:20
%     plot((bin_i*60)*ones(1,2),[0.5,1 + 0.5],'k', 'color',[0,0,0]./255, 'LineWidth', 0.3);
% end
% plot([0.5,1200 + 0.5], [0.5, 0.5],'k', 'color',[0,0,0]./255, 'LineWidth', 0.3);
% plot([0.5,1200 + 0.5], [1.5, 1.5],'k', 'color',[0,0,0]./255, 'LineWidth', 0.3);
% ax=gca;hold on;
% ax.XColor='none';
% ax.YColor='none';
% print(gcf, '-dpng', '-r600', './new_sketch/ets_rss_global_reorganized_imagesc');
% print(gcf, '-dmeta', './new_sketch/ets_rss_global_reorganized_imagesc');
% close;

figure("Position", [100 100 700 40]);
ets_reorganized_averaged = nan(1, 1200);
ets_rss_sorted = ets_rss(idxsort);
ets_reorganized_averaged(1:60) = mean(ets_rss_sorted(1:60));
for bin_i=2:20
    bin_index = (bin_i-1)*60:bin_i*60;
    ets_reorganized_averaged(bin_index) = mean(ets_rss_sorted(bin_index));
end
imagesc(ets_reorganized_averaged);
colormap(meep_red);
hold on;
plot((1)*ones(1,2),[0.5,1 + 0.5],'k', 'color',[0,0,0]./255, 'LineWidth', 0.3);
for bin_i=1:20
    plot((bin_i*60)*ones(1,2),[0.5,1 + 0.5],'k', 'color',[0,0,0]./255, 'LineWidth', 0.3);
end
plot([0.5,1200 + 0.5], [0.5, 0.5],'k', 'color',[0,0,0]./255, 'LineWidth', 0.3);
plot([0.5,1200 + 0.5], [1.5, 1.5],'k', 'color',[0,0,0]./255, 'LineWidth', 0.3);
ax=gca;hold on;
ax.XColor='none';
ax.YColor='none';
print(gcf, '-dpng', '-r600', './new_sketch/ets_rss_global_reorganized_averaged_imagesc');
print(gcf, '-dmeta', './new_sketch/ets_rss_global_reorganized_averaged_imagesc');
close;


%% dan roi 1.
figure("Position", [100 100 700 40]);
imagesc(ets_dan_roi_1_reltaed_rss.');
colormap(meep_blue_reverse);
hold on;
plot((1)*ones(1,2),[0.5,1 + 0.5],'k', 'color',[0,0,0]./255, 'LineWidth', 0.3);
plot((1200)*ones(1,2),[0.5,1 + 0.5],'k', 'color',[0,0,0]./255, 'LineWidth', 0.3);
plot([0.5,1200 + 0.5], [0.5, 0.5],'k', 'color',[0,0,0]./255, 'LineWidth', 0.3);
plot([0.5,1200 + 0.5], [1.5, 1.5],'k', 'color',[0,0,0]./255, 'LineWidth', 0.3);
ax=gca;hold on;
ax.XColor='none';
ax.YColor='none';
print(gcf, '-dpng', '-r600', './new_sketch/ets_rss_dan_roi_1_imagesc');
print(gcf, '-dmeta', './new_sketch/ets_rss_dan_roi_1_imagesc');
exportgraphics(gca, './new_sketch/ets_rss_dan_roi_1_imagesc_new.emf', 'ContentType','image');
close;

% figure("Position", [100 100 700 40]);
% imagesc(ets_dan_roi_1_reltaed_rss(idxsort_dan_roi_1).');
% colormap(meep_blue_reverse);
% hold on;
% plot((1)*ones(1,2),[0.5,1 + 0.5],'k', 'color',[0,0,0]./255, 'LineWidth', 0.3);
% for bin_i=1:20
%     plot((bin_i*60)*ones(1,2),[0.5,1 + 0.5],'k', 'color',[0,0,0]./255, 'LineWidth', 0.3);
% end
% plot([0.5,1200 + 0.5], [0.5, 0.5],'k', 'color',[0,0,0]./255, 'LineWidth', 0.3);
% plot([0.5,1200 + 0.5], [1.5, 1.5],'k', 'color',[0,0,0]./255, 'LineWidth', 0.3);
% ax=gca;hold on;
% ax.XColor='none';
% ax.YColor='none';
% print(gcf, '-dpng', '-r600', './new_sketch/ets_rss_dan_roi_1_reorganized_imagesc');
% print(gcf, '-dmeta', './new_sketch/ets_rss_dan_roi_1_reorganized_imagesc');
% close;

figure("Position", [100 100 700 40]);
ets_reorganized_averaged = nan(1, 1200);
ets_rss_sorted = ets_dan_roi_1_reltaed_rss(idxsort_dan_roi_1);
ets_reorganized_averaged(1:60) = mean(ets_rss_sorted(1:60));
for bin_i=2:20
    bin_index = (bin_i-1)*60:bin_i*60;
    ets_reorganized_averaged(bin_index) = mean(ets_rss_sorted(bin_index));
end
imagesc(ets_reorganized_averaged);
colormap(meep_blue_reverse);
hold on;
plot((1)*ones(1,2),[0.5,1 + 0.5],'k', 'color',[0,0,0]./255, 'LineWidth', 0.3);
for bin_i=1:20
    plot((bin_i*60)*ones(1,2),[0.5,1 + 0.5],'k', 'color',[0,0,0]./255, 'LineWidth', 0.3);
end
plot([0.5,1200 + 0.5], [0.5, 0.5],'k', 'color',[0,0,0]./255, 'LineWidth', 0.3);
plot([0.5,1200 + 0.5], [1.5, 1.5],'k', 'color',[0,0,0]./255, 'LineWidth', 0.3);
ax=gca;hold on;
ax.XColor='none';
ax.YColor='none';
print(gcf, '-dpng', '-r600', './new_sketch/ets_rss_dan_roi_1_reorganized_averaged_imagesc');
print(gcf, '-dmeta', './new_sketch/ets_rss_dan_roi_1_reorganized_averaged_imagesc');
close;


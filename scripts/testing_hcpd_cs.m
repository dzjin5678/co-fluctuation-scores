%% plot roi-level contribution by age.
site = "HCPD";
num_of_bins = 20;
load(strcat("../../data/hcp-d/hcpd_subj_scan_ages.mat"));
fd=2;
load(strcat("../outputs_bins/", site ,"/bins_", num2str(num_of_bins) ,"_bin_schaefer200x17_default_fd", num2str(fd), "mm.mat"), "ets_region_level_contribution_all");
brain_map_sa_axis = readmatrix("../results_bins/mechanism/SA_Axis/s_a_axis_shaefer200x17.csv");


%%
% final_labels = ages>=6 & ages<19;
final_labels = ages>=1;
final_subj_ids = subj_ids(final_labels, :);
final_subj_ids_unique = unique(final_subj_ids);
final_sex = sexs(final_labels, :);
final_age = ages(final_labels, :);

num_of_females = 0;
age_all = [];
for s_i=1:530
    sex = sexs(final_subj_ids==final_subj_ids_unique(s_i));
    % disp(sex);
    if sex(1)==0
        num_of_females = num_of_females + 1;
    end
    age = ages(final_subj_ids==final_subj_ids_unique(s_i));
    age_all = [age_all age(1)];
end


%% coloobar, corrs.
maxcolor    = [247/255,126/255,105/255]; %   
mediancolor = [255/255 255/255 255/255]; %    
mincolor    = [52/255,96/255,141/255]; %  
ColorMapSize = 50;
int1 = zeros(ColorMapSize,3); 
int2 = zeros(ColorMapSize,3);
for k=1:3
    int1(:,k) = linspace(mincolor(k), mediancolor(k), ColorMapSize);
    int2(:,k) = linspace(mediancolor(k), maxcolor(k), ColorMapSize);
end
colormap_corrs = [int1(1:end-1,:); int2];


%% bins auto-corr.
[u_bins,v_bins] = find(triu(ones(num_of_bins),1));
idx_bins = (v_bins - 1)*num_of_bins + u_bins;

std_age = nan(18, 1);
diff_age = nan(18, 1);
bins_corr_all_age = nan(20, size(u_bins, 1));
for age=6:18
    scans_idx = (ages >= age) & (ages < (age+1));
    ets_region_level_contribution_age_i = ets_region_level_contribution_all(scans_idx, :, :);
    roi_contri_group = squeeze(mean(ets_region_level_contribution_age_i, 1));

    % figure("Position", [100 100 390 320]);
    bins_corr = ones(num_of_bins, num_of_bins);
    for bin_i=1:num_of_bins
        for bin_j=(bin_i+1):num_of_bins
            [coef, ~] = corr(roi_contri_group(:, bin_i), roi_contri_group(:, bin_j));
            bins_corr(bin_i, bin_j) = coef;
            bins_corr(bin_j, bin_i) = coef;
        end
    end

    figure("Position", [100, 100, 200, 185]);
    imagesc(bins_corr);
    colormap(colormap_corrs);
    clim([-1 1]);  
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    % title(strcat(num2str(age), " year"), 'fontsize', 22);
    print(gcf, '-dpng', '-r300', ['../results_bins/DRIVER/' char(site) '/region_level_fd' char(num2str(fd)) '/bins_', char(num2str(num_of_bins)) ,'_corr_' char(site) '_age_' , num2str(age), '_', 'schaefer200x17', '_dynamic_sa_axis_raw_bins20']);
    close;

    corrs = bins_corr(idx_bins);
    bins_corr_all_age(age, :) = corrs;

    pos_mean = mean(corrs(corrs>0));
    neg_mean = mean(corrs(corrs<0));
    diff = pos_mean - neg_mean;
    diff_age(age) = diff;

    std_age(age) = std(corrs);
end


%% distribution of auto-corr.
mincolor    = [189/255 189/255 189/255]; % gray
% lgLabel = {'6 year', '7 year', '8 year', '9 year', '10 year', '11 year', ...
%     '12 year', '13 year', '14 year', '15 year', '16 year', '17 year', '18 year'};
lgLabel = {' ', ' ', ' ', ' ', ' ', ' ', ...
    ' ', ' ', ' ', ' ', ' ', ' ', ' '};
p = 1;
% lineLength = 0.04; % 竖线的长度
% figure("Position", [100, 100, 260, 400]);
figure;
set(gcf, 'PaperUnits', 'inches');             % 设置单位为英寸
set(gcf, 'PaperPosition', [0 0 1.5 5]);         % 设置 6x4 英寸大小
hold on;
for age = 6:1:18
    [f, x] = ksdensity(bins_corr_all_age(age, :));
    fShifted = f + age * p;
    pHandle = plot(x, fShifted, 'color', mincolor, 'LineWidth', 0.3, 'HandleVisibility', 'off');
    % yline(age * p, '-', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    Xfill = [x, fliplr(x)];
    Yfill = [fShifted, ones(1, length(x)) * age * p];
    fill(Xfill, Yfill, pHandle.Color, 'EdgeColor', 'none', 'FaceAlpha', 1);
    % for j = 1:length(subnet_i_cg)
    %     line([subnet_i_cg(j), subnet_i_cg(j)], [i * p, i * p + lineLength], 'Color', 'k', 'LineWidth', 1);
    % end
end
yTick = (6:18);
xTick = [-2, 0, 2];
set(gca, ...
    'YTick', yTick, ...
    'xTick', xTick, ...
    'YTickLabel', lgLabel, ...
    'XTickLabel', {' ', ' ', ' '}, ...
    'fontname', 'Airal', 'fontsize', 15);
% set(gca, 'YTick', yTick, ...
%     'YTickLabel', lgLabel, ...
%     'XMinorTick', 'on', 'TickDir', 'none', ...
%     'fontname', 'Airal', 'fontsize', 15);
% xlabel(['Auto Correlation ' newline 'between 20 Bins']);
% ylabel("Age");
ylim([5.8 19]);
ax=gca;hold on;
% ax.FontSize=17;
ax.XAxis.FontSize = 18;  % 设置 X 轴刻
ax.YAxis.FontSize = 17;  % 设置 X 轴刻
save_path_filename = ['../results_bins/DRIVER/' char(site) '/region_level_fd' char(num2str(fd)) '/bins_', char(num2str(num_of_bins)), '_auto_corr_' char(site) '_all_', 'schaefer200x17'];
print(gcf, '-dpng', '-r600', save_path_filename);
print(gcf, '-dmeta', save_path_filename);
close;


%% similarity between SA axis and cofluctuation score changed with development.
maxcolor    = [227,74,51]./255; % purple  
mediancolor = [255/255 255/255 255/255]; % white   
mincolor    = [8,81,156]./255; % blue    

ColorMapSize = 9;
int1 = zeros(ColorMapSize,3); 
int2 = zeros(ColorMapSize,3);
for k=1:3
    int1(:,k) = linspace(mincolor(k), mediancolor(k), ColorMapSize);
    int2(:,k) = linspace(mediancolor(k), maxcolor(k), ColorMapSize);
end
meep = [int1(1:end-1,:); int2];
meep(8:11, :)=[];

figure("Position", [100 100 550 400]);
plot([0, num_of_bins], [0, 0], '-o', 'Color', [235/255 235/255 235/255], 'LineWidth', 0.5, 'MarkerFaceColor',[218,165,32]./255, 'MarkerSize', 0.001); hold on;
x = linspace(1, num_of_bins, num_of_bins);
hcpd_cs_sa_corr = nan(19, num_of_bins);
h_all = [];
for age=6:18
    scans_idx = (ages >= age) & (ages < (age+1));
    ets_region_level_contribution_age_i = ets_region_level_contribution_all(scans_idx, :, :);
    roi_contri_group = squeeze(mean(ets_region_level_contribution_age_i, 1));
    coefs = nan(num_of_bins, 1);
    for bin_i=1:num_of_bins
        [coef, pval] = corr(brain_map_sa_axis, roi_contri_group(:, bin_i), "type", "Spearman");
        coefs(bin_i) = coef;
    end
    hcpd_cs_sa_corr(age, :) = coefs;
    h_age_i = plot(x, coefs, '-o', 'Color', meep(age-5, :), 'LineWidth', 1.2, 'MarkerSize', 2); hold on;
    h_all = [h_all h_age_i];
    scatter(x, coefs, 8, meep(age-5, :), 'filled'); hold on;
end
xlabel("Bins");
ylabel(['Similarity between ' newline ' Cofluctuation Score and SA Rank (r)']);
ylim([-1 1]);
xlim([0, num_of_bins]);
xticks([1,5,10,15,20]);
xticklabels(["95-100%", "75-80%", "50-55%", "25-30%", "0-5%"]);
% xticks([1,4,7,10]);
% xticklabels(["90-100%", "60-70%", "30-40%", "0-10%"]);
yticks(-1:0.2:1);
ax=gca;hold on;
ax.FontSize=12;
ax.FontSize=12;
ax.XTickLabelRotation=35;
set(gca,'Box','off','TickDir','none');
set(gca, 'FontName', 'Arial');
my_legend = legend(h_all, '6 year', '7 year', '8 year', '9 year', '10 year', '11 year', ...
    '12 year', '13 year', '14 year', '15 year', '16 year', '17 year', '18 year', 'NumColumns', 2);
legend('boxoff');
set(my_legend, "position", [0.5 0.2 0.4 0.33]);
print(gcf, '-dpng', '-r800', ['../results_bins/DRIVER/' char(site) '/region_level_fd' char(num2str(fd)) '/bins_', char(num2str(num_of_bins)), '_' char(site) '_corr_with_sa_axis_spearman']);
close;


%% age effect of cofluctuation scores.
site = "HCPD";
num_of_bins = 20;
load(strcat("../../../data/hcp-d/hcpd_subj_scan_ages.mat"));
load(strcat("../../outputs_bins/", site ,"/bins_", num2str(num_of_bins) ,"_bin_schaefer200x17_default_fd3mm.mat"), "ets_region_level_contribution_all");

num_of_rois = 200;
load('../../../atlas/Schaefer200x17.mat');

C1=[162, 81, 172]./255;
C2=[120, 154, 192]./255;
C3=[64, 152, 50]./255;
C4=[224, 102, 254]./255;
C5=[246, 253, 201]./255;
C6=[238, 185, 67]./255;
C7=[217, 113, 125]./255;
C8=[192, 0, 0]./255;
colorList=[C1; C2; C3; C4; C5; C6; C7; C8];

corr_r_bin_roi = nan(num_of_bins, num_of_rois);
corr_p_bin_roi = nan(num_of_bins, num_of_rois);
corr_t_bin_roi = nan(num_of_bins, num_of_rois);
corr_lme_p_bin_roi = nan(num_of_bins, num_of_rois);
for bin_i=1:num_of_bins
    for roi_i=1:num_of_rois
        ets_region_level_contribution_age_i = ets_region_level_contribution_all(:, roi_i, bin_i);
        ets_region_level_contribution_age_i = squeeze(ets_region_level_contribution_age_i);

        % 换成lme.
        column_names = {'cofluc_score', 'age'};
        data_table = table(ets_region_level_contribution_age_i, ages, 'VariableNames', column_names);
        % plot(data_table.age,data_table.cofluc_score,'ro');

        lme = fitlme(data_table,'cofluc_score ~ 1 + age');
        corr_t_bin_roi(bin_i, roi_i) = lme.Coefficients.tStat(2);
        corr_lme_p_bin_roi(bin_i, roi_i) = lme.Coefficients.pValue(2);

        % [RHO,PVAL] = corr(ets_region_level_contribution_age_i, ages);
        % corr_r_bin_roi(bin_i, roi_i) = RHO;
        % corr_p_bin_roi(bin_i, roi_i) = PVAL;
    end
end

save(['../../results_bins/DRIVER/' char(site) '/region_level/bins_' char(num2str(num_of_bins)) '_age_effect_lme/cofluctuation_score_age_effect'], "corr_t_bin_roi", "corr_lme_p_bin_roi");


%% session-level.
% bin1-2.
% for roi_i=1:num_of_rois
%     ets_region_level_contribution_age_i = ets_region_level_contribution_all(:,:, 1:2);
%     ets_region_level_contribution_age_i = squeeze(ets_region_level_contribution_age_i);
%     ets_region_level_contribution_age_i = mean(ets_region_level_contribution_age_i, 2);
%     % 换成lme.
%     column_names = {'cofluc_score', 'age'};
%     data_table = table(ets_region_level_contribution_age_i, ages, 'VariableNames', column_names);
%     lme = fitlme(data_table,'cofluc_score ~ 1 + age');
%     corr_t_combin_bin_roi(1, roi_i) = lme.Coefficients.tStat(2);
%     corr_lme_p_combin_bin_roi(1, roi_i) = lme.Coefficients.pValue(2);
% end

% bin_1_2_ages_all = nan(size(ets_region_level_contribution_all, 1));
% bin_1_2_sims_all = nan(size(ets_region_level_contribution_all, 1));
bin_1_2_sims_all = [];
bin_9_12_sims_all = [];
bin_19_20_sims_all = [];

age_6_18 = [];
sex_6_18 = [];

for scan_i=1:size(ets_region_level_contribution_all, 1)

    if ages(scan_i) < 6
        continue;
    end
    if ages(scan_i) >= 19
        continue;
    end

    age_6_18 = [age_6_18 ages(scan_i)];
    sex_6_18 = [sex_6_18 sexs(scan_i)];

    ets_region_level_contribution_scan_i = ets_region_level_contribution_all(scan_i, :, 1:2);
    ets_region_level_contribution_scan_i = squeeze(ets_region_level_contribution_scan_i);
    ets_region_level_contribution_scan_i = mean(ets_region_level_contribution_scan_i, 2);
    [coef, pval] = corr(brain_map_sa_axis, ets_region_level_contribution_scan_i, "type", "Spearman");
    bin_1_2_sims_all = [bin_1_2_sims_all coef];

    ets_region_level_contribution_scan_i = ets_region_level_contribution_all(scan_i, :, 9:12);
    ets_region_level_contribution_scan_i = squeeze(ets_region_level_contribution_scan_i);
    ets_region_level_contribution_scan_i = mean(ets_region_level_contribution_scan_i, 2);
    [coef, pval] = corr(brain_map_sa_axis, ets_region_level_contribution_scan_i, "type", "Spearman");
    bin_9_12_sims_all = [bin_9_12_sims_all coef];

    ets_region_level_contribution_scan_i = ets_region_level_contribution_all(scan_i, :, 19:20);
    ets_region_level_contribution_scan_i = squeeze(ets_region_level_contribution_scan_i);
    ets_region_level_contribution_scan_i = mean(ets_region_level_contribution_scan_i, 2);
    [coef, pval] = corr(brain_map_sa_axis, ets_region_level_contribution_scan_i, "type", "Spearman");
    bin_19_20_sims_all = [bin_19_20_sims_all coef];
end

subplot(3, 1, 1);
scatter(age_6_18, bin_1_2_sims_all);ylim([-1 1]);
subplot(3, 1, 2);
scatter(age_6_18, bin_9_12_sims_all);ylim([-1 1]);
subplot(3, 1, 3);
scatter(age_6_18, bin_19_20_sims_all);ylim([-1 1]);

save(strcat("../../outputs_bins/HCPD/bins_", num2str(num_of_bins), "_top_middle_bottom_sa_sim_individual_spearman"), "age_6_18", "sex_6_18", "bin_1_2_sims_all", "bin_9_12_sims_all", "bin_19_20_sims_all");


%% subject-level.
ets_region_level_contribution_final = ets_region_level_contribution_all(final_labels, :, :);
subj_ids_final = subj_ids(final_labels);
ages_final = ages(final_labels);
sexs_final = sexs(final_labels);
subj_ids_final_unique = unique(subj_ids_final);

bin_1_2_sims_all = [];
bin_9_12_sims_all = [];
bin_19_20_sims_all = [];

age_6_18 = [];
sex_6_18 = [];

for subj_i=1:size(subj_ids_final_unique, 1)

    sess_subj_i = subj_ids_final==subj_ids_final_unique(subj_i);
    disp(sum(sess_subj_i));

    disp(ages_final(sess_subj_i));
    disp(sexs_final(sess_subj_i));

    age_6_18 = [age_6_18 mean(ages_final(sess_subj_i))];
    sex_6_18 = [sex_6_18 mean(sexs_final(sess_subj_i))];

    ets_region_level_contribution_scan_i = ets_region_level_contribution_final(sess_subj_i, :, 1:2);
    ets_region_level_contribution_scan_i = mean(ets_region_level_contribution_scan_i, 1);
    ets_region_level_contribution_scan_i = squeeze(ets_region_level_contribution_scan_i);
    ets_region_level_contribution_scan_i = mean(ets_region_level_contribution_scan_i, 2);
    [coef, ~] = corr(brain_map_sa_axis, ets_region_level_contribution_scan_i, "type", "Spearman");
    bin_1_2_sims_all = [bin_1_2_sims_all coef];

    ets_region_level_contribution_scan_i = ets_region_level_contribution_final(sess_subj_i, :, 9:12);
    ets_region_level_contribution_scan_i = mean(ets_region_level_contribution_scan_i, 1);
    ets_region_level_contribution_scan_i = squeeze(ets_region_level_contribution_scan_i);
    ets_region_level_contribution_scan_i = mean(ets_region_level_contribution_scan_i, 2);
    [coef, pval] = corr(brain_map_sa_axis, ets_region_level_contribution_scan_i, "type", "Spearman");
    bin_9_12_sims_all = [bin_9_12_sims_all coef];

    ets_region_level_contribution_scan_i = ets_region_level_contribution_final(sess_subj_i, :, 19:20);
    ets_region_level_contribution_scan_i = mean(ets_region_level_contribution_scan_i, 1);
    ets_region_level_contribution_scan_i = squeeze(ets_region_level_contribution_scan_i);
    ets_region_level_contribution_scan_i = mean(ets_region_level_contribution_scan_i, 2);
    [coef, pval] = corr(brain_map_sa_axis, ets_region_level_contribution_scan_i, "type", "Spearman");
    bin_19_20_sims_all = [bin_19_20_sims_all coef];
end

subplot(3, 1, 1);
scatter(age_6_18, bin_1_2_sims_all);ylim([-1 1]);
subplot(3, 1, 2);
scatter(age_6_18, bin_9_12_sims_all);ylim([-1 1]);
subplot(3, 1, 3);
scatter(age_6_18, bin_19_20_sims_all);ylim([-1 1]);

save(strcat("../../outputs_bins/HCPD/bins_", num2str(num_of_bins), "_top_middle_bottom_sa_sim_individual_spearman"), "age_6_18", "sex_6_18", "bin_1_2_sims_all", "bin_9_12_sims_all", "bin_19_20_sims_all");


%%
load(strcat("../../outputs_bins/HCPD/bins_", num2str(20), "_top_middle_bottom_sa_sim_individual_spearman"), "age_6_18", "sex_6_18", "bin_1_2_sims_all", "bin_9_12_sims_all", "bin_19_20_sims_all");

column_names = {'sim', 'age'};
data_table = table(bin_1_2_sims_all.', age_6_18.', 'VariableNames', column_names);
lme = fitlme(data_table,'sim ~ 1 + age');

[coef, pval] = corr(bin_1_2_sims_all.', age_6_18.', "type","Spearman");
[coef, pval] = corr(bin_9_12_sims_all.', age_6_18.', "type","Spearman");

%% 
corr_t_combin_bin_roi = nan(3, num_of_rois);
corr_lme_p_combin_bin_roi = nan(3, num_of_rois);

% bin1-2.
for roi_i=1:num_of_rois
    ets_region_level_contribution_age_i = ets_region_level_contribution_all(:, roi_i, 1:2);
    ets_region_level_contribution_age_i = squeeze(ets_region_level_contribution_age_i);
    ets_region_level_contribution_age_i = mean(ets_region_level_contribution_age_i, 2);
    % 换成lme.
    column_names = {'cofluc_score', 'age'};
    data_table = table(ets_region_level_contribution_age_i, ages, 'VariableNames', column_names);
    lme = fitlme(data_table,'cofluc_score ~ 1 + age');
    corr_t_combin_bin_roi(1, roi_i) = lme.Coefficients.tStat(2);
    corr_lme_p_combin_bin_roi(1, roi_i) = lme.Coefficients.pValue(2);
end

% bin9-12.
for roi_i=1:num_of_rois
    ets_region_level_contribution_age_i = ets_region_level_contribution_all(:, roi_i, 9:12);
    ets_region_level_contribution_age_i = squeeze(ets_region_level_contribution_age_i);
    ets_region_level_contribution_age_i = mean(ets_region_level_contribution_age_i, 2);
    % 换成lme.
    column_names = {'cofluc_score', 'age'};
    data_table = table(ets_region_level_contribution_age_i, ages, 'VariableNames', column_names);
    lme = fitlme(data_table,'cofluc_score ~ 1 + age');
    corr_t_combin_bin_roi(2, roi_i) = lme.Coefficients.tStat(2);
    corr_lme_p_combin_bin_roi(2, roi_i) = lme.Coefficients.pValue(2);
end

% bin 19-20.
for roi_i=1:num_of_rois
    ets_region_level_contribution_age_i = ets_region_level_contribution_all(:, roi_i, 19:20);
    ets_region_level_contribution_age_i = squeeze(ets_region_level_contribution_age_i);
    ets_region_level_contribution_age_i = mean(ets_region_level_contribution_age_i, 2);
    % 换成lme.
    column_names = {'cofluc_score', 'age'};
    data_table = table(ets_region_level_contribution_age_i, ages, 'VariableNames', column_names);
    lme = fitlme(data_table,'cofluc_score ~ 1 + age');
    corr_t_combin_bin_roi(3, roi_i) = lme.Coefficients.tStat(2);
    corr_lme_p_combin_bin_roi(3, roi_i) = lme.Coefficients.pValue(2);
end

save(['../../results_bins/DRIVER/' char(site) '/region_level/bins_' char(num2str(num_of_bins)) '_age_effect_lme/combins_cofluctuation_score_age_effect'], "corr_t_combin_bin_roi", "corr_lme_p_combin_bin_roi");


%%
schaefer200_names = readtable("G:datasets\atlas\Schaefer\Schaefer2018_200Parcels_17Networks_order_FSLMNI152_1mm.Centroid_RAS.csv");

site = "HCPD";
num_of_bins = 20;
load(['../../results_bins/DRIVER/' char(site) '/region_level/bins_' char(num2str(num_of_bins)) '_age_effect_lme/combins_cofluctuation_score_age_effect'], "corr_t_combin_bin_roi", "corr_lme_p_combin_bin_roi");

% top
corr_t_combin_bin_roi = corr_t_combin_bin_roi(2, :).';
corr_lme_p_combin_bin_roi = corr_lme_p_combin_bin_roi(2, :).';

pos_t_idx = corr_lme_p_combin_bin_roi<0.005 & corr_t_combin_bin_roi>0;
neg_t_idx = corr_lme_p_combin_bin_roi<0.005 & corr_t_combin_bin_roi<0;

pos_names = schaefer200_names;
neg_names = schaefer200_names;

pos_names(pos_t_idx==0, :) = [];
neg_names(neg_t_idx==0, :) = [];

pos_names.T = corr_t_combin_bin_roi(pos_t_idx==1);
pos_names.p = corr_lme_p_combin_bin_roi(pos_t_idx==1);

neg_names.T = corr_t_combin_bin_roi(neg_t_idx==1);
neg_names.p = corr_lme_p_combin_bin_roi(neg_t_idx==1);

writetable(pos_names, "D:\matlab_proj\MyMatlab\ets\results_bins\DRIVER\HCPD\region_level\bins_20_age_effect_lme/combins_cofluctuation_score_age_effect_middle_pos.csv");
writetable(neg_names, "D:\matlab_proj\MyMatlab\ets\results_bins\DRIVER\HCPD\region_level\bins_20_age_effect_lme/combins_cofluctuation_score_age_effect_middle_neg.csv");


%% association between age effect and SA rank.
num_of_bins=20;
load(['../../results_bins/DRIVER/HCPD/region_level/bins_' char(num2str(num_of_bins)) '_age_effect_lme/combins_cofluctuation_score_age_effect'], "corr_t_combin_bin_roi");
brain_map_sa_axis = readmatrix("../mechanism/s_a_axis_shaefer200x17.csv");

[coef_1, pval_1] = corr(brain_map_sa_axis, squeeze(corr_t_combin_bin_roi(1, :).'));
[coef_2, pval_2] = corr(brain_map_sa_axis, squeeze(corr_t_combin_bin_roi(2, :).'));
[coef_3, pval_3] = corr(brain_map_sa_axis, squeeze(corr_t_combin_bin_roi(3, :).'));


%% age effect of co-fluc score.
site = "HCPD";
load(['../../results_bins/DRIVER/' char(site) '/region_level/bins_' char(num2str(num_of_bins)) '_age_effect_lme/combins_cofluctuation_score_age_effect'], "corr_t_combin_bin_roi", "corr_lme_p_combin_bin_roi");

parcellation = fetch_parcellation('fslr32k', 'schaefer', 200);
[surf_lh, surf_rh] = load_conte69();

labels = ["Bin 1-2", "Bin 9-12", "Bin 19-20"];

for bin_i=1:3
    lme_t_bin_i = corr_t_combin_bin_roi(bin_i, :);
    lme_t_bin_i_old = lme_t_bin_i;
    lme_p_bin_i = corr_lme_p_combin_bin_roi(bin_i, :);
    lme_t_bin_i(lme_p_bin_i > (0.001/200)) = NaN;
    plot_hemispheres([lme_t_bin_i.'], {surf_lh,surf_rh},'parcellation',parcellation, 'labeltext',{labels(bin_i)});
    print(gcf, '-dpng', '-r300', ['../../results_bins/DRIVER/' char(site) '/region_level/bins_' char(num2str(num_of_bins)) '_age_effect_lme/combins_cofluctuation_score_age_effect_bin_' num2str(bin_i)]);
    close;
end


%% bins 20, auto corr of co-fluc score.
num_of_bins = 20;
load(['../../results_bins/DRIVER/HCPD/region_level/bins_' char(num2str(num_of_bins)) '_age_effect_lme/combins_cofluctuation_score_age_effect']);
age_effect_corrs = corr(corr_t_combin_bin_roi.');
figure("Position", [100 100 300 250]);
heatmap(round(age_effect_corrs, 2), "FontSize", 15);
ax=gca;
ax.FontSize=11;
ax.XData = ["Bin 1-2" "Bin 9-12" "Bin 19-20"];
ax.YData = ["Bin 1-2" "Bin 9-12" "Bin 19-20"];
set(gca, 'FontName', 'Arial'); % 设置坐标轴字体为Arial
colormap(colormap_corrs);
colorbar('off');
print(gcf, '-dpng', '-r600', ['../../results_bins/DRIVER/HCPD/region_level/bins_' char(num2str(num_of_bins)) '_age_effect_lme/age_effects_bins_3_corrs']);
close;


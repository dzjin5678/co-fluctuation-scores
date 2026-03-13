%% settings.
% subject ids.
load('../../data/hcp/subj_ids.mat');
% phenotype data.
phenotype_data = readtable("G:\datasets\HCP\Phenotypes\unrestricted_dzjin_7_30_2023_20_45_23.csv");
pipeline = 'proc_regress_fix_ts_mat'; % proc_regress_fixWglob_ts_mat, proc_regress_fix_ts_mat
atlas = 'schaefer200'; % schaefer200, schaefer400, hcpmmp
num_of_rois = 200;
% load atlas.
load('../../atlas/Schaefer200x17.mat');
[~,idxsort] = sort(lab16);
num_of_bins = 20;
session_names = ['REST1_LR'; 'REST1_RL'; 'REST2_LR'; 'REST2_RL'];
high_fd=1;

cal_and_plot_bin_corr(num_of_bins, pipeline, atlas, "ets_ps_region_level_contribution_ratioOfMeans_all", "phase_shift_1_region_level_ratioOfMeans", 1);
cal_and_plot_bin_corr(num_of_bins, pipeline, atlas, "ets_ps_region_level_contribution_ratioOfMeans_all", "phase_shift_2_region_level_ratioOfMeans", 2);
cal_and_plot_cs_sa_association(num_of_bins, pipeline, atlas, "ets_ps_region_level_contribution_ratioOfMeans_all", "phase_shift_1_region_level_ratioOfMeans", 1);
cal_and_plot_cs_sa_association(num_of_bins, pipeline, atlas, "ets_ps_region_level_contribution_ratioOfMeans_all", "phase_shift_2_region_level_ratioOfMeans", 2);


%% functions.
function cal_and_plot_bin_corr(num_of_bins, pipeline, atlas, cs_name, save_dir_name, shift_i)
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
    rest1_lr = load(['../outputs_bins\HCP\phase_random/phase_shift_' num2str(shift_i) '_bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(1, :) '_' atlas '_hfd_1.mat'], cs_name);
    rest1_rl = load(['../outputs_bins\HCP\phase_random/phase_shift_' num2str(shift_i) '_bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(2, :) '_' atlas '_hfd_1.mat'], cs_name);
    rest2_lr = load(['../outputs_bins\HCP\phase_random/phase_shift_' num2str(shift_i) '_bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(3, :) '_' atlas '_hfd_1.mat'], cs_name);
    rest2_rl = load(['../outputs_bins\HCP\phase_random/phase_shift_' num2str(shift_i) '_bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(4, :) '_' atlas '_hfd_1.mat'], cs_name);
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
    figure("Position", [100, 100, 260, 250]);
    imagesc(bins_corr);
    colormap(colormap_corrs);
    % colorbar('fontsize',10, 'TickDirection', 'none');
    clim([-1 1]);
    set(gca, 'XTick', [1, 5, 10, 15, 20], 'TickDir', 'None', ...
        'XTickLabel', ["95-100%", "75-80%", "50-55%", "25-30%", "0-5%"], "XTickLabelRotation", 45);
    set(gca, 'YTick', [1, 5, 10, 15, 20], 'TickDir', 'None', ...
        'YTickLabel', ["95-100%", "75-80%", "50-55%", "25-30%", "0-5%"], "YTickLabelRotation", 0);
    % set(gca, 'XTick', [1, 5, 10, 15, 20], 'TickDir', 'None', ...
    %     'XTickLabel', ["", "", "", "", ""]);
    % set(gca, 'YTick', [1, 5, 10, 15, 20], 'TickDir', 'None', ...
    %     'YTickLabel', ["", "", "", "", ""]);
    ax=gca;hold on;
    ax.FontSize=9;
    save_path_filename=strcat("../results_bins/DRIVER/HCP/",save_dir_name,"/bins_",num2str(num_of_bins),"_",pipeline,"_",atlas,"_corr_HCP_schaefer200x17");
    print(gcf, '-dpng', '-r1000', save_path_filename);
    print(gcf, '-dmeta', save_path_filename);
    close;
end


function cal_and_plot_cs_sa_association(num_of_bins, pipeline, atlas, cs_name, save_dir_name, shift_i)
    % cal group-level cs.
    session_names = ['REST1_LR'; 'REST1_RL'; 'REST2_LR'; 'REST2_RL'];
    rest1_lr = load(['../outputs_bins\HCP\phase_random/phase_shift_' num2str(shift_i) '_bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(1, :) '_' atlas '_hfd_1.mat'], cs_name);
    rest1_rl = load(['../outputs_bins\HCP\phase_random/phase_shift_' num2str(shift_i) '_bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(2, :) '_' atlas '_hfd_1.mat'], cs_name);
    rest2_lr = load(['../outputs_bins\HCP\phase_random/phase_shift_' num2str(shift_i) '_bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(3, :) '_' atlas '_hfd_1.mat'], cs_name);
    rest2_rl = load(['../outputs_bins\HCP\phase_random/phase_shift_' num2str(shift_i) '_bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(4, :) '_' atlas '_hfd_1.mat'], cs_name);
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
    brain_map_sa_axis = readmatrix("../results_bins/mechanism/SA_Axis/s_a_axis_shaefer200x17.csv");
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
    x_1 = linspace(1, num_of_bins, num_of_bins);
    plot(x_1, coefs, '-o', 'Color', [189,189,189]./255, 'LineWidth', 1.5, ...
        'MarkerFaceColor',[99,99,99]./255, 'MarkerSize',6, 'MarkerEdgeColor','none'); hold on;
    ylim([-1 1]);
    xlim([0, num_of_bins+1]);
    xticks([1,5,10,15,20]);
    % xticklabels([" ", " ", " ", " ", " "]);
    xticklabels(["95-100%", "75-80%", "50-55%", "25-30%", "0-5%"]);
    yticks(-1:0.2:1);
    ax=gca;hold on;
    ax.FontSize=13;
    ax.XTickLabelRotation=45;
    set(gca,'Box','off');
    save_filename=strcat("bins_cs_sa");
    save_path=strcat("../results_bins/DRIVER/HCP/",save_dir_name);
    print(gcf, '-dpng', '-r1000', strcat(save_path, "/", save_filename));
    print(gcf, '-dmeta', strcat(save_path, "/", save_filename));
    close;
end

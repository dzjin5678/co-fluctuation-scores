% subject ids.
load('../../data/hcp/subj_ids.mat');
num_of_bins = 20;
pipeline = 'proc_regress_fix_ts_mat'; % proc_regress_fixWglob_ts_mat, proc_regress_fix_ts_mat
atlas = 'schaefer200'; % schaefer200, schaefer400, hcpmmp
num_of_rois = 200;
regressCofounds = 1; % regress global signal, fd, respiration, heart rate.
session_names = ['REST1_LR'; 'REST1_RL'; 'REST2_LR'; 'REST2_RL'];
load("../../atlas/Schaefer200x17.mat");


% Yeo7子网络标准配色列表(+1 TP)
C1=[162, 81, 172]./255;
C2=[120, 154, 192]./255;
C3=[64, 152, 50]./255;
C4=[224, 102, 254]./255;
C5=[169, 169, 169]./255;
C6=[238, 185, 67]./255;
C7=[217, 113, 125]./255;
C8=[0, 0, 128]./255;
colorList=[C1; C2; C3; C4; C5; C6; C7; C8];


%% load eta_out_all.
var_name="rss_eta_out_all";
rest1_lr = load(['../outputs_bins/HCP/bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(1, :) '_' atlas '_hfd_0_regressCofounds_' num2str(regressCofounds) '.mat'], var_name);
rest1_rl = load(['../outputs_bins/HCP/bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(2, :) '_' atlas '_hfd_0_regressCofounds_' num2str(regressCofounds) '.mat'], var_name);
rest2_lr = load(['../outputs_bins/HCP/bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(3, :) '_' atlas '_hfd_0_regressCofounds_' num2str(regressCofounds) '.mat'], var_name);
rest2_rl = load(['../outputs_bins/HCP/bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(4, :) '_' atlas '_hfd_0_regressCofounds_' num2str(regressCofounds) '.mat'], var_name);

% cs_eta_out_all plot.
rest1_lr_out = rest1_lr.rss_eta_out_all;
rest1_rl_out = rest1_rl.rss_eta_out_all;
rest2_lr_out = rest2_lr.rss_eta_out_all;
rest2_rl_out = rest2_rl.rss_eta_out_all;
% Plot ETA for one network k
figure("Position", [100 100 500 360]);
for net_i=1:8
    plot_data = zeros(41, 1);
    for subj_i=1:100
        if size(rest1_lr_out(subj_i).tau)>0
            subj_i_etaMean_data = rest1_lr_out(subj_i).etaMean;
            plot_data = plot_data + mean(subj_i_etaMean_data(:,lab17to8==net_i), 2);
        end
        if size(rest1_rl_out(subj_i).tau)>0
            subj_i_etaMean_data = rest1_rl_out(subj_i).etaMean;
            plot_data = plot_data + mean(subj_i_etaMean_data(:,lab17to8==net_i), 2);
        end
        if size(rest2_lr_out(subj_i).tau)>0
            subj_i_etaMean_data = rest2_lr_out(subj_i).etaMean;
            plot_data = plot_data + mean(subj_i_etaMean_data(:,lab17to8==net_i), 2);
        end
        if size(rest2_rl_out(subj_i).tau)>0
            subj_i_etaMean_data = rest2_rl_out(subj_i).etaMean;
            plot_data = plot_data + mean(subj_i_etaMean_data(:,lab17to8==net_i), 2);
        end

    end
    plot(linspace(-20, 20, 41).', plot_data, 'LineWidth', 1, 'Color', colorList(net_i, :)); hold on;
end
% xline(0,'k-');
% xlabel('Lag (TR)'); ylabel('ETA (z-scored RSS)');
% title(sprintf('Network %d ETA aligned to high RSS_{GLOBAL} events (n=%d)', k, out.nEvents));
ax=gca;hold on;
ax.FontSize=13;
% ax.XTickLabelRotation=45;
set(gca,'Box','off');
save_filename=strcat("rss_eta_out_all");
save_path=strcat("../results_bins/local_global_eta");
print(gcf, '-dpng', '-r1000', strcat(save_path, "/", save_filename));
print(gcf, '-dmeta', strcat(save_path, "/", save_filename));
close;


%% load eta_out_all.
var_name="cs_eta_out_all";
rest1_lr = load(['../outputs_bins/HCP/bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(1, :) '_' atlas '_hfd_0_regressCofounds_' num2str(regressCofounds) '.mat'], var_name);
rest1_rl = load(['../outputs_bins/HCP/bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(2, :) '_' atlas '_hfd_0_regressCofounds_' num2str(regressCofounds) '.mat'], var_name);
rest2_lr = load(['../outputs_bins/HCP/bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(3, :) '_' atlas '_hfd_0_regressCofounds_' num2str(regressCofounds) '.mat'], var_name);
rest2_rl = load(['../outputs_bins/HCP/bins_', num2str(num_of_bins), '_' pipeline '_scan_' session_names(4, :) '_' atlas '_hfd_0_regressCofounds_' num2str(regressCofounds) '.mat'], var_name);

% cs_eta_out_all plot.
rest1_lr_out = rest1_lr.cs_eta_out_all;
rest1_rl_out = rest1_rl.cs_eta_out_all;
rest2_lr_out = rest2_lr.cs_eta_out_all;
rest2_rl_out = rest2_rl.cs_eta_out_all;
% Plot ETA for one network k
figure("Position", [100 100 500 360]);
for net_i=1:8
    plot_data = zeros(41, 1);
    for subj_i=1:100
        if size(rest1_lr_out(subj_i).tau)>0
            subj_i_etaMean_data = rest1_lr_out(subj_i).etaMean;
            plot_data = plot_data + mean(subj_i_etaMean_data(:,lab17to8==net_i), 2);
        end
        if size(rest1_rl_out(subj_i).tau)>0
            subj_i_etaMean_data = rest1_rl_out(subj_i).etaMean;
            plot_data = plot_data + mean(subj_i_etaMean_data(:,lab17to8==net_i), 2);
        end
        if size(rest2_lr_out(subj_i).tau)>0
            subj_i_etaMean_data = rest2_lr_out(subj_i).etaMean;
            plot_data = plot_data + mean(subj_i_etaMean_data(:,lab17to8==net_i), 2);
        end
        if size(rest2_rl_out(subj_i).tau)>0
            subj_i_etaMean_data = rest2_rl_out(subj_i).etaMean;
            plot_data = plot_data + mean(subj_i_etaMean_data(:,lab17to8==net_i), 2);
        end

    end
    plot(linspace(-20, 20, 41).', plot_data, 'LineWidth', 1.5, 'Color', colorList(net_i, :)); hold on;
end
% xline(0,'k-');
% xlabel('Lag (TR)'); ylabel('ETA (z-scored RSS)');
% title(sprintf('Network %d ETA aligned to high RSS_{GLOBAL} events (n=%d)', k, out.nEvents));
ax=gca;hold on;
ax.FontSize=13;
% ax.XTickLabelRotation=45;
set(gca,'Box','off');
save_filename=strcat("cs_eta_out_all");
save_path=strcat("../results_bins/local_global_eta");
print(gcf, '-dpng', '-r1000', strcat(save_path, "/", save_filename));
print(gcf, '-dmeta', strcat(save_path, "/", save_filename));
close;
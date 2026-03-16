%% Setup
addpath(genpath("D:\matlab_proj\arousal-waves-main\Dependencies\"));


%% Import files
load("./outputs_bins/HCP/physio_phase/RV_PLVs_fix_schaefer400.mat", "physios");
physios = reshape(physios,1191,[]);

% load func data. fixWglob fix
func_data_rest1lr = load("./outputs_bins/HCP/bins_20_proc_regress_fixWglob_ts_mat__scan_REST1_LR_schaefer400_hfd_1.mat", "ets_region_rss_all", "ets_global_rss_all");
func_data_rest1rl = load("./outputs_bins/HCP/bins_20_proc_regress_fixWglob_ts_mat__scan_REST1_RL_schaefer400_hfd_1.mat", "ets_region_rss_all", "ets_global_rss_all");
func_data_rest2lr = load("./outputs_bins/HCP/bins_20_proc_regress_fixWglob_ts_mat__scan_REST2_LR_schaefer400_hfd_1.mat", "ets_region_rss_all", "ets_global_rss_all");
func_data_rest2rl = load("./outputs_bins/HCP/bins_20_proc_regress_fixWglob_ts_mat__scan_REST2_RL_schaefer400_hfd_1.mat", "ets_region_rss_all", "ets_global_rss_all");

ets_global_rss_all = nan(1200, 100, 4);
ets_global_rss_all(:, :, 1) = permute(func_data_rest1lr.ets_global_rss_all, [2 1]);
ets_global_rss_all(:, :, 2) = permute(func_data_rest1rl.ets_global_rss_all, [2 1]);
ets_global_rss_all(:, :, 3) = permute(func_data_rest2lr.ets_global_rss_all, [2 1]);
ets_global_rss_all(:, :, 4) = permute(func_data_rest2rl.ets_global_rss_all, [2 1]);
ets_global_rss_all = ets_global_rss_all(5:end-5, :, :);
ets_region_rss_all = nan(1200, 400, 100, 4);
ets_region_rss_all(:, :, :, 1) = permute(func_data_rest1lr.ets_region_rss_all, [3 2 1]);
ets_region_rss_all(:, :, :, 2) = permute(func_data_rest1rl.ets_region_rss_all, [3 2 1]);
ets_region_rss_all(:, :, :, 3) = permute(func_data_rest2lr.ets_region_rss_all, [3 2 1]);
ets_region_rss_all(:, :, :, 4) = permute(func_data_rest2rl.ets_region_rss_all, [3 2 1]);
ets_region_rss_all = ets_region_rss_all(5:end-5, :, :, :);

nanmask = sum(isnan(physios));
nanmask = nanmask | sum(physios)==0;
physios = detrend(physios(:,~nanmask));


%% Define params
clear struct
tr = .72; % sampling interval (s)
Fs = 1/tr; % sampling rate (Hz)
struct.Fs = Fs;
struct.fpass = [0 0.12];
struct.trialave = 1;
struct.tapers = [6 6];


%% Plot coherence spectra
cs_rv_corr_all=nan(1,400);
figure;
for n = 1:400

    net_sigs=squeeze(ets_region_rss_all(:,n,:,:));
    net_global_sigs=ets_global_rss_all(:,:,:);
    net_sigs = net_sigs./net_global_sigs;
    net_sigs = zscore(net_sigs);
    disp(mean(net_sigs(:, 1, 1)));
    disp(mean(net_sigs(:, 1, 2)));
    disp(mean(net_sigs(:, 1, 3)));
    disp(std(net_sigs(:, 1, 1)));
    disp(std(net_sigs(:, 1, 2)));
    disp(std(net_sigs(:, 1, 3)));
    net_sigs=reshape(net_sigs,1191,[]);
    net_sigs = net_sigs(:,~nanmask);

    d = round(7/tr);
    net_sigs = [net_sigs(d+1:end,:);zeros(d,size(net_sigs,2));]; % 7 sec time shift
    [C,phi,S12,S1,S2,f] = coherencyc(net_sigs,physios,struct);
    plot(f,C,'linewidth',0.5) % for coherence plot
    hold on;
    % plot(f,unwrap(phi-phi_global),'linewidth',1) % for phase plot

    r_all=[];
    for i=1:size(net_sigs, 2)
        r_all = [r_all corr(net_sigs(:, i), physios(:, i))];
    end
    cs_rv_corr_all(n)=mean(r_all);
end
ylim([0 0.3]);
set(gca,'fontsize',20,'fontweight','bold')
ax = gca;
ax.LineWidth = 2;
print(gcf, '-dpng', '-r600', "./outputs_bins/HCP/physio_phase/cs_RV_PLVs_fixWglob_schaefer400_coherence_spectra.png");
close;


%% plot corr.
addpath(genpath("D:\Program Files\MATLAB\R2023a\toolbox\dzjin/fieldtrip-master/external/freesurfer"));
addpath(genpath("D:\Program Files\MATLAB\R2023a\toolbox\dzjin\spm12"));
addpath(genpath("D:\Program Files\MATLAB\R2023a\toolbox\dzjin/brainstat_matlab/io"));
addpath(genpath("D:\Program Files\MATLAB\R2023a\toolbox\dzjin/BrainSpace-0.1.10"));
addpath(genpath("D:\Program Files\MATLAB\R2023a\toolbox\dzjin\gifti-main"));
addpath(genpath("../fcn"));
X = ft_read_cifti(char(strcat("../atlas/hcp_fslr32k_cifti/Schaefer2018_400Parcels_17Networks_order.dlabel.nii")), 'mpname', 'array');
vertex_to_roi = X.parcels;
vertex_to_roi(vertex_to_roi==0)=NaN;
% plot phase.
plot_data = vertex_to_roi;
for roi_i=1:400
    plot_data(plot_data==roi_i)=cs_rv_corr_all(roi_i);
end
save_dir=strcat("./outputs_bins/HCP/physio_phase/cs_RV_PLVs_fixWglob_schaefer400");
save_filename = "corr";
plot_my_fslr32k_v_parcel(plot_data, turbo, strcat(save_dir, "/", save_filename));


%% ets rss data.
pipeline = "proc_regress_fix_ts_mat"; % proc_regress_fix_ts_mat proc_regress_fixWglob_ts_mat

func_data_rest1lr = load(strcat("./outputs_bins/HCP/bins_20_proc_regress_fix_ts_mat_scan_REST1_LR_schaefer200_hfd_1_regressCofounds_0.mat"), "ets_region_rss_all", "ets_global_rss_all");
func_data_rest1rl = load(strcat("./outputs_bins/HCP/bins_20_proc_regress_fix_ts_mat_scan_REST1_RL_schaefer200_hfd_1_regressCofounds_0.mat"), "ets_region_rss_all", "ets_global_rss_all");
func_data_rest2lr = load(strcat("./outputs_bins/HCP/bins_20_proc_regress_fix_ts_mat_scan_REST2_LR_schaefer200_hfd_1_regressCofounds_0.mat"), "ets_region_rss_all", "ets_global_rss_all");
func_data_rest2rl = load(strcat("./outputs_bins/HCP/bins_20_proc_regress_fix_ts_mat_scan_REST2_RL_schaefer200_hfd_1_regressCofounds_0.mat"), "ets_region_rss_all", "ets_global_rss_all");

ets_global_rss_all = nan(1191, 100, 4);
ets_global_rss_all(:, :, 1) = permute(func_data_rest1lr.ets_global_rss_all, [2 1]);
ets_global_rss_all(:, :, 2) = permute(func_data_rest1rl.ets_global_rss_all, [2 1]);
ets_global_rss_all(:, :, 3) = permute(func_data_rest2lr.ets_global_rss_all, [2 1]);
ets_global_rss_all(:, :, 4) = permute(func_data_rest2rl.ets_global_rss_all, [2 1]);
% ets_global_rss_all = ets_global_rss_all(5:end-5, :, :, :);

ets_region_rss_all = nan(1191, 200, 100, 4);
ets_region_rss_all(:, :, :, 1) = permute(func_data_rest1lr.ets_region_rss_all, [3 2 1]);
ets_region_rss_all(:, :, :, 2) = permute(func_data_rest1rl.ets_region_rss_all, [3 2 1]);
ets_region_rss_all(:, :, :, 3) = permute(func_data_rest2lr.ets_region_rss_all, [3 2 1]);
ets_region_rss_all(:, :, :, 4) = permute(func_data_rest2rl.ets_region_rss_all, [3 2 1]);
% ets_region_rss_all = ets_region_rss_all(5:end-5, :, :, :);

% Set parameters
num_nodes = 200; % schaefer200x17
tr = .72; % sampling interval (s)
Fs = 1/tr; % sampling rate (Hz)
num_frames = 1191;

% HCP subjects
root_dir = '/mnt/system_v2022/my_data/dzjin_data/hcp/'; % HCP data directory
tasks = {'rfMRI_REST1_LR','rfMRI_REST1_RL','rfMRI_REST2_LR','rfMRI_REST2_RL'};
subjects = load("../data/hcp/subj_ids.mat");
subjects = subjects.subj_ids;
subjects = string(subjects);

% Initialize group matrices for phase-locking values, mean network signals
% for each brain structure (cortex, striatum, thalamus, and cerebellum, and
% physiological time series
%
% (Network signals not required for computing PLVs; to save memory,
% recommended to run script separately to obtain PLVs and to obtain average
% network signals (which are used in later scripts)

plvs_ds = single(nan(num_nodes,numel(subjects),numel(tasks)));
plvs_avg1 = single(nan(num_nodes,numel(subjects),numel(tasks)));
physios_ds = single(nan(num_frames,numel(subjects),numel(tasks)));
physios_avg1 = single(nan(num_frames,numel(subjects),numel(tasks)));

for s = 1:numel(subjects)
    tic
    % load RSS_REGION time series
    subj = num2str(subjects(s));
    disp(['Processing subject ' num2str(s) ' out of ' num2str(numel(subjects))]);
    for t = 1:numel(tasks)
        disp(['Processing task ' num2str(t) ' out of ' num2str(numel(tasks))]);
        task = tasks{t};
        % ets_region_rss_all = nan(1200, 200, 100, 4);
        RSS_REGION = squeeze(ets_region_rss_all(:, :, s, t));
        RSS_GLOBAL = squeeze(ets_global_rss_all(:, s, t));
        RSS_REGION = RSS_REGION./RSS_GLOBAL;
        % Physio
        try
            physio = importdata([root_dir '/physio/' subj '_3T_' task '_Physio_log.txt']);
        catch error
            disp([subj ' ' task 'physio missing:']);
            disp(error)
            continue
        end
        % Get 6 sec windows
        if length(physio)/400<860
            continue
        end
        time_vec_bold = tr*(1:1200)';
        time_vec_phys = (0:length(physio)-1)'/400;
        physio_ds = zeros(size(time_vec_bold));
        physio_avg1 = zeros(size(time_vec_bold));
        physio(:,3) = zscore(physio(:,3));
        for i = 5:length(physio_ds)-4
            % For RV
            [~,phys_start] = min(abs(time_vec_phys-(time_vec_bold(i)-3)));
            [~,phys_end] = min(abs(time_vec_phys-(time_vec_bold(i)+3)));
            physio_ds(i) = std(physio(phys_start:phys_end,2));
            % For HRV
            [pks,locs] = findpeaks(physio(phys_start:phys_end,3),'minpeakdistance',round(400/(180/60)));%,'minpeakwidth',400/(1/(200/60))); % max heart rate = 180 bpm; at 400 Hz, minimum of 100 samples apart
            locs = locs(pks>prctile(physio(phys_start:phys_end,3),60));
            % Avg1(i) = mean(diff(locs))/400;
            physio_avg1(i) = mean(diff(locs))/400;
        end
        physio_ds = physio_ds(5:end-5);
        physio_avg1 = physio_avg1(5:end-5);
        % Get mean network signals
        if size(RSS_REGION,1)~=1191, continue, end
        physios_ds(:,s,t) = physio_ds;
        physios_avg1(:,s,t) = physio_avg1;

        disp(size(physio_ds));
        disp(size(physio_avg1));
        disp(size(RSS_REGION));

        % Filter
        disp('Filtering')
        hp_thresh = .01; % lower bound
        lp_thresh = .05; % higher bound
        [b,a] = butter(2,[hp_thresh,lp_thresh]/(Fs/2));
        physio_ds = single(filtfilt(b,a,double(physio_ds)));
        RSS_REGION = single(filtfilt(b,a,double(RSS_REGION)));
        % Phase-locking values
        h1 = hilbert(physio_ds);
        h2 = hilbert(RSS_REGION);
        plvs_ds(:,s,t) = nanmean(exp(1i*(bsxfun(@minus,unwrap(angle(h1)),unwrap(angle(h2))))));

        % Filter
        disp('Filtering')
        hp_thresh = .01; % lower bound
        lp_thresh = .05; % higher bound
        [b,a] = butter(2,[hp_thresh,lp_thresh]/(Fs/2));
        physio_avg1 = single(filtfilt(b,a,double(physio_avg1)));
        RSS_REGION = single(filtfilt(b,a,double(RSS_REGION)));
        % Phase-locking values
        h1 = hilbert(physio_avg1);
        h2 = hilbert(RSS_REGION);
        plvs_avg1(:,s,t) = nanmean(exp(1i*(bsxfun(@minus,unwrap(angle(h1)),unwrap(angle(h2))))));

        toc
    end
end

save(string(['./outputs_bins/HCP/physio_phase/ETS_REGION_RSS_RV_PLVs_' char(pipeline) '.mat']), ...
    'plvs_ds', 'plvs_avg1', 'physios_ds', 'physios_avg1', '-v7.3');

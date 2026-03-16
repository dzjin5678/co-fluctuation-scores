%% Setup
outdir = './outputs_bins/HCP/physio_phase'; % set directory for saving out images here
mkdir(outdir);

% Set parameters
num_nodes = 200; % schaefer200x17
tr = .72; % sampling interval (s)
Fs = 1/tr; % sampling rate (Hz)
num_frames = 1191;
pipeline = 'proc_regress_fixWglob_ts_mat'; % proc_regress_fix_ts_mat proc_regress_fixWglob_ts_mat

% HCP subjects
root_dir = '/mnt/system_v2022/my_data/dzjin_data/hcp/';
tasks = {'rfMRI_REST1_LR','rfMRI_REST1_RL','rfMRI_REST2_LR','rfMRI_REST2_RL'};
subjects = load("../data/hcp/subj_ids.mat");
subjects = subjects.subj_ids;
subjects = string(subjects);
plvs_ds = single(nan(num_nodes,numel(subjects),numel(tasks)));
plvs_avg1 = single(nan(num_nodes,numel(subjects),numel(tasks)));
physios_ds = single(nan(num_frames,numel(subjects),numel(tasks)));
physios_avg1 = single(nan(num_frames,numel(subjects),numel(tasks)));
bold_all = single(nan(num_frames,num_nodes,numel(subjects),numel(tasks)));
for s = 1:numel(subjects)
    tic
    % load BOLD time series
    subj = num2str(subjects(s));
    disp(['Processing subject ' num2str(s) ' out of ' num2str(numel(subjects))]);
    for t = 1:numel(tasks)
        disp(['Processing task ' num2str(t) ' out of ' num2str(numel(tasks))]);
        task = tasks{t};
        data_dir = [root_dir pipeline '_schaefer200/hcp_' subj '_' task '_schaefer200.mat'];
        try
            roi_ts_data = load(string(data_dir));
        catch error
            disp([subj ' ' task ' missing:']);
            disp(error)
            continue
        end
        BOLD = roi_ts_data.ts_roi';
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
        time_vec_bold = tr*(1:size(BOLD,1))';
        time_vec_phys = (0:length(physio)-1)'/400;
        physio_ds = zeros(size(time_vec_bold));
        physio_avg1 = zeros(size(time_vec_bold));
        for i = 5:length(physio_ds)-4
            % For RV
            [~,phys_start] = min(abs(time_vec_phys-(time_vec_bold(i)-3)));
            [~,phys_end] = min(abs(time_vec_phys-(time_vec_bold(i)+3)));
            physio_ds(i) = std(physio(phys_start:phys_end,2));
            % For HRV
            [pks,locs] = findpeaks(physio(phys_start:phys_end,3),'minpeakdistance',round(400/(180/60)));%,'minpeakwidth',400/(1/(200/60))); % max heart rate = 180 bpm; at 400 Hz, minimum of 100 samples apart
            locs = locs(pks>prctile(physio(phys_start:phys_end,3),60));
            physio_avg1(i) = mean(diff(locs))/400;
        end
        physio_ds = physio_ds(5:end-5);
        physio_avg1 = physio_avg1(5:end-5);
        BOLD = BOLD(5:end-5,:);
        % Get mean network signals
        if size(BOLD,1) ~= 1191, continue, end
        physios_ds(:,s,t) = physio_ds;
        physios_avg1(:,s,t) = physio_avg1;
        bold_all(:, :, s, t) = BOLD;

        % Filter
        disp('Filtering')
        hp_thresh = .01; % lower bound
        lp_thresh = .05; % higher bound
        [b,a] = butter(2,[hp_thresh,lp_thresh]/(Fs/2));
        physio_ds = single(filtfilt(b,a,double(physio_ds)));
        BOLD = single(filtfilt(b,a,double(BOLD)));
        % Phase-locking values
        h1 = hilbert(physio_ds);
        h2 = hilbert(BOLD);
        plvs_ds(:,s,t) = nanmean(exp(1i*(bsxfun(@minus,unwrap(angle(h1)),unwrap(angle(h2))))));

        % Filter
        disp('Filtering')
        hp_thresh = .01; % lower bound
        lp_thresh = .05; % higher bound
        [b,a] = butter(2,[hp_thresh,lp_thresh]/(Fs/2));
        physio_avg1 = single(filtfilt(b,a,double(physio_avg1)));
        BOLD = single(filtfilt(b,a,double(BOLD)));
        % Phase-locking values
        h1 = hilbert(physio_avg1);
        h2 = hilbert(BOLD);
        plvs_avg1(:,s,t) = nanmean(exp(1i*(bsxfun(@minus,unwrap(angle(h1)),unwrap(angle(h2))))));

        toc
    end
end
save(string([outdir '/BOLD_PLVs_' pipeline '_schaefer200.mat']), 'plvs_ds', 'plvs_avg1', 'physios_ds', 'physios_avg1', 'bold_all', '-v7.3');


%% Setup
addpath(genpath("D:\matlab_proj\arousal-waves-main\Dependencies\"));

%% load data.
physio_fix = load("./outputs_bins/HCP/physio_phase/BOLD_PLVs_proc_regress_fix_ts_mat_schaefer200.mat");
physio_fixWglob = load("./outputs_bins/HCP/physio_phase/BOLD_PLVs_proc_regress_fixWglob_ts_mat_schaefer200.mat");


%% Import files
load("./outputs_bins/HCP/physio_phase/RV_PLVs_fixWglob_schaefer400.mat");
physios = reshape(physios,1191,[]);

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
figure;
for n = 1:400
    disp(mean(bold_all(:,n, 1, 1)));
    disp(mean(bold_all(:,n, 1, 2)));
    disp(mean(bold_all(:,n, 1, 3)));
    disp(std(bold_all(:,n, 1, 1)));
    disp(std(bold_all(:,n, 1, 2)));
    disp(std(bold_all(:,n, 1, 3)));
    % tic
    net_sigs=reshape(bold_all(:,n,:,:),1191,[]);

    net_sigs = net_sigs(:,~nanmask);
    d = round(7/tr);
    net_sigs = [net_sigs(d+1:end,:);zeros(d,size(net_sigs,2));]; % 7 sec time shift
    [C,phi,S12,S1,S2,f] = coherencyc(net_sigs,physios,struct);
    plot(f,C,'linewidth',0.5) % for coherence plot
    hold on;
    % plot(f,unwrap(phi-phi_global),'linewidth',1) % for phase plot
    % toc
end
set(gca,'fontsize',20,'fontweight','bold')
ax = gca;
ax.LineWidth = 2;
print(gcf, '-dpng', '-r600', "./outputs_bins/HCP/physio_phase/bold_RV_PLVs_fixWglob_schaefer400_coherence_spectra.png");
close;


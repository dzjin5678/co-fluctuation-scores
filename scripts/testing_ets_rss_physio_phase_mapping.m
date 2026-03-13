addpath(genpath("D:\Program Files\MATLAB\R2023a\toolbox\dzjin/fieldtrip-master/external/freesurfer"));
addpath(genpath("D:\Program Files\MATLAB\R2023a\toolbox\dzjin\spm12"));
addpath(genpath("D:\Program Files\MATLAB\R2023a\toolbox\dzjin/brainstat_matlab/io"));
addpath(genpath("D:\Program Files\MATLAB\R2023a\toolbox\dzjin/BrainSpace-0.1.10"));
addpath(genpath("D:\Program Files\MATLAB\R2023a\toolbox\dzjin\gifti-main"));
addpath(genpath("../../fcn"));

X = ft_read_cifti(char(strcat("../../atlas/hcp_fslr32k_cifti/Schaefer2018_200Parcels_17Networks_order.dlabel.nii")), 'mpname', 'array');
vertex_to_roi = X.parcels;
vertex_to_roi(vertex_to_roi==0)=NaN;


%% load data.
load("../outputs_bins/HCP/physio_phase/BOLD_RV_PLVs_proc_regress_fix_ts_mat_schaefer200.mat");
plv_mean = mean(plvs, [2, 3], "omitmissing");
plv_strength = abs(plv_mean);
plv_phase    = angle(plv_mean);

% plot phase.
plot_data = vertex_to_roi;
for roi_i=1:200
    plot_data(plot_data==roi_i)=plv_phase(roi_i);
end
save_dir=strcat("../outputs_bins/HCP/physio_phase/BOLD_RV_PLVs_proc_regress_fix_ts_mat_schaefer200");
save_filename = "phase";
plot_my_fslr32k_v_parcel(plot_data, turbo, strcat(save_dir, "/", save_filename));

% plot strength.
plot_data = vertex_to_roi;
for roi_i=1:200
    plot_data(plot_data==roi_i)=plv_strength(roi_i);
end
save_dir=strcat("../outputs_bins/HCP/physio_phase/BOLD_RV_PLVs_proc_regress_fix_ts_mat_schaefer200");
save_filename = "strength";
plot_my_fslr32k_v_parcel(plot_data, hot, strcat(save_dir, "/", save_filename));



%% load data.
load("../outputs_bins/HCP/physio_phase/ETS_REGION_RSS_RV_PLVs_proc_regress_fix_ts_mat.mat");

plv_mean = mean(plvs_ds, [2, 3], "omitmissing");
plv_strength = abs(plv_mean);
plv_phase    = angle(plv_mean);
% plot phase.
plot_data = vertex_to_roi;
for roi_i=1:200
    plot_data(plot_data==roi_i)=plv_phase(roi_i);
end
save_dir=strcat("../outputs_bins/HCP/physio_phase/ETS_REGION_RSS_RV_PLVs_proc_regress_fix_ts_mat");
save_filename = "RV_PLVs_ds_phase";
plot_my_fslr32k_v_parcel(plot_data, turbo, strcat(save_dir, "/", save_filename));
% plot strength.
plot_data = vertex_to_roi;
for roi_i=1:200
    plot_data(plot_data==roi_i)=plv_strength(roi_i);
end
save_dir=strcat("../outputs_bins/HCP/physio_phase/ETS_REGION_RSS_RV_PLVs_proc_regress_fix_ts_mat");
save_filename = "RV_PLVs_ds_strength";
plot_my_fslr32k_v_parcel(plot_data, hot, strcat(save_dir, "/", save_filename));

plv_mean = mean(plvs_avg1, [2, 3], "omitmissing");
plv_strength = abs(plv_mean);
plv_phase    = angle(plv_mean);
% plot phase.
plot_data = vertex_to_roi;
for roi_i=1:200
    plot_data(plot_data==roi_i)=plv_phase(roi_i);
end
save_dir=strcat("../outputs_bins/HCP/physio_phase/ETS_REGION_RSS_RV_PLVs_proc_regress_fix_ts_mat");
save_filename = "RV_PLVs_avg1_phase";
plot_my_fslr32k_v_parcel(plot_data, turbo, strcat(save_dir, "/", save_filename));
% plot strength.
plot_data = vertex_to_roi;
for roi_i=1:200
    plot_data(plot_data==roi_i)=plv_strength(roi_i);
end
save_dir=strcat("../outputs_bins/HCP/physio_phase/ETS_REGION_RSS_RV_PLVs_proc_regress_fix_ts_mat");
save_filename = "RV_PLVs_avg1_strength";
plot_my_fslr32k_v_parcel(plot_data, hot, strcat(save_dir, "/", save_filename));


%% load data.
load("../outputs_bins/HCP/physio_phase/ETS_REGION_RSS_GLOBAL_RSS_PLVs_proc_regress_fixWglob_ts_mat.mat");
plv_mean = mean(plvs, [2, 3], "omitmissing");
plv_strength = abs(plv_mean);
plv_phase    = angle(plv_mean);

% plot phase.
plot_data = vertex_to_roi;
for roi_i=1:200
    plot_data(plot_data==roi_i)=plv_phase(roi_i);
end
save_dir=strcat("../outputs_bins/HCP/physio_phase/ETS_REGION_RSS_GLOBAL_RSS_PLVs_proc_regress_fixWglob_ts_mat");
save_filename = "RV_PLVs_phase";
plot_my_fslr32k_v_parcel(plot_data, turbo, strcat(save_dir, "/", save_filename));

% plot strength.
plot_data = vertex_to_roi;
for roi_i=1:200
    plot_data(plot_data==roi_i)=plv_strength(roi_i);
end
save_dir=strcat("../outputs_bins/HCP/physio_phase/ETS_REGION_RSS_GLOBAL_RSS_PLVs_proc_regress_fixWglob_ts_mat");
save_filename = "RV_PLVs_strength";
plot_my_fslr32k_v_parcel(plot_data, hot, strcat(save_dir, "/", save_filename));


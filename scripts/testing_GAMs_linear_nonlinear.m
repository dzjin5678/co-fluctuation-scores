addpath(genpath("D:\Program Files\MATLAB\R2023a\toolbox\dzjin/fieldtrip-master/external/freesurfer"));
addpath(genpath("D:\Program Files\MATLAB\R2023a\toolbox\dzjin\spm12"));
addpath(genpath("D:\Program Files\MATLAB\R2023a\toolbox\dzjin/brainstat_matlab/io"));
addpath(genpath("D:\Program Files\MATLAB\R2023a\toolbox\dzjin/BrainSpace-0.1.10"));
addpath(genpath("D:\Program Files\MATLAB\R2023a\toolbox\dzjin\gifti-main"));
addpath(genpath("../../fcn"));


%% load data. 
anova_results = readmatrix("../results_bins/DRIVER/HCP/1_region_level_ratioOfMeans/" + ...
    "gam_nonlinear_linear_compare_high_fd-0_num_of_bins-20_atlas-schaefer200_regress_cofounds-0_gsr-proc_regress_fix_ts_mat.csv");
dev_values = anova_results(:, 1);
p_values = anova_results(:, 2);
p_values_corrected = p_values*200;

sig_brain_map = dev_values;
sig_brain_map(p_values_corrected>=0.05)=NaN;

not_sig_brain_map = dev_values;
not_sig_brain_map(p_values_corrected<0.05)=NaN;


%% plot results.
X = ft_read_cifti(char(strcat("../../atlas/hcp_fslr32k_cifti/Schaefer2018_200Parcels_17Networks_order.dlabel.nii")), 'mpname', 'array');
vertex_to_roi = X.parcels;
vertex_to_roi(vertex_to_roi==0)=NaN;


%% colorbar.
cb = hot;
cb = cb(20:256, :);


%% plot (sig).
plot_data = vertex_to_roi;
for roi_i=1:200
    plot_data(plot_data==roi_i)=sig_brain_map(roi_i);
end
save_dir=strcat("../results_bins/DRIVER/HCP/1_region_level_ratioOfMeans/");
save_filename = "gam_nonlinear_linear_compare_high_fd-0_num_of_bins-20_atlas-schaefer200_regress_cofounds-0_gsr-proc_regress_fix_ts_mat_sig";
plot_my_fslr32k_v_parcel(plot_data, cb, strcat(save_dir, "/", save_filename));


%% plot (nosig).
plot_data = vertex_to_roi;
for roi_i=1:200
    plot_data(plot_data==roi_i)=not_sig_brain_map(roi_i);
end
save_dir=strcat("../results_bins/DRIVER/HCP/1_region_level_ratioOfMeans/");
save_filename = "gam_nonlinear_linear_compare_high_fd-0_num_of_bins-20_atlas-schaefer200_regress_cofounds-0_gsr-proc_regress_fix_ts_mat_not_sig";
plot_my_fslr32k_v_parcel(plot_data, cb, strcat(save_dir, "/", save_filename));


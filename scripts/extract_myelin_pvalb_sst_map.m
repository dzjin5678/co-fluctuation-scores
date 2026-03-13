addpath(genpath("D:\Program Files\MATLAB\R2023a\toolbox\dzjin/fieldtrip-master/external/freesurfer"));
addpath(genpath("D:\Program Files\MATLAB\R2023a\toolbox\dzjin\spm12"));
addpath(genpath("D:\Program Files\MATLAB\R2023a\toolbox\dzjin/brainstat_matlab/io"));
addpath(genpath("D:\Program Files\MATLAB\R2023a\toolbox\dzjin/BrainSpace-0.1.10"));
addpath(genpath("D:\Program Files\MATLAB\R2023a\toolbox\dzjin\gifti-main"));
% load the conte69 hemisphere surfaces
[surf_lh, surf_rh] = load_conte69();


%% myelination.
myelin_map = ft_read_cifti(['G:\datasets\atlas\HCPMMP\fromBALSA\Q1-Q6_RelatedParcellation210\MNINonLinear\' ...
    'fsaverage_LR32k\Q1-Q6_RelatedParcellation210.MyelinMap_BC_MSMAll_2_d41_WRN_DeDrift.32k_fs_LR.dscalar.nii'], ...
    'mapname', 'field');
myelin_map = myelin_map.myelinmap_bc_msmall_2_d41_wrn_dedrift;

myelin_map_plot = myelin_map;
myelin_map_plot(isnan(myelin_map))=0;

plot_hemispheres([myelin_map_plot], {surf_lh,surf_rh}, 'labeltext',{'myelin'});
print(gcf, '../results_bins/mechanism/myelin/myelin_fslr32k.png', '-dpng', '-r600');
close;

parcellation_glasser = fetch_parcellation('fslr32k', 'glasser', 360);
myelin_map_glasser = nan(360, 1);
for roi_i=1:360
    myelin_map_glasser(roi_i, 1) = mean(myelin_map(parcellation_glasser==roi_i, 1), "omitmissing");
end
plot_hemispheres(myelin_map_glasser, {surf_lh,surf_rh},'parcellation',parcellation_glasser, 'labeltext',{'myelin'});
print(gcf, '../results_bins/mechanism/myelin/myelin_fslr32k_glaaser.png', '-dpng', '-r600');
close;

parcellation_schaefer200 = fetch_parcellation('fslr32k', 'schaefer', 200);
myelin_map_schaefer200x17 = nan(200, 1);
for roi_i=1:200
    myelin_map_schaefer200x17(roi_i, 1) = mean(myelin_map(parcellation_schaefer200==roi_i, 1), "omitmissing");
end
plot_hemispheres(myelin_map_schaefer200x17, {surf_lh,surf_rh},'parcellation',parcellation_schaefer200, 'labeltext',{'myelin'});
print(gcf, '../results_bins/mechanism/myelin/myelin_fslr32k_schaefer200x17.png', '-dpng', '-r600');
close;

parcellation_schaefer400 = fetch_parcellation('fslr32k', 'schaefer', 400);
myelin_map_schaefer400x17 = nan(400, 1);
for roi_i=1:400
    myelin_map_schaefer400x17(roi_i, 1) = mean(myelin_map(parcellation_schaefer400==roi_i, 1), "omitmissing");
end
plot_hemispheres(myelin_map_schaefer400x17, {surf_lh,surf_rh},'parcellation',parcellation_schaefer400, 'labeltext',{'myelin'});
print(gcf, '../results_bins/mechanism/myelin/myelin_fslr32k_schaefer400x17.png', '-dpng', '-r600');
close;

save('../results_bins/mechanism/myelin/myelin_data.mat', "myelin_map_glasser", "myelin_map_schaefer200x17", "myelin_map_schaefer400x17", "myelin_map");


%% sst, pvalm, pvalm-sst.
atlas="glasser360";
abha_expression_data = readtable(strcat("G:\datasets\AHBA\expression_ahba_", atlas, ".csv"));
abha_expression_data_pvalb_sst_glasser360 = abha_expression_data.PVALB - abha_expression_data.SST;
parcellation_glasser = fetch_parcellation('fslr32k', 'glasser', 360);
plot_hemispheres(abha_expression_data_pvalb_sst_glasser360, {surf_lh,surf_rh},'parcellation',parcellation_glasser, 'labeltext',{'pvalb-sst'});
print(gcf, strcat("../results_bins/mechanism/pvalb_sst/pvalb_", atlas, ".png"), '-dpng', '-r600');
close;

atlas="schaefer200x17";
abha_expression_data = readtable(strcat("G:\datasets\AHBA\expression_ahba_", atlas, ".csv"));
abha_expression_data_pvalb_sst_schaefer200x17 = abha_expression_data.PVALB - abha_expression_data.SST;
parcellation_schaefer200 = fetch_parcellation('fslr32k', 'schaefer', 200);
plot_hemispheres(abha_expression_data_pvalb_sst_schaefer200x17, {surf_lh,surf_rh},'parcellation',parcellation_schaefer200, 'labeltext',{'pvalb-sst'});
print(gcf, strcat("../results_bins/mechanism/pvalb_sst/pvalb_", atlas, ".png"), '-dpng', '-r600');
close;

atlas="schaefer400x17";
abha_expression_data = readtable(strcat("G:\datasets\AHBA\expression_ahba_", atlas, ".csv"));
abha_expression_data_pvalb_sst_schaefer400x17 = abha_expression_data.PVALB - abha_expression_data.SST;
parcellation_schaefer400 = fetch_parcellation('fslr32k', 'schaefer', 400);
plot_hemispheres(abha_expression_data_pvalb_sst_schaefer400x17, {surf_lh,surf_rh},'parcellation',parcellation_schaefer400, 'labeltext',{'pvalb-sst'});
print(gcf, strcat("../results_bins/mechanism/pvalb_sst/pvalb_", atlas, ".png"), '-dpng', '-r600');
close;

save('../results_bins/mechanism/pvalb_sst/pvalb_sst_data.mat', ...
    "abha_expression_data_pvalb_sst_glasser360", ...
    "abha_expression_data_pvalb_sst_schaefer200x17", ...
    "abha_expression_data_pvalb_sst_schaefer400x17");


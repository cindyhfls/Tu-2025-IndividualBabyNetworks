% main_PFM_Lynch_postprocessing
clear;close all;clc;
addpath(genpath('/data/wheelock/data1/people/Cindy/PublicRepo/NLA_toolbox_070319'))
addpath(genpath('/data/wheelock/data1/people/Cindy/PublicRepo/PFM')); % Repository from Lynch 2024 Nature
addpath('./subjlists');
addpath('./CIFTI_read_save');

WorkbenchBinary = 'wb_command';
load('priors.mat');% load the priors (https://github.com/cindyhfls/PFM/blob/main/PFM-Tutorial/Utilities/priors.mat)
empty_cifti = load('empty_cifti.mat');

here = pwd;
curr_fs = 'bestk23'
IM = smartload(['./template_matching_results/IM_consensus_eLABE_Y2_Y3_assn_',curr_fs,'.mat']);
%% Group average assignment
input_dir ='./template_matching_results/cifti/bestk23_groupsummary'
cd(input_dir);
group_Spatial = ft_read_cifti_mod('eLABE_Y2_Y3_templates_network_probability_maps.dtseries.nii');
group_FC = ft_read_cifti_mod('eLABE_Y2_Y3_templates_network_FC_maps.dtseries.nii');
Output = 'eLABE_Y2_Y3_templates_network'
Group = struct();
Group.Spatial = group_Spatial.data;
Group.FC = group_FC.data;
Group.NetworkLabels = IM.Nets;
Group.NetworkColors = IM.cMap;

pfm_identify_networks_precalculated(Group,Priors,Output,pwd,WorkbenchBinary,empty_cifti); 
% N.B. the code to plot legend is in plot_hierarchical_similarity(Group);

% Also plot the FC similarity heatmap organized into hierarchical clusters
plot_hierarchical_similarity(Group);

return

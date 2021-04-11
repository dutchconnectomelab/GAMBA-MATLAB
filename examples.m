clc, clear, close all;
filepath = fileparts(mfilename('fullpath'));

% Select an example to run
exampleID = 1; % 1, 2, 3, 4, 5


%% Example 1
% I have an imaging map (.nii file) and a gene set. I want to test if the 
% imaging pattern correlates to the expression pattern.

if exampleID == 1
    disp(strcat('Example 1. Examining association between Alzheimer', ...
        "'", 's VBM map and APOE, APP, PSEN2 expression patterns.'));

    % 1.1 Co-register the imaging file to MNI152 space
    input_img_file = fullfile(filepath, 'src', 'examples', ...
        'alzheimers_ALE.nii.gz'); % meta-analysis of alzheimer's VBM studies
    input_img_anat_file = fullfile(filepath, 'src', 'examples', ...
        'Colin27_T1_seg_MNI_2x2x2.nii.gz'); % anatomical file in the same space
    ref_img_file = fullfile(filepath, 'src', 'atlas', ...
        'brain.nii.gz'); % reference file in MNI152 space
    reg_file = fullfile(filepath, 'output', 'registration.mat');
    output_img_file = fullfile(filepath, 'output', 'coreg_alzheimers_ALE.nii.gz');

    % here using FSL flirt
    system(['flirt -in ', input_img_anat_file, ' -ref ', ref_img_file, ...
        ' -omat ', reg_file]);

    system(['flirt -in ', input_img_file, ' -ref ', ref_img_file, ...
        ' -applyxfm -init ', reg_file, ' -out ', output_img_file]);

    % 1.2 Group the imaging map into brain regions (here regional mean)
    res_Y = group_regions(output_img_file, 'DK114');

    % 1.3 Association between the imaging profile and the expression profile
    img_data = res_Y.data(contains(res_Y.regionDescriptions, 'ctx-lh-'));

    % 1.3.1 Null-coexpression model
    res_nullcoexp = permutation_null_coexp(img_data, {'APOE', 'APP', 'PSEN2'});

    % 1.3.2 Null-brain model
    res_nullbrain = permutation_null_brain(img_data, {'APOE', 'APP', 'PSEN2'});

    % 1.3.3 Null-spin model
    res_nullspin = permutation_null_spin(img_data, {'APOE', 'APP', 'PSEN2'});
end


%% Example 2
% I have an imaging data matrix (region by feature), a gene expression data 
% matrix (region by gene), and a gene set. I want to test if the imaging 
% pattern correlates to the expression pattern.

% UNDER CONSTRUCTION

if exampleID == 2

    load('data/example_HSE.mat');
    imgDescriptions = {'NOS','FA','SD','Degree','FC'};

    % null-coexpression model
    res_nullcoexp = permutation_null_coexp(img_data, geneset);

    % null-brain model
    res_nullbrain = permutation_null_brain(img_data, geneset);

    % null-spin model
end


%% Example 3
% I have a gene-set. I want to test in which brain regions the gene-set is 
% over-expressed.

% UNDER CONSTRUCTION


%% Example 4
% I have a gene expression matrix and a gene-set. I want to test in which 
% brain regions the gene-set is over-expressed.

% UNDER CONSTRUCTION


%% Example 5
% I have an imaging map (.nii file) and I want to look for the most
% correlated genes

% UNDER CONSTRUCTION








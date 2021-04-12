clc, clear, close all;
filepath = fileparts(mfilename('fullpath'));
addpath(genpath(filepath));

% Select an example to run
exampleID = 1; % 1, 2, 3, 4, 5


%% Example 1
% I have an imaging map (a nifti file) and a gene set. I want to test if 
% the imaging pattern correlates to the expression pattern.'
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
if exampleID == 2
    disp(strcat('Example 2. Examining association between connectome ', ...
    'metrics and the expression pattern of supragranular-enriched genes.'));

    % 2.1 load example data
    load('src/examples/example_conn_5k_genes.mat');
    % expression data matrix: 57 regions by 5000 genes
    % imaging data matrix: 57 regions by 5 phenotypes
    % gene set: 19 Human-supragranular genes
    
    % 2.2 Association between the imaging profile and the expression profile
    % 2.2.1 null-coexpression model
    res_nullcoexp = permutation_null_coexp(img_data, geneset, ...
        gene_expression, gene_symbols);

    % 2.2.2 null-brain model
    res_nullbrain = permutation_null_brain(img_data, geneset, ...);
        gene_expression, gene_symbols);
    
    % 2.2.3 null-spin model
    % ATTENTION: null-spin model ONLY works for DK114 atlas. If you use
    % other atlas, please refer to Alexander-Bloch et al., 2018 to first
    % generate gene expression data for 'spinned' atlas.
    res_nullspin = permutation_null_spin(img_data, geneset);
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








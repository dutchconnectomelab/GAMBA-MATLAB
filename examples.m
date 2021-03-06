clc, clear, close all;
filepath = fileparts(mfilename('fullpath'));
addpath(genpath(filepath));

% Select an example to run
exampleID = 3; % 1, 2, 3, 4, 5


%% Example 1
% Research question:
% - I have an imaging map (a nifti file) and a gene set. I want to test if 
%   the imaging pattern correlates to the expression pattern.'
if exampleID == 1
    disp(strcat('Example 1. Examining association between Alzheimer', ...
        "'", 's VBM map and APOE, APP, PSEN2 expression patterns.'));

    % 1.1 Co-register the imaging file to MNI152
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
% Research question:
% - I have an imaging data matrix (region by feature), a gene expression data 
%   matrix (region by gene), and a GOI. I want to test if the imaging 
%   pattern correlates to the expression pattern.
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
% Research question:
% - I have a GOI. I want to test in which brain regions the GOI is 
%   over-expressed.

if exampleID == 3
    disp(strcat('Example 3. Testiassociation between connectome ', ...
    'metrics and the expression pattern of supragranular-enriched genes.'));

    % 3.1 load example data
    load('src/examples/example_conn_5k_genes.mat', 'geneset');
    % gene set: 19 Human-supragranular genes
    
    % 3.2.1 null-coexp model
    res_nullcoexp = permutation_expression_null_coexp(geneset);
    
    % 3.2.2 null-brain model
    res_nullbrain = permutation_expression_null_brain(geneset);
    
    % 3.2.3 null-spin model
    res_nullspin = permutation_expression_null_spin(geneset);
end

%% Example 4
% Research question:
% - I have a gene expression matrix and a GOI. I want to test in which 
%   brain regions the gene-set is over-expressed.

if exampleID == 4
    disp(strcat('Example 4. Testiassociation between connectome ', ...
    'metrics and the expression pattern of supragranular-enriched genes.'));

    % 4.1 load example data
    load('src/examples/example_conn_5k_genes.mat', 'geneset', ...
        'gene_expression', 'gene_symbols');
    % gene set: 19 Human-supragranular genes
    
    % 4.2.1 null-coexp model
    res_nullcoexp = permutation_expression_null_coexp(geneset, ...
        gene_expression, gene_symbols);
    
    % 4.2.2 null-brain model
    res_nullbrain = permutation_expression_null_brain(geneset, ...
        gene_expression, gene_symbols);
    
    % 4.2.3 null-spin model
    res_nullspin = permutation_expression_null_spin(geneset, ...
        gene_expression, gene_symbols);
end

%% Example 5
% Research question:
% - I have an imaging map (.nii file) and I want to look for the most
%   correlated genes

if exampleID == 5
    % 5.1 Co-register the imaging file to MNI152 space
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

    % 5.2 Group the imaging map into brain regions (here regional mean)
    res_Y = group_regions(output_img_file, 'DK114');
    
    % 5.3 Using the null-spin model to look for the most correlated genes    
    res_nullspin = permutation_null_spin_correlated_genes(res_Y.data(1:57));
end

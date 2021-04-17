# GAMBA-MATLAB

A Matlab toolbox to study whether the expression of the gene(s) of interest (GOI) and neuroimaging-derived brain phenotypes show overlapped spatial patterns. Different statistical **null models** are available to examine both ***gene specificity*** and ***spatial specificity***. This toolbox is an extension of the web application [GAMBA](www.dutchconnectomelab.nl/GAMBA). 

For details, please see:

> Wei Y. et al. (2021), Statistical testing and annotation of gene transcriptomic-neuroimaging associations, bioRxiv

## Installation

### Requirements

Before you start, make sure you have **Matlab** on your machine. 

No other software required, but in our examples we use [FSL - FMRIB Software Library](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki) to perform coregistration. You can use other equivalent tools for the same purpose. Just modify the `example.m` scripts if needed. 

### Download

Through the command line in terminal:

`git clone https://github.com/yongbin-wei/GAMBA_MATLAB.git`.

Alternatively, you can click 'Code -- Download ZIP' and unzip the downloaded file.

### Extract required data

Run `install.m` in Matlab to retrieve all required data.

Alternatively, you can download files here: https://www.dropbox.com/sh/psfudnzktyd0860/AABtx7ESvEphO60dcV_xbQ4qa?dl=0. Unzip downloaded files to the folder `src/`.

## Examples

We use examples to show the utility of this toolbox. Examples cover the following usage of the toolbox:

 "I have an imaging map (e.g., a nifti file) and a gene set. I want to test if the imaging pattern correlates to the gene expression pattern."
"I have an imaging data matrix (region by feature), a gene expression data matrix (region by gene), and a gene set. I want to test if the imaging pattern correlates to the gene expression pattern."
"I have a gene-set. I want to test in which brain regions the gene-set is differentially expressed."
"I have a gene expression data matrix and a gene-set. I want to test in which brain regions the gene-set is differentially expressed."
"I have an imaging map (e.g., a nifti file) and I want to look for the most correlated genes."

Scripts concerning each question are included in `examples.m`. A detailed tutorial is can be found below:

#### 1. "I have an imaging map (a nifti file) and a gene set. I want to test if the imaging pattern correlates to the gene expression pattern."

This example uses VBM meta-analysis result -- a brain map showing the vulnerability of brain volume -- for Alzheimer's disease (AD) and examines whether the VBM pattern is associated with the pattern of brain gene expression of three AD risk genes ('APOE', 'APP', 'PSEN2') through GAMBA null-models. Preprocessed gene expression data from the Allen Human Brain Atlas (https://human.brain-map.org) will be used.

In the first step, you need to co-register the imaging file to '/src/atlas/brain.nii.gz' (which is the MNI152 atlas used for brain parcellation and segmentation). The following codes will do the co-registration using FSL FLIRT. If you do not use FSL, please adjust the codes to *coregister the input imaging file to /src/atlas/brain.nii.gz* and then continue.

```matlab    
    input_img_file = fullfile(filepath, 'src', 'examples', 'alzheimers_ALE.nii.gz');
    input_img_anat_file = fullfile(filepath, 'src', 'examples', 'Colin27_T1_seg_MNI_2x2x2.nii.gz'); % anatomical file in the same space
    ref_img_file = fullfile(filepath, 'src', 'atlas', 'brain.nii.gz'); % reference file in MNI152 space
    reg_file = fullfile(filepath, 'output', 'registration.mat');
    output_img_file = fullfile(filepath, 'output', 'coreg_alzheimers_ALE.nii.gz');

    system(['flirt -in ', input_img_anat_file, ' -ref ', ref_img_file, ...
        ' -omat ', reg_file]);
    system(['flirt -in ', input_img_file, ' -ref ', ref_img_file, ...
        ' -applyxfm -init ', reg_file, ' -out ', output_img_file]);
```

Next, you need to compute region-wise measurements based on the registered imaging map. Here, voxels within each brain region in the DK-114 atlas are averaged.

    res_Y = group_regions(output_img_file, 'DK114');

After this you will have a 114 x 1 brain-phenotypic data array, together with region descriptions of the 114 regions. Due to that only AHBA transcriptomic data from the left hemisphere is used, you need to extract brain-phenotypic data for the left hemisphere.

    img_data = res_Y.data(contains(res_Y.regionDescriptions, 'ctx-lh-'));

Now, you will be able to examine transcriptome-neuroimaging associations. It is highly recommended to first examine whether there is a **spatially specific** association between the brain-phenotypic pattern and the mean expression pattern of the input gene set of interest (here 'APOE', 'APP', 'PSEN2').

    res_nullspin = permutation_null_spin(img_data, {'APOE', 'APP', 'PSEN2'});

Next, you're recommended to examine the level of **gene specificity**, namely, the observed association goes beyond what you can expect by any other gene sets with similar co-expression levels.

    res_nullcoexp = permutation_null_coexp(img_data, {'APOE', 'APP', 'PSEN2'});

Afterward, it would be also more informative to examine whether these genes are **brain-enriched** genes and whether the observed association goes beyond what you can expect by any **brain-enriched** genes.

    res_nullbrain = permutation_null_brain(img_data, {'APOE', 'APP', 'PSEN2'});

By doing the above analysis, you will know whether there is a significant association between your brain-phenotypic pattern and the gene set of your interest, and, more importantly, whether the association is **spatial specific** and **gene-specific**.




#### 2. "I have an imaging data matrix (region by feature), a gene expression data matrix (region by gene), and a gene set. I want to test if the imaging pattern correlates to the expression pattern."

This example examines whether the gene expression pattern of 19 human-supragranular genes is associated with regional connectome metrics. The input data include: 1) a gene expression data matrix (57 regions x 5000 genes, a subset of the entire AHBA data); 2) an imaging data matrix (57 regions by 5 phenotypes, such as NOS-, FA-, SD-weighted nodal strength, nodal degree, and FC strength); 3) a gene set (19 genes).

NOTE: null-spin model ONLY works for DK114 atlas. If you use another atlas, please refer to Alexander-Bloch et al. (2018) to first generate gene expression data matrices for the 'spun' atlases.


#### 3. "I have a gene-set. I want to test in which brain regions the gene-set is differentially expressed."

In the third example, we show how to use this toolbox to examine in which brain regions the input GOI is differentially expressed considering different types of null models. The human-supragranular-enriched (HSE) genes are used here.

We load the gene set first:

`load('src/examples/example_conn_5k_genes.mat', 'geneset');`
 
Now we have a gene set of 19 HSE genes. We then use the null-coexpression model to examine whether the mean expression of HSE genes is significantly higher/lower than random genes with the similar coexpression conserved for each brain region.

`res_nullcoexp = permutation_expression_null_coexp(geneset);`
    
Next, we use the null-brain model to examine whether the mean expression of HSE genes is significantly higher/lower than random brain-enriched genes in general for each brain region.

`res_nullbrain = permutation_expression_null_brain(geneset);`

Using the above null models we will be able to examine the level of *gene specificity*. Another question we may ask is whether the mean expression of the GOI is specifically higher in one region in contrast to random brain regions. We can answer this “spatial specificity” question using the null-spin model.

`res_nullspin = permutation_expression_null_spin(geneset);`

Now, we have results for the questions we ask.


#### 4. "I have a gene expression data matrix and a gene-set. I want to test in which brain regions the gene-set is over-expressed."

In the fourth example, we show how to use this toolbox to examine in which brain regions the input GOI is differentially expressed considering different types of null models, given a customized gene expression matrix. The human-supragranular-enriched (HSE) genes are used as an example here.

You need to load the example data first. Gene symbols of the gene set, the gene expression matrix of all genes, and symbols of all genes are needed.

    load('src/examples/example_conn_5k_genes.mat', 'geneset', 'gene_expression', 'gene_symbols');

Then you can examine in which brain region the input GOI is differentially expressed compared to random genes, using the null-coexpressed-gene model where the null distribution of gene expression of random genes with similar coexpression level is generated and is used in permutation testing.

    res_nullcoexp = permutation_expression_null_coexp(geneset, gene_expression, gene_symbols);

You can also examine the same question using the null-brain-gene model, where the null distribution of gene expression of brain-expressed genes is generated.

    res_nullbrain = permutation_expression_null_brain(geneset, gene_expression, gene_symbols);

For the other research question -- "in which brain region the input GOI is differentially expressed compared to random brain regions" -- you can use the following code only if you work on the Desikan-Killiany atlas with 57 left-hemisphere regions. Otherwise please use the spin model to generate random gene expression matrices first.

    res_nullspin = permutation_expression_null_spin(geneset, gene_expression, gene_symbols);


#### 5. "I have an imaging map (.nii file) and I want to look for the most correlated genes"



In the last example, we show how to use our toolbox to look for the most correlated genes given an imaging map. This example uses VBM meta-analysis results for Alzheimer's disease (AD).

You need to first co-register the imaging file to MNI152 space.

    input_img_file = fullfile(filepath, 'src', 'examples', 'alzheimers_ALE.nii.gz');
    input_img_anat_file = fullfile(filepath, 'src', 'examples', 'Colin27_T1_seg_MNI_2x2x2.nii.gz'); % anatomical file in the same space
    ref_img_file = fullfile(filepath, 'src', 'atlas', 'brain.nii.gz'); % reference file in MNI152 space
    reg_file = fullfile(filepath, 'output', 'registration.mat');
    output_img_file = fullfile(filepath, 'output', 'coreg_alzheimers_ALE.nii.gz');

    system(['flirt -in ', input_img_anat_file, ' -ref ', ref_img_file, ' -omat ', reg_file]);
    system(['flirt -in ', input_img_file, ' -ref ', ref_img_file, ' -applyxfm -init ', reg_file, ' -out ', output_img_file]);

Next, you need to compute region-wise measurements based on the registered imaging map. Here, voxels within each brain region in the DK-114 atlas are averaged.

    res_Y = group_regions(output_img_file, 'DK114');

Using the resulted phenotypic data you can now do correlations between single-gene expression profiles and the phenotypic profile, and then do permutations per gene using the null-spatial model.

    res_nullspin = permutation_null_spin_correlated_genes(img_data);


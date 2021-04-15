# GAMBA-MATLAB
A Matlab toolbox to perform statistical analyses examining whether the gene expression profile(s) of the input gene(s) and neuroimaging-derived brain phenotypes show overlapped spatial patterns. Different statistical null models are available to examine both gene specificity and spatial specificity. This toolbox is an extension of the web application GAMBA (www.dutchconnectomelab.nl/GAMBA).  

For details, please see:

> Wei Y. et al. (2021), Statistical testing and annotation of gene transcriptomic-neuroimaging associations, bioRxiv

## Installation
### Requirements
Before you start, make sure you have **Matlab** on your machine. There is no other requirement, but in our examples we use **FSL - FMRIB Software Library** to perform coregistration. It's not necessary of course, you can choose other equivalent tools for the same purpose, just modify the `examples.m` scripts if needed. 

### Download
Download the Github repository through the command line 

`git clone https://github.com/yongbin-wei/GAMBA_MATLAB.git`.

Alternatively, you can click 'Code -- Download ZIP' and unzip the downloaded file.


### Extract required data
Run `install.m` in Matlab to retrieve all required data.

Alternatively, you can download files here: https://www.dropbox.com/sh/psfudnzktyd0860/AABtx7ESvEphO60dcV_xbQ4qa?dl=0. Unzip downloaded files to the folder `src/`.


## Examples
We use examples to show the utility of this toolbox. Condition with respect to each example is described below and also in `examples.m`. 

### 1. "I have an imaging map (a nifti file) and a gene set. I want to test if the imaging pattern correlates to the expression pattern."

This example uses VBM meta-analysis result -- a brain map showing the vulnerability of brain volume -- for Alzheimer's disease (AD), and examines whether the VBM pattern is associated with the pattern of brain gene expression of three AD risk genes ('APOE', 'APP', 'PSEN2') through GAMBA null-models. Preprocessed gene expression data from the Allen Human Brain Atlas (https://human.brain-map.org) will be used.

NOTE: The example requires FSL (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/) installed. For users do not use FSL, please adjust the script to *coregister your brain map to /src/atlas/brain.nii.gz* and then run the rest of the codes. 


### 2. "I have an imaging data matrix (region by feature), a gene expression data matrix (region by gene), and a gene set. I want to test if the imaging pattern correlates to the expression pattern."

This example examines whether the gene expression pattern of 19 human-supragranular genes is associated with regional connectome metrics. The input data include: 1) a gene expression data matrix (57 regions x 5000 genes, a subset of the entire AHBA data); 2) an imaging data matrix (57 regions by 5 phenotypes, such as NOS-, FA-, SD-weighted nodal strength, nodal degree, and FC strength); 3) a gene set (19 genes).

NOTE: null-spin model ONLY works for DK114 atlas. If you use other atlas, please refer to Alexander-Bloch et al. (2018) to first generate gene expression data matrices for the 'spinned' atlases.

### 3. "I have a gene-set. I want to test in which brain regions the gene-set is differentially expressed."

In the third example, we show how to use this toolbox to examine in which brain regions the input GOI is differentially expressed considering different types of null models. The human-supragranular-enriched (HSE) genes are used here.

We load the gene set first:

`load('src/examples/example_conn_5k_genes.mat', 'geneset');`
 
Now we have a gene set of 19 HSE genes. We then use the null-coexpression model to examine whether the mean expression of HSE genes is significantly higher/lower than random genes with the similar coexpression conserved for each brain region.

`res_nullcoexp = permutation_expression_null_coexp(geneset);`
    
Next, we use the null-brain model to examine whether the mean expression of HSE genes is significantly higher/lower than random brain-enriched genes in general for each brain region.

`res_nullbrain = permutation_expression_null_brain(geneset);`

Using the above null models we will be able to examine the level of *gene specificity*. Another question we may ask is whether the mean expression of the GOI is specifically higher in one region in contrast to random brain regions. We are able to answer this “spatial specificity” question using the null-spin model.

`res_nullspin = permutation_expression_null_spin(geneset);`

Now, we have results for the questions we ask.


### 4. "I have a gene expression data matrix and a gene-set. I want to test in which brain regions the gene-set is over-expressed."

### 5. "I have an imaging map (.nii file) and I want to look for the most correlated genes"


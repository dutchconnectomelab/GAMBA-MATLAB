# GAMBA-MATLAB

A Matlab toolbox to study whether the expression of the gene(s) of interest (GOI) and neuroimaging-derived brain phenotypes show overlapped spatial patterns. Different statistical **null models** are available to examine both ***gene specificity*** and ***spatial specificity***. This toolbox is an extension of the web application [GAMBA](http://www.dutchconnectomelab.nl/GAMBA/). 

For details, please see:

> Wei Y. et al. (2021), Statistical testing and annotation of gene transcriptomic-neuroimaging associations, bioRxiv

## Installation

### Requirements

Before you start, make sure you have **Matlab** on your machine. 

No other software required, but in our examples we use [FSL - FMRIB Software Library](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki) to perform image coregistration. You can use other equivalent tools for the same purpose. See for details in `example.m`. 

### Download

Through the command line in terminal:

`git clone https://github.com/dutchconnectomelab/GAMBA-MATLAB.git`.

Alternatively, you can click 'Code -- Download ZIP' and unzip the downloaded file.

### Add to path

Open Matlab, click "Set Path" -- "Add with Subfolder..." and add the `GAMBA-MATLAB` folder to the search path.

### Extract required data

Run `install.m` in Matlab to retrieve all required data.

Alternatively, you can download files here: https://www.dropbox.com/sh/psfudnzktyd0860/AABtx7ESvEphO60dcV_xbQ4qa?dl=0. Unzip downloaded files to the folder `src/`.

## Examples

We use examples to show the utility of this toolbox. Examples cover the following usage of the toolbox:

- "I have an imaging map (e.g., a nifti file) and a gene set. I want to test if the imaging pattern correlates to the gene expression pattern."
- "I have an imaging data matrix (region by feature), a gene expression data matrix (region by gene), and a gene set. I want to test if the imaging pattern correlates to the gene expression pattern."
- "I have a gene-set. I want to test in which brain regions the gene-set is differentially expressed."
- "I have a gene expression data matrix and a gene-set. I want to test in which brain regions the gene-set is differentially expressed."
- "I have an imaging map (e.g., a nifti file) and I want to look for the most correlated genes."

Scripts concerning each question are included in `examples.m`. Here is a detailed tutorial:

#### 1. "I have an imaging map (a nifti file) and a gene set. I want to test if the imaging pattern correlates to the gene expression pattern."

This example uses VBM meta-analysis result -- a brain map showing the vulnerability of brain volume -- for Alzheimer's disease (AD) and examines whether the VBM pattern is associated with the pattern of brain gene expression of three AD risk genes ('APOE', 'APP', 'PSEN2') through GAMBA null-models. 

In the first step, you need to co-register the imaging file to *'/src/atlas/brain.nii.gz'* (which is the MNI152 atlas used for brain parcellation and segmentation). The following codes will do the co-registration using FSL FLIRT.

```matlab    
input_img_file = fullfile(filepath, 'src', 'examples', 'alzheimers_ALE.nii.gz');
input_img_anat_file = fullfile(filepath, 'src', 'examples', 'Colin27_T1_seg_MNI_2x2x2.nii.gz'); % anatomical file in the same space
ref_img_file = fullfile(filepath, 'src', 'atlas', 'brain.nii.gz'); % reference file in MNI152 space
reg_file = fullfile(filepath, 'output', 'registration.mat');
output_img_file = fullfile(filepath, 'output', 'coreg_alzheimers_ALE.nii.gz');

system(['flirt -in ', input_img_anat_file, ' -ref ', ref_img_file, ' -omat ', reg_file]);
system(['flirt -in ', input_img_file, ' -ref ', ref_img_file, ' -applyxfm -init ', reg_file, ' -out ', output_img_file]);
```

Note 1: if `flirt` doesn’t work you may need to use the absolute path to FSL (which can be viewed in terminal through `echo $FSLDIR`). Then the Matlab script should be something like `system([‘/usr/local/fsl/bin/flirt -in input_img_anat_file … … `
 
Note 2: If you do not use FSL, please adjust the codes to coregister the input imaging file to */src/atlas/brain.nii.gz* and then continue.

Next, you need to compute region-wise measures based on the registered imaging map. Here, voxels within each brain region in the Desikan-Killiany-114-region (DK114) atlas are averaged.

```matlab    
res_Y = group_regions(output_img_file, 'DK114');
```

After this you will get a 114 x 1 brain-phenotypic data array, together with region descriptions of the 114 regions. 

Now it’s time to compute the correlation between the phenotypic data array and gene expression data matrix. 

Preprocessed gene expression data matrix (57 x 17879) from the [Allen Human Brain Atlas](https://human.brain-map.org) will be used. Due to that only AHBA transcriptomic data from the left hemisphere is used, you also need to extract brain-phenotypic data for the left hemisphere.

```matlab 
img_data = res_Y.data(contains(res_Y.regionDescriptions, 'ctx-lh-'));
```

Then, you can examine for instance whether there is a **spatially specific** association between the brain-phenotypic pattern and the mean expression pattern of the input GOI (here 'APOE', 'APP', 'PSEN2') using the **null-spatial** model.

```matlab 
res_nullspatial = permutation_null_spin(img_data, {'APOE', 'APP', 'PSEN2'});
```

As we showed in our paper, **spatial specificity** is *however* not enough, you need to examine the level of **gene specificity**, namely, whether the observed association exceeds what one can expect by any other gene sets with similar co-expression levels (i.e., the **null-coexpressed-gene model**).

```matlab 
res_nullcoexpGene = permutation_null_coexp(img_data, {'APOE', 'APP', 'PSEN2'});
```

A more stringent null model can be used to examine how many genes of the GOI are **brain-expressed** genes and whether the observed association exceeds what one can expect by any brain-expressed genes (i.e., the **null-brain-gene model**).

```matlab 
res_nullbraingene = permutation_null_brain(img_data, {'APOE', 'APP', 'PSEN2'});
```

After doing the above analysis, you will know whether there is a significant association between your brain-phenotypic pattern and the gene set of your interest, and, more importantly, whether the association is **spatial specific** and **gene-specific**.




#### 2. "I have an imaging data matrix (region by feature), a gene expression data matrix (region by gene), and a gene set. I want to test if the imaging pattern correlates to the expression pattern."

This example examines whether the expression pattern of 19 human-supragranular-enriched (HSE) genes is associated with the pattern of regional connectome metrics. The input data for this example include: 1) a gene expression data matrix (57 regions x 5000 genes, a subset of the entire AHBA data); 2) an imaging data matrix (57 regions by 5 phenotypes, such as NOS-, FA-, SD-weighted nodal strength, nodal degree, and FC strength); 3) a gene set (19 genes).

You can load the example data first:

```matlab 
load('src/examples/example_conn_5k_genes.mat', 'geneset');
```

Then you can use the **null-spatial model** to examine whether there is a spatially specific association between the brain-phenotypic pattern and the mean expression pattern of the input GOI. **ATTENTION** the null-spatial model ONLY works for the DK114 atlas here. If you use other atlas, please refer to [Alexander-Bloch et al., 2018](https://github.com/spin-test/spin-test) to first generate random gene expression matrices using the spin-model.

```matlab 
res_nullspatial = permutation_null_spin(img_data, geneset);
```

Going beyond the null-spatial model, you need to examine the level of **gene specificity** using the **null-coexpressed-gene model**.

```matlab 
res_nullcoexpGene = permutation_null_coexp(img_data, geneset, gene_expression, gene_symbols);
```

A more stringent **null-brain-gene model** is recommended to examine how many genes of the GOI are **brain-expressed** genes and whether the observed association exceeds what one can expect by any brain-expressed genes

```matlab 
res_nullbraingene = permutation_null_brain(img_data, geneset, gene_expression, gene_symbols);
```

The above codes will help you get to know whether there is a significant association between your brain-phenotypic pattern and the gene set of your interest, and, more importantly, whether the association is **spatial specific** and **gene-specific**.

#### 3. "I have a gene-set. I want to test in which brain regions the gene-set is differentially expressed."

In the third example, we show how to use this toolbox to examine in which brain regions the input GOI is differentially expressed considering different types of null models. The human-supragranular-enriched (HSE) genes are used here.

You can load the example data first:

```matlab 
load('src/examples/example_conn_5k_genes.mat', 'geneset');
```

Now you get a gene set of 19 HSE genes. Then you can use the **null-coexpressed-gene model** to examine whether the mean expression of HSE genes is significantly higher/lower than random genes with the similar coexpression conserved for each brain region.

```matlab 
res_nullcoexpGene = permutation_expression_null_coexp(geneset);
```
    
You can also use the **null-brain-gene** model to examine whether the mean expression of HSE genes is significantly higher/lower than random brain-expressed genes for each brain region.

```matlab 
res_nullbraingene = permutation_expression_null_brain(geneset);
```

Using the above null models we will be able to examine the level of *gene specificity*.

Another question we may ask is *whether the mean expression of the GOI is specifically higher in one region in contrast to random brain regions*. You can get the answer using the **null-spatial** model.

```matlab 
res_nullspatial = permutation_expression_null_spin(geneset);
```

#### 4. "I have a gene expression data matrix and a gene-set. I want to test in which brain regions the gene-set is over-expressed."

In the fourth example, we show how to use this toolbox to examine in which brain regions the input GOI is differentially expressed considering different types of null models, given a customized gene expression matrix. The human-supragranular-enriched (HSE) genes are used as an example here.

You can load the example data first. Gene symbols of the GOI, the gene expression matrix of all genes, and symbols of all genes are needed.

```matlab 
load('src/examples/example_conn_5k_genes.mat', 'geneset', 'gene_expression', 'gene_symbols');
```

Then you can examine in which brain region the input GOI is differentially expressed compared to random genes, using the **null-coexpressed-gene model* where the null distribution of gene expression of random genes with similar coexpression level is generated and is used for permutation testing.

```matlab 
res_nullcoexpGene = permutation_expression_null_coexp(geneset, gene_expression, gene_symbols);
```

You can also examine the same question using the more stringent **null-brain-gene** model, where the null distribution of gene expression of brain-expressed genes is generated.

```matlab 
res_nullbraingene = permutation_expression_null_brain(geneset, gene_expression, gene_symbols);
```

For the other research question -- "in which brain region the input GOI is differentially expressed compared to random brain regions" -- you can use the following code only if you work on the Desikan-Killiany atlas with 57 left-hemisphere regions. Otherwise please use the spin model [Alexander-Bloch et al., 2018](https://github.com/spin-test/spin-test) or other equivalents to generate random gene expression matrices first.

```matlab 
res_nullspatial = permutation_expression_null_spin(geneset, gene_expression, gene_symbols);
```

#### 5. "I have an imaging map (.nii file) and I want to look for the most correlated genes"

In the last example, we show how to use our toolbox to look for the most correlated genes given an imaging map. This example uses VBM meta-analysis results for Alzheimer's disease (AD).

You need to first coregister the imaging file to MNI152 space.

```matlab
input_img_file = fullfile(filepath, 'src', 'examples', 'alzheimers_ALE.nii.gz');
input_img_anat_file = fullfile(filepath, 'src', 'examples', 'Colin27_T1_seg_MNI_2x2x2.nii.gz'); ref_img_file = fullfile(filepath, 'src', 'atlas', 'brain.nii.gz'); % reference file in MNI152 space
reg_file = fullfile(filepath, 'output', 'registration.mat');
output_img_file = fullfile(filepath, 'output', 'coreg_alzheimers_ALE.nii.gz');

system(['flirt -in ', input_img_anat_file, ' -ref ', ref_img_file, ' -omat ', reg_file]);
system(['flirt -in ', input_img_file, ' -ref ', ref_img_file, ' -applyxfm -init ', reg_file, ' -out ', output_img_file]);
```

Note 1: if `flirt` doesn’t work you may need to use the absolute path to FSL (which can be viewed in terminal through `echo $FSLDIR`). Then the Matlab script should be something like `system([‘/usr/local/fsl/bin/flirt -in input_img_anat_file … … `
 
Note 2: If you do not use FSL, please adjust the codes to coregister the input imaging file to */src/atlas/brain.nii.gz* and then continue.

Next, you need to compute region-wise measurements based on the registered imaging map. Here, voxels within each brain region in the DK-114 atlas are averaged.

```matlab 
res_Y = group_regions(output_img_file, 'DK114');
```

Using the resulted phenotypic data you can now do correlations between single-gene expression profiles and the phenotypic profile, and then do permutations per gene using the null-spatial model.

```matlab 
res_nullspatial = permutation_null_spin_correlated_genes(img_data);
```

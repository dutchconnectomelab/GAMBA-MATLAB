# GAMBA-MATLAB
A Matlab toolbox to perform statistical analyses included in GAMBA (www.dutchconnectomelab.nl/GAMBA), which can be used to test whether the gene expression profile(s) of the input gene(s) and neuroimaging-derived brain phenotypes show overlapped spatial patterns. Different statistical null models are availale to examine both gene specificity and spatial specificity.

For details, please see:

> Wei Y. et al. (2021), Statistical testing and annotation of gene transcriptomic-neuroimaging associations, bioRxiv

> Wei, Y. et al. (2019), Genetic mapping and evolutionary analysis of human-expanded cognitive networks. Nat Commun. https://doi.org/10.1038/s41467-019-12764-8

## Installation
#### Download matlab toolbox
- Download the github repository through the command line `git clone https://github.com/yongbin-wei/GAMBA_MATLAB.git`.

Alternative:

- Click 'Code -- Download ZIP' and unzip the downloaded file.

#### Extract required data
- Run `install.m` in Matlab to retrieve all required data.

Alternative:

- Manually download files under https://www.dropbox.com/sh/psfudnzktyd0860/AABtx7ESvEphO60dcV_xbQ4qa?dl=0 and unzip downloaded files under the folder `src/`.

## Examples
We use examples to show the utility of this toolbox. Condition with respect to each example is described below and also in `examples.m`. 

#### 1. 'I have an imaging map (a nifti file) and a gene set. I want to test if the imaging pattern correlates to the expression pattern.'

This example uses VBM meta-analysis result -- a brain map showing the vulnerability of brain volume -- for Alzheimer's disease (AD), and examines whether the VBM pattern is associated with the pattern of brain gene expression of three AD risk genes ('APOE', 'APP', 'PSEN2') through GAMBA null-models. Preprocessed gene expression data from the Allen Human Brain Atlas (https://human.brain-map.org) will be used.

NOTE: The example requires FSL (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/) installed. For users do not use FSL, please adjust the script to *coregister your brain map to /src/atlas/brain.nii.gz* and then run the rest of the codes. 


#### 2. "I have an imaging data matrix (region by feature), a gene expression data matrix (region by gene), and a gene set. I want to test if the imaging pattern correlates to the expression pattern."

#### 3. "I have a gene-set. I want to test in which brain regions the gene-set is over-expressed."

#### 4. "I have a gene expression data matrix and a gene-set. I want to test in which brain regions the gene-set is over-expressed."

#### 5. "I have an imaging map (.nii file) and I want to look for the most correlated genes"



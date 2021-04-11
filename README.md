# GAMBA-MATLAB
A Matlab toolbox to perform statistical analyses included in GAMBA (www.dutchconnectomelab.nl/GAMBA), which can be used to test whether the gene expression profile(s) of the input gene(s) and neuroimaging-derived brain phenotypes show overlapped spatial patterns. Different statistical null models are availale to examine both gene specificity and spatial specificity.

For details, please see:

> Wei Y. et al. (2021), Statistical testing and annotation of gene transcriptomic-neuroimaging associations, bioRxiv

> Wei, Y. et al. (2019), Genetic mapping and evolutionary analysis of human-expanded cognitive networks. Nat Commun. https://doi.org/10.1038/s41467-019-12764-8

## Install
### Download matlab toolbox
- Download the github repository through the command line `git clone https://github.com/yongbin-wei/GAMBA_MATLAB.git`.
OR
- Click 'Code -- Download ZIP' and unzip the downloaded file.

### Extract required data
- Run `install.m` in Matlab to retrieve all required data
OR
- Manually download data from https://www.dropbox.com/sh/psfudnzktyd0860/AABtx7ESvEphO60dcV_xbQ4qa?dl=0 and unzip downloaded files under the folder `src/`

## Examples
We use three simple examples that show analyses commonly performed in literature to illustrate the usage of different statistical null models. Examples include human-supragranular-enriched (HSE) genes, APOE gene, and risk genes of autism spectrum disorder (ASD). To get started, please see `/examples/README.txt` for details.

- Example of associations between the spatial pattern of HSE gene expression and the connectome metrics:

`/examples/scripts_example1_HSE.m`

- Example of associations between the spatial pattern of APOE gene expression and the pattern of brain atrophy in diseases:

`/examples/scripts_example2_apoe.m`

- Example of associations between the spatial expression pattern of ASD risk genes and the pattern of functional changes in diseases:

`/examples/scripts_example3_ASD.m`

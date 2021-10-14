function res = permutation_expression_null_brain(geneset, expressions, gene_symbols, background)
% PERMUTATION_EXPRESSION_NULL_BRAIN(GENESET, EXPRESSIONS, GENE_SYMBOLS) 
% performs permutation testing to examine in which regions the input gene 
% set is overexpressed.
%
% INPUT
%   geneset -- a cell array of gene symbols of the genes of interest
%   
%   OPTIONAL:
%   expressions -- a NxK matrix of gene expressions of all genes. N is the 
%       number of genes, K is the number of genes. if not available, 
%       default gene expression data in the DK114 atlas will be loaded.
%   gene_symbols -- a cell array of gene symbols of all genes. This must be
%       provided if EXPRESSIONS is provided. If not available, gene symbols
%       will be loaded from the default gene expression data.
%   background -- a string indicating the selection of brackground genes
%       Options: "brain" -- genes over-expressed in brain compared to other 
%                body sites (N=2655), default
%                "body" -- genes over-/similarly expressed in brain tissue
%                compared to other body sites (N=8296)
%                "general" -- genes expressed in brain without contrast to
%                other body sites (N=16778)
%
% OUTPUT
%   res.p -- two-tailed p-value in permutation testing
%   res.mean_expressions -- regional mean expressions of the input GOI
%   res.null_expressions -- regional mean expressions of random BRAIN genes
%   res.difference -- difference between mean_expressions and the mean of
%       null_expressions, indicating the effect direction.
%
% REFERENCE
%   Wei Y. et al., (2021) Statistical testing and annotation of gene 
%   transcriptomic-neuroimaging associations, bioRxiv


% ========================== Check input data =============================
disp('Runing null-brain model');
filepath = fileparts(mfilename('fullpath'));

if nargin == 1
    data_ge = load(fullfile(filepath, 'gene_expression.mat'));
    disp('## Loading default gene expression data in DK114 atlas ...');
    expressions = data_ge.mDataGEctx;
    gene_symbols = data_ge.gene_symbols;
    regionDescriptions = data_ge.regionDescriptionCtx;
    background = 'brain';
elseif nargin == 2
    error('Please provide gene symbols of all genes included in the expression data.');
elseif nargin == 3
    background = 'brain';
end

if isempty(expressions)
    warning('"expressions" is empty. Loading default expression matrix.');
    data_ge = load(fullfile(filepath, 'gene_expression.mat'));
    expressions = data_ge.mDataGEctx;
end

if isempty(gene_symbols)
    warning('"gene_symbols" is empty. Loading default expression matrix.');
    data_ge = load(fullfile(filepath, 'gene_expression.mat'));
    gene_symbols = data_ge.gene_symbols;
end

[N, K] = size(expressions);
disp(['## ', num2str(K), ' genes detected totally.']);
disp(['## ', num2str(N), ' brain regions detected.']);

NG = numel(geneset);
disp(['## ', num2str(NG), ' gene(s) of the GOI detected.']);

NGA = numel(gene_symbols);
if NGA~=K
    error('The number of gene symbols is different from the number of genes in the expression data');
end

II = ismember(gene_symbols, geneset);
if nnz(II) == 0
    error('None of the genes in the input gene set found in gene data.');
end
disp(['## ', num2str(nnz(II)), '/', num2str(NG), ' genes with gene expression data available.']);

% load brain genes
data_ge = load(fullfile(filepath, 'gene_expression.mat'), 'BRAINgene_idx',...
    'BRAINandBODYgene_idx', 'BRAIN_expressed_gene_idx', 'gene_symbols');

switch background
    case 'brain'
        BRAINgenes = data_ge.gene_symbols(data_ge.BRAINgene_idx);
        disp(['## Background is selected as genes over-expressed in brain, N=',num2str(numel(BRAINgenes))]);
    case 'body'
        BRAINgenes = data_ge.gene_symbols(data_ge.BRAINandBODYgene_idx);
        disp(['## Background is selected as genes over- or similarly expressed in brain compared to other body sites, N=',num2str(numel(BRAINgenes))]);
    case 'general'
        BRAINgenes = data_ge.gene_symbols(data_ge.BRAIN_expressed_gene_idx);
        disp(['## Background is selected as genes expressed in brain, without contrast to other body sites, N=',num2str(numel(BRAINgenes))]);
    otherwise
        warning('Background genes are not properly selected. Setting to "brain" by default.')
        BRAINgenes = data_ge.gene_symbols(data_ge.BRAINgene_idx);
        disp(['## Background is selected as genes over-expressed in brain, N=',num2str(numel(BRAINgenes))]);
end
geneset = intersect(geneset, BRAINgenes);
II = ismember(gene_symbols, geneset);
disp(['## ', num2str(nnz(II)), ' genes of the input GOI are brain-enriched genes.']);


% ========================= Perform permutation ===========================
nPerm = 10000;

% raw mean expressions
mGE = nanmean(expressions(:, ismember(gene_symbols, geneset)), 2);
res.mean_expressions = mGE;

% initialize permutation
fprintf('%s', '## Progress:     ');
idx_rand_genes = nan(nPerm, nnz(II));
[~, idx_background] = ismember(BRAINgenes, gene_symbols);
idx_background(idx_background==0) = [];
tmpGE = nan(N, nPerm);

% permutation
for kk = 1:nPerm
    fprintf('\b\b\b\b%.3d%%', round(kk/nPerm*100));
    rid = idx_background(randperm(numel(idx_background), nnz(II)));    
    idx_rand_genes(kk, :) = rid;
    tmpGE(:, kk) = nanmean(expressions(:, rid), 2);
end
res.null_expressions = tmpGE;
res.difference = res.mean_expressions - nanmean(res.null_expressions, 2);

% compute p-value
for ii = 1: N
    P = nnz(double(tmpGE(ii, :) > mGE(ii))) ./ nPerm;
    if P > 0.5
        res.p(ii, 1) = (1 - P) * 2;
    else
        res.p(ii, 1) = P * 2;
    end
end

if exist('regionDescriptions', 'var')
    res.regionDescriptions = regionDescriptions;
end

disp(' >> finished without errors');

end




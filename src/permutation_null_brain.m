function res = permutation_null_brain(img_data, geneset, expressions, gene_symbols, background)
% PERMUTATION_NULL_BRAIN(...) performs permutation testing to examine
% whether imaging profiles associate with gene expression profiles, based
% on the null-brain model (where random genes were selected from 
% brain-expressed genes).
%
% INPUT
%   img_data -- a NxM matrix of imaging profiles. N is the number of
%       brain regions, M is the number of imaging traits.
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
%   res.lr.beta -- standardized beta in linear regression
%   res.lr.p -- two-tailed p-value in linear regression
%   res.permut_beta -- standardized beta in the null model
%   res.permut_gene_idx -- indexes of random genes in each permutation
%
% REFERENCE
%   Wei Y. et al., (2021) Statistical testing and annotation of gene 
%   transcriptomic-neuroimaging associations, bioRxiv


% ========================== Check input data =============================
disp('Runing null-brain model');
filepath = fileparts(mfilename('fullpath'));

if nargin == 2
    data_ge = load(fullfile(filepath, 'gene_expression.mat'));
    disp('## Loading default gene expression data in DK114 atlas ...');
    expressions = data_ge.mDataGEctx;
    gene_symbols = data_ge.gene_symbols;
    background = 'brain';
elseif nargin == 3
    error('Please provide gene symbols of all genes included in the expression data.');
elseif nargin == 4
    background = 'brain';
elseif (nargin > 5) || (nargin < 2)
    error('Input error. Please check input data.')
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

[N, M] = size(img_data);
disp(['## ', num2str(N), ' brain regions detected.']);
disp(['## ', num2str(M), ' imaging trait(s) detected.']);

[NE, K] = size(expressions);
if N~=NE
    error('Different amount of regions in imaging data and gene data. Please check input data.');
end
disp(['## ', num2str(K), ' genes detected totally.']);

NG = numel(geneset);
disp(['## ', num2str(NG), ' gene(s) of the GOI detected.']);

NGA = numel(gene_symbols);
if NGA~=K
    error('The number of gene symbols is different from the number of genes in the expression data');
end

II = ismember(gene_symbols, geneset);
if nnz(II) == 0
    error('None of the genes in the input gene set found in gene data');
end
disp(['## ', num2str(nnz(II)), '/', num2str(NG), ' genes with gene expression data available']);

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
disp(['## ', num2str(nnz(II)), ' gene(s) of the input GOI found in background.']);


% ====================== Perform linear regression ========================
beta = nan(M, 1);
pval = nan(M, 1);

% Standardize
X = nanmean(expressions(:, II), 2);
X = (X - nanmean(X)) ./ nanstd(X);
Y = img_data;
Y = (Y - repmat(nanmean(Y, 1), N, 1)) ./ repmat(nanstd(Y, '', 1), N, 1);

for ii = 1: M
    stats = regstats(Y(:, ii),  X, 'linear', {'beta','tstat'});            
    beta(ii, 1) = stats.beta(2);
    pval(ii, 1) = stats.tstat.pval(2);
end
res.lr.beta = beta;
res.lr.p = pval;


% ========================= Perform permutation ===========================
nPerm = 10000;

idx_rand_genes = nan(nPerm, nnz(II));
beta_null = nan(nPerm, M);
fprintf('%s', '## Progress:     ');

[~, idx_background] = ismember(BRAINgenes, gene_symbols);
idx_background(idx_background==0) = [];

for kk = 1:nPerm
    fprintf('\b\b\b\b%.3d%%', round(kk/nPerm*100));

    rid = idx_background(randperm(numel(idx_background), nnz(II)));    
    idx_rand_genes(kk, :) = rid;

    % gene expressions of random genes
    X = nanmean(expressions(:, rid), 2);
    X = (X - nanmean(X)) ./ nanstd(X);

    for ii = 1: M
        stats = regstats(Y(:, ii), X, 'linear', {'beta'});            
        beta_null(kk, ii) = stats.beta(2);
    end
end
res.permut_gene_idx = idx_rand_genes;
res.permut_beta = beta_null;

% compute p-value
for ii = 1: M
    P = nnz(double(beta_null(:, ii) > beta(ii))) ./ nPerm;
    if P > 0.5
        res.p(ii, 1) = (1 - P) * 2;
    else
        res.p(ii, 1) = P * 2;
    end
end

disp(' >> finished without errors');

end

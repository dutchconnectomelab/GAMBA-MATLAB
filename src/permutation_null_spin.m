function res = permutation_null_spin(img_data, geneset, expressions, gene_symbols)
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
disp('Runing null-spin model');

if nargin == 2
    filepath = fileparts(mfilename('fullpath'));
    data_ge = load(fullfile(filepath, 'gene_expression.mat'));
    disp('## Loading default gene expression data in DK114 atlas ...');
    expressions = data_ge.mDataGEctx;
    gene_symbols = data_ge.gene_symbols;
elseif nargin == 3
    error('Please provide gene symbols of all genes included in the expression data.');
elseif (nargin > 4) || (nargin < 2)
    error('Input error. Please check input data.')
end

[N, M] = size(img_data);
disp(['## ', num2str(N), ' brain regions detected.']);
disp(['## ', num2str(M), ' imaging traits detected.']);

[NE, K] = size(expressions);
if N~=NE
    error('Different amount of regions in imaging data and gene data. Please check input data.');
end
disp(['## ', num2str(K), ' genes detected totally.']);

NG = numel(geneset);
disp(['## ', num2str(NG), ' gene(s) of the GOI detected.']);
if NG == 1
    error('Only 1 gene included in the GOI. Coexpression cannot be computed.')
end

NGA = numel(gene_symbols);
if NGA~=K
    error('The number of gene symbols is different from the number of genes in the expression data');
end

II = ismember(gene_symbols, geneset);
if nnz(II) == 0
    error('None of the genes in the input gene set found in gene data');
end
disp(['## ', num2str(nnz(II)), '/', num2str(NG), ' genes with gene expression data available']);

% load available gene symbols
filepath = fileparts(mfilename('fullpath'));
data_ge = load(fullfile(filepath, 'gene_expression.mat'), 'gene_symbols');

geneset = intersect(geneset, data_ge.gene_symbols);
NG = numel(geneset);
II = ismember(gene_symbols, geneset);
disp(['## ', num2str(NG), ' genes of the input GOI are available in Null-spin model.']);


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
% Load data for the null-spin model
null_spin_expression = nan(1000, 57, numel(geneset));
for ii = 1:NG
   filename = fullfile(filepath, 'gene_expression_spin', ...
       ['GE_spin_', geneset{ii}, '.txt']);
   null_spin_expression(:, :, ii) = dlmread(filename); 
end

% Average across genes
null_spin_exp_mean = nanmean(null_spin_expression, 3);

nPerm = 1000; % maximum 1000 here
fprintf('%s', '## Progress:     ');
beta_null = nan(nPerm, M);
for kk = 1:nPerm
    fprintf('\b\b\b\b%.3d%%', round(kk/nPerm*100));

    % randomized gene expressions
    X = null_spin_exp_mean(kk, :);
    X = (X - nanmean(X)) ./ nanstd(X);

    for ii = 1: M
        stats = regstats(Y(:, ii), X, 'linear', {'beta'});            
        beta_null(kk, ii) = stats.beta(2);
    end
end
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

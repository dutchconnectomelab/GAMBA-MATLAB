function res = permutation_expression_null_coexp(geneset, expressions, gene_symbols)
% PERMUTATION_EXPRESSION_NULL_COEXP(GENESET, EXPRESSIONS, GENE_SYMBOLS) 
% performs permutation testing to examine in which brain regions the input
% gene set is differentially expressed, in comparison to random genes with
% similar level of coexpression conserved.
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
%
% OUTPUT
%   res.p -- two-tailed p-value in permutation testing
%   res.mean_expressions -- regional mean expressions of the input GOI
%   res.null_expressions -- regional mean expressions of random BRAIN genes
%   res.difference -- difference between mean_expressions and the mean of
%       null_expressions, indicating the effect direction.
%   res.coexp_mean -- mean coexpression level among the input GOI
%   res.permut_coexp_mean -- mean coexpression level among random genes
%
% REFERENCE
%   Wei Y. et al., (2021) Statistical testing and annotation of gene 
%   transcriptomic-neuroimaging associations, bioRxiv


% ========================== Check input data =============================
disp('Runing null-coexpression model');

if nargin == 1
    filepath = fileparts(mfilename('fullpath'));
    data_ge = load(fullfile(filepath, 'gene_expression.mat'));
    disp('## Loading default gene expression data in DK114 atlas ...');
    expressions = data_ge.mDataGEctx;
    gene_symbols = data_ge.gene_symbols;
    regionDescriptions = data_ge.regionDescriptionCtx;
elseif nargin == 2
    error('Please provide gene symbols of all genes included in the expression data.');
end

[N, K] = size(expressions);
disp(['## ', num2str(K), ' genes detected totally.']);
disp(['## ', num2str(N), ' brain regions detected.']);

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


% ========================= Perform permutation ===========================
nPerm = 1000;

% raw mean expressions
mGE = nanmean(expressions(:, ismember(gene_symbols, geneset)), 2);
res.mean_expressions = mGE;

% compute coexpression of the input GOI
G = expressions(:, II);
coexp_mat = corr(G, 'rows', 'pairwise');
mask_tril = tril(ones(size(coexp_mat)), -1);
coexp = nanmean(coexp_mat(mask_tril == 1));
disp(['## Mean coexpression: ', num2str(coexp)]); 
res.coexp_mean = coexp;

coexp_null = nan(nPerm, 1);
idx_rand_genes = nan(nPerm, nnz(II));
tmpGE = nan(N, nPerm);

fprintf('%s', '## Progress:     ');
for kk = 1:nPerm
    tmp_status = false;
    while tmp_status ~= true
        [rid, coexp_null(kk), tmp_status] = y_rand_gs_coexp(...
            expressions, coexp, nnz(II));
    end
    fprintf('\b\b\b\b%.3d%%', round(kk/nPerm*100));
    idx_rand_genes(kk, :) = rid;
    
    % gene expressions of random genes
    tmpGE(:, kk) = nanmean(expressions(:, rid), 2);
end
res.permut_gene_idx = idx_rand_genes;
res.permut_coexp_mean = coexp_null;
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
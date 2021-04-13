function res = permutation_null_spin_correlated_genes(img_data)
% PERMUTATION_NULL_SPIN_CORRELATED_GENES(IMG_DATA) performs permutation testing to examine
% whether imaging profiles associate with gene expression profiles, based
% on the null-brain model (where random genes were selected from 
% brain-expressed genes).
%
% INPUT
%   img_data -- a NxM matrix of imaging profiles. N is the number of
%       brain regions, M is the number of imaging traits.
%
% OUTPUT
%   res.p -- two-tailed p-value in permutation testing
%   res.lr.beta -- standardized beta in linear regression
%   res.lr.p -- two-tailed p-value in linear regression
%   res.gene_symbols -- two-tailed p-value in permutation testing
%
% REFERENCE
%   Wei Y. et al., (2021) Statistical testing and annotation of gene 
%   transcriptomic-neuroimaging associations, bioRxiv

% ========================== Check input data =============================
disp('Looking for top correlated genes');

filepath = fileparts(mfilename('fullpath'));
data_ge = load(fullfile(filepath, 'gene_expression.mat'));
disp('## Loading default gene expression data in DK114 atlas ...');
expressions = data_ge.mDataGEctx;
gene_symbols = data_ge.gene_symbols;

[N, M] = size(img_data);
disp(['## ', num2str(N), ' brain regions detected.']);
disp(['## ', num2str(M), ' imaging traits detected.']);

[NE, K] = size(expressions);
if N~=NE
    error('Different amount of regions in imaging data and gene data. Please check input data.');
end
disp(['## ', num2str(K), ' genes totally.']);


% =========================== Loop all genes ==============================
beta = nan(K, M);
pval = nan(K, M);

% Standardize
Y = img_data;
Y = (Y - repmat(nanmean(Y, 1), N, 1)) ./ repmat(nanstd(Y, '', 1), N, 1);

fprintf('%s', '## Progress:     ');
for ii = 1:K
    fprintf('\b\b\b\b%.3d%%', round(ii/K*100));

    X = expressions(:, ii);
    X = (X - nanmean(X)) ./ nanstd(X);

    beta_null = nan(1000, M);
    for jj = 1:M
        % linear regression
        stats = regstats(Y(:, jj),  X, 'linear', {'beta', 'tstat'});            
        beta(ii, jj) = stats.beta(2);
        pval(ii, jj) = stats.tstat.pval(2);
        
        % null-spin
        filename = fullfile(filepath, 'gene_expression_spin', ...
           ['GE_spin_', gene_symbols{ii}, '.txt']);
        if exist(filename, 'file')
            null_spin_expression = dlmread(filename); 
        else
            null_spin_expression = nan(1000, N);
        end
        null_spin_expression = null_spin_expression';
        for kk = 1:1000
            X = null_spin_expression(:, kk);
            X = (X - nanmean(X)) ./ nanstd(X);
           
            % linear regression
            stats = regstats(Y(:, jj),  X, 'linear', {'beta','tstat'});            
            beta_null(kk, jj) = stats.beta(2);
        end
        
        % p value
        P = nnz(double(beta_null(:, jj) > beta(ii, jj))) ./ 1000;
        if P > 0.5
            res.p(ii, jj) = (1 - P) * 2;
        else
            res.p(ii, jj) = P * 2;
        end
    end
end
res.lr.beta = beta;
res.lr.p = pval;
res.gene_symbols = gene_symbols;

disp(' >> finished without errors');
end


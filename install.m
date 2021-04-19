% Run this script for the first time to retrieve necessary data 
clc, clear, close all
filepath = fileparts(mfilename('fullpath'));

if ~exist(fullfile(filepath, 'output'), 'dir')
    mkdir(fullfile(filepath, 'output'));
end

% Download atlas data
fprintf('%s', '## Download atlas data ...');
outfilename = websave(fullfile(filepath, 'src', 'atlas.zip'), ...
    'https://www.dropbox.com/s/a3ov1ztgaen5u2m/atlas.zip?dl=1');
unzip(outfilename, fullfile(filepath, 'src'));
delete(outfilename);
fprintf('%s\n', 'finished');

% Download expression data
fprintf('%s', '## Download expression data ...');
outfilename = websave(fullfile(filepath, 'src', 'gene_expression.mat.zip'), ...
    'https://www.dropbox.com/s/fujw4i8600xrc9k/gene_expression.mat.zip?dl=1');
unzip(outfilename, fullfile(filepath, 'src'));
delete(outfilename);
fprintf('%s\n', 'finished');

% Download examples
fprintf('%s', '## Download examples ...');
outfilename = websave(fullfile(filepath, 'src', 'examples.zip'), ...
    'https://www.dropbox.com/s/3iteb4btcbvi39p/examples.zip?dl=1');
unzip(outfilename, fullfile(filepath, 'src'));
delete(outfilename);
fprintf('%s\n', 'finished');

% Download data for the null-spatial model
fprintf('%s', '## Download expression data for the null-spatial model (~1.5G) ...');
outfilename = websave(fullfile(filepath, 'src', 'gene_expression_spin.zip'), ...
    'https://www.dropbox.com/s/nwmtqro3dkma3u1/gene_expression_spin.zip?dl=1');
unzip(outfilename, fullfile(filepath, 'src'));
delete(outfilename);
fprintf('%s\n', 'finished');

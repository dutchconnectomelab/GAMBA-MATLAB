function res = group_regions(input_file, atlas)
% GROUP_REGIONS(...) computes region-wise mean values from a voxel-wise
% brain map.  
%
% INPUT
%   input_file -- input brain map (.nii file) that has been co-registered 
%       to the same space as './data/brain.nii.gz' (i.e., MNI152 brain)
%
%   OPTIONAL:
%   atlas -- atlas name. Default: 'DK114'. Other options: 'aparc','DK250'
%
% OUTPUT
%   res.data -- N x 1 array of regional mean value extracted from the
%       input imaging data. N is the number of regions.
%   res.regionDescriptions -- N x 1 cell str of region descriptions.
%   res.regionIndexes -- N x 1 array of region indexes.
%
% REFERENCE
%   Wei Y. et al., (2021) Statistical testing and annotation of gene 
%   transcriptomic-neuroimaging associations, bioRxiv
% 
%   NOTE: This function uses FreeSurfer's MRIread(...) to read nifti file.
%   Fischl, B. 2012. “FreeSurfer.” NeuroImage 62 (2): 774–81.


if nargin == 1
    atlas = 'DK114';
end

filepath = fileparts(mfilename('fullpath'));

switch atlas
    case 'DK114'
        lookuptable = fullfile(filepath, 'src', 'atlas', 'lausanne120.txt');
        ref_file = fullfile(filepath, 'src', 'atlas', 'lausanne120+aseg.nii.gz');
    case 'aparc'
        lookuptable = fullfile(filepath, 'src', 'atlas', 'aparc.txt');
        ref_file = fullfile(filepath, 'src', 'atlas', 'aparc+aseg.nii.gz');
    case 'DK219'
        lookuptable = fullfile(filepath, 'src', 'atlas', 'lausanne250.txt');
        ref_file = fullfile(filepath, 'src', 'atlas', 'lausanne250+aseg.nii.gz');
    otherwise
        error('Atlas not supported.');
end

% load data
hdr = load_nifti(input_file);
vol = hdr.vol;

ref = load_nifti(ref_file);

% read color table
tbl = readtable(lookuptable);
res.regionDescriptions = tbl.Var2;
res.regionIndexes = tbl.Var1;

% compute regional mean
res.data = nan(numel(res.regionDescriptions), 1);
for ii = 1:numel(res.regionIndexes)
    tmp = vol(ref.vol == res.regionIndexes(ii));
    res.data(ii, 1) = mean(tmp);
end

end
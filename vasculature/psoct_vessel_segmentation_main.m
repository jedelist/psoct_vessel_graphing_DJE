%% Main file for calling segmentation functions
% Author: Mack Hyman
% Date Created: March 16, 2023
%
% Detailed Description
%{
This script performs the following:
- segment the original volume
- apply a mask to the segmentation
- convert segmentation to graph
- Remove loops from graph

To Do:
- find optimal range for remove_mask_islands
- prune spurs from segments in graph
%}
clear; clc; close all;

%% Add top-level directory of code repository to path
% This allows Matlab to find the functions in the project folders

% Start in current directory
mydir  = pwd;
% Find indices of slashes separating directories
if ispc
    idcs = strfind(mydir,'\');
elseif isunix
    idcs = strfind(mydir,'/');
end
% Truncate path to reach top-level directory (psoct_vessel_graphing)
topdir = mydir(1:idcs(end));
addpath(genpath(topdir));

%% Initialize data path for linux or personal machine (debugging)

%%% Local machine
if ispc
    dpath = 'C:\Users\mack\Documents\BU\Boas_Lab\psoct_data_and_figures\test_data\Ann_Mckee_samples_10T\';
    % Subject IDs
    subid = 'NC_6047';
    subdir = '\dist_corrected\volume\';
    % Filename to parse (this is test data)
    fname = 'ref_4ds_norm_inv_crop';
    % filename extension
    ext = '.tif';
    % sigma for Gaussian smoothing
    gsigma = [2,3,4];

%%% Computing cluster (SCC)
elseif isunix
    % Path to top-level directory
    dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/';
    % Subfolder containing data
    subdir = '/dist_corrected/volume/';
    % Filename to parse (this will be the same for each subject)
    fname = 'ref_4ds_norm_inv';
    % filename extension
    ext = '.tif';
    
    %%% Complete subject ID list for Ann_Mckee_samples_10T
    subid = {'AD_10382', 'AD_20832', 'AD_20969',...
             'AD_21354', 'AD_21424',...
             'CTE_6489','CTE_6912',...
             'CTE_7019','CTE_8572','CTE_7126',...
             'NC_6047', 'NC_6839',...
             'NC_6974', 'NC_7597',...
             'NC_8095', 'NC_8653',...
             'NC_21499','NC_301181'};
    
    %%% Gaussian sigma arrays:
    % Small vessel sigma array = [1, 3, 5]
    % Medium vessel sigma array = [5, 7, 9]
    % Large vessel sigma array = [7, 9, 11]
    sigmas = [3,5,7; 5,7,9; 7,9,11];
    sigmas = [2,3,4];    
    
    %%% Create cell array of subject ID and sigma for job array on the SCC 
    nrow = length(subid)*size(sigmas,1);
    nsigma = size(sigmas,1);
    sub_sigma = cell(length(subid).*size(sigmas,1), 2);
    idx = 1;
    if size(sigmas,1) > 1
        % Fill sub_sigma cell array with each sigma array for each subject
        for i = 1:size(sigmas,1):nrow
            sub_sigma{i,1} = subid{idx};
            sub_sigma{(i+1),1} = subid{idx};
            sub_sigma{(i+2),1} = subid{idx};
            idx = idx + 1;
            for j = 1:nsigma
                sub_sigma{(i+j-1),2} = sigmas(j,:);
            end
        end
    else
        for i = 1:nrow
            sub_sigma{i,1} = subid{i};
            sub_sigma{i,2} = sigmas;
        end
    end
    %%% Reassign subid and sigma based on job array counter
    % Retrieve SGE_TASK_ID from system (job array index)
    batch_idx = getenv('SGE_TASK_ID');
    
    % If this is a job array, then batch_idx will not be empty.
    if ~strcmp(batch_idx,'undefined')
        % Convert from ASCII to double
        batch_idx = str2double(batch_idx);
        % Retrieve corresponding row from sub_sigma
        [subid, gsigma] = sub_sigma{batch_idx, :};
    % Otherwise, set the Gaussian sigma manually
    elseif strcmp(batch_idx,'undefined')
        subid = 'NC_7597';
        gsigma = [5,7,9];
    end
end

%% Initialization parameters (same for both 

%%% Assign PS-OCT voxel dimension [x, y, z] according to downsample factor
% Downasample factor = 4 --> Voxel = [12, 12, 15] micron
% Downasample factor = 10 --> Voxel = [30, 30, 35] micron
% 2P microscopy pixel will always be [2, 2] micron
if regexp(fname, '4ds')
    vox_dim = [12, 12, 15];
elseif regexp(fname, '10ds')
    vox_dim = [30, 30, 35];
else
    vox_dim = [30, 30, 35];
end

%%% Size of the Gaussian kernel. This should be a 3-element array of
% positive, odd integers. Default size is 2*ceil(2*gsigma)+1
gsize = 2.*ceil(2.*gsigma)+1;

%%% Minimum fringi filter probability to classify voxel as vessel
min_prob = 0.18:0.04:0.26;
min_prob = 0.23;

%%% A segment with < "min_conn" voxels will be removed
min_conn = 5;

%%% Boolean for converting segment to graph (0 = do not convert. 1 = convert)
graph_boolean = 0;

%%% Boolean for visualizing the graph debugging plots
viz = false;

%% Load raw volume (TIF) and convert to MAT
% Define entire filepath 
fullpath = fullfile(dpath, subid, subdir);
filename = strcat(fullpath, strcat(fname, ext));
% Convert .tif to .MAT
vol_uint16 = TIFF2MAT(filename);

%%% Create subfolder for Gaussian sigma and kernel size
% Create string of Gaussian sigmas
gsigma_str = num2str(gsigma);
% Replace spaces with hyphens
gsigma_str = strrep(gsigma_str, '  ', '-');
gsigma_subfolder = strcat('gsigma_',gsigma_str);

% Create string of Gaussian kernel sizes
gsize_str = num2str(gsize);
% Replace spaces with hyphens
gsize_str = strrep(gsize_str, '  ', '-');
gsize_subfolder = strcat('_gsize_',gsize_str);

% concatenate sigma and kernel into single directory
subfolder = strcat(gsigma_subfolder, gsize_subfolder);

% Create string for entire directory path to subfolder
fullpath = fullfile(fullpath, subfolder);

% Create subfolder with Gaussian sigma and kernel size
if ~exist(fullpath, 'dir')
   mkdir(fullpath)
   % Add metadata text file
   metadata = {strcat('gaussian sigma =  ', num2str(gsigma)),...
       strcat('gaussian kernel =  ', num2str(gsize)),...
       strcat('minimum probability =  ', num2str(min_prob)),...
       };
   writelines(metadata, fullfile(fullpath, 'metadata.txt'));
end


%% Segment volume
% convert volume to double matrix
vol = double(vol_uint16);
% Segment volume. Threshold with first element of probability matrix.
[pmat, seg] = vesSegment(vol, gsigma, gsize, min_prob(1), min_conn);
% Save probability map for posterity
fout = strcat(fullfile(fullpath, 'probability_map'), '.mat');
save(fout, 'pmat', '-v7.3');

for j = 1:length(min_prob)
    %%% Threshold probability matrix with min_prob array
    I_seg = pmat;
    I_seg(pmat < min_prob(j)) = 0;
    I_seg(pmat >= min_prob(j)) = 1;
    % Convert binary matrix to unsigned 8-bit to save memory
    I_seg = uint8(I_seg);
    % Remove segments with fewer than voxmin connected voxels
    if min_conn > 0
        I_seg = rm_short_vessels(I_seg, min_conn);
    end

    %%% Save unmasked & thresholded segmentation to TIF
    % Create filename for probability
    fname_seg = strcat(fname,'_segment_pmin_',num2str(min_prob(j)));
    fout = strcat(fullfile(fullpath, fname_seg), '.tif');
    segmat2tif(I_seg, fout);

    %%% Overlay volume (grayscale) and unmasked segmentation (green)
    % Create output filename
    overlay_name = strcat(fname_seg, '_overlay.tif');
    overlay_fout = fullfile(fullpath, overlay_name);
    % Call function to overlay mask and segmentation
    overlay_vol_seg(vol_uint16, I_seg, 'green', overlay_fout);
    
    %%% Create a graph of the segmentation
    if graph_boolean
        seg_graph_init(I_seg, vox_dim, fullpath, fname_seg, viz);
    end    
end
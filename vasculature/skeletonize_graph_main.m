%% Main file for skeletonizing/graphing the vasculature
%
%{
This script performs the following:
- Calculate volume of masked PSOCT volume
- skeletonize/graph the segmentation
- calculate vascular metrics from graph
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
% This section is for setting the directory path to your datasets.
% The code is written to create a filepath with the following structure:
% [dpath]\[subid]\[subdir]\[fname].[ext]

%James' notes: All dir info commented out, will hardcode my own dir paths

% Your local directory to the segmentation data
dpath = '/projectnb/npbssmic/ns/Hui_Wang_samples/'
% Subject IDs
subid = 'caa_17_DJE';    %caa_17
subdir = '/caa17_fullfiles/';     %
% Filename to parse (this is test data)
fname = 'I_mosaic_1_0_0';      %I_mosaic_1_0_0 only b/c this is OCT volume
% filename extension
ext = '.mgz';
% sigma for Gaussian smoothing
gsigma = [7, 9, 11];
% Voxel dimensions (verify this with Dylan or another RA at Martinos)
vox_dim = [10, 10, 10];
% Volume of a single voxel (microns cubed)
vox_vol = vox_dim(1) .* vox_dim(2) .* vox_dim(3);

%% Load the masked PSOCT volume

% CHANGE THIS: This is the filepath to the original PSOCT volume with the
% background (non-tissue voxels) removed (set to white).
fullpath = fullfile(dpath, subid, subdir);          %FIND WHICH ONE TO HARDCODE AND COMMENT OUT
filename = strcat(fullpath, strcat(fname, ext));  %COMMENT THIS OUT AND HARDCODE PATH

% Convert .tif to .MAT
%tissue = TIFF2MAT(filename);

% Read nii and convert to .MAT
nii = MRIread(filename);
tissue = nii.vol;

% Invert so that non-tissue voxels are zeros
tissue_inv = imcomplement(tissue);

% Calculate total number of non-zero voxels in tissue sample
tissue_logical = logical(tissue_inv);
voxels = sum(tissue_logical(:));

% Convert voxels to metric volume (cubic microns)
vol = voxels .* vox_vol;

%% Load segmentation volume and convert to MAT
% We usually use .TIF for the segmentation volume. I am unsure if you use a
% different file format though. If so, you may need to use a different
% function to load the segmentation volume. This section uses "TIFF2MAT" to
% load a TIF file into a Matlab matrix

% Filename to parse (this is test data)
fname = 'seg_mosaic_1_0_0';      %fin.seg_mosaic_1_0_0_DJE OR seg_mosaic_1_0_0
% filename extension
ext = '.mgz';

% Define entire filepath 
fullpath = fullfile(dpath, subid, subdir);
filename = strcat(fullpath, strcat(fname, ext));
% CHECK THIS: Convert .tif to .MAT
%segment = TIFF2MAT(filename);  % comment out for mgz file

nii = MRIread(filename);
segment = nii.vol;

%% Skeletonize and convert to Graph
% This section calls the function at the bottom of the script. This
% function then calls another function "seg_to_graph" which actually
% performs the skeletonization and converts the skeleton to the three
% dimensional graph (nodes + edges). I wrote an example below, where the
% output "segment_graph_data" will be a data structure containing the graph
% data.

segment_graph_data = seg_graph_init(segment, vox_dim, fullpath, filename);

%% Calculate the vessel metrics from the graph data
% Metrics: length of each vessel, length density of volume , tortuosity
% NOTE: This section may require modification, depending on how you save
% your data in the previous section. 

met = struct();

% Load graph and segmentation (angio)
graph = segment_graph_data.Graph; % might be segment_graph_data.Data.Graph
seg = segment_graph_data.angio;

% Load end nodes from graph
endnodes = graph.segInfo.segEndNodes;

% Load position of all nodes from graph
nodepos = graph.segInfo.segPos;

%%% Calculate total length (microns) from graph
len = graph.segInfo.segLen_um;
len_tot = sum(len(:));
met.total_length = len_tot;

%%% Calculate mean length (microns)
len_avg = mean(len);
met.avg_length = len_avg;

%%% Calculate length density (length of vess/unit volume)
len_density = len_tot ./ vol;
met.length_density = len_density;

%%% Total number of vessels
nves = length(len);
met.total_vessels = nves;

%%% tortuosity arc-chord ratio (curve length / euclidean)
% Initalize matrix for storing tortuosity
tort = zeros(nves, 1);
for j=1:nves
    % convert segment end nodes to cartesian coordinate
    node1 = graph.nodes(endnodes(j,1), :);
    node2 = graph.nodes(endnodes(j,2), :);
    % Calcualte euclidean distance of segment
    euc = sqrt((node1(1) - node2(1)).^2 +...
                (node1(2) - node2(2)).^2 +...
                (node1(3) - node2(3)).^2);
    % Calculate tortuosity (single arc-chord ratio)
    tort(j) = len(j) ./ euc;
end
% Add to metrics structures
met.tortuosity = mean(tort);   

% Save the output
save(ad_cte_fout, 'met', '-v7.3');

%% Function to skeletonize and convert to Graph
function [Data] = seg_graph_init(seg, vox_dim, fullpath, fname_seg) % PRBLM fname_seg not initialized
% Initialize graph from segmentation
% INPUTS:
%   seg (mat): segmentation matrix
%   vox_dim (array): 3-element array of voxel dimensions (microns)
% OUTPUTS:
%   Data (struct): graph data for the segmentation volume

%%% Convert segmentation to graph (just the nodes and segments)
graph_nodes_segs = seg_to_graph(seg, vox_dim);

%%% Initialize graph metadata (Graph.Data)
[Data] = init_graph(graph_nodes_segs);

%%% Append "angio" data (segmentation matrix)
% Rearrange [x,y,z] --> [z,x,y]. This is the standard in
% the graph validation GUI.
angio = permute(seg, [3,1,2]);
Data.angio = angio;

% Create new filename for graph and add .MAT extension
fname_graph = strcat(fname_seg, '_graph_data.mat');
fout = fullfile(fullpath, fname_graph);
save(fout,'Data', '-v7.3');
end

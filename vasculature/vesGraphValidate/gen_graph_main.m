%% Main script for creating graph from segmentation
% This script will call the function for generating a graph. This will also
% remove loops from the graph.
clear; clc; close all;
%% Add top-level directory of code repository to path
% Print current working directory
mydir  = pwd;
% Find indices of slashes separating directories
if ispc
    idcs = strfind(mydir,'\');
elseif isunix
    idcs = strfind(mydir,'/');
end
% Remove the two sub folders to reach parent
% (psoct_human_brain\vasculature\vesSegment)
topdir = mydir(1:idcs(end-1));
addpath(genpath(topdir));

%% Lookup table of pia pixel intensity reference values
subid = {'AD_10382', 'AD_20832', 'AD_20969',...
             'AD_21354', 'AD_21424',...
             'CTE_6489','CTE_6912',...
             'CTE_7019','CTE_8572','CTE_7126',...
             'NC_6047', 'NC_6839',...
             'NC_6974', 'NC_7597',...
             'NC_8095', 'NC_8653',...
             'NC_21499','NC_301181'};

% Subset of subjects with good segmentation
subid =     {'AD_10382', 'AD_20832', 'AD_20969',...
             'AD_21354', 'AD_21424',...
             'CTE_6489','CTE_6912',...
             'CTE_7019','CTE_7126',...
             'NC_6839','NC_6974',...
             'NC_7597','NC_8653',...
             'NC_21499','NC_301181'};

%% Import files
if ispc
    % Laptop directory structure
    dpath = ['C:\Users\mack\Documents\BU\Boas_Lab\psoct_data_and_figures\' ...
            'test_data\Ann_Mckee_samples_10T\'];
    % Subfolder containing data
    subdir1 = '\dist_corrected\volume\';
elseif isunix
    % Upper level directory
    dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_55T';    
    % Subfolder containing data
    subdir1 = '/dist_corrected/volume/';
end

% Filename with PSOCT scattering tissue volume
vol_name = 'ref_4ds_norm_inv'; 
% Combined segmentations subfolder
subdir2 = 'combined_segs';
% Combined segmentation subfolder
segdir = 'gsigma_3-5-7_5-7-9_7-9-11';
% Filename with combined segmentations
seg_name = 'seg_masked';

%%% Assign subid based on job array counter
if isunix
    % Retrieve SGE_TASK_ID from system (job array index)
    batch_idx = getenv('SGE_TASK_ID');
    % Convert from ASCII to double
    batch_idx = str2double(batch_idx);
    % Assign local subject ID
    sub = subid{batch_idx};
    % Disable visualization of plots
    viz = false;
    % Set the remove loop boolean to zero
    rmloop_bool = 0;
else
    sub = subid{1};
    % Enable visualization of plots
    viz = false;
    % Set the remove loop boolean to zero
    rmloop_bool = 0;
end

%% Iterate through subjects. Generate graph

%%% Import combined segmentation file
fpath = fullfile(dpath, sub, subdir1, subdir2, segdir, strcat(seg_name,'.mat'));
tmp = load(fpath,'segm');
seg = tmp.segm;

%%% Pre-process combined segmentation prior to graphing
% Remove small segments (fewer than voxmin)
voxmin = 30;
seg = rm_short_vessels(seg, voxmin);
% Remove large islands (likley pia boundary false positives)
% Some of the segmentations contain large blobs in the pia. I need to
% write a function to identify + remove these.

%%% Convert to graph
% Set the voxel size
if regexp(vol_name, '4ds')
    vox_dim = [12, 12, 15];
elseif regexp(vol_name, '10ds')
    vox_dim = [30, 30, 35];
else
    vox_dim = [30, 30, 35];
end
%%% Initialize graph + remove loops + save output
graph_path = fullfile(dpath, sub, subdir1, subdir2, segdir);
seg_graph_init(seg, vox_dim, graph_path, seg_name, viz, rmloop_bool);


% Radial Density Eccentricity
% AOIP
% Created by: Jenna Grieshop
% Date created: 4/21/2025
%
% 
%
%

clear all;
close all;
clc;

scaleFactor = 0.25; %umpp

basepath = which('Radial_Density_Eccentricity.m');
[basepath] = fileparts(basepath);
path(path,fullfile(basepath,'lib')); % Add our support library to the path.

% select matrix
[filename, pathname] = uigetfile('*.csv', 'MultiSelect', 'off');

% load in the data
data = readmatrix(fullfile(pathname,filename));

% select and load in filename of the LUT with CDC
[LUTfilename, LUTpathname] = uigetfile('*.csv', 'Select file with CDC coords', pathname);

% Load in the LUT file
[~, lutData] = load_LUT_file(fullfile(LUTpathname,LUTfilename));

% Find the information from the LUT file for the data                                
LUTindex=find( cellfun(@(s) contains(filename,s ), lutData{1} ) );

CDC_x = lutData{2}(LUTindex);
CDC_y = lutData{3}(LUTindex);

if CDC_x ~= CDC_y
    warning('CDC is not centered in matrix');
end

dimension = size(data, 1);

px_num_25um = (1/scaleFactor)*25;

[Tics,Average] = radial_profile(data, px_num_25um);

% combine ring averages
first = 1;
for i=1:6
    if first ==1
        Avg(i) = Average(i);
        first = 0;
    else
        Avg(i) = (Average(i) + Avg(i-1))/2;
    end
end











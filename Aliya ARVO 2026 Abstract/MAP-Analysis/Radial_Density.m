% Radial Density
% AOIP
% Created by: Jenna Grieshop
% Date created: 4/22/2025
% Last edited by Aliya Siddiqui 11/24/25
%
% Description: Script that takes density matrix (i.e. averaged bound
% density map) and calculates average density radially from the CDC.
%
% Input: The density matrix, CDC LUT file (file name, CDC coordinates)
%
% Output: csv with average radial densities going out from the center.
%


clear all
close all
clc

scaleFactor = 0.25; %umpp
spacing = 1; %5 microns

basepath = which('Radial_Density.m');
[basepath] = fileparts(basepath);
path(path,fullfile(basepath,'lib')); % Add our support library to the path.

% Select matrix
[filename, pathname] = uigetfile('*.csv', 'MultiSelect', 'off');

% Load in the matrix
data = readmatrix(fullfile(pathname,filename));

% Select and load in filename of the LUT with CDC
[LUTfilename, LUTpathname] = uigetfile('*.csv', 'Select file with CDC coords', pathname);

% Load in the LUT file
[~, lutData] = load_LUT_file(fullfile(LUTpathname,LUTfilename));

% Find the information from the LUT file for the data                                
LUTindex=find( cellfun(@(s) contains(filename,s ), lutData{1} ) );

% Extract cdc from the lut file
CDC_x = lutData{2}(LUTindex);
CDC_y = lutData{3}(LUTindex);
CDC_density = data(CDC_y, CDC_x);  % row = y, col = x
CDC = [CDC_y, CDC_x];

% Warn user if the matrix is not square
if CDC_x ~= CDC_y
    warning('CDC is not centered in matrix');
end

% Figure out how many pixels are in the selected spacing (microns)
px_num = (1/scaleFactor)*spacing;
dimension = size(data, 1);

% Determine how many eccentricities to take a density reading at
points = floor((dimension/px_num)/2);

% go through getting 
for i=1:points

    % Create Meshgrid of the matrix
    [X,Y] = meshgrid(1:size(data,1),1:size(data,2));
    % Distance formula to get all data in meshgrid within radius of
    % eccentricity (5 microns * iteration)
    disk_locations = sqrt((X-CDC(1)).^2+(Y-CDC(2)).^2) <= (px_num*i);
    
    % Get the outer boundary of the circle (1 pixel)
    outerboundary = imdilate(disk_locations,strel('disk',1))&~disk_locations;
    outerboundary = outerboundary';
    imshow(outerboundary)
    
    % Get values from the matrix in the locations of the circle boundary
    density_ring = data(outerboundary==1);
    % Get average of the density ring
    Avg(i) = mean(density_ring);
    eccentricity(i)=i


end

plot(Avg);

% Format data to be saved
Avg = Avg';
header = {"Eccentricity", filename};
densityveccentricity = [num2cell(eccentricity)', num2cell(Avg)]
header2 = {0, CDC_density}
final_results = [header; header2; densityveccentricity];
result_fname =  fullfile(LUTpathname, [filename(1:end-4), '_' num2str(spacing) '_um_spacing_radial_density.csv']);

% Save csv
writecell(final_results, result_fname);





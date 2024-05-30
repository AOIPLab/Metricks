% Cells within Radius from PCD
%
% Created by: Jenna Grieshop
% Date created:  11/22/23
%
% Description: Finds the nujmber of cells within a circle with a user given radius from
% the PCD in a density matrix. User enters the radius in um (microns) from 
% the PCD. 
%
% NOTE: Currently not able to handle both eyes from the same subject.
%
% Input: Denisty matrices, cone coordinate files, a LUT file ("LUT" in the
% name of the file) (subject ID, axial length, and ppd) all in the same 
% data folder.
%
% Output: .csv saved to the data folder.


clear all
close all
clc

% Add our support library to the path.
basePath = which('cells_within_radius_from_PCD.m');
[basePath] = fileparts(basePath);
path(path,fullfile(basePath,'lib')); % Add our support library to the path.

% User selects desired folder:
dataPath = uigetdir('.','Select directory containing analyses');
dataDir = dir(dataPath);
dataDir = struct2cell(dataDir)';
addpath(genpath(dataPath));

% Make output folder
impath = fullfile(dataPath,'Matlab_Outputs');
mkdir(impath);

% Looks for the Batch Info spreadsheet that contains the ID initial, ID #, eye, axial length, and device #:
batchDir = dataDir(...
    ~cellfun(@isempty, strfind(dataDir(:,1), 'LUT')),:);
% Loads in the Batch Info file:
batchInfo = table2cell(readtable(fullfile(batchDir{1,2},batchDir{1,1})));

% Looks for all bound density matrices
boundDensityDir = dataDir(...
    ~cellfun(@isempty, strfind(dataDir(:,1), 'bounddensity_matrix')),:);

% Looks for all bound density matrices
coordsDir = dataDir(...
    ~cellfun(@isempty, strfind(dataDir(:,1), 'coords')),:);


% Get user input radius in microns
input = NaN;
while isnan(input)

    input = inputdlg('Input the radius from the PCD in um:', 'Input the radius from the PCD in um:');
    
    input = str2double(input);

    if isempty(input)
            error('radius input cancelled by user.');
    end
end


count = 1;

for i = 1:size(batchInfo, 1)
    
    %% Pulling information from batch info:
    subjectId = num2str(batchInfo{i,1});
    axialLength = batchInfo{i,2};
    ppd = batchInfo{i,3};
    
    subjMatrix = boundDensityDir(...
    ~cellfun(@isempty, strfind(boundDensityDir(:,1), subjectId)),:);

    matrixName = char(subjMatrix(...
    ~cellfun(@isempty, strfind(subjMatrix(:,1), subjectId)), 1));

    subjCoords = coordsDir(...
    ~cellfun(@isempty, strfind(coordsDir(:,1), subjectId)),:);

    coordName = char(subjCoords(...
    ~cellfun(@isempty, strfind(subjCoords(:,1), subjectId)), 1)); 

    %% Load in data
    densityMap = csvread(fullfile(boundDensityDir{1,2}, matrixName));
    coords = dlmread(coordName);
    x = coords(:,1);
    y = coords(:,2);

    %% 
    % Find PCD
    [maxVal, maxInd] = max(densityMap(:));

    % Finding the weighted (mean) centoid if multiple max locations found
    [maxYCoords, maxXCoords] = find(densityMap == maxVal);
    centroidX = mean(maxXCoords);
    centroidY = mean(maxYCoords);

    %% Calculate scale and find radius in pixels
    micronsPerDegree = (291*axialLength)/24; % These numbers are from a book in Joe's office, Joe has confirmed it is correct
    scaleVal = ppd / micronsPerDegree; % Pixel/um
    radiusPx = round(input * scaleVal); % Radius in pixels

    %% Find cells in radius

    % Creating a circle mask of 0s on the density map 
    [xMesh,yMesh] = meshgrid(1:size(densityMap,1));
    circleMask = densityMap;
    circleMask((xMesh - centroidX).^2 + (yMesh - centroidY).^2 < radiusPx^2) = 0;
    circleMask2 = (circleMask>0);

    contour = edge(circleMask2); % Binary contour from thresholded density map
    [newIndeces] = find_close_indeces(contour); % Using function from matlab file exchange
    in = inpolygon(x,y,newIndeces(:,2),newIndeces(:,1)); % Finds which cell coordinates are inside/on the contour
    
   
    % % plot for sanity check
    % plot(x,y, '.');
    % hold on
    % plot(x(in),y(in),'r+') % points inside
    % set(gca,'XAxisLocation','bottom','YAxisLocation','left','ydir','reverse');
    % hold off
    
    cellsInContour = numel(x(in)); % Gets the number of cells in and on the contour

    %% Storing data
    if (i ==1)
        data = [cellsInContour];
    else
        data = [data; cellsInContour];
    end
    fnameList{i} = matrixName;

end

%% Compiling and writing data
data = num2cell(data);
header = {'File Name', '# of Cells'};
compiled = cat(2, fnameList', data);
compiled2 = cat(1, header, compiled);
writecell(compiled2, fullfile(batchDir{1,2}, ['Cells_in_', num2str(input),'um_radius_from_PCD_', datestr(now, 'dd_mmm_yyyy'), '.csv']))


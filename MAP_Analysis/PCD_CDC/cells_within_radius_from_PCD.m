% Finds the nujmber of cells within a circle with a user given radius from
% the PCD in a density matrix. User enters the radius in um (microns) from 
% the PCD. Currently not able to handle both eyes from the same subject.

% Inputs: Denisty matrices, cone coordinate files, a LUT file
% (subject ID, axial length, and ppd) all in the same data folder.

% Outputs: .csv saved to the data folder.

% Created by Jenna Grieshop 11/22/23 - most code was recycled from
% PCD_CDC_Analysis.m code.


clear all
close all
clc

% User selects desired folder:
root_path = uigetdir('.','Select directory containing analyses');
root_dir = dir(root_path);
root_dir = struct2cell(root_dir)';
addpath(genpath(root_path));
impath = fullfile(root_path,'Matlab_Outputs');
mkdir(impath);

basepath = which('cells_within_radius_from_PCD.m');
[basepath] = fileparts(basepath);
path(path,fullfile(basepath,'lib')); % Add our support library to the path.

% looks for the Batch Info spreadsheet that contains the ID initial, ID #, eye, axial length, and device #:
batch_dir = root_dir(...
    ~cellfun(@isempty, strfind(root_dir(:,1), 'LUT')),:);
% Loads in the Batch Info file:
batch_info = table2cell(readtable(fullfile(batch_dir{1,2},batch_dir{1,1})));

% looks for all bound density matrices
bounddensity_dir = root_dir(...
    ~cellfun(@isempty, strfind(root_dir(:,1), 'bounddensity_matrix')),:);

% looks for all bound density matrices
coords_dir = root_dir(...
    ~cellfun(@isempty, strfind(root_dir(:,1), 'coords')),:);


% get user input radius in microns
input = NaN;
while isnan(input)

    input = inputdlg('Input the radius from the PCD in um:', 'Input the radius from the PCD in um:');
    
    input = str2double(input);

    if isempty(input)
            error('radius input cancelled by user.');
    end
end


count = 1;

for i = 1:size(batch_info, 1)
    
    %% Pulling information from batch info:
    subject_ID = num2str(batch_info{i,1});
    axiallength = batch_info{i,2};
    ppd = batch_info{i,3};
    
    subj_matrix = bounddensity_dir(...
    ~cellfun(@isempty, strfind(bounddensity_dir(:,1), subject_ID)),:);

    matrixname = char(subj_matrix(...
    ~cellfun(@isempty, strfind(subj_matrix(:,1), subject_ID)), 1));

    subj_coords = coords_dir(...
    ~cellfun(@isempty, strfind(coords_dir(:,1), subject_ID)),:);

    coordname = char(subj_coords(...
    ~cellfun(@isempty, strfind(subj_coords(:,1), subject_ID)), 1)); 

    %% load in data
    densitymap = csvread(fullfile(bounddensity_dir{1,2}, matrixname));
    coords = dlmread(coordname);
    x = coords(:,1);
    y = coords(:,2);

    %% find PCD
    [maxval, maxind] = max(densitymap(:));

    % finding the weighted (mean) centoid if multiple max locations found
    [max_y_coords, max_x_coords] = find(densitymap == maxval);
    centroid_x = mean(max_x_coords);
    centroid_y = mean(max_y_coords);

    %% Calculate scale and find radius in pixels
    micronsperdegree = (291*axiallength)/24; % these numbers are from a book in Joe's office, Joe has confirmed it is correct
    scaleval = ppd / micronsperdegree; %pixel/um
    radius_px = round(input * scaleval); % radius in pixels

    %% Find cells in radius

    % creating a circle mask of 0s on the density map 
    [x_mesh,y_mesh] = meshgrid(1:size(densitymap,1));
    circle_mask = densitymap;
    circle_mask((x_mesh - centroid_x).^2 + (y_mesh - centroid_y).^2 < radius_px^2) = 0;
    circle_mask2 = (circle_mask>0);

    contour = edge(circle_mask2); % binary contour from thresholded density map
    [new_indeces] = find_close_indeces(contour); % using function from matlab file exchange
    in = inpolygon(x,y,new_indeces(:,2),new_indeces(:,1)); % finds which cell coordinates are inside/on the contour
    
    % 
    % % plot for sanity check
    % plot(x,y, '.');
    % hold on
    % plot(x(in),y(in),'r+') % points inside
    % set(gca,'XAxisLocation','bottom','YAxisLocation','left','ydir','reverse');
    % hold off
    
    cellsInContour = numel(x(in)); % gets the number of cells in and on the contour

    %% storing data
    if (i ==1)
        data = [cellsInContour];
    else
        data = [data; cellsInContour];
    end
    fnamelist{i} = matrixname;

end

%% compiling and writing data
data = num2cell(data);
header = {'File Name', '# of Cells'};
compiled = cat(2, fnamelist', data);
compiled2 = cat(1, header, compiled);
writecell(compiled2, fullfile(batch_dir{1,2}, ['Cells_in_', num2str(input),'um_radius_from_PCD_', datestr(now, 'dd-mmm-yyyy'), '.csv']))


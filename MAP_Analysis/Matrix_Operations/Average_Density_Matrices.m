% Author: Jenna Grieshop
% Date of creation: 2/22/22, v0

% Description: Script that performs element-wise mean and standard devation
% of multiple density matrices. The resulting averaged matrix is displayed
% and saved as a matlab figure. Additionally the averaged matrix and 
% standard deviation matrices are saved as csv files.
% Note: The viridis.m colormap is needed to run this script.

% Description of edits (include date & author):
% Added image saving (2/23/22, Joe Carroll), v1

clear all;
close all;
clc;


% select and load in multiple filenames
[filenames, pathnames] = uigetfile('*.csv', 'MultiSelect', 'on');
filenames = filenames';

% load in the data
for i=1:length(filenames)
    data{i} = readmatrix(filenames{i});
end

% stack the matrices
stackedData = cat(3, data{:});

% take the element-wise mean
meanMatrix = mean(stackedData,3);

% take the element-wise standard deviation
stdevMatrix = std(stackedData,[],3);

% taken from Coordinate_Mosaic_Metrics_MAP_jc.m
vmap=viridis; % calls viridis colormap function
clims = [50000 200000]; % added to set limits of color scale, so all images use the same scale
imagesc(meanMatrix,clims); % added to use limits of color scale
axis image;
colormap(vmap); 
colorbar; 
[minval, minind] = min(meanMatrix(:));
[maxval, maxind] = max(meanMatrix(:));

[minrow,mincol]=ind2sub(size(meanMatrix),minind);
[maxrow,maxcol]=ind2sub(size(meanMatrix),maxind);

max_x_vals = maxcol;
max_y_vals = maxrow;
title(['Minimum value: ' num2str(minval) '(' num2str(mincol) ',' num2str(minrow) ') Maximum value: ' num2str(maxval) '(' num2str(maxcol) ',' num2str(maxrow) ')'])
% end of section taken from Coordinate_Mosaic_Metrics_MAP_jc.m
figname = ['avg_density_map_' date];
saveas(gcf,fullfile([figname '_fig.png']));

% save averaged matrix
filename1 = ['avg_density_matrix_' date '_data.csv'];
csvwrite(filename1,meanMatrix);


% STDEV
[minval2, minind2] = min(stdevMatrix(:));
[maxval2, maxind2] = max(stdevMatrix(:));

[minrow2,mincol2]=ind2sub(size(stdevMatrix),minind2);
[maxrow2,maxcol2]=ind2sub(size(stdevMatrix),maxind2);

figure 

max_x_vals2 = maxcol2;
max_y_vals2 = maxrow2;
vmap=viridis; % calls viridis colormap function
clims = [0 10000]; % added to set limits of color scale, so all images use the same scale
imagesc(stdevMatrix,clims); % added to use limits of color scale
axis image;
colormap(vmap); 
colorbar; 
title(['Minimum value: ' num2str(minval2) '(' num2str(mincol2) ',' num2str(minrow2) ') Maximum value: ' num2str(maxval2) '(' num2str(maxcol2) ',' num2str(maxrow2) ')'])
figname2 = ['stdev_density_map_' date];
saveas(gcf,fullfile([figname2 '_fig.png']));


% save standard deviation matrix
filename3 = ['stdev_density_matrix_' date '.csv'];
csvwrite(filename3,stdevMatrix);
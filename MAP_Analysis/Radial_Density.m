



clear all
close all
clc

scaleFactor = 0.25; %umpp

basepath = which('Radial_Density.m');
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
CDC = [CDC_y, CDC_x];

if CDC_x ~= CDC_y
    warning('CDC is not centered in matrix');
end

px_num = (1/scaleFactor)*5;
dimension = size(data, 1);

points = floor((dimension/px_num)/2);

for i=1:points

    [X,Y]=meshgrid(1:size(data,1),1:size(data,2));
    disk_locations=sqrt((X-CDC(1)).^2+(Y-CDC(2)).^2) <= (px_num*i);
    
    outerboundary = imdilate(disk_locations,strel('disk',1))&~disk_locations;
    outerboundary = outerboundary';
    imshow(outerboundary)
    
    res1 = data(outerboundary==1);
    Avg(i) = mean(res1);


end

plot(Avg);
Avg = Avg';


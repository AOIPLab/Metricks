% Average Density Maps
% AOIP
% Created by: Jenna Grieshop
% Date created: 7/12/2024
%
% Description: Script takes multiple metrics maps (e.g. density maps) and
% aligns then accordingly to the map CDC points to then create average and
% standard deviation maps for the data. This script creates different
% variations of the average and standard deviation maps. Including all
% possible data even in spots where no overlapping occurs, data that has
% at least N (50) maps overlapping, and data that has all matrices
% overlapping.
%
% Input: Directory containing the matrices, a LUT file that contains the
% names of the matrices and the CDC coordinates. The LUT can be in the same
% folder as the matrices.
%
% Output: 
% Outputs are saved in the folder containing the LUT file
% 1. 3 versions of the average and standard deviation matricies. 
% Including all possible data even in spots where no overlapping occurs (ALL), 
% data that has at least N (50) maps overlapping, and data that has all 
% matrices overlapping. These are saved as figures, raw images, and csv
% files. 
% 2. Image and csv of the amount of overlapping datapoints.
% 3. CDC csv of the master CDC for each of the resulting types of map.
%
%



clear all;
close all;
clc;

N=50; % Minimum amount of subjects overlapping (make sure you have at least more that N input maps)

basepath = which('Density_Matrix_Averaging.m');
[basepath] = fileparts(basepath);
path(path,fullfile(basepath,'lib')); % Add our support library to the path.

% User selects folder with data
dataPath = uigetdir('.','Select directory containing analyses');

% Read in csv names and then have user select the LUT
[fnameList] = read_folder_contents(dataPath,'csv');
% select and load in filename of the LUT with CDC
[LUTfilename, LUTpathname] = uigetfile('*.csv', 'Select file with CDC coords', dataPath);

% Remove LUT file from fnameList
fnameList(ismember(fnameList,LUTfilename))=[];

% Load in the LUT file
[~, lutData] = load_LUT_file(fullfile(LUTpathname,LUTfilename));

numFiles = size(fnameList,1);
numFilesDim = size(fnameList);

% Initialize variables
CDC_x = cell(numFilesDim);
CDC_y = cell(numFilesDim);
data = cell(numFilesDim);

tl = cell(numFilesDim);
tr = cell(numFilesDim);
bl = cell(numFilesDim);
br = cell(numFilesDim);

tla = cell(numFilesDim);
tra = cell(numFilesDim);
bla = cell(numFilesDim);
bra = cell(numFilesDim);

tlao = cell(numFilesDim);
trao = cell(numFilesDim);
blao =cell(numFilesDim);
brao = cell(numFilesDim);


for i=1:numFiles % Go through all files in list 

    % Find the information from the LUT file for the data                           
    LUTindex=find( cellfun(@(s) ~isempty(strfind(fnameList{i},s )), lutData{1} ) );
    
    % Extract info from lut file and read in data
    CDC_x{i} = lutData{2}(LUTindex);
    CDC_y{i} = lutData{3}(LUTindex);
    data{i} = readmatrix(fullfile(dataPath,fnameList{i}));

    % Flip OS maps to match the OD maps so that the orientation matches
    if contains(fnameList{i}, 'OS')
        data{i} = fliplr(data{i});
        CDC_x{i} = length(data{i})-(CDC_x{i}-1);
    end

    % Get the orignal coordinates for the corners of the matrices
    l = size(data{i});
    tl{i} = [2,2];
    tr{i} = [l(2)+1,2];
    bl{i} = [2,l(1)+1];
    br{i} = [l(2)+1,l(1)+1];

    % Get the adjusted coordinates for the corners of the matrix. Offset
    % by the CDC coords
    tla{i} = tl{i}-[CDC_x{i},CDC_y{i}];
    tra{i} = tr{i}-[CDC_x{i},CDC_y{i}];
    bla{i} = bl{i}-[CDC_x{i},CDC_y{i}];
    bra{i} = br{i}-[CDC_x{i},CDC_y{i}];

    % Subtract from cdc coords too
    CDC_x{i} = CDC_x{i} - CDC_x{i}+1;
    CDC_y{i} = CDC_y{i} - CDC_y{i}+1;

    % Organize data so it can be used later
    array{i,1} = tla{i}(1);
    array{i,2} = tla{i}(2);
    array{i,3} = tra{i}(1);
    array{i,4} = tra{i}(2);
    array{i,5} = bla{i}(1);
    array{i,6} = bla{i}(2);
    array{i,7} = bra{i}(1);
    array{i,8} = bra{i}(2);

end

% Reformat the to a matrix
array = cell2mat(array);

% Find the minimum of all the coordinates
[minimumx, minIx] = min(array(:,1));
[minimumy, minIy] = min(array(:,2));

% Find the maximum of all the coordinates
[maximumx, maxIx] = max(array(:,7));
[maximumy, maxIy] = max(array(:,8));

% Determine the offset
offset = abs([minimumx-1, minimumy-1]);

for j=1:numFiles

    % Adjust all coordinates by the minimum by adding the offset
    tlao{j} = tla{j} + offset;
    trao{j} = tra{j} + offset;
    blao{j} = bla{j} + offset;
    brao{j} = bra{j} + offset;
    
    % Adjust the cdc by the offset
    CDC_x{j} = CDC_x{j} + offset(1);
    CDC_y{j} = CDC_y{j} + offset(2);

end

for j=1:numFiles

    % Figure out how much padding is needed on each size

    left = abs(tlao{minIx}(1)-tlao{j}(1));
    top = abs(tlao{minIy}(2)-tlao{j}(2));
    right = abs(brao{maxIx}(1)-brao{j}(1));
    bottom = abs(brao{maxIy}(2)-brao{j}(2));

    % Make a copy of the data so that it can be padded without corrupting
    % the original
    paddedData = data{j};

    leftPad = zeros(size(paddedData,1), left);
    paddedData = horzcat(leftPad, paddedData);

    topPad = zeros(top, size(paddedData,2));
    paddedData = vertcat(topPad, paddedData);

    rightPad = zeros(size(paddedData,1), right(1));
    paddedData = horzcat(paddedData, rightPad);

    bottomPad = zeros(bottom, size(paddedData, 2));
    paddedData = vertcat(paddedData, bottomPad);
  
    % Change the padded portion to be Nan
    paddedData(paddedData==0) = NaN;

    % Overwrite data with the padded version
    data{j} = paddedData;
    % Combine all data into 3D matrix so they are stacked
    A(:,:,j) = data{j};

end
    
    
AN = A;
Atotal = A;
% Figure out how many maps contributed to each location
NanCount = NaN(size(A,1),size(A,2));
dim = size(A,3);
for m=1:size(A,1)
    for n = 1:size(A,2)
        NanCount(m,n) = dim - sum(isnan(A(m,n,:)));
        if NanCount(m,n) < N
            % If less than N get rid of the spot for the N map one
            AN(m,n,:) = NaN;
        end
        if NanCount(m,n) < numFiles
            % If less than total number of maps get rid of the spot for the N map one
            Atotal(m,n,:) = NaN;
        end
    end
end

% These are the x and y coordinates in order
CDC_og = [CDC_x{1}, CDC_y{1}];
CDC_N = CDC_og;
CDC_tot = CDC_og;


%%
% Get the average of the stacked maps
average_all = mean(A,3, "omitnan");
standard_dev_all = std(A, [], 3, "omitmissing");

og_size = size(average_all);

average_tot = mean(Atotal,3, "omitnan");
standard_dev_tot = std(Atotal, [], 3, "omitmissing");

average_N = mean(AN,3, "omitnan");
standard_dev_N = std(AN, [], 3, "omitmissing");


% Figure out how much will need to adjust CDC coords for new maps once full
% NaN rows and columns deleted

% Total  
for q=1:size(average_tot,2)
    if sum(average_tot(:,q), "omitnan") > 0
        break;
    else
        CDC_tot(1) = CDC_tot(1)- 1;
    end
end

for r=1:size(average_tot,1)
    if sum(average_tot(r,:), "omitnan") > 0
        break;
    else
        CDC_tot(2) = CDC_tot(2)- 1;
    end
end

% N
for s=1:size(average_N,2)
    if sum(average_N(:,s), "omitnan") > 0
        break;
    else
        CDC_N(1) = CDC_N(1)- 1;
    end
end

for t=1:size(average_N,1)
    if sum(average_N(t,:), "omitnan") > 0
        break;
    else
        CDC_N(2) = CDC_N(2)- 1;
    end
end


% Get rid of Rows and Columns with only Nans
average_tot = average_tot(:,~all(isnan(average_tot))); 
average_tot = average_tot(~all(isnan(average_tot),2), :); 
average_N = average_N(:,~all(isnan(average_N))); 
average_N = average_N(~all(isnan(average_N),2), :); 
standard_dev_tot = standard_dev_tot(:,~all(isnan(standard_dev_tot))); 
standard_dev_tot = standard_dev_tot(~all(isnan(standard_dev_tot),2), :); 
standard_dev_N = standard_dev_N(:,~all(isnan(standard_dev_N))); 
standard_dev_N = standard_dev_N(~all(isnan(standard_dev_N),2), :); 


% Figure out the max square area without NaNs for the N matrix
new_dim = min(CDC_N(1)) - 1;

for j = 1:size(average_N,1)
    if isnan(average_N(CDC_N(1) - new_dim,j)) || isnan(average_N(CDC_N(1) + new_dim,j)) || isnan(average_N(j,CDC_N(2)- new_dim)) || isnan(average_N(j,CDC_N(2) + new_dim))
        new_dim = new_dim - 1;
    else
        break;
    end
end

% Adjust the CDC coordinates
x_diff = CDC_N(1) - new_dim;
CDC_N(1) = CDC_N(1) - x_diff;

y_diff = CDC_N(2) - new_dim;
CDC_N(2) = CDC_N(2) - y_diff;

% Delete the rows and columns that do not fit in the largest square
try
    average_N(1:y_diff,:) = [];
    average_N(:,1:x_diff) = [];
    
    average_N((new_dim*2):end,:) = [];
    average_N(:,(new_dim*2):end) = [];

catch
    error("N is too large, please reduce on line 35.")
end
standard_dev_N(1:y_diff,:) = [];
standard_dev_N(:,1:x_diff) = [];

standard_dev_N((new_dim*2):end,:) = [];
standard_dev_N(:,(new_dim*2):end) = [];

% Added to set limits of color scale
clims = [50000 225000]; %Density%[1.5 4.5]; %NND %
clims2 = [5000 45000]; %Density%[0.1 0.3]; %NND %

vmap=viridis; %calls viridis colormap function

% Display average of everything map
dispfig=figure(1); 
imagesc(average_all,clims); % added to use limits of color scale
axis image;
colormap(vmap); 
colorbar; 

saveas(gcf,fullfile(LUTpathname, ['Averaged_ALL_bound_density_map_' datestr(now, 'yyyymmdd') '_fig.png']));

% Display average of just overlapping map
dispfig=figure(2); 
imagesc(average_tot,clims); % added to use limits of color scale
axis image;
colormap(vmap); 
colorbar; 

saveas(gcf,fullfile(LUTpathname, ['Averaged_'  num2str(numFiles) '_bound_density_map_' datestr(now, 'yyyymmdd') '_fig.png']));

% Display average of overlapping N subjects
dispfig=figure(3); 
imagesc(average_N,clims); % added to use limits of color scale
axis image;
colormap(vmap); 
colorbar; 

saveas(gcf,fullfile(LUTpathname, ['Averaged_' num2str(N) '_bound_density_map_' datestr(now, 'yyyymmdd') '_fig.png']));


% Display stdev of everything map
dispfig=figure(4); 
imagesc(standard_dev_all,clims2); % added to use limits of color scale
axis image;
colormap(vmap); 
colorbar; 

saveas(gcf,fullfile(LUTpathname, ['Stdev_ALL_bound_density_map_' datestr(now, 'yyyymmdd') '_fig.png']));

% Display stdev of just overlapping map
dispfig=figure(5); 
imagesc(standard_dev_tot,clims2); % added to use limits of color scale
axis image;
colormap(vmap); 
colorbar; 

saveas(gcf,fullfile(LUTpathname, ['Stdev_'  num2str(numFiles) '_bound_density_map_' datestr(now, 'yyyymmdd') '_fig.png']));

% Display stdev of just overlapping map of N
dispfig=figure(6); 
imagesc(standard_dev_N,clims2); % added to use limits of color scale
axis image;
colormap(vmap); 
colorbar; 

saveas(gcf,fullfile(LUTpathname, ['Stdev_' num2str(N) '_bound_density_map_' datestr(now, 'yyyymmdd') '_fig.png']));


% Display NanCount map
dispfig=figure(7); 
imagesc(NanCount);
axis image;
colormap(vmap); 
colorbar; 

saveas(gcf,fullfile(LUTpathname, ['Included_datapoint_map_' datestr(now, 'yyyymmdd') '_fig.png']));
imwrite(NanCount, vmap, fullfile(LUTpathname,['Included_datapoint_map_' datestr(now, 'yyyymmdd') '_raw.tif'])); 

% Average
scaled_map_all = average_all-min(clims);
scaled_map_all(scaled_map_all <0) =0; %in case there are min values below this
scaled_map_all = uint8(255*scaled_map_all./(max(clims)-min(clims)));
scaled_map_all(scaled_map_all  >255) = 255; %in case there are values above this
imwrite(scaled_map_all, vmap, fullfile(LUTpathname,['Averaged_ALL_bound_density_map_' datestr(now, 'yyyymmdd') '_raw.tif'])); 


scaled_map = average_tot-min(clims);
scaled_map(scaled_map <0) =0; %in case there are min values below this
scaled_map = uint8(255*scaled_map./(max(clims)-min(clims)));
scaled_map(scaled_map  >255) = 255; %in case there are values above this
imwrite(scaled_map, vmap, fullfile(LUTpathname,['Averaged_' num2str(numFiles) '_bound_density_map_' datestr(now, 'yyyymmdd') '_raw.tif'])); 

scaled_mapN = average_N-min(clims);
scaled_mapN(scaled_mapN <0) =0; %in case there are min values below this
scaled_mapN = uint8(255*scaled_mapN./(max(clims)-min(clims)));
scaled_mapN(scaled_mapN  >255) = 255; %in case there are values above this
imwrite(scaled_mapN, vmap, fullfile(LUTpathname,['Averaged_' num2str(N) '_bound_density_map_' datestr(now, 'yyyymmdd') '_raw.tif'])); 


%stdev
scaled_map_stdev_all = standard_dev_all-min(clims2);
scaled_map_stdev_all(scaled_map_stdev_all <0) =0; %in case there are min values below this
scaled_map_stdev_all = uint8(255*scaled_map_stdev_all./(max(clims2)-min(clims2)));
scaled_map_stdev_all(scaled_map_stdev_all  >255) = 255; %in case there are values above this
imwrite(scaled_map_stdev_all, vmap, fullfile(LUTpathname,['Stdev_ALL_bound_density_map_' datestr(now, 'yyyymmdd') '_raw.tif'])); 


scaled_map_stdev = standard_dev_tot-min(clims2);
scaled_map_stdev(scaled_map_stdev <0) =0; %in case there are min values below this
scaled_map_stdev = uint8(255*scaled_map_stdev./(max(clims2)-min(clims2)));
scaled_map_stdev(scaled_map_stdev  >255) = 255; %in case there are values above this
imwrite(scaled_map_stdev, vmap, fullfile(LUTpathname,['Stdev_'  num2str(numFiles) '_bound_density_map_' datestr(now, 'yyyymmdd') '_raw.tif'])); 

scaled_map_stdevN = standard_dev_N-min(clims2);
scaled_map_stdevN(scaled_map_stdevN <0) =0; %in case there are min values below this
scaled_map_stdevN = uint8(255*scaled_map_stdevN./(max(clims2)-min(clims2)));
scaled_map_stdevN(scaled_map_stdevN  >255) = 255; %in case there are values above this
imwrite(scaled_map_stdevN, vmap, fullfile(LUTpathname,['Stdev_' num2str(N) '_bound_density_map_' datestr(now, 'yyyymmdd') '_raw.tif'])); 


% Save matrices
writematrix(average_all, fullfile(LUTpathname,['Averaged_ALL_bound_density_map_' datestr(now, 'yyyymmdd') '_raw.csv']));
writematrix(average_tot,  fullfile(LUTpathname,['Averaged_' num2str(numFiles) '_bound_density_map_' datestr(now, 'yyyymmdd') '_raw.csv']));
writematrix(average_N,  fullfile(LUTpathname,['Averaged_' num2str(N) '_bound_density_map_' datestr(now, 'yyyymmdd') '_raw.csv']));
writematrix(standard_dev_all, fullfile(LUTpathname,['Stdev_ALL_bound_density_map_' datestr(now, 'yyyymmdd') '_raw.csv']));
writematrix(standard_dev_tot,  fullfile(LUTpathname,['Stdev_' num2str(numFiles) '_bound_density_map_' datestr(now, 'yyyymmdd') '_raw.csv']));
writematrix(standard_dev_N,  fullfile(LUTpathname,['Stdev_' num2str(N) '_bound_density_map_' datestr(now, 'yyyymmdd') '_raw.csv']));
writematrix(NanCount, fullfile(LUTpathname,['Included_datapoint_map_' datestr(now, 'yyyymmdd') '_raw.csv']));
header = {['ALL_bound_density_map_' datestr(now, 'yyyymmdd') '_raw']; [ num2str(numFiles) '_bound_density_map_' datestr(now, 'yyyymmdd') '_raw']; [num2str(N) '_bound_density_map_' datestr(now, 'yyyymmdd') '_raw']};
result = [header, num2cell([CDC_og; CDC_tot; CDC_N])];
writecell(result, fullfile(LUTpathname,['CDC_' datestr(now, 'yyyymmdd') '.csv']));



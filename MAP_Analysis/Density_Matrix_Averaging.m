% Average Density Maps
% 
% Created by: Jenna Grieshop
% Date created: 7/12/2024
%
%

clear all;
close all;
clc;

basepath = which('Density_Matrix_Averaging.m');
[basepath] = fileparts(basepath);
path(path,fullfile(basepath,'lib')); % Add our support library to the path.

% User selects folder with data
dataPath = uigetdir('.','Select directory containing analyses');

% Read in csv names and then have user select the LUT
[fnameList] = read_folder_contents(dataPath,'csv');
% select and load in filename of the LUT with CDC
[LUTfilename, LUTpathname] = uigetfile('*.csv', 'Select file with CDC coords');

% Remove LUT file from fnameList
fnameList(ismember(fnameList,LUTfilename))=[];

% load in the LUT file
[~, lutData] = load_LUT_file(fullfile(LUTpathname,LUTfilename));

numFiles = size(fnameList,1);
numFilesDim = size(fnameList);

CDC_x = cell(numFilesDim);
CDC_y = cell(numFilesDim);
data = cell(numFilesDim);
sz = cell(numFilesDim);

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

    % Calculate the scale for this identifier.                                
    LUTindex=find( cellfun(@(s) ~isempty(strfind(fnameList{i},s )), lutData{1} ) );
    
    CDC_x{i} = lutData{2}(LUTindex);
    CDC_y{i} = lutData{3}(LUTindex);
    data{i} = readmatrix(fullfile(dataPath,fnameList{i}));

    % figure out the sizes of the matricies
    sz{i} = size(data{i},1);

    % get the orignal coordinates for the corners of the matrices
    l = size(data{i});
    tl{i} = [2,2];
    tr{i} = [l(2)+1,2];
    bl{i} = [2,l(1)+1];
    br{i} = [l(2)+1,l(1)+1];

    % get the adjusted coordinates for the corners of the matrix. Offset
    % by the CDC coords
    tla{i} = tl{i}-[CDC_x{i},CDC_y{i}];
    tra{i} = tr{i}-[CDC_x{i},CDC_y{i}];
    bla{i} = bl{i}-[CDC_x{i},CDC_y{i}];
    bra{i} = br{i}-[CDC_x{i},CDC_y{i}];

    %subtract from cdc coords too
    CDC_x{i} = CDC_x{i} - CDC_x{i}+1;
    CDC_y{i} = CDC_y{i} - CDC_y{i}+1;

    array{i,1} = tla{i}(1);
    array{i,2} = tla{i}(2);
    array{i,3} = tra{i}(1);
    array{i,4} = tra{i}(2);
    array{i,5} = bla{i}(1);
    array{i,6} = bla{i}(2);
    array{i,7} = bra{i}(1);
    array{i,8} = bra{i}(2);

end


array = cell2mat(array);

% find the minimum of all the coordinates
[minimumx, minIx] = min(array(:,1));
[minimumy, minIy] = min(array(:,2));

[maximumx, maxIx] = max(array(:,7));
[maximumy, maxIy] = max(array(:,8));


offset = abs([minimumx-1, minimumy-1]);

for j=1:numFiles

    % adjust all coordinates by the minimum by adding the offset
    tlao{j} = tla{j} + offset;
    trao{j} = tra{j} + offset;
    blao{j} = bla{j} + offset;
    brao{j} = bra{j} + offset;
    
    
    CDC_x{j} = CDC_x{j} + offset(1);
    CDC_y{j} = CDC_y{j} + offset(2);

end

for j=1:numFiles

    % figure out how much padding is needed

    left = abs(tlao{minIx}(1)-tlao{j}(1));
    top = abs(tlao{minIy}(2)-tlao{j}(2));
    right = abs(brao{maxIx}(1)-brao{j}(1));
    bottom = abs(brao{maxIy}(2)-brao{j}(2));

    paddedData = data{j};

    leftPad = zeros(size(paddedData,1), left);
    paddedData = horzcat(leftPad, paddedData);

    topPad = zeros(top, size(paddedData,2));
    paddedData = vertcat(topPad, paddedData);

    rightPad = zeros(size(paddedData,1), right(1));
    paddedData = horzcat(paddedData, rightPad);

    bottomPad = zeros(bottom, size(paddedData, 2));
    paddedData = vertcat(paddedData, bottomPad);
  
    % change the padded portion to be Nan
    paddedData(paddedData==0) = NaN;

    data{j} = paddedData;
    A(:,:,j) = data{j};

end
    
    
    A50 = A;
    % figure out how many maps contributed to each location
    NanCount = NaN(size(A,1),size(A,2));
    dim = size(A,3);
    for m=1:size(A,1)
        for n = 1:size(A,2)
            NanCount(m,n) = dim - sum(isnan(A(m,n,:)));
            if NanCount(m,n) < 30
                % if less than 50 get rid of the spot for the 50 map one
                A50(m,n,:) = NaN;
            end
        end
    end

%%
    % get the average of the stacked maps
    average_all = mean(A,3, "omitnan");
    standard_dev_all = std(A, [], 3, "omitmissing");

    average = mean(A,3);
    standard_dev = std(A, [], 3);

    average50 = mean(A50,3, "omitnan");
    standard_dev50 = std(A50, [], 3, "omitmissing");

    % writematrix(average_all,fullfile(LUTpathname, ['Averaged_bound_density_map_wNans' datestr(now, 'yyyymmdd') '.csv']));

    clims = [50000 225000]; % added to set limits of color scale
    clims2 = [5000 45000];

    %average
    % average_all = average_all(:,~all(isnan(average_all))); % get rid of columns with only nans
    % average_all = average_all(~all(isnan(average_all),2), :); % get rid of row with only nans
    average = average(:,~all(isnan(average))); % get rid of columns with only nans
    average = average(~all(isnan(average),2), :); % get rid of row with only nans
    average50 = average50(:,~all(isnan(average50))); % get rid of columns with only nans
    average50 = average50(~all(isnan(average50),2), :); % get rid of row with only nans
    % get rid of any NaNs
    average50(any(isnan(average50), 2), :) = [];
    average50(any(isnan(average50), 1), :) = [];
    % average_all(any(isnan(average_all), 2), :) = [];
    % average_all(any(isnan(average_all), 1), :) = [];
    % average(any(isnan(average), 2), :) = [];
    % average(any(isnan(average), 1), :) = [];

    %standard deviation
    % standard_dev_all = standard_dev_all(:,~all(isnan(standard_dev_all))); % get rid of columns with only nans
    % standard_dev_all = standard_dev_all(~all(isnan(standard_dev_all),2), :); % get rid of row with only nans
    standard_dev = standard_dev(:,~all(isnan(standard_dev))); % get rid of columns with only nans
    standard_dev = standard_dev(~all(isnan(standard_dev),2), :); % get rid of row with only nans
    standard_dev50 = standard_dev50(:,~all(isnan(standard_dev50))); % get rid of columns with only nans
    standard_dev50 = standard_dev50(~all(isnan(standard_dev50),2), :); % get rid of row with only nans
    % get rid of any NaNs
    standard_dev50(any(isnan(standard_dev50), 2), :) = [];
    standard_dev50(any(isnan(standard_dev50), 1), :) = [];
    % standard_dev_all(any(isnan(standard_dev_all), 2), :) = [];
    % standard_dev_all(any(isnan(standard_dev_all), 1), :) = [];
    % standard_dev(any(isnan(standard_dev), 2), :) = [];
    % standard_dev(any(isnan(standard_dev), 1), :) = [];

    vmap=viridis; %calls viridis colormap function

    % display average of everything map
    dispfig=figure(1); 
    imagesc(average_all,clims); % added to use limits of color scale
    axis image;
    colormap(vmap); 
    colorbar; 

    saveas(gcf,fullfile(LUTpathname, ['Averaged_ALL_bound_density_map_' datestr(now, 'yyyymmdd') '_fig.png']));

    % display average of just overlapping map
    dispfig=figure(2); 
    imagesc(average,clims); % added to use limits of color scale
    axis image;
    colormap(vmap); 
    colorbar; 

    saveas(gcf,fullfile(LUTpathname, ['Averaged_bound_density_map_' datestr(now, 'yyyymmdd') '_fig.png']));

    % display average of overlapping 50 subjects
    dispfig=figure(3); 
    imagesc(average50,clims); % added to use limits of color scale
    axis image;
    colormap(vmap); 
    colorbar; 

    saveas(gcf,fullfile(LUTpathname, ['Averaged_50_bound_density_map_' datestr(now, 'yyyymmdd') '_fig.png']));


    % display stdev of everything map
    dispfig=figure(4); 
    imagesc(standard_dev_all,clims2); % added to use limits of color scale
    axis image;
    colormap(vmap); 
    colorbar; 

    saveas(gcf,fullfile(LUTpathname, ['Stdev_ALL_bound_density_map_' datestr(now, 'yyyymmdd') '_fig.png']));

    % display stdev of just overlapping map
    dispfig=figure(5); 
    imagesc(standard_dev,clims2); % added to use limits of color scale
    axis image;
    colormap(vmap); 
    colorbar; 

    saveas(gcf,fullfile(LUTpathname, ['Stdev_bound_density_map_' datestr(now, 'yyyymmdd') '_fig.png']));

    % display stdev of just overlapping map
    dispfig=figure(6); 
    imagesc(standard_dev50,clims2); % added to use limits of color scale
    axis image;
    colormap(vmap); 
    colorbar; 

    saveas(gcf,fullfile(LUTpathname, ['Stdev_50_bound_density_map_' datestr(now, 'yyyymmdd') '_fig.png']));


    % display NanCount map
    dispfig=figure(7); 
    imagesc(NanCount);
    axis image;
    colormap(vmap); 
    colorbar; 

    saveas(gcf,fullfile(LUTpathname, ['Included_datapoint_map_' datestr(now, 'yyyymmdd') '_fig.png']));
    imwrite(NanCount, vmap, fullfile(LUTpathname,['Included_datapoint_map_' datestr(now, 'yyyymmdd') '_raw.tif'])); 
  
    % average
    scaled_map_all = average_all-min(clims);
    scaled_map_all(scaled_map_all <0) =0; %in case there are min values below this
    scaled_map_all = uint8(255*scaled_map_all./(max(clims)-min(clims)));
    scaled_map_all(scaled_map_all  >255) = 255; %in case there are values above this
    imwrite(scaled_map_all, vmap, fullfile(LUTpathname,['Averaged_ALL_bound_density_map_' datestr(now, 'yyyymmdd') '_raw.tif'])); 


    scaled_map = average-min(clims);
    scaled_map(scaled_map <0) =0; %in case there are min values below this
    scaled_map = uint8(255*scaled_map./(max(clims)-min(clims)));
    scaled_map(scaled_map  >255) = 255; %in case there are values above this
    imwrite(scaled_map, vmap, fullfile(LUTpathname,['Averaged_bound_density_map_' datestr(now, 'yyyymmdd') '_raw.tif'])); 

    scaled_map50 = average50-min(clims);
    scaled_map50(scaled_map50 <0) =0; %in case there are min values below this
    scaled_map50 = uint8(255*scaled_map50./(max(clims)-min(clims)));
    scaled_map50(scaled_map50  >255) = 255; %in case there are values above this
    imwrite(scaled_map50, vmap, fullfile(LUTpathname,['Averaged_50_bound_density_map_' datestr(now, 'yyyymmdd') '_raw.tif'])); 


    %stdev
    scaled_map_stdev_all = standard_dev_all-min(clims2);
    scaled_map_stdev_all(scaled_map_stdev_all <0) =0; %in case there are min values below this
    scaled_map_stdev_all = uint8(255*scaled_map_stdev_all./(max(clims2)-min(clims2)));
    scaled_map_stdev_all(scaled_map_stdev_all  >255) = 255; %in case there are values above this
    imwrite(scaled_map_stdev_all, vmap, fullfile(LUTpathname,['Stdev_ALL_bound_density_map_' datestr(now, 'yyyymmdd') '_raw.tif'])); 


    scaled_map_stdev = standard_dev-min(clims2);
    scaled_map_stdev(scaled_map_stdev <0) =0; %in case there are min values below this
    scaled_map_stdev = uint8(255*scaled_map_stdev./(max(clims2)-min(clims2)));
    scaled_map_stdev(scaled_map_stdev  >255) = 255; %in case there are values above this
    imwrite(scaled_map_stdev, vmap, fullfile(LUTpathname,['Stdev_bound_density_map_' datestr(now, 'yyyymmdd') '_raw.tif'])); 

    scaled_map_stdev50 = standard_dev50-min(clims2);
    scaled_map_stdev50(scaled_map_stdev50 <0) =0; %in case there are min values below this
    scaled_map_stdev50 = uint8(255*scaled_map_stdev50./(max(clims2)-min(clims2)));
    scaled_map_stdev50(scaled_map_stdev50  >255) = 255; %in case there are values above this
    imwrite(scaled_map_stdev50, vmap, fullfile(LUTpathname,['Stdev_50_bound_density_map_' datestr(now, 'yyyymmdd') '_raw.tif'])); 


    %save matrices
    writematrix(average_all, fullfile(LUTpathname,['Averaged_ALL_bound_density_map_' datestr(now, 'yyyymmdd') '_raw.csv']));
    writematrix(average,  fullfile(LUTpathname,['Averaged_bound_density_map_' datestr(now, 'yyyymmdd') '_raw.csv']));
    writematrix(average50,  fullfile(LUTpathname,['Averaged_50_bound_density_map_' datestr(now, 'yyyymmdd') '_raw.csv']));
    writematrix(standard_dev_all, fullfile(LUTpathname,['Stdev_ALL_bound_density_map_' datestr(now, 'yyyymmdd') '_raw.csv']));
    writematrix(standard_dev,  fullfile(LUTpathname,['Stdev_bound_density_map_' datestr(now, 'yyyymmdd') '_raw.csv']));
    writematrix(standard_dev50,  fullfile(LUTpathname,['Stdev_50_bound_density_map_' datestr(now, 'yyyymmdd') '_raw.csv']));
    writematrix(NanCount, fullfile(LUTpathname,['Included_datapoint_map_' datestr(now, 'yyyymmdd') '_raw.csv']));
    writecell([CDC_x(1), CDC_y(1)], fullfile(LUTpathname,['CDC_' datestr(now, 'yyyymmdd') '.csv']));
   


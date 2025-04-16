% Average Density Maps
% AOIP
% Created by: Jenna Grieshop
% Date created: 7/12/2024
%
% Description: This is a modified script of the Density Matrix Averaging
% script that aligns density maps to the common area centered on the CDC
% and then extracts the row of density values from the row centered on the
% CDC for each subject only in the overlapping columns. It then averages
% the density values in that row for each subject.
%
% Input: Directory containing the matrices, a LUT file that contains the
% names of the matrices and the CDC coordinates. The LUT can be in the same
% folder as the matrices.
%
% Output: 
% Outputs are saved in the folder containing the LUT file
% 1. A csv containing subject name and averaged density through the CDC
% point (horizontally)
%
%



clear all;
close all;
clc;

N=50; % Minimum amount of subjects overlapping (make sure you have at least more that N input maps)

basepath = which('Density_Matrix_Averaging_Common_Area_Extraction.m');
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
            template(m,n,:) = NaN;
        else
            template(m,n,:) = 0;
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


% Extract data and ensure CDC point is centered in available square area

for w=1:size(data,1)
    combo = template + data{w,1};
    common_area_1side = size(average_tot,2) - CDC_tot(1);
    left = CDC_og(1) - common_area_1side;
    right = CDC_og(1) + common_area_1side;
    common(:,:,w) = combo(CDC_og(2), left:right);
    common_avg(w) = mean(common(:,:,w));

end
CDC_tot_old = CDC_tot(1);
CDC_tot_extracted = common_area_1side+1; % plus 1 because there are common area 1size on either side of the CDC

% average the common_row_extraction to see if I get the same matrix as the
% tot matrix

target = average_tot(CDC_tot(2),:);
test = mean(common, 3); % used as a test to ensure that the averaged CDC from our extracted rows match what our original averaged was - to ensure the correct columns were taken

% sanity check
if(target(1,CDC_tot_old)) ~= test(1,CDC_tot(1))
    disp('warning!! Something is wrong')
end

common_avg = common_avg';
avg_result = [fnameList, num2cell(common_avg)];

writecell(avg_result, fullfile(LUTpathname,['Avg_Horizontal_row_CDC_' num2str(size(common,2)) '_pxWidth_' datestr(now, 'yyyymmdd') '.csv']));


%% graph individual plots
% basic plot of the individual results
% 
% figure(1)
% plot(common(1,:,i));
% 
% figure(2)
% plot(test);


   




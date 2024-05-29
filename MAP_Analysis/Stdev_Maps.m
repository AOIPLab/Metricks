% Standard Deviation MAPs (& Average and CoV Maps)
%
% 1/26/2024
% Jenna Grieshop
%
% Description: This scripts determines standard deviation of density from
% previously made density maps. Originally used to compare density maps
% made using different window sizes. This script uses subject ID to know
% which ones to do standard deviation on. The data must be the same scale
% and the same size, and across the same subjects.
%
% Input: User to select the directory with the bound density analyses and
% the directory with the cdc analyses. The script uses the bound density
% matrix matfiles (with MATFILE in the name) and the PCD_CDC analysis
% summary file.
%
% Output: raw tif of the resulting standard deviation map, marked tif of
% the standard deviation map with the master CDC location, raw .csvs for
% standard deviation, coefficient of variance, and average. A .csv of the
% master CDC list.
%
%

clear all
clc
addpath('lib');

clims = [0 25000]; % added to set limits of color scale, so all images use the same scale

sub_id = {''};

% user to select the directories
root_path_bd = uigetdir('.','Select directory containing bound density analyses');

root_path_cdc = uigetdir('.','Select directory containing CDC analyses');

% get all the file paths that we are interested in in each of the folders
[file_paths_bd] = read_folder_contents_rec(root_path_bd, 'mat', 'MATFILE');

% get all the file paths that we are interested in in each of the folders
[file_paths_cdc] = read_folder_contents_rec(root_path_cdc, 'csv', 'PCD_CDC_Analysis_Summary');


% finding all the subject Ids within the folders
for i=1:size(file_paths_bd,1)
    spl = split(file_paths_bd{i,1}, "_bound");
    spl = split(spl{1,1}, "\");
    spl = spl(end);
    ispresent = cellfun(@(s) ~isempty(strfind(spl{:}, s)), sub_id);
    if any(ispresent)
        continue
    else
        sub_id{i} = spl{1};
    end
    
end

% initialize variables to have the correct dimensions
sum_x = cell(size(sub_id, 2),1);
sum_x(:,1) = {0};
sum_y = cell(size(sub_id, 2),1);
sum_y(:,1) = {0};

% load in the cdc files and sum the cdc coordinates for each subject ID
for m=1:size(file_paths_cdc,1)
    cdc_data = readtable(file_paths_cdc{m});
    for n=1:size(cdc_data, 1)
        sum_x{n} = sum_x{n} + cdc_data{n,8};
        sum_y{n} = sum_y{n} + cdc_data{n,9};

    end
end

sum_x = cell2mat(sum_x(:)); % convert so can be divided
sum_y = cell2mat(sum_y(:)); % convert so can be divided

% get the average CDC for each subject
avg_x(:) = sum_x(:)/size(file_paths_cdc,1);
avg_y(:) = sum_y(:)./size(file_paths_cdc,1);

% go through loop for each subject
for j=1:size(sub_id,2)

    clear maps;
    clear A;
    
    % find indices in the file name list that belong to the subject ID
    index = find(contains(file_paths_bd,sub_id(j)));
    for m=1:size(index,1)
        % load in and stack the maps
        data = load(file_paths_bd{index(m)});
        maps{m} = round(data.interped_map,2);
        A(:,:,m) = maps{m}; % maps is an unnecesary middle step - but good for trouble shooting
    end

    % take the standard deviation across the stacked maps
    standard_dev = std(A, [], 3);
    % get the average of the stacked maps
    average = mean(A,3);
    % calculate the coefficient of variance
    coeffofvar = standard_dev/average;
    vmap=viridis; %calls viridis colormap function
   
    
    % display the standard deviation map
    dispfig=figure(1); 
    imagesc(standard_dev); % added to use limits of color scale
    axis image;
    colormap(vmap); 
    colorbar; 

    % find the min and max of the standard deviation
    [minval, minind] = min(standard_dev(:));
    [maxval, maxind] = max(standard_dev(:));
    
    [minrow,mincol]=ind2sub(size(standard_dev),minind);
    [maxrow,maxcol]=ind2sub(size(standard_dev),maxind);
    
    max_x_vals = maxcol;
    max_y_vals = maxrow;

    %updated to scale to the max of clims 10/11/23     
    scaled_map = standard_dev-min(clims);
    scaled_map(scaled_map <0) =0; %in case there are min values below this
    scaled_map = uint8(255*scaled_map./(max(clims)-min(clims)));
    scaled_map(scaled_map  >255) = 255; %in case there are values above this

    % format the subject ID
    subjectID = sub_id(j); 
    subjectID = subjectID{1};

    % save results
    result_fname = [subjectID '_stdev_' datestr(now, 'dd_mmm_yyyy') '_raw.tif'];
    imwrite(scaled_map, vmap, fullfile(root_path_bd,result_fname));

    result_fname2 = [subjectID '_stdev_' datestr(now, 'dd_mmm_yyyy') '_marked.tif'];
    scaled_map_mark = uint8(255*standard_dev./max(clims));

    MARK = insertShape(scaled_map_mark,'circle',[avg_x(j) avg_y(j) 2], 'LineWidth' ,3, 'Color' , 'red');
    imwrite(MARK, vmap, fullfile(root_path_bd,result_fname2));

    result_fname3 = [subjectID '_stdev_' datestr(now, 'dd_mmm_yyyy') '_raw.csv'];
    csvwrite(fullfile(root_path_bd, result_fname3), standard_dev);

    result_fname5 = [subjectID '_coeffvar_' datestr(now, 'dd_mmm_yyyy') '_raw.csv'];
    csvwrite(fullfile(root_path_bd, result_fname5), coeffofvar);

    result_fname6 = [subjectID '_average_' datestr(now, 'dd_mmm_yyyy') '_raw.csv'];
    csvwrite(fullfile(root_path_bd, result_fname6), average);
    
    % record the master CDC for each subject
    master_cdc(j,1) = {sub_id{j}};
    master_cdc(j,2) = {avg_x(j)};
    master_cdc(j,3) = {avg_y(j)};
  
end

% save the master CDC list
header = {'Subject ID', 'X CDC master', 'Y CDC master'};
finaloutput = [header; master_cdc];
result_fname4 = ['stdev_' datestr(now, 'dd_mmm_yyyy') '_master_cdc.xlsx'];
xlswrite(fullfile(root_path_bd, result_fname4), master_cdc);



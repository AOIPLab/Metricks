% Row and Column Extraction Simplified
%
% Created by: Jenna Grieshop
% Date created: 6/10/2026
%
% Description: Extracts straight or maltese (radial wedges) meridians 
% originating at the at the CDC coordinates for each subject.
%
% Inputs: User selects directory containing matrices in mat format
% (could be results from Stdev_Maps, matrix subtraction, etc.). NOTE: if
% the matrices are not bound density mm matrices, change the file
% identifier on line 25. User is asked if straight cross or maltese. Then
% how many rows/columns to average if straight ran or angle for the
% maltese.
% 
% Then user selects the PCD CDC analysis summary file with OS/OD eye column 
% added after the file name column and PPD column added as the last column.
%
% Outputs: .csv files with N, T, S, I strips extracted from the density
% matrices through the CDC for straight, maltese, or both.
% 

close all
clear all
clc

maltese = 0;
straight = 0;

identifier = 'bounddensity_mm_matrix';

Meridian = {'Temporal', 'Nasal', 'Superior', 'Inferior'};

basePath = which('Row_Column_Extraction_Simplified.m');

[basePath ] = fileparts(basePath);
path(path,fullfile(basePath,'lib')); % Add our support library to the path.

% user to select directory with data
root_path_raw = uigetdir('.','Select directory containing raw mat analyses');
mkdir(root_path_raw,'Results');

output_root = fullfile(root_path_raw,'Results');


% get all the file paths that we are interested in in each of the folders
[file_names] = read_folder_contents(root_path_raw, 'mat'); % used for stdev map

% select cdc analysis file
[filename_cdc, pathname_cdc] = uigetfile('*.csv', 'Please select the pcd cdc analysis summary file', 'MultiSelect', 'off');


% load in cdc LUT data and extract all the subject IDs
LUT_data = readtable(fullfile(pathname_cdc, filename_cdc));

% Unit selection by the user
list2 = {'Microns', 'Degrees'};
[indx2, tf2] = listdlg('PromptString', 'Select the x-axis units.', 'SelectionMode', 'single', 'ListString', list2);
if tf2 == 1
    unitStr = list2{indx2};
else
    % Canceled dialog box - end the program
    return
end

% Have user select what kind of row/column analysis they'd like to do
list = {'Straight Cross', 'Maltese Cross'};
[indx, tf] = listdlg('PromptString', 'Select the following row/column analyses to perform.', 'SelectionMode', 'multiple', 'ListString', list);
if tf == 1
    selections = size(indx,2);
    if selections == 2
        maltese = 1;
        straight = 1;
    else
        indx_max = max(indx);
        if indx_max == 1
            straight = 1;
        else
            maltese = 1;
        end
    end
    
else
    % Canceled dialog box - end the program
    return
end

% Have the user select the thickness of the row/column extraction
if straight == 1
    m1 = 'Please enter the desired straight extraction thickness (# of rows/columns) MUST be an odd number:';
    thickness = inputdlg(m1);
    thickness = str2double(thickness{1});
    while mod(thickness, 2) ~= 1 
        warndlg('WARNING: Must enter odd number for thickness.')
        thickness = inputdlg(m1);
        thickness = str2double(thickness{1});
    end
end

if maltese ==1
    % Have the user select the thickness of the row/column extraction
    m2 = 'Please enter the desired maltese cross angle (degrees):';
    angle = inputdlg(m2);
    angle = str2double(angle{1});
end



%% INDIVIDUAL SUBJECTS

% go through all the subjects
for i=1:size(file_names,1)

    % need to find the index in the lut to match the item in the folder
    for index = 1:size(LUT_data, 1)
        trigger = find(strcmp(LUT_data{index,1},(file_names{i, end})));
        if length(trigger) == 1
            break
        end
    end

     if unitStr == "Microns"
        % get the scale factor
        scale = LUT_data{index,11}; % umpp
     else
         % get the scale factor for deg
         pixelsPerDegree = LUT_data{index,12};
         scale = 1/pixelsPerDegree;
    end
    % get the eye
    eye = LUT_data{index,2};
    

    % Load in current map
    CurrentMap = load(fullfile(root_path_raw,file_names{i}));
    dummyName = fieldnames(CurrentMap);
    mat_data = CurrentMap.(dummyName{1});

    cdc_x = LUT_data{index,9};
    cdc_y = LUT_data{index,10};

    if straight

        straight_results = straight_cross(mat_data, cdc_x, cdc_y, thickness, scale, eye);
    
        %% ---- Save output ----
    
        % organize all the data into a format that can be saved to the csv
        for k = 1:numel(straight_results)/2
            % create the output file name
            outname = strcat(file_names{i},'_', eye, '_',  Meridian{k}, '_', unitStr, '_Straight_', string(datetime('now','TimeZone','local','Format','yyyyMMdd_hhmmss')), '.csv');
            outfile = fullfile(output_root, outname);
            output = [straight_results{k*2-1}, straight_results{k*2}];
            writematrix(output, outfile);
        end

    end
    if maltese
        [distance_density_mean, degree_center_points_full] = Extract_radial_densities_one(mat_data, [cdc_x, cdc_y], angle);
        % now need to save the density curves from 0, 90, 180, 270

        right = distance_density_mean{degree_center_points_full == 0};
        top = distance_density_mean{degree_center_points_full == 90};
        left = distance_density_mean{degree_center_points_full == 180};
        bottom = distance_density_mean{degree_center_points_full == 270};

        ecc_l = (0:size(left,2)-1) * scale;
        ecc_r = (0:size(right,2)-1) * scale;
        ecc_t = (0:size(top,2)-1) * scale;
        ecc_b = (0:size(bottom,2)-1) * scale;

         v_dens{1} = top;
         v_dens{2} = bottom;
         h_dens{1} = left;
         h_dens{2} = right;

        % compile results and create headers for output sheet
        % OD - left is temporal, OS - right is temporal
        % order if temporal, nasal, superior, inferior
        if strcmp(eye,'OD')
            results = {ecc_l', h_dens{1}', ecc_r', h_dens{2}', ecc_t', v_dens{1}', ecc_b', v_dens{2}'};
        else
            results = {ecc_r', h_dens{2}', ecc_l', h_dens{1}', ecc_t', v_dens{1}', ecc_b', v_dens{2}'};
        end

        
        %% ---- Save output ----
    
        % organize all the data into a format that can be saved to the csv
        for k = 1:numel(results)/2
            % create the output file name
            outname = strcat(file_names{i},'_', eye, '_',  Meridian{k}, '_', unitStr, '_Maltese_', string(datetime('now','TimeZone','local','Format','yyyyMMdd_hhmmss')), '.csv');
            outfile = fullfile(output_root, outname);
            output = [results{k*2-1}, results{k*2}];
            writematrix(output, outfile);
        end
    end

end







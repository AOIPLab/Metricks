% Row and Column Extraction Simplified
%
% Created by: Jenna Grieshop
% Date created: 6/10/2026
%
% Description: Extracts individual rows and columns at the CDC
% coordinates for each subject. 
%
% Inputs: User selects directory containing matrices in mat format
% (could be results from Stdev_Maps, matrix subtraction, etc.). NOTE: if
% the matrices are not bound density mm matrices, change the file
% identifier on line 25.
% 
% Then user selects the PCD CDC analysis summary file with OS/OD eye column 
% added after the file name column and PPD column added as the last column.
%
% Outputs: .csv files with N, T, S, I strips extracted from the density
% matrices through the CDC
% 

close all
clear all
clc

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
    % get the scale factor
    eye = LUT_data{index,2};
    

    % Load in current map
    CurrentMap = load(fullfile(root_path_raw,file_names{i}));
    dummyName = fieldnames(CurrentMap);
    mat_data = CurrentMap.(dummyName{1});

    cdc_x = LUT_data{index,9};
    cdc_y = LUT_data{index,10};


    % get row and column of the data

    h_strip = mat_data(:,cdc_x);
    v_strip = mat_data(cdc_y, :);

    % create x data arrays, reformat y data
    data.x_h = (0:length(h_strip)-1)';
    data.x_v = (0:length(v_strip)-1)';
    data.y_h = h_strip';
    data.y_v = v_strip';
    
    % initialize matrix for the converted data
    xy_h_converted = zeros(length(data.x_h), 2);
    xy_v_converted = zeros(length(data.x_v), 2);

    % scale data and store
    for j=1:length(data.x_h)
        xy_h_converted(j,1) = (data.x_h(j) - cdc_y + 1) * scale; % x values for h now in um
        xy_v_converted(j,1) = (data.x_v(j) - cdc_x + 1) * scale; %  x values for v now in um
        % y values are not scaled
        xy_h_converted(j,2) = data.y_h(j);
        xy_v_converted(j,2) = data.y_v(j);
    end

    % store data for each subject
    all_subjects_raw_h_rc{i} = xy_h_converted;
    all_subjects_raw_v_rc{i} = xy_v_converted;


%% graph individual plots
    % basic plot of the individual results
    figure(1)
    plot(xy_h_converted(:,1), xy_h_converted(:,2));
    title("Horizontal Density Through CDC Point");
    xlabel(unitStr);
    ylabel("Density");
    hold off

    figure(2)
    plot(xy_v_converted(:,1), xy_v_converted(:,2));
    title("Vertical Density Through CDC Point");
    xlabel(unitStr);
    ylabel("Density");
    hold off


    %% split data at the cdc point

    % find the 0 point in the eccentricity
    v_0_index = find(all_subjects_raw_v_rc{1,i} == 0); 
    h_0_index = find(all_subjects_raw_h_rc{1,i} == 0); 

    % split the vertical data into top and bottom
    v_ecc_t = -all_subjects_raw_v_rc{1,i}(1:v_0_index,1);
    v_dens{1} = all_subjects_raw_v_rc{1,i}(1:v_0_index,2)'; % top

    v_ecc_b = all_subjects_raw_v_rc{1,i}(v_0_index:end,1);
    v_dens{2} = all_subjects_raw_v_rc{1,i}(v_0_index:end,2)'; % bottom

    % split the horizontal data into left and right
    h_ecc_l = -all_subjects_raw_h_rc{1,i}(1:h_0_index,1);
    h_dens{1} = all_subjects_raw_h_rc{1,i}(1:h_0_index,2)'; % left
    
    h_ecc_r = all_subjects_raw_h_rc{1,i}(h_0_index:end,1);
    h_dens{2} = all_subjects_raw_h_rc{1,i}(h_0_index:end,2)'; % right



    % compile results and create headers for output sheet
    % OD - left is temporal, OS - right is temporal
    % order if temporal, nasal, superior, inferior
    if strcmp(eye,'OD')
        results = {flipud(h_ecc_l), flipud(h_dens{1}'), h_ecc_r, h_dens{2}', flipud(v_ecc_t), flipud(v_dens{1}'), v_ecc_b, v_dens{2}'};
    else
        results = {h_ecc_r, h_dens{2}', flipud(h_ecc_l), flipud(h_dens{1}'), flipud(v_ecc_t), flipud(v_dens{1}'), v_ecc_b, v_dens{2}'};
    end
   
   
    
    %% Plotting for sanity check
    figure(3)
    plot(v_ecc_t, v_dens{1});
    hold on
    plot(v_ecc_b, v_dens{2});
    plot(h_ecc_l, h_dens{1});
    plot(h_ecc_r, h_dens{2});
    title("Density Through CDC Point");
    xlabel(unitStr);
    ylabel("Density");
    hold off


    %% ---- Save output ----

    % organize all the data into a format that can be saved to the csv
    for k = 1:numel(results)/2
        % D(1:numel(results{k}),k)=num2cell(results{k});
        % create the output file name
        outname = strcat(file_names{i},'_', eye, '_',  Meridian{k}, '_', unitStr, '_', string(datetime('now','TimeZone','local','Format','yyyyMMdd_hhmmss')), '.csv');
        outfile = fullfile(output_root, outname);
        output = [results{k*2-1}, results{k*2}];
        % results = {h_ecc_r, h_dens{2}', h_ecc_l, h_dens{1}', v_ecc_t, v_dens{1}', v_ecc_b, v_dens{2}'};
        writematrix(output, outfile);
    end



end

%% save raw data for all subjects

% fname_h_all = fullfile(pathname_cdc, 'all_h_rc_data.mat');
% fname_v_all = fullfile(pathname_cdc, 'all_v_rc_data.mat');
% save(fname_h_all, 'all_subjects_raw_h_rc');
% save(fname_v_all, 'all_subjects_raw_v_rc');






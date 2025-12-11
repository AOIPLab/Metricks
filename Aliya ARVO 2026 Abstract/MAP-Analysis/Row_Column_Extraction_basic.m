% Row and Column Extraction 
%
% Created by: Jenna Grieshop
% Date created: 1/26/2024
%
% Description: Extracts individual rows and columns at the master CDC
% coordinates for each subject. If the master CDC location is not an integer
% the extracted information is a weighted average of the two columns (or rows).
% Extracted row and column is binned and averaged for each subject.
% Additionally, this data is then binned and averaged again across all
% subjects. Only overlapping data across all subjects is included.
% 
% Note: The files in the folder must match exacly 
% what is in the master CDC list and be in the correct order.
%
% Inputs: User selects directory containing matrices in csv format
% (could be results from Stdev_Maps, matrix subtraction, etc.). User enters
% a unique identifier of the files to be used (MINUS, coeffvar, etc).
% Then
% user selects the PCD CDC analysis summary file with OS/OD eye column added after the subject ID column. 
%
% Outputs: .mat files with the individual data
% 

clear all
clc

basePath = which('Row_Column_Extraction_basic.m');

[basePath ] = fileparts(basePath);
path(path,fullfile(basePath,'lib')); % Add our support library to the path.

% user to select directory with data
root_path_raw = uigetdir('.','Select directory containing raw csv analyses');

identifier = inputdlg("please enter the unique identifier of the files to use");

% get all the file paths that we are interested in in each of the folders
[file_paths_raw] = read_folder_contents_rec(root_path_raw, 'csv', identifier{1}); % used for stdev map
% [file_paths_raw] = read_folder_contents_rec(root_path_raw, 'csv','coeffvar'); % used for CoV

% select cdc analysis file
[filename_cdc, pathname_cdc] = uigetfile('*.csv', 'Please select the pcd cdc analysis summary file', 'MultiSelect', 'off');


% load in cdc LUT data and extract all the subject IDs
LUT_data = readtable(fullfile(pathname_cdc, filename_cdc));
file_names = split(file_paths_raw, '\');
% sub_id = cdc_data(:,1);
% spl = split(sub_id{:,1}, "_");
% subjectID = spl(:,1);


%% INDIVIDUAL SUBJECTS

% go through all the subjects
for i=1:size(file_paths_raw,1)

    % need to find the index in the lut to match the item in the folder
    for index = 1:size(LUT_data, 1)
        trigger = find(strcmp(LUT_data{index,1},(file_names{i, end})));
        if length(trigger) == 1
            break
        end
    end

    % get the scale factor
    scale = LUT_data{index,11}; % umpp
    eye = LUT_data{index,2};

    % load in the csv data and master cdc coords
    csv_data = load(file_paths_raw{i});
    cdc_x = LUT_data{index,9};
    cdc_y = LUT_data{index,10};


    % get row and column of the data

    h_strip = csv_data(:,cdc_x);
    v_strip = csv_data(cdc_y, :);

    % create x data arrays, reformat y data
    data.x_h = (0:length(h_strip)-1)';
    data.x_v = (0:length(v_strip)-1)';
    data.y_h = h_strip';
    data.y_v = v_strip';
    
    % initialize matrix for the converted data
    xy_h_converted = zeros(length(data.x_h), 2);
    xy_v_converted = zeros(length(data.x_v), 2);

    % convert the data to be all in the same orientation.
    for j=1:length(data.x_h)
        % check what eye - if left flip sign so that temporal is negative -
        % need eye from the lut file to do this
        if strcmp(eye,'OD')
            xy_h_converted(j,1) = (data.x_h(j) - cdc_y) * scale; % x values for h now in um
        else
            xy_h_converted(j,1) = -(data.x_h(j) - cdc_y) * scale; % x values for h now in um
        end
        xy_v_converted(j,1) = (data.x_v(j) - cdc_x) * scale; %  x values for v now in um
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
    title("Horizontal Stdev Through CDC Point");
    xlabel("Microns");
    ylabel("Density");
    hold on

    figure(2)
    plot(xy_v_converted(:,1), xy_v_converted(:,2));
    title("Vertical Stdev Through CDC Point");
    xlabel("Microns");
    ylabel("Density");
    hold on


end

%% save raw data for all subjects

fname_h_all = fullfile(pathname_cdc, 'all_h_rc_data.mat');
fname_v_all = fullfile(pathname_cdc, 'all_v_rc_data.mat');
save(fname_h_all, 'all_subjects_raw_h_rc');
save(fname_v_all, 'all_subjects_raw_v_rc');






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

identifier = 'bounddensity_matrix';

% get all the file paths that we are interested in in each of the folders
[file_paths_raw] = read_folder_contents_rec(root_path_raw, 'csv', identifier); % used for stdev map
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
    % figure(1)
    % plot(xy_h_converted(:,1), xy_h_converted(:,2));
    % title("Horizontal Stdev Through CDC Point");
    % xlabel("Microns");
    % ylabel("Density");
    % hold off
    % 
    % figure(2)
    % plot(xy_v_converted(:,1), xy_v_converted(:,2));
    % title("Vertical Stdev Through CDC Point");
    % xlabel("Microns");
    % ylabel("Density");
    % hold off


    %% split data at the cdc point

    % find the 0 point in the eccentricity
    v_0_index = find(all_subjects_raw_v_rc{1,1} == 0); 
    h_0_index = find(all_subjects_raw_h_rc{1,1} == 0); 

    % split the vertical data into top and bottom
    v_ecc_t = -flipud(all_subjects_raw_v_rc{1,1}(1:v_0_index,1));
    v_dens{1} = flipud(all_subjects_raw_v_rc{1,1}(1:v_0_index,2))'; % top
    v_ecc_b = all_subjects_raw_v_rc{1,1}(v_0_index:end,1);
    v_dens{2} = all_subjects_raw_v_rc{1,1}(v_0_index:end,2)'; % bottom
    % find minimum eccentricity for the average
    if size(v_ecc_b, 1) > size(v_ecc_t,1)
        v_ecc_avg = v_ecc_t;
    else
        v_ecc_avg = v_ecc_b;
    end

    % split the horizontal data into left and right
    h_ecc_l = flipud(all_subjects_raw_h_rc{1,1}(1:h_0_index,1));
    h_dens{1} = all_subjects_raw_h_rc{1,1}(1:h_0_index,2)'; % left
    h_ecc_r = -all_subjects_raw_h_rc{1,1}(h_0_index:end,1);
    h_dens{2} = all_subjects_raw_h_rc{1,1}(h_0_index:end,2)'; % right
    % find minimum eccentricity for the average
    if size(h_ecc_r, 1) > size(h_ecc_l,1)
        h_ecc_avg = h_ecc_l;
    else
        h_ecc_avg = h_ecc_r;
    end


    % prepare for vertical averaging
    maxNumCol_v = max(cellfun(@(c) size(c,2), v_dens));  % max number of columns
    combined_v = cell2mat(cellfun(@(c){padarray(c,[0,maxNumCol_v-size(c,2)],NaN,'Post')}, v_dens)');

    % prepare for horizontal averaging
    maxNumCol_h = max(cellfun(@(c) size(c,2), h_dens));  % max number of columns
    combined_h = cell2mat(cellfun(@(c){padarray(c,[0,maxNumCol_h-size(c,2)],NaN,'Post')}, h_dens)');

    % agerage the vertical and horizontal densities
    avg_v_dens = rmmissing(mean(combined_v, 1));
    avg_h_dens = rmmissing(mean(combined_h, 1));


    % compile results and create headers for output sheet
    results = {h_ecc_l, h_dens{1}', h_ecc_r, h_dens{2}', v_ecc_t, v_dens{1}', v_ecc_b, v_dens{2}', h_ecc_avg, avg_h_dens', v_ecc_avg, avg_v_dens'};
    headers = {'Ecc_Horz_L_um', 'Density_Horz_L', 'Ecc_Horz_R_um', 'Density_Horz_R', 'Ecc_Vert_T_um,', 'Density_Vert_T', 'Ecc_Vert_B_um', 'Density_Vert_B', 'Ecc_Avg_Horz_um', 'Density_Avg_Horz', 'Ecc_Avg_Vert_um', 'Density_Avg_Vert'};

    %% Plotting for sanity check
    plot(v_ecc_t, v_dens{1});
    hold on
    plot(v_ecc_b, v_dens{2});
    plot(h_ecc_l, h_dens{1});
    plot(h_ecc_r, h_dens{2});
    
    % 
    % plot(v_ecc_avg, avg_v_dens);
    % plot(h_ecc_avg, avg_h_dens);

    %% ---- Save output ----

    % create the output file name
    outname = strcat(subj,'_', eye, '_HorizVertDensity_', string(datetime('now','TimeZone','local','Format','yyyyMMdd_hhmmss')), '.csv');
    outfile = fullfile(output_root, outname);

    % organize all the data into a format that can be saved to the csv
    D = {};
    for k = 1:numel(results)
      D(1:numel(results{k}),k)=num2cell(results{k});
    end
    combined_output = [headers; D];
    writecell(combined_output, outfile);


end

%% save raw data for all subjects

fname_h_all = fullfile(pathname_cdc, 'all_h_rc_data.mat');
fname_v_all = fullfile(pathname_cdc, 'all_v_rc_data.mat');
save(fname_h_all, 'all_subjects_raw_h_rc');
save(fname_v_all, 'all_subjects_raw_v_rc');






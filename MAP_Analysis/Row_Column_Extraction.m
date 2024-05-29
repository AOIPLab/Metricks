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
% Then user selects the master cdc list (result from stdev_Maps.m), then
% user selects the PCD CDC analysis summary file. The user is then asked to
% input the spacing and window in microns (required for binning).
%
% Outputs: Horizontal and vertical bin analysis csvs for each subject, all 
% raw data for horizontal and vertical mat file for all subjects,
% Horizontal and Vertical graphs as svgs for data between subjects
% 

clear all
clc
addpath('lib');

% user to select directory with data
root_path_raw = uigetdir('.','Select directory containing raw csv analyses');

identifier = inputdlg("please enter the unique identifier of the files to use");

% get all the file paths that we are interested in in each of the folders
[file_paths_raw] = read_folder_contents_rec(root_path_raw, 'csv', identifier{1}); % used for stdev map
% [file_paths_raw] = read_folder_contents_rec(root_path_raw, 'csv','coeffvar'); % used for CoV

% select master cdc file
[filename_master_cdc, pathname_master_cdc] = uigetfile('*.xlsx', 'Please select the master cdc list', 'MultiSelect', 'off');

% select cdc analysis file
[filename_cdc, pathname_cdc] = uigetfile('*.csv', 'Please select the pcd cdc analysis summary file', 'MultiSelect', 'off');


% load in cdc LUT data and extract all the subject IDs
LUT_data = readtable(fullfile(pathname_cdc, filename_cdc));
% sub_id = cdc_data(:,1);
% spl = split(sub_id{:,1}, "_");
% subjectID = spl(:,1);

master_cdc = readtable(fullfile(pathname_master_cdc, filename_master_cdc));


% set up messages to be displayed to user to set spacing and window in the desired lateral unit
m1 = 'Please enter desired SPACING (um)';
m2 = 'Please enter the desired WINDOW (um)';

% user sets spacing and window
spacing = inputdlg(m1);
spacing = str2double(spacing{1});
window = inputdlg(m2);
window = str2double(window{1});

% spacing must be >= window for no overlap
% if not loop until it is or the user can override and use the overlap
while spacing < (window)
    contin = questdlg('WARNING: Window overlap with current spacing and window selections Continue with current selection?', ...
      '', ...
      'YES', 'NO', 'NO');
  if strcmpi(contin, 'YES')
    break;
  else
    % user sets spacing and window
    spacing = inputdlg(m1);
    spacing = str2double(spacing{1});
    window = inputdlg(m2);
    window = str2double(window{1});
  end
end

%% INDIVIDUAL SUBJECTS

% go through all the subjects
for i=1:size(file_paths_raw,1)

    % get the scale factor
    scale = LUT_data{i,10}; % umpp
    eye = LUT_data{i,2};

    % load in the csv data and master cdc coords
    csv_data = load(file_paths_raw{i});
    master_cdc_x = master_cdc{i,2};
    master_cdc_y = master_cdc{i,3};

    % check if the cdc coords are integers
    integerTest_x =~ mod(master_cdc_x,1);
    integerTest_y =~ mod(master_cdc_y,1);

    % if integers just get that row, if not do weighted average the two rows it is in
    % between.
    if integerTest_x
        h_strip = csv_data(:,master_cdc_x);
    else
        down_weight_h = ceil(master_cdc_x) - master_cdc_x;
        up_weight_h = 1-down_weight_h;
        down_h = csv_data(:,floor(master_cdc_x)) * down_weight_h;
        up_h = csv_data(:,ceil(master_cdc_x)) * up_weight_h;
        h_strip = up_h + down_h;
    end

    % if integers just get that row, if not do weighted average the two rows it is in
    % between.
    if integerTest_y
        v_strip = csv_data(master_cdc_y, :);
    else
        down_weight_v = ceil(master_cdc_y) - master_cdc_y;
        up_weight_v = 1-down_weight_v;
        down_v = csv_data(:,floor(master_cdc_y)) * down_weight_v;
        up_v = csv_data(:,ceil(master_cdc_y)) * up_weight_v;
        v_strip = up_v + down_v;
    end

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
            xy_h_converted(j,1) = (data.x_h(j) - master_cdc_y) * scale; % x values for h now in um
        else
            xy_h_converted(j,1) = -(data.x_h(j) - master_cdc_y) * scale; % x values for h now in um
        end
        xy_v_converted(j,1) = (data.x_v(j) - master_cdc_x) * scale; %  x values for v now in um
        % y values are not scaled
        xy_h_converted(j,2) = data.y_h(j);
        xy_v_converted(j,2) = data.y_v(j);
    end

    % store data for each subject
    all_subjects_raw_h_rc{i} = xy_h_converted;
    all_subjects_raw_v_rc{i} = xy_v_converted;


    % find max x values of the row and column
    x_h_max = max(xy_h_converted(:,1));
    x_v_max = max(xy_v_converted(:,1));
    x_h_min = min(xy_h_converted(:,1));
    x_v_min = min(xy_v_converted(:,1));

    % find all the bin center values based on maximums and spacing
    bin_centers_left_h = (0:-spacing:x_h_min);
    bin_centers_right_h = (0:spacing: x_h_max);
    bin_centers_left_v = (0:-spacing:x_v_min);
    bin_centers_right_v = (0:spacing: x_v_max);

    bin_centers_right_h(1) = []; % get rid of second zero before combining
    bin_centers_right_v(1) = []; % get rid of second zero before combining
    bin_centers_h = [flip(bin_centers_left_h), bin_centers_right_h]; % combine
    bin_centers_v = [flip(bin_centers_left_v), bin_centers_right_v]; % combine

%% horizontal
    % initialize values and arrays
    sum_h = 0;
    count_h = 0;
    averages_h = zeros(length(bin_centers_h), 1); 
    items_in_bin_h = zeros(length(bin_centers_h),1); 

    % loops to go through bin centers and x values to determine if the data
    % point belongs in each bin
    for bin = bin_centers_h 
        for data_point = xy_h_converted(:,1)' 
            if (data_point >= (bin - (window/2))) && (data_point < (bin + (window/2))) % check if the x value is within the bin range
                value_index = find(xy_h_converted(:,1)==data_point);
                sum_h = sum_h + xy_h_converted(value_index,2); % add y value to the sum
                count_h = count_h + 1;
            end
        end
        bin_index = find(bin_centers_h==bin);
        averages_h(bin_index) = sum_h/count_h;
        items_in_bin_h(bin_index) = count_h;
        sum_h = 0; % reset value
        count_h = 0; % reset value
    end

    
%% vertical
    % initialize values and arrays
    sum_v = 0;
    count_v = 0;
    averages_v = zeros(length(bin_centers_v), 1); 
    items_in_bin_v = zeros(length(bin_centers_v),1); 


    % loops to go through bin centers and x values to determine if the data
    % point belongs in each bin
    for bin = bin_centers_v 
        for data_point = xy_v_converted(:,1)' 
            if (data_point >= (bin - (window/2))) && (data_point < (bin + (window/2))) % check if the x value is within the bin range
                value_index = find(xy_v_converted(:,1)==data_point);
                sum_v = sum_v + xy_v_converted(value_index,2); % add y value to the sum
                count_v = count_v + 1;
            end
        end
        bin_index = find(bin_centers_v==bin);
        averages_v(bin_index) = sum_v/count_v;
        items_in_bin_v(bin_index) = count_v;
        sum_v = 0; % reset value
        count_v = 0; % reset value
    end

%% graph individual plots
    % basic plot of the individual results
    figure(1)
    plot(xy_h_converted(:,1), xy_h_converted(:,2));
    title("Horizontal Stdev Through CDC Point");
    xlabel("Microns");
    ylabel("Standard Deviation");
    hold on

    figure(2)
    plot(xy_v_converted(:,1), xy_v_converted(:,2));
    title("Vertical Stdev Through CDC Point");
    xlabel("Microns");
    ylabel("Standard Deviation");
    hold on

%% format output and save to file
    % horizontal
    output_fname_h = strcat(master_cdc{i,1}, '_Horizontal_Bin_Analysis_', string(datetime('now','TimeZone','local','Format','yyyyMMdd')), '.csv');

    % setting up for table creation
    BinCenter_h = num2cell(bin_centers_h');
    Average_h = num2cell(averages_h);
    ItemsInBin_h = num2cell(items_in_bin_h);
    % Table creation
    T_h = table(BinCenter_h, Average_h, ItemsInBin_h);
    
    % write output file
    writetable(T_h, fullfile(pathname_cdc,output_fname_h));
    all_h_data{i} = {bin_centers_h', averages_h}; 


    % vertical
    output_fname_v = strcat(master_cdc{i,1}, '_Vertical_Bin_Analysis_', string(datetime('now','TimeZone','local','Format','yyyyMMdd')), '.csv');

    % setting up for table creation
    BinCenter_v = num2cell(bin_centers_v');
    Average_v = num2cell(averages_v);
    ItemsInBin_v = num2cell(items_in_bin_v);
    % Table creation
    T_v = table(BinCenter_v, Average_v, ItemsInBin_v);
    
    % write output file
    writetable(T_v, fullfile(pathname_cdc,output_fname_v));
    all_v_data{i} = {bin_centers_v', averages_v};
 

end

%% save raw data for all subjects

fname_h_all = fullfile(pathname_cdc, 'all_h_rc_data.mat');
fname_v_all = fullfile(pathname_cdc, 'all_v_rc_data.mat');
save(fname_h_all, 'all_subjects_raw_h_rc');
save(fname_v_all, 'all_subjects_raw_v_rc');



%% ACROSS SUBJECTS

% Average Values in bins 

% Figure out the minimum and maximum bin location that overlaps for each
% subject - these are the new min and max bins so only overlapping data
% ends up in the graph/data comparing all subjects
curr_min_h_bin_center = -1000000;
curr_min_v_bin_center = -1000000;
curr_max_h_bin_center = 1000000;
curr_max_v_bin_center = 1000000;

for j=1:size(all_h_data,2)
    new_min_h_bin_center = min(min(all_h_data{1,j}{1,1}));    
    if new_min_h_bin_center > curr_min_h_bin_center
        curr_min_h_bin_center = new_min_h_bin_center;
    end
    new_min_v_bin_center = min(min(all_v_data{1,j}{1,1}));    
    if new_min_v_bin_center > curr_min_v_bin_center
        curr_min_v_bin_center = new_min_v_bin_center;
    end
    new_max_h_bin_center = max(max(all_h_data{1,j}{1,1}));    
    if new_max_h_bin_center < curr_max_h_bin_center
        curr_max_h_bin_center = new_max_h_bin_center;
    end
    new_max_v_bin_center = max(max(all_v_data{1,j}{1,1}));    
    if new_max_v_bin_center < curr_max_v_bin_center
        curr_max_v_bin_center = new_max_v_bin_center;
    end

end

count_h = 1;
sum_h_bin = zeros(((curr_max_h_bin_center-curr_min_h_bin_center)/window) + 1,1);
clear all_h_data_in_bin

%% horizontal
% go through all the subjects
for subj=1:size(all_h_data,2)
    % go through all the bin centers
    for bin=1:size((all_h_data{1,subj}{1,1}))

        % if the bin for the subject is not within the minimum and maximum of
        % for all subjects skip it, otherwise add it to the bin sum
        if all_h_data{1,subj}{1,1}(bin) < curr_min_h_bin_center || all_h_data{1,subj}{1,1}(bin) > curr_max_h_bin_center
            continue
        else
            sum_h_bin(count_h) = all_h_data{1,subj}{1,2}(bin) + sum_h_bin(count_h);
            all_h_data_in_bin(subj,count_h) = all_h_data{1,subj}{1,2}(bin);
            count_h = count_h + 1;
        end
    end
    count_h = 1;
    
end

% get average in the bins and standard deviation
avg_h_bin = sum_h_bin(:)/size(all_h_data,2);
stdev_h_bin = std(all_h_data_in_bin,0,1); % standard deviation of columns (bins)

% get +- 2 standard deviations to add to the graph
plus_stdev_h_bin = avg_h_bin + (stdev_h_bin' * 2);
minus_stdev_h_bin = avg_h_bin - (stdev_h_bin' * 2);

% graph the binned results across subjects
x_h_bin = (curr_min_h_bin_center:window:curr_max_h_bin_center)';
f = figure(3);
plot(x_h_bin, avg_h_bin);
hold on 
plot(x_h_bin, plus_stdev_h_bin, ':b');
plot(x_h_bin, minus_stdev_h_bin, ':b');
hold off
title("Average Horizontal Result Through CDC Point");
xlabel("Microns");
ylabel("Result");

% save graph
output_fname_horz_graph = strcat('Horizontal_Graph_', string(datetime('now','TimeZone','local','Format','yyyyMMdd')), '.svg');
print(f, '-dsvg', fullfile(pathname_cdc,output_fname_horz_graph));

%% vertical

count_v = 1;
sum_v_bin = zeros(((curr_max_v_bin_center-curr_min_v_bin_center)/window) + 1,1);
clear all_v_data_in_bin

% go through all the subjects
for subj=1:size(all_v_data,2)
    % all_v_data_in_bin = zeros(size((all_v_data{1,k}{1,1}),1),1);
    % go through all the bin centers
    for bin=1:size((all_v_data{1,subj}{1,1}))

        % if the bin for the subject is not within the minimum and maximum of
        % for all subjects skip it, otherwise add it to the bin sum
        if all_v_data{1,subj}{1,1}(bin) < curr_min_v_bin_center || all_v_data{1,subj}{1,1}(bin) > curr_max_v_bin_center
            continue
        else
            sum_v_bin(count_v) = all_v_data{1,subj}{1,2}(bin) + sum_v_bin(count_v);
            all_v_data_in_bin(subj,count_v) = all_v_data{1,subj}{1,2}(bin);
            count_v = count_v + 1;
        end
    end
    count_v = 1;
    
end

% get average in the bins and standard deviation
avg_v_bin = sum_v_bin(:)/size(all_v_data,2);
stdev_v_bin = std(all_v_data_in_bin,0,1);

% get +- 2 standard deviations to add to the graph
plus_stdev_v_bin = avg_v_bin + (stdev_v_bin' * 2);
minus_stdev_v_bin = avg_v_bin - (stdev_v_bin' * 2);

% graph the binned results across subjects
x_v_bin = (curr_min_v_bin_center:window:curr_max_v_bin_center)';
g = figure(4);
plot(x_v_bin, avg_v_bin);
hold on 
plot(x_v_bin, plus_stdev_v_bin, ':b')
plot(x_v_bin, minus_stdev_v_bin, ':b')
title("Average Vertical Result Through CDC Point");
xlabel("Microns");
ylabel("Result");

% save graph
output_fname_vert_graph = strcat('Vertical_Graph_', string(datetime('now','TimeZone','local','Format','yyyyMMdd')), '.svg');
print(g, '-dsvg', fullfile(pathname_cdc,output_fname_vert_graph));








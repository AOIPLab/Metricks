% Loads in the raw extracted rows from avg and stdev and gets the cov
% The files in the folder must match exactly what is in the
% master cdc list in the correct order.
% 2/27/2024
% Jenna Grieshop

clear all
clc
addpath('lib');

% get path that the data in is
root_path_stdev = uigetdir('.','Select directory containing raw stdev csv results');
root_path_avg = uigetdir('.','Select directory containing raw avg csv results');

% find files for average h and v
[root_path_avg_h] = read_folder_contents_rec(root_path_avg, 'mat', '_h');
[root_path_avg_v] = read_folder_contents_rec(root_path_avg, 'mat', '_v');

% fine files for stdev h and v
[root_path_stdev_h] = read_folder_contents_rec(root_path_stdev, 'mat', '_h');
[root_path_stdev_v] = read_folder_contents_rec(root_path_stdev, 'mat', '_v');


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



% load in data
avg_data_h = load(root_path_avg_h{1});
avg_data_v = load(root_path_avg_v{1});

stdev_data_h = load(root_path_stdev_h{1});
stdev_data_v = load(root_path_stdev_v{1});

for i=1:44
    xy_h_converted = zeros(length(stdev_data_h.all_subjects_raw_h_rc_stdev{1,i}(:,1)), 2);
    xy_v_converted = zeros(length(stdev_data_v.all_subjects_raw_v_rc_stdev{1,i}(:,1)), 2);

    xy_h_converted(:,1) = stdev_data_h.all_subjects_raw_h_rc_stdev{1,i}(:,1);
    xy_v_converted(:,1) = stdev_data_v.all_subjects_raw_v_rc_stdev{1,i}(:,1);
    xy_h_converted(:,2) = stdev_data_h.all_subjects_raw_h_rc_stdev{1,i}(:,2)./avg_data_h.all_subjects_raw_h_rc_avg{1,i}(:,2);
    xy_v_converted(:,2) = stdev_data_v.all_subjects_raw_v_rc_stdev{1,i}(:,2)./avg_data_v.all_subjects_raw_v_rc_avg{1,i}(:,2);

    CoV_h{i} = stdev_data_h.all_subjects_raw_h_rc_stdev{1,i}(:,2)./avg_data_h.all_subjects_raw_h_rc_avg{1,i}(:,2);
    CoV_v{i} = stdev_data_v.all_subjects_raw_v_rc_stdev{1,i}(:,2)./avg_data_v.all_subjects_raw_v_rc_avg{1,i}(:,2);

    %  % basic plot of the individual results
    % figure(1)
    % plot(stdev_data_h.all_subjects_raw_h_rc_stdev{1,i}(:,1),  CoV_h{i});
    % title("Horizontal CoV Through CDC Point");
    % xlabel("Microns");
    % ylabel("CoV");
    % hold on
    % 
    % figure(2)
    % plot(stdev_data_h.all_subjects_raw_h_rc_stdev{1,i}(:,1),CoV_v{i});
    % title("Vertical CoV Through CDC Point");
    % xlabel("Microns");
    % ylabel("CoV");
    % hold on

    % get list of bin centers
    x_h_max = max(xy_h_converted(:,1));
    x_v_max = max(xy_v_converted(:,1));
    x_h_min = min(xy_h_converted(:,1));
    x_v_min = min(xy_v_converted(:,1));

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

    % loops to go through bin centers and x values
    for m = bin_centers_h 
        for n = xy_h_converted(:,1)' 
            if (n >= (m - (window/2))) && (n < (m + (window/2))) % check if the x value is within the bin range
                value_index = find(xy_h_converted(:,1)==n);
                sum_h = sum_h + xy_h_converted(value_index,2); % add y value to the sum
                count_h = count_h + 1;
            end
        end
        bin_index = find(bin_centers_h==m);
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


    % loops to go through bin centers and x values
    for m = bin_centers_v 
        for n = xy_v_converted(:,1)' 
            if (n >= (m - (window/2))) && (n < (m + (window/2))) % check if the x value is within the bin range
                value_index = find(xy_v_converted(:,1)==n);
                sum_v = sum_v + xy_v_converted(value_index,2); % add y value to the sum
                count_v = count_v + 1;
            end
        end
        bin_index = find(bin_centers_v==m);
        averages_v(bin_index) = sum_v/count_v;
        items_in_bin_v(bin_index) = count_v;
        sum_v = 0; % reset value
        count_v = 0; % reset value
    end

%% format output and save to file
    % output_fname_h = strcat(num2str(master_cdc{i,1}), '_Horizontal_Bin_Analysis_', string(datetime('now','TimeZone','local','Format','yyyyMMdd')), '.csv');

    % setting up for table creation
    BinCenter_h = num2cell(bin_centers_h');
    Average_h = num2cell(averages_h);
    ItemsInBin_h = num2cell(items_in_bin_h);
    % Table creation
    T_h = table(BinCenter_h, Average_h, ItemsInBin_h);
    
    
    % write output file
    % writetable(T_h, fullfile(pathname_cdc,output_fname_h));


    all_h_data{i} = {bin_centers_h', averages_h}; 


    % output_fname_v = strcat(num2str(master_cdc{i,1}), '_Vertical_Bin_Analysis_', string(datetime('now','TimeZone','local','Format','yyyyMMdd')), '.csv');

    % setting up for table creation
    BinCenter_v = num2cell(bin_centers_v');
    Average_v = num2cell(averages_v);
    ItemsInBin_v = num2cell(items_in_bin_v);
    % Table creation
    T_v = table(BinCenter_v, Average_v, ItemsInBin_v);
    
    % write output file
    % writetable(T_v, fullfile(pathname_cdc,output_fname_v));

    all_v_data{i} = {bin_centers_v', averages_v};

end


%% Average Values in bins across subjects

% get list of min and max for bin locations across all subjects
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

%% horizontal
% go through all the subjects
for k=1:size(all_h_data,2)
    % go through all the bin centers
    for m=1:size((all_h_data{1,k}{1,1}))

        if all_h_data{1,k}{1,1}(m) < curr_min_h_bin_center || all_h_data{1,k}{1,1}(m) > curr_max_h_bin_center
            continue
        else
            sum_h_bin(count_h) = all_h_data{1,k}{1,2}(m) + sum_h_bin(count_h);
            count_h = count_h + 1;
        end
    end
    count_h = 1;
    
end

avg_h_bin = sum_h_bin(:)/size(all_h_data,2);
stdev_h_bin = std(avg_h_bin);

plus_stdev_h_bin = avg_h_bin + (stdev_h_bin * 2);
minus_stdev_h_bin = avg_h_bin - (stdev_h_bin * 2);

x_h_bin = (curr_min_h_bin_center:window:curr_max_h_bin_center)';
f = figure(3);
plot(x_h_bin, avg_h_bin);
hold on 
plot(x_h_bin, plus_stdev_h_bin, ':b');
plot(x_h_bin, minus_stdev_h_bin, ':b');
hold off
title("Average Horizontal CoV Through CDC Point");
xlabel("Microns");
ylabel("CoV");

output_fname_horz_graph = strcat('Horizontal_CoV_Graph_', string(datetime('now','TimeZone','local','Format','yyyyMMdd')), '.svg');
print(f, '-dsvg', fullfile(root_path_stdev,output_fname_horz_graph));

%% vertical

count_v = 1;
sum_v_bin = zeros(((curr_max_v_bin_center-curr_min_v_bin_center)/window) + 1,1);


% go through all the subjects
for k=1:size(all_v_data,2)
    % go through all the bin centers
    for m=1:size((all_v_data{1,k}{1,1}))

        if all_v_data{1,k}{1,1}(m) < curr_min_v_bin_center || all_v_data{1,k}{1,1}(m) > curr_max_v_bin_center
            continue
        else
            sum_v_bin(count_v) = all_v_data{1,k}{1,2}(m) + sum_v_bin(count_v);
            count_v = count_v + 1;
        end
    end
    count_v = 1;
    
end

avg_v_bin = sum_v_bin(:)/size(all_v_data,2);
stdev_v_bin = std(avg_v_bin);

plus_stdev_v_bin = avg_v_bin + (stdev_v_bin * 2);
minus_stdev_v_bin = avg_v_bin - (stdev_v_bin * 2);

x_v_bin = (curr_min_v_bin_center:window:curr_max_v_bin_center)';
g = figure(4);
plot(x_v_bin, avg_v_bin);
hold on 
plot(x_v_bin, plus_stdev_v_bin, ':b')
plot(x_v_bin, minus_stdev_v_bin, ':b')
title("Average Vertical CoV Through CDC Point");
xlabel("Microns");
ylabel("CoV");

output_fname_vert_graph = strcat('Vertical_CoV_Graph_', string(datetime('now','TimeZone','local','Format','yyyyMMdd')), '.svg');
print(g, '-dsvg', fullfile(root_path_stdev,output_fname_vert_graph));






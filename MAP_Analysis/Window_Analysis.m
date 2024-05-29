% Window results analysis
%
% Created by: Jenna Grieshop
% Date created: 1/24/2024
% 
% Description: Script to be used as a sanity check againsts Metrics MAP
% script. This takes in the raw data window analysis matfiles and 
% independently finds the min, max, average, and range of the selected
% metric across the windows. Used originally to find the actual ammount of cones that were
% being used for each sampling window.
%
% Input: User selects the directory containing the window results matfiles
% from Metricks. User chooses which metric to analyze.
% 
% Output: A csv that contains the min, max, mean, and range of the metric
% for each subject ID/window results file 


% user prompted to select directory with the window analysis matfiles
root_path = uigetdir('.','Select directory containing analyses');
root_dir = dir(root_path);
root_dir = struct2cell(root_dir)';

% user prompted to select the metric they want to check
liststr = {'bound_area','unbound_area','bound_num_cells', 'unbound_num_cells'};
[selectedmap, oked] = listdlg('PromptString','Select map type:',...
                              'SelectionMode','single',...
                              'ListString',liststr);
if oked == 0
    error('Cancelled by user.');
end

selectedmap = liststr{selectedmap};    

% looks for all the window results
win_results_dir = root_dir(...
    ~cellfun(@isempty, strfind(root_dir(:,1), 'window_results')),:);

comptable = [];

% loop goes through all of the window results files
for i=1:size(win_results_dir,1)

    % load in data and extract subject ID
    data = load(fullfile(win_results_dir{i,2}, win_results_dir{i,1}));
    subject_id = win_results_dir{i,1}(1:8);

    % separate the data the user is interested in
     if selectedmap == "bound_area"
        selected_data = data.win_res.bound_area;
    elseif selectedmap == "unbound_area"
        selected_data = data.win_res.unbound_area;
    elseif selectedmap == "bound_num_cells"
        selected_data = data.win_res.bound_num_cells;
    elseif selectedmap == "unbound_num_cells"
        selected_data = data.win_res.unbound_num_cells;
    else
        disp("something is wrong");
     end

    % compile outputs into a table
    output(:,1) = {subject_id};
    output(:,2) = {max(selected_data)};
    output(:,3) = {min(selected_data)};
    output(:,4) = {mean(selected_data)};
    output(:,5) = {max(selected_data)-min(selected_data)};

    comptable = [comptable;output];

end

% save the results
header = {'Subject_ID','Max', 'Min', 'Mean', 'Range'};
finaloutput = [header;comptable];
newname = [selectedmap, '_Window_Analysis_', datestr(now, 'dd_mmm_yyyy'), '.csv'];
writecell(finaloutput,fullfile(root_path,newname));


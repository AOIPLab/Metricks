% Window results analysis
% 1/24/2024
% Jenna Grieshop

root_path = uigetdir('.','Select directory containing analyses');
root_dir = dir(root_path);
root_dir = struct2cell(root_dir)';

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

for i=1:size(win_results_dir,1)

    data = load(fullfile(win_results_dir{i,2}, win_results_dir{i,1}));
    subject_id = win_results_dir{i,1}(1:8);

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

    output(:,1) = {subject_id};;
    output(:,2) = {max(selected_data)};
    output(:,3) = {min(selected_data)};
    output(:,4) = {mean(selected_data)};
    output(:,5) = {max(selected_data)-min(selected_data)};

    comptable = [comptable;output];

end

header = {'Subject_ID','Max', 'Min', 'Mean', 'Range'};
finaloutput = [header;comptable];
newname = [selectedmap, '_Window_Analysis_', date, '.csv'];
writecell(finaloutput,fullfile(root_path,newname));


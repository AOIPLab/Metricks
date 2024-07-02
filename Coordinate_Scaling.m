

clear all
close all
clc


basePath = which('Coordinate_Scaling.m');
[basePath] = fileparts(basePath);
path(path,fullfile(basePath,'lib')); 

% User selects folder with data
dataPath = uigetdir('.','Select directory containing analyses');

%User selects results folder
resultPath = uigetdir('.','Select directory to save new coordinates in');

% Read in csv names and then have user select the LUT
[fnameList] = read_folder_contents(dataPath,'csv');
[scalingFname, scalingPath] = uigetfile(fullfile(dataPath,'*.csv'),'Select scaling LUT.');

% Remove LUT file from fnameList
fnameList(ismember(fnameList,scalingFname))=[];

% load in the LUT
[~, lutData] = load_scaling_file_2(fullfile(scalingPath,scalingFname));


for i=1:size(fnameList,1) % Go through all files in list 
    
    % Match LUT entry with file                                
    LUTindex=find( cellfun(@(s) ~isempty(strfind(fnameList{i},s )), lutData{1} ) );

    % Extract info from LUT
    og_center = lutData{2}(LUTindex);
    scaling_factor = lutData{3}(LUTindex);
    new_center = lutData{4}(LUTindex);
    mpp = lutData{5}(LUTindex);
    ppd = lutData{6}(LUTindex);

    % Load in the coordinates
    og_coords = csvread(fullfile(dataPath,fnameList{i}));

    % subtract og center
    sub_center(:,1) = og_coords(:,1)-og_center;
    sub_center(:,2) = og_coords(:,2)-og_center;

    % multiply by scaling factor
    mult_scale(:,1) = sub_center(:,1) * scaling_factor;
    mult_scale(:,2) = sub_center(:,2) * scaling_factor;

    % add new center to get new coords
    new_coords(:,1) = mult_scale(:,1) + new_center;
    new_coords(:,2) = mult_scale(:,2) + new_center;

    % format new scales to be correct for file name
    str_mpp = num2str(mpp);
    new_mpp = strrep(str_mpp, ".", "p");

    str_ppd = num2str(ppd);
    new_ppd = strrep(str_ppd, ".", "p");

    % create new file name
    parts = split(fnameList{i}, "_");
    
    new_name = [parts{1}, "_" parts{2} "_" parts{3} "_" parts{4} "_" new_mpp "mpp_" new_ppd "ppd_" parts{7} "_" parts{8} "_" parts{9} "_scaled_" parts{10} ".csv"];
    new_name = strjoin(new_name);
    new_name = strrep(new_name,' ','');

    writecell(num2cell(new_coords), fullfile(resultPath, new_name) )

    clear parts
    clear sub_center
    clear mult_scale
    clear new_coords






end




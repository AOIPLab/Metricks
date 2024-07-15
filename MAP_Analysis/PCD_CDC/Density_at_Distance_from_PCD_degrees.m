% Density at Distance from PCD in degrees (with LUT file)
%
% Created by: Jenna Grieshop
% Date created: 7/9/2024
%
% Description: Finds density at a specific distance away from the PCD in a 
% density matrix. User enters the distance in degrees from the 
% PCD. Then selects the direction (Nasal, Temporal, Superior, Inferior)
%
% Input: Denisty matrices in an input folder, a LUT file
% (file name, axial length, ppd, x PCD coordinates, y PCD coordinates) (LUT file can be in the same folder as
% the data or in a different folder), user input about distance in degrees
% and direction
%
% Output: .csv saved to folder containing the LUT file.


clear all
close all
clc

basepath = which('Density_at_Distance_from_PCD_degrees.m');
[basepath] = fileparts(basepath);
path(path,fullfile(basepath,'lib')); % Add our support library to the path.


% User selects folder with data
dataPath = uigetdir('.','Select directory containing analyses');

% Read in csv names and then have user select the LUT
[fnameList] = read_folder_contents(dataPath,'csv');
[scalingFname, scalingPath] = uigetfile(fullfile(dataPath,'*.csv'),'Select scaling LUT.');

% Remove LUT file from fnameList
fnameList(ismember(fnameList,scalingFname))=[];

% If LUT wasn't selected the user has the chance to enter the scale
% manually
scaleInput = NaN;
if scalingFname == 0        
    
    while isnan(scaleInput)                
        
        scaleInput = inputdlg('Input the scale in ppd:','Input the scale in ppd:');
        
        scaleInput = str2double(scaleInput);
        
        if isempty(scaleInput)
            error('Cancelled by user.');
        end
    end
else
    [~, lutData] = load_scaling_file_ddd(fullfile(scalingPath,scalingFname));
end

% prompts user to enter degrees out from PCD
DegreeInput = NaN;
while isnan(DegreeInput)

    DegreeInput = inputdlg('Input the distance from the PCD in degrees:', 'Input the distance from the PCD in degrees:');
    
    DegreeInput = str2double(DegreeInput);

    if isempty(DegreeInput)
            error('Degree input Cancelled by user.');
    end
end

% Direction selected by user
list = {'Superior', 'Temporal', 'Inferior', 'Nasal'};
[indx, tf] = listdlg('PromptString', 'Select the desired direction.', 'SelectionMode', 'single', 'ListString', list);
if tf == 1
    Direction = list{indx};
else
    % Canceled dialog box - end the program
    return
end

count = 1;
for i=1:size(fnameList,1) % Go through all files in list 
    
     if isnan(scaleInput)
                % Calculate the scale for this identifier.                                
                LUTindex=find( cellfun(@(s) ~isempty(strfind(fnameList{i},s )), lutData{1} ) );

                % Use whichever scale is most similar to our filename.
                sim = 1000*ones(length(LUTindex),1);
                for l=1:length(LUTindex)
                    sim(l) = lev(fnameList{i}, lutData{1}{LUTindex(l)});
                end
                [~,simind]=min(sim);
                LUTindex = LUTindex(simind);
                
                % get information from LUT file
                axialLength = lutData{2}(LUTindex);
                pixelsPerDegree = lutData{3}(LUTindex);

                % Get PCD coordinates from LUT file
                PCD_x = lutData{4}(LUTindex);
                PCD_y = lutData{5}(LUTindex);
                
            else
                pixelsPerDegree = scaleinput;
    end

    % Load in the density map
    densityMap = csvread(fullfile(dataPath,fnameList{i}));
    if contains(fnameList{i}, 'OS')
        densityMap = fliplr(densityMap);
        PCD_x = length(densityMap)-(PCD_x-1);
    end

    % All for OD orientation
    if strcmp(Direction,'Nasal')
        xDist = round(DegreeInput * pixelsPerDegree);
        yDist = 0;
    elseif strcmp(Direction,'Temporal')
        xDist = -round(DegreeInput * pixelsPerDegree);
        yDist = 0;
    elseif strcmp(Direction,'Superior')
        xDist = 0;
        yDist = round(DegreeInput * pixelsPerDegree);
    elseif strcmp(Direction,'Inferior')
        xDist = 0;
        yDist = -round(DegreeInput * pixelsPerDegree);
    else
        error('No valid direction.');
    end

    % Get the location that we need to get the density value from
    newX = PCD_x + xDist;
    newY = PCD_y + yDist;

    % Get the density at the point
    % if the new point is negative it does not exist
    if newX < 0 || newY <0 
        d_at_point = NaN;
    else
        d_at_point = densityMap(newY, newX); % x and y appear flipped bc of row and column rules in matlab
    end

    % Compile data for file
    if (i ==1)
        data = [d_at_point];
    else
        data = [data; d_at_point];
    end


end

% Write the output csv
data = num2cell(data);
header = {'File Name', 'Density at Point'};
compiled = cat(2, fnameList, data);
compiled2 = cat(1, header, compiled);
writecell(compiled2, fullfile(dataPath, ['Density_at_point_', num2str(DegreeInput), 'deg_', num2str(Direction), '_from_PCD_', datestr(now, 'dd_mmm_yyyy'), '.csv']))






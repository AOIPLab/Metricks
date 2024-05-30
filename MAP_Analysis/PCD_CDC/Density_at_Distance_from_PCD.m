% Density at Distance from PCD
%
% Created by: Jenna Grieshop
% Date created: 11/22/23
%
% Description: Finds density at a specific distance away from the PCD in a 
% density matrix. User enters the x and y distance in um (microns) from the 
% PCD. Negative distances are to the left and up.
%
% Input: Denisty matrices in the folder the code is running from, a LUT file
% (file name, axial length, and ppd) (LUT file can be in the same folder as
% the data or in a different folder), user input about distance in x and y 
% from the PCD
%
% Output: .csv saved to folder containing the LUT file.



clear all
close all
clc

basepath = which('Density_at_Distance_from_PCD.m');
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
        
        scaleInput = inputdlg('Input the scale in UNITS/PIXEL:','Input the scale in UNITS/PIXEL:');
        
        scaleInput = str2double(scaleInput);
        
        if isempty(scaleInput)
            error('Cancelled by user.');
        end
    end
else
    [~, lutData] = load_scaling_file(fullfile(scalingPath,scalingFname));
end

% prompts user to enter x distance from PCD
xInput = NaN;
while isnan(xInput)

    xInput = inputdlg('Input the X distance from the PCD in um:', 'Input the X distance from the PCD in um:');
    
    xInput = str2double(xInput);

    if isempty(xInput)
            error('X input Cancelled by user.');
    end
end

% prompts user to enter y distance from PCD
yInput = NaN;
while isnan(yInput)
    
    yInput = inputdlg('Input the Y distance from the PCD in um:', 'Input the Y distance from the PCD in um:');
    
    yInput = str2double(yInput);
    if isempty(yInput)
            error('Y input cancelled by user.');
    end
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

                % calculate microns per degree
                micronsPerDegree = (291*axialLength)/24;
                
                % get the scale value
                scaleVal = 1 / (pixelsPerDegree / micronsPerDegree); %pixels/um
            else
                scaleVal = scaleinput;
    end

    
    % Load in the density map and find the max value (peak)
    densityMap = csvread(fullfile(dataPath,fnameList{i}));
    peak = max(densityMap(:));

    % Check that the max value/peak is unique         
    [maxYCoords, maxXCoords] = find(densityMap == peak);
    allMaxCoords = [maxXCoords, maxYCoords];
    for j = 1:length(maxXCoords)
        allMaxes{count,1} = {fnameList{i}, allMaxCoords(j,1), allMaxCoords(j,2)}; % store all maxes
        count = count +1;
    end

    % Finding the weighted (mean) centoid if multiple max locations found
    centroidX = mean(maxXCoords);
    centroidY = mean(maxYCoords);

    % Finds the distance in pixels
    xDist = round(xInput * scaleVal);
    yDist = round(yInput * scaleVal);

    % Get the location that we need to get the density value from
    newX = centroidX + xDist;
    newY = centroidY + yDist;

    % Get the density at the point
    d_at_point = densityMap(newY, newX); % x and y appear flipped bc of row and column rules in matlab
    
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
writecell(compiled2, fullfile(dataPath, ['Density_at_point_', num2str(xInput), 'um_', num2str(yInput), 'um_from_PCD_', datestr(now, 'dd_mmm_yyyy'), '.csv']))


% Stand Alone Map Creator
% Purpose: To make maps separately from the Mosaic Metricks Map script
% Input: Folder with window_results .mat files
% - The script can br run on a batch of files
% Output: Map tif
% Date created: 1/24/24
% Jenna Grieshop
% Updated by MG on 4/13/2026 To display maps for updated Map script output

clear all;
close all;
clc;

basepath = which('Stand_Alone_Map_Creator.m');
[basepath] = fileparts(basepath);

path(path,fullfile(basepath,'lib')); % Add our support library to the path.

% Get the file and path 
filepath = uigetdir('.','Select a directory of .mat files (one type) to generate a map figure');
[filenames, isadir ] = read_folder_contents(filepath,'mat');


% Clims can either be found automatically or the user can set their own
% bounds to ensure all figures of a given folder will be on the same clims
climType = questdlg('How would you like to set your clims?', ...
	'Specify clims', ...
	'Auto','Manual','Auto');

Map_min = inf;
Map_max = -inf;

for i=1:length(filenames)

    % Load in current map
    CurrentMap = load(fullfile(filepath,filenames{i}));
    dummyName = fieldnames(CurrentMap);

    temp_min = min(CurrentMap.(dummyName{1}),[], 'all');
    temp_max = max(CurrentMap.(dummyName{1}),[], 'all');

    if temp_min < Map_min
        Map_min = temp_min;
    end
    if temp_max > Map_max
        Map_max = temp_max;
    end

end
    
    
    switch climType
        case 'Auto'
            % Find the min and max value of the map and +/- 10% to pad the
            % clims slightly - per JC 4/13/2026

            Clim_min = Map_min - (0.10*Map_min);
            Clim_max = Map_max + (0.10*Map_max);

        case 'Manual'
            % Launch another popup window so the user can type in the min and
            % max clim values of their choice... handy for plotting multiple
            % maps on the same set of clims
            manual_prompt = {'Clim Min', 'Clim Max'};
            dlgtitle = 'Please define Clim values';
            fieldsize = [1 20; 1 20];
            definput = {'0','10000'};
            manual_clims = inputdlg(manual_prompt,dlgtitle,fieldsize,definput);
    
            Clim_min = str2double(manual_clims{1});
            Clim_max = str2double(manual_clims{2});
    end

for i=1:length(filenames)

    % Load in current map
    CurrentMap = load(fullfile(filepath,filenames{i}));
    dummyName = fieldnames(CurrentMap);
    
    % breaking apart the map filenames to find the subject ID and the map type
    % + unit. Note for VCAR it will say matrix for unit...(it's unitless)
    fileparts = split(filenames{i}, '_');
    subjectID = [fileparts{1}, '_', fileparts{2}];
    selectedMap = [fileparts{5}, '_', fileparts{6}];
    
    interpedMap = CurrentMap.(dummyName{1});
    smoothedInterpedMap = imgaussfilt(interpedMap,20); % Filters interpedMap with a 2-D Gaussian smoothing kernel with standard deviation 20
    
    % Set nans to 0
    interpedMap(isnan(interpedMap)) =0;
    smoothedInterpedMap(isnan(smoothedInterpedMap)) =0;
    
    vmap=viridis; %calls viridis colormap function from library
    
    clims = [Clim_min Clim_max];
    
    % Plot the map
    dispfig=figure(1); 
    imagesc(interpedMap,clims); % use limits of color scale set above
    axis image;
    colormap(vmap); 
    colorbar; 
    
    % Get the min and max values
    [minVal, minInd] = min(interpedMap(:));
    [maxVal, maxInd] = max(interpedMap(:));
    
    [minRow,minCol]=ind2sub(size(interpedMap),minInd);
    [maxRow,maxCol]=ind2sub(size(interpedMap),maxInd);
    
    % Set title and file name
    title(['Minimum value: ' num2str(minVal) '(' num2str(minCol) ',' num2str(minRow) ') Maximum value: ' num2str(maxVal) '(' num2str(maxCol) ',' num2str(maxRow) ')']);
    resultFname = [selectedMap, '_map_' datestr(now, 'dd_mmm_yyyy')];
    
    % Scale to the max of clims 
    scaledMap = interpedMap-min(clims);
    scaledMap(scaledMap <0) =0; % In case there are min values below this
    scaledMap = uint8(255*scaledMap./(max(clims)-min(clims)));
    scaledMap(scaledMap  >255) = 255; % In case there are values above this
    imwrite(scaledMap, vmap, fullfile(filepath,[subjectID, '_', resultFname '_clims_' num2str(Clim_min) '_' num2str(Clim_max) '_raw.tif'])); % Save map image



end





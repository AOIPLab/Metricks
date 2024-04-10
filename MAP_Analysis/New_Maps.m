% Stand Alone Map Creator
% Purpose: To make maps separately from the Mosaic Metricks Map script
% Input: Folder with window_results .mat files, coordinate files, and
% tif images for each subject.
% - The script can run multiple subjects at a time or just one
% Output: Map or Map and csv if bound_density_deg selected
% Date created: 1/24/24
% Jenna Grieshop

clear all
clc

% Update clims based on min and max of your data
clims = [500 2500]; % Set limits of color scale so all images use the same scale

% Options of maps to create for user to select
liststr = {'bound_area','unbound_area','bound_num_cells', 'unbound_num_cells', 'bound_density_deg', 'bound_density'};
[selectedMap, oked] = listdlg('PromptString','Select map type:',...
                              'SelectionMode','single',...
                              'ListString',liststr);
if oked == 0
    error('Cancelled by user.');
end

selectedMap = liststr{selectedMap};  

% Get directory of data
rootPath = uigetdir('.','Select directory containing analyses');
rootDir = dir(rootPath);
rootDir = struct2cell(rootDir)';

% Find all csv and txt files
[fnameList, isadir ] = read_folder_contents(rootPath,'csv');
[fnameListTxt, isDirTxt ] = read_folder_contents(rootPath,'txt');

fnameList = [fnameList; fnameListTxt];
isadir = [isadir;isDirTxt];


% Looks for all the window results
winResultsDir = rootDir(...
    ~cellfun(@isempty, strfind(rootDir(:,1), 'window_results')),:);


for i=1:size(fnameList,1)

    subjectID = fnameList{i}(1:8);

    % Read in coordinates - assumes x,y
    coords=dlmread(fullfile(rootPath,fnameList{i}));
    
    % It should ONLY be a coordinate list, that means x,y, and
    % nothing else.
    if size(coords,2) ~= 2
        warning('Coordinate list contains more than 2 columns! Skipping...');
        continue;
    end

    % If the image exists use it for the dimensions
    if exist(fullfile(rootPath, [fnameList{i}(1:end-length('_coords.csv')) '.tif']), 'file')

        im = imread( fullfile(rootPath, [fnameList{i}(1:end-length('_coords.csv')) '.tif']));

        width = size(im,2);
        height = size(im,1);
        maxRowVal = height;
        maxColVal = width;
    % If it doesn't exist, warn the user and then use coords for the
    % dimensions
    else
        warning(['No matching image file found for ' fnameList{i}]);
        coords = coords-min(coords)+1;
        width  = ceil(max(coords(:,1)));
        height = ceil(max(coords(:,2)));
        maxRowVal = max(coords(:,2));
        maxColVal = max(coords(:,1));
    end


    % Load in the data
    data = load(fullfile(winResultsDir{i,2}, winResultsDir{i,1}));

    % Initialize the map and meshgrid
    interpedMap=zeros([height width]);
    [Xq, Yq] = meshgrid(1:size(im,2), 1:size(im,1));
    
    
    % Based on the map the user selected, create the interpolation with the
    % correct data in the struct
    if selectedMap == "bound_density_deg"
        scattah = scatteredInterpolant(coords(:,1), coords(:,2), data.win_res.bound_density_DEG);  
    elseif selectedMap == "bound_density"
        scattah = scatteredInterpolant(coords(:,1), coords(:,2), data.win_res.bound_density); 
    elseif selectedMap == "bound_area"
        scattah = scatteredInterpolant(coords(:,1), coords(:,2), data.win_res.bound_area);       
    elseif selectedMap == "unbound_area"
        scattah = scatteredInterpolant(coords(:,1), coords(:,2), data.win_res.unbound_area);
    elseif selectedMap == "bound_num_cells"
        scattah = scatteredInterpolant(coords(:,1), coords(:,2), data.win_res.bound_num_cells);
    elseif selectedMap == "unbound_num_cells"
        scattah = scatteredInterpolant(coords(:,1), coords(:,2), data.win_res.unbound_num_cells);
    else
        disp("something is wrong");
    end

    interpedMap = scattah(Xq,Yq);
    smoothedInterpedMap = imgaussfilt(interpedMap,20); % Filters interpedMap with a 2-D Gaussian smoothing kernel with standard deviation 20
    
    % Set nans to 0
    interpedMap(isnan(interpedMap)) =0;
    smoothedInterpedMap(isnan(smoothedInterpedMap)) =0;
    
    vmap=viridis; %calls viridis colormap function from library
    
   
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
    imwrite(scaledMap, vmap, fullfile(rootPath,[subjectID, '_', resultFname '_raw.tif'])); % Save map image

    % If bound density deg selected save csv and mat file as well
    if selectedMap == "bound_density_deg"
        fileName = fullfile(rootPath,'Results',[subjectID '_bounddensity_matrix_DEG_' datestr(now, 'dd_mmm_yyyy') '.csv']);
        writematrix(interpedMap, fileName);
        % Save matrix as matfile
        save(fullfile(rootPath,'Results',[subjectID '_bounddensity_matrix_DEG_MATFILE_' datestr(now, 'dd_mmm_yyyy') '.mat']), "interpedMap");
    end

end





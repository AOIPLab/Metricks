% Standard deviation of maps
% 1/26/2024
% Jenna Grieshop

clear all
clc
addpath('lib');

sub_id = {''};


root_path = uigetdir('.','Select directory containing analyses');
root_dir = dir(root_path);
root_dir = root_dir(~ismember({root_dir.name}, {'.', '..'}));
root_dir = struct2cell(root_dir)';

% get all the file paths that we are interested in in each of the folders
[file_paths] = read_folder_contents_rec(root_path, 'mat', 'MATFILE');


% finding all the subject Ids within the folders
for i=1:size(file_paths,1)
    spl = split(file_paths{i,1}, "_");
    spl = split(spl{1,1}, "\");
    spl = spl(end);
    if contains(sub_id{:}, spl)
        continue
    else
        sub_id{i} = spl{1};
    end
    
end

for j=1:size(sub_id,2)

    index = find(contains(file_paths,sub_id(j)));
    for m=1:size(index,1)
        data = load(file_paths{m});
        maps{m} = round(data.interped_map,2);
        A(:,:,m) = maps{m};
    end

    standard_dev = std(A, [], 3);
    vmap=viridis; %calls viridis colormap function
    
    clims = [0 25000]; % added to set limits of color scale, so all images use the same scale
    
    dispfig=figure(1); 
    imagesc(standard_dev); % added to use limits of color scale
    axis image;
    colormap(vmap); 
    colorbar; 

    [minval, minind] = min(standard_dev(:));
    [maxval, maxind] = max(standard_dev(:));
    
    [minrow,mincol]=ind2sub(size(standard_dev),minind);
    [maxrow,maxcol]=ind2sub(size(standard_dev),maxind);
    
    max_x_vals = maxcol;
    max_y_vals = maxrow;

    %updated to scale to the max of clims 10/11/23     
    scaled_map = standard_dev-min(clims);
    scaled_map(scaled_map <0) =0; %in case there are min values below this
    scaled_map = uint8(255*scaled_map./(max(clims)-min(clims)));
    scaled_map(scaled_map  >255) = 255; %in case there are values above this


    subjectID = sub_id(j); 
    subjectID = subjectID{1};
    result_fname = [subjectID '_stdev_' date '_raw.tif'];
    imwrite(scaled_map, vmap, fullfile(root_path,result_fname));

    
end





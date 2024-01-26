% Standard deviation of maps
% 1/26/2024
% Jenna Grieshop

clear all
clc


root_path = uigetdir('.','Select directory containing analyses');
root_dir = dir(root_path);
root_dir = root_dir(~ismember({root_dir.name}, {'.', '..'}));
root_dir = struct2cell(root_dir)';


temp_dir = dir(fullfile(root_dir{1,2}, root_dir{1,1}));
temp_dir = struct2cell(temp_dir)';
% looks for all the MATFILE matrices
MATFILE_dir = temp_dir(...
    ~cellfun(@isempty, strfind(temp_dir(:,1), 'MATFILE')),:);

% finding all the subject Ids within the folders
for i=1:size(MATFILE_dir,1)
    spl = split(MATFILE_dir{i,1}, "_");
    spl = spl(1);
    SubID(i) = spl;
end

N = size(root_dir,1);
for j=1:size(SubID,2)
    for k=1:size(root_dir,1)
        fol_dir = dir(fullfile(root_dir{k,2}, root_dir{k,1}));
        fol_dir = struct2cell(fol_dir)';
        % looks for all the MATFILE matrices
        MATFILE_dir = temp_dir(...
            ~cellfun(@isempty, strfind(temp_dir(:,1), 'MATFILE')),:);
        index = strfind(MATFILE_dir(:,1), SubID(j));
        data = load(fullfile(MATFILE_dir{index{1,1},2}, MATFILE_dir{index{1,1},1}));
        maps{j,k} = round(data.interped_map,2);
        A(:,:,k) = maps{j,k};
    end
    standard_dev = std(A, [], 3);
    vmap=viridis; %calls viridis colormap function, added by Joe 2/19/22
    
    clims = [0 7810^-11]; % added to set limits of color scale, so all images use the same scale by Joe 2/19/22
    
    dispfig=figure(1); 
    imagesc(standard_dev); % added to use limits of color scale, by Joe 2/19/22
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


    subjectID = SubID(j); 
    subjectID = subjectID{1,1};
    result_fname = [subjectID '_stdev_' date '_raw.tif'];
    imwrite(scaled_map, vmap, fullfile(root_path,result_fname));

    
end





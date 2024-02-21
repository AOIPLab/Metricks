% Standard deviation of maps of different window sizes
% 1/26/2024
% Jenna Grieshop

clear all
clc
addpath('lib');

sub_id = {''};


root_path_bd = uigetdir('.','Select directory containing bound density analyses');

root_path_cdc = uigetdir('.','Select directory containing CDC analyses');

% get all the file paths that we are interested in in each of the folders
[file_paths_bd] = read_folder_contents_rec(root_path_bd, 'mat', 'MATFILE');

%get all the file paths that we are interested in in each of the folders
[file_paths_cdc] = read_folder_contents_rec(root_path_cdc, 'csv', 'PCD_CDC_Analysis_Summary');


% finding all the subject Ids within the folders
for i=1:size(file_paths_bd,1)
    spl = split(file_paths_bd{i,1}, "_bound");
    spl = split(spl{1,1}, "\");
    spl = spl(end);
    ispresent = cellfun(@(s) ~isempty(strfind(spl{:}, s)), sub_id);
    if any(ispresent)
        continue
    else
        sub_id{i} = spl{1};
    end
    
end

sum_x = cell(size(sub_id, 2),1);
sum_x(:,1) = {0};
sum_y = cell(size(sub_id, 2),1);
sum_y(:,1) = {0};

for m=1:size(file_paths_cdc,1)
    cdc_data = readtable(file_paths_cdc{m});
    for n=1:size(cdc_data, 1)
        sum_x{n} = sum_x{n} + cdc_data{n,8};
        sum_y{n} = sum_y{n} + cdc_data{n,9};

    end
end

sum_x = cell2mat(sum_x(:)); % convert so can be divided
sum_y = cell2mat(sum_y(:)); % convert so can be divided


avg_x(:) = sum_x(:)/size(file_paths_cdc,1);
avg_y(:) = sum_y(:)./size(file_paths_cdc,1);


for j=1:size(sub_id,2)

    clear maps;
    clear A;
    
    index = find(contains(file_paths_bd,sub_id(j)));
    for m=1:size(index,1)
        data = load(file_paths_bd{index(m)});
        maps{m} = round(data.interped_map,2);
        A(:,:,m) = maps{m}; % maps is an unnecesary middle step - but good for trouble shooting
    end

    standard_dev = std(A, [], 3);
    average = mean(A,3);
    coeffofvar = standard_dev/average;
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
    imwrite(scaled_map, vmap, fullfile(root_path_bd,result_fname));

    result_fname2 = [subjectID '_stdev_' date '_marked.tif'];
    scaled_map_mark = uint8(255*standard_dev./max(clims));
    MARK = insertShape(scaled_map_mark,'circle',[avg_x(j) avg_y(j) 2], 'LineWidth' ,3, 'Color' , 'red');
    imwrite(MARK, vmap, fullfile(root_path_bd,result_fname2));

    result_fname3 = [subjectID '_stdev_' date '_raw.csv'];
    csvwrite(fullfile(root_path_bd, result_fname3), standard_dev);

    result_fname5 = [subjectID '_coeffvar_' date '_raw.csv'];
    csvwrite(fullfile(root_path_bd, result_fname5), coeffofvar);
    
    master_cdc(j,1) = {sub_id{j}};
    master_cdc(j,2) = {avg_x(j)};
    master_cdc(j,3) = {avg_y(j)};
  
end

header = {'Subject ID', 'X CDC master', 'Y CDC master'};
finaloutput = [header; master_cdc];
result_fname4 = ['stdev_' date '_master_cdc.xlsx'];
% csvwrite(fullfile(root_path_bd, result_fname4), cell2mat(master_cdc));
% writecell(finaloutput, fullfile(root_path_bd, result_fname4));
xlswrite(fullfile(root_path_bd, result_fname4), master_cdc);



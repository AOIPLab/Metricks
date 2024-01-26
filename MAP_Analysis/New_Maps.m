% Script to make maps separately from the Mosaic map script
% 1/24/24
% Jenna Grieshop

clear all
clc


liststr = {'bound_area','unbound_area','bound_num_cells', 'unbound_num_cells', 'bound_density_deg', 'bound_density'};
[selectedmap, oked] = listdlg('PromptString','Select map type:',...
                              'SelectionMode','single',...
                              'ListString',liststr);
if oked == 0
    error('Cancelled by user.');
end

root_path = uigetdir('.','Select directory containing analyses');
root_dir = dir(root_path);
root_dir = struct2cell(root_dir)';

selectedmap = liststr{selectedmap};   

[fnamelist, isadir ] = read_folder_contents(root_path,'csv');
[fnamelisttxt, isdirtxt ] = read_folder_contents(root_path,'txt');

fnamelist = [fnamelist; fnamelisttxt];
isadir = [isadir;isdirtxt];


% looks for all the window results
win_results_dir = root_dir(...
    ~cellfun(@isempty, strfind(root_dir(:,1), 'window_results')),:);


for i=1:size(fnamelist,1)

    subject_ID = fnamelist{i}(1:8);

    %Read in coordinates - assumes x,y
    coords=dlmread(fullfile(root_path,fnamelist{i}));
    
    % It should ONLY be a coordinate list, that means x,y, and
    % nothing else.
    if size(coords,2) ~= 2
        warning('Coordinate list contains more than 2 columns! Skipping...');
        continue;
    end

    if exist(fullfile(root_path, [fnamelist{i}(1:end-length('_coords.csv')) '.tif']), 'file')

        im = imread( fullfile(root_path, [fnamelist{i}(1:end-length('_coords.csv')) '.tif']));

        width = size(im,2);
        height = size(im,1);
        maxrowval = height;
        maxcolval = width;
    else
        warning(['No matching image file found for ' fnamelist{i}]);
        coords = coords-min(coords)+1;
        width  = ceil(max(coords(:,1)));
        height = ceil(max(coords(:,2)));
        maxrowval = max(coords(:,2));
        maxcolval = max(coords(:,1));
    end


    data = load(fullfile(win_results_dir{i,2}, win_results_dir{i,1}));

    interped_map=zeros([height width]);
    [Xq, Yq] = meshgrid(1:size(im,2), 1:size(im,1));
    
    
    if selectedmap == "bound_density_deg"
        scattah = scatteredInterpolant(coords(:,1), coords(:,2), data.win_res.bound_density_DEG);  
    elseif selectedmap == "bound_density"
        scattah = scatteredInterpolant(coords(:,1), coords(:,2), data.win_res.bound_density); 
    elseif selectedmap == "bound_area"
        scattah = scatteredInterpolant(coords(:,1), coords(:,2), data.win_res.bound_area);       
    elseif selectedmap == "unbound_area"
        scattah = scatteredInterpolant(coords(:,1), coords(:,2), data.win_res.unbound_area);
    elseif selectedmap == "bound_num_cells"
        scattah = scatteredInterpolant(coords(:,1), coords(:,2), data.win_res.bound_num_cells);
    elseif selectedmap == "unbound_num_cells"
        scattah = scatteredInterpolant(coords(:,1), coords(:,2), data.win_res.unbound_num_cells);
    else
        disp("something is wrong");
    end

    interped_map = scattah(Xq,Yq);
    smoothed_interped_map = imgaussfilt(interped_map,20);
    
    interped_map(isnan(interped_map)) =0;
    smoothed_interped_map(isnan(smoothed_interped_map)) =0;
    
    vmap=viridis; %calls viridis colormap function, added by Joe 2/19/22
    
    clims = [500 2500]; % added to set limits of color scale, so all images use the same scale by Joe 2/19/22
    
    dispfig=figure(1); 
    imagesc(interped_map,clims); % added to use limits of color scale, by Joe 2/19/22
    axis image;
    colormap(vmap); 
    colorbar; 
    [minval, minind] = min(interped_map(:));
    [maxval, maxind] = max(interped_map(:));
    
    [minrow,mincol]=ind2sub(size(interped_map),minind);
    [maxrow,maxcol]=ind2sub(size(interped_map),maxind);
    
    max_x_vals = maxcol;
    max_y_vals = maxrow;
    
   % disp([subject_ID ' Maximum value: ' num2str(round(maxval)) '(' num2str(maxcol) ',' num2str(maxrow) ')' ]) % display added by Katie Litts in 2019
               
    title(['Minimum value: ' num2str(minval) '(' num2str(mincol) ',' num2str(minrow) ') Maximum value: ' num2str(maxval) '(' num2str(maxcol) ',' num2str(maxrow) ')']);
    
    result_fname = [selectedmap, '_map_' date];
    
    %updated to scale to the max of clims 10/11/23     
    scaled_map = interped_map-min(clims);
    scaled_map(scaled_map <0) =0; %in case there are min values below this
    scaled_map = uint8(255*scaled_map./(max(clims)-min(clims)));
    scaled_map(scaled_map  >255) = 255; %in case there are values above this
    imwrite(scaled_map, vmap, fullfile(root_path,[subject_ID, '_', result_fname '_raw.tif'])); %added by Joe Carroll 

    if selectedmap == "bound_density_deg"
        filename = fullfile(root_path,'Results',[subject_ID '_bounddensity_matrix_DEG_' date '.csv']);
        writematrix(interped_map, filename);
        %save matrix as matfile
        save(fullfile(root_path,'Results',[subject_ID '_bounddensity_matrix_DEG_MATFILE_' date '.mat']), "interped_map");
    end

end





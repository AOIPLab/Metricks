% Copyright (C) 2019 Robert F Cooper
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
% Metricks - A MATLAB package for analyzing the cone photoreceptor mosaic.
% 
% Coordinate_Mosiac_Metrics calculates the metrics for every
% image/coordinate pair in a given folder.
%
% When run, the script will prompt the user to select a folder with image/coordinate pairs.
% 
% **At present, images must be 8-bit grayscale tifs, coordinates must be formatted 
%   as a 2 column matrix (x,y), and must be named using the following convention,
%   where [imagename] can be any valid filename:**
% * Image File: [imagename].tif
% * Coordinate File: [imagename]\_coords.csv
% 
% It will then prompt the user to select what the output unit should be. At present,
% the options are:
% * Microns (using millimeters^2 for density)
% * Degrees
% * Arcminutes
% 
% Once the output unit is select, it will give the user the option to pick a 
% lookup table. The lookup table allows the software to analyze a folder of 
% images from different subjects/timepoints/conditions. The lookup table itself
% **must** be a 3 column 'csv' file, where the **first column** is a common 
% identifier for image/coordinate pairs, the **second column** is the axial 
% length (or '24' if the axial length is unknown) of the image/coordinate pairs,
% and the **third column** is the pixels per degree of the image/coordinate pairs.
% Each row must contain a different identifier/axial length/pixels per degree tuple.
% 
% An example common identifier could be a subject number, e.g, when working with the files:
% - 1235_dateoftheyear_OD_0004.tif
% - 1235_dateoftheyear_OD_0005.tif
% 
% Common identifiers could be "1235", "1235_dateoftheyear", "1235_dateoftheyear_OD".
% If all three were placed in a LUT, then the one that matches the most (as determined
% via levenshtein distance) will be used. In this example, we would use "1235_dateoftheyear_OD".
% 
% If we had another date, say: 1235_differentdateoftheyear_OD_0005.tif, then 
% _only_ the identifier "1235" would match between all images. However, say the
% two dates have different scales, then you would want to create two rows in the
% look up table for each date, with identifiers like: "1235_dateoftheyear" and
% "1235_differentdateoftheyear". To have the best chance at sucess, it is a
% good practice to just use the file name as the identifier.
% 
% **If you do not wish to use a lookup table, then press "cancel", and the software
% will allow you put in your own scale in UNITS/pixel.**
% 
% **This software has the ability to pre-crop the input data (if, for example,
% you have 80 pixels of coordinates and you only want to analyze the middle 50).**
% 
% To specify a cropping window, input the size (in the units you are going to 
% use) in to the brackets on line ~129 (of the variable 
% windowsize) of Coordinate_Mosaic_Metrics.m.
% 
% **Cropping is governed by the following rules:**
% 
% 1) If the tif is present and windowsize is not specified, the analysis will 
% be done on everything within the dimensions of the image.
% 2) If the tif is present and windowsize is specified, the assumed center of 
% the image is calculated according to the borders of the tif. **In either case,
% it doesn’t “care” how many (or even if there are any) cells in the image.**
% 3) If the tif is not present and windowsize is not specified, the analysis will
% be done on everything within the min and max coordinates in both x and y directions.
% So if you have an image in which there is an absence of cells on one side, 
% for example, you might end up with a clipped area that is not a square.
% 4) If the tif is not present and windowsize is specified, the assumed center 
% of the image is calculated according to the min and max coordinates in both 
% x and y directions. So if you have an image in which there is an absence of 
% cells on one side, the center will shift towards the other side of the image.
% 
% 
% The software will then run, and calculate every metric currently validated.
% 
% At present, it calculates the following metrics from each image and coordinate pair:
% 
% - Number of Unbound Cells
% - Number of Bound Cells
% - Total Area
% - Total Bounded Area
% - Mean Voronoi Area
% - Percent Six-Sided Voronoi
% - Density (uncorrected/corrected)
% - Nearest Neighbor Distance (uncorrected/corrected)
% - Inter-Cell Distance (uncorrected/corrected)
% - Furthest Neighbor Distance (uncorrected/corrected)
% - Density Recovery Profile Distance
% - Voronoi Area Regularity Index
% - Voronoi Number of Sides Regularity Index
% - Nearest Neighbor Regularity Index
% - Inter-Cell Regularity Index
% 
% The results will then be placed in to a datestamped file within a "Results" 
% folder as a subfolder of the one selected for analysis.
% 
%
% Don't thank me; cite me:
% 
% Every metric that is run via the main "Coordinate_Mosaic_Metrics.m" script
% has been validated and used in the following manuscript: 
% 
% Cooper RF, Wilk MA, Tarima S, Dubra A, Carroll J. 
% “Evaluating descriptive metrics of the human cone mosaic.”
% Invest Ophthalmol Vis Sci. 2016 57(7):2993.
% 
% You can also find formal definitions of each metric calculated here in that paper.
% 
% **This package is free for use under GPL v3, but I ask that you please cite 
% the above paper if you use this package.**
%
% This version (Taylor's Version) is specifically a batch mode for
% metrics. You will need to direct it to a folder of the folders you want
% to run metrics on.
% 


clear;
close all force;


% User to select directory with data
folder_path = uigetdir('.','Select directory containing folders of data to run');
all_subfolders = genpath(folder_path);
subfolders = regexp(all_subfolders, ';', 'split');

windowsize = [];
%% Crop the coordinates/image to this size in [scale], and calculate the area from it.
% If left empty, it uses the size of the image.

if length(windowsize) > 1
   error('Window size can only be empty ([]), or a single value!');
end

basePath = which('Coordinate_Mosaic_Metrics.m');

[basePath ] = fileparts(basePath);
path(path,fullfile(basePath,'lib')); % Add our support library to the path.

[basepath] = folder_path;


first = true;
                          

[scalingfname, scalingpath] = uigetfile(fullfile(folder_path,'*.csv'),'Select scaling LUT, OR cancel if you want to input the scale directly.');

scaleinput = NaN;
if scalingfname == 0        
    
    while isnan(scaleinput)                
        
        scaleinput = inputdlg('Input the scale in UNITS/PIXEL:','Input the scale in UNITS/PIXEL:');
        
        scaleinput = str2double(scaleinput);
        
        if isempty(scaleinput)
            error('Cancelled by user.');
        end
    end
else
    lutData = readtable(fullfile(scalingpath,scalingfname));
    lut_identifier = table2array(lutData(:,"Var1"));
    ALs = table2cell(lutData(:,"Var2"));
    ppds = table2cell(lutData(:,"Var3"));
end


for q = 2 :  length(subfolders)-1

    [fnamelist_1, isadir_1 ] = read_folder_contents(subfolders{q},'csv');
    [fnamelisttxt, isadirtxt ] = read_folder_contents(subfolders{q},'txt');
    
    fnamelist = [fnamelist_1; fnamelisttxt];
    isadir = [isadir_1;isadirtxt];


    %% Process the data.
    proghand = waitbar(0,'Processing...');
    first = true;
    
    for i=1:size(fnamelist,1)
    
        try
            if ~isadir{i}
    
                
                if length(fnamelist{i})>42
                    
                    waitbar(i/size(fnamelist,1), proghand, strrep(fnamelist{i}(1:42),'_','\_') );
                else
                    waitbar(i/size(fnamelist,1), proghand, strrep(fnamelist{i},'_','\_') );
                end
    
                if isnan(scaleinput)
                    % Calculate the scale for this identifier.                                
                    LUTindex=find( cellfun(@(s) ~isempty(strfind(fnamelist{i},s )), lut_identifier ) );
                    axiallength = ALs{LUTindex};
                    pixelsperdegree = ppds{LUTindex};

                    micronsperdegree = (291*axiallength)/24;
    
                    % microns or cones/mm^2 for density
                    scaleval_um = 1 / (pixelsperdegree / micronsperdegree);

                    % degrees
                    scaleval_deg = 1/pixelsperdegree;
    
                    % arcmin
                    scaleval_arcmin = 60/pixelsperdegree;
                    
                    
                else
                    % TODO: MG deal with this if we still want to have the
                    % possibility for user input scale
                    scaleval = scaleinput;
                end
    
                
                
                %Read in coordinates - assumes x,y
                coords=dlmread(fullfile(subfolders{q},fnamelist{i}));
                
                % It should ONLY be a coordinate list, that means x,y, and
                % nothing else.
                if size(coords,2) ~= 2
                    warning('Coordinate list contains more than 2 columns! Skipping...');
                    continue;
                end
    
                % If the corresponding image exists in the folder, use the image bounds to calculate our sizes
                if exist(fullfile(subfolders{q}, [fnamelist{i}(1:end-length('_coords.csv')) '.tif']), 'file')
    
                    im = imread( fullfile(subfolders{q}, [fnamelist{i}(1:end-length('_coords.csv')) '.tif']));
    
                    width = size(im,2);
                    height = size(im,1);
    
                    if ~isempty(windowsize)
                        pixelwindowsize = windowsize/scaleval;
    
                        diffwidth  = (width-pixelwindowsize)/2;
                        diffheight = (height-pixelwindowsize)/2;
                        
                        if diffwidth<0
                            diffwidth=0;
                        end                    
                        if diffheight<0
                            diffheight=0;
                        end
                    else
    
                        pixelwindowsize = [height width];
                        diffwidth=0;
                        diffheight=0;
                    end
    
                    clipped_coords =coordclip(coords,[diffwidth  width-diffwidth],...
                                                     [diffheight height-diffheight],'i');
    
                    clip_start_end = [diffwidth  width-diffwidth diffheight height-diffheight];
                else
                    warning(['No matching image file found for ' fnamelist{i}]);
                    width  = max(coords(:,1)) - min(coords(:,1));
                    height = max(coords(:,2)) - min(coords(:,2));
    
                    if ~isempty(windowsize)
                        pixelwindowsize = windowsize/scaleval;
    
                        diffwidth  = (width-pixelwindowsize)/2;
                        diffheight = (height-pixelwindowsize)/2;
                    else
                        pixelwindowsize = [height width];
                        diffwidth=0;
                        diffheight=0;
                    end
    
                    clipped_coords =coordclip(coords,[min(coords(:,1))-0.01+diffwidth  max(coords(:,1))-diffwidth+0.01],...
                                                     [min(coords(:,2))-0.01+diffheight max(coords(:,2))-diffheight+0.01],'i');
    
                    clip_start_end = [min(coords(:,1))+diffwidth-0.01  max(coords(:,1))-diffwidth+0.01 min(coords(:,2))+diffheight-0.01 max(coords(:,2))-diffheight+0.01];
                end
    
    
                [statistics_um, statistics_deg, statistics_arcmin] = determine_mosaic_stats(clipped_coords, pixelsperdegree, micronsperdegree, clip_start_end ,NaN, 4);
    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Determine FFT Power Spectra %%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % TODO update to have all units not just um
                if (exist('fit_fourier_spacing') == 2) && exist(fullfile(subfolders{q}, [fnamelist{i}(1:end-length('_coords.csv')) '.tif']), 'file')==2
                    
                    clipped_im = im(round(clip_start_end(3)+1:clip_start_end(4)), round(clip_start_end(1)+1:clip_start_end(2)) );
                    
                    [pixel_spac, ~, quality] = fit_fourier_spacing(clipped_im, min(size(clipped_im)), false,'row');
                    statistics.DFT_Row_Spacing = pixel_spac*scaleval_um;
                    statistics.DFT_Row_Quality = quality;
                    
                    [pixel_spac, ~, quality] = fit_fourier_spacing(clipped_im, min(size(clipped_im)), false,'cell');
                    statistics.DFT_Cell_Spacing = pixel_spac*scaleval_um;
                    statistics.DFT_Cell_Quality = quality;
                    
                    
                end
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Write Results %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if isempty(windowsize)
                result_fname_um = [getparent(subfolders{q},'short') '_coordstats_' date, '_um', '.csv'];
                result_fname_deg = [getparent(subfolders{q},'short') '_coordstats_' date, '_deg', '.csv'];
                result_fname_arcmin = [getparent(subfolders{q},'short') '_coordstats_' date, '_arcmin', '.csv'];
            else
                result_fname_um = [getparent(subfolders{q},'short') '_coordstats_' date '_' num2str(windowsize), '_um',  '.csv'];
                result_fname_deg = [getparent(subfolders{q},'short') '_coordstats_' date '_' num2str(windowsize) , '_deg', '.csv'];
                result_fname_arcmin = [getparent(subfolders{q},'short') '_coordstats_' date '_' num2str(windowsize) , '_arcmin', '.csv'];
            end

            % call write_metrics function
            write_metrics_results(subfolders{q}, result_fname_um, fnamelist{i}, statistics_um, first);
            write_metrics_results(subfolders{q}, result_fname_deg, fnamelist{i}, statistics_deg, first);
            write_metrics_results(subfolders{q}, result_fname_arcmin, fnamelist{i}, statistics_arcmin, first);
            first = false;
             
    
            end
        catch ex
            warning(['Unable to analyze ' fnamelist{i} ':']);
            warning([ex.message ', In file: ' ex.stack(1).file '  Line: ' num2str(ex.stack(1).line)]);
            warning('If warning on line 199, change {} to [] and vice versa on that line.');
        end
    end
    close(proghand);
end
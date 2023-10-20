% Copyright (C) 2019 Robert F Cooper, created 2017-09-29
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
% Coordinate_Mosiac_Metrics_MAP creates a map of an image/coordinate set
% across a set of image/coordinates.
% 
% To run this script, the script will prompt the user to select an folder containing
% image/coordinate pair.
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
% "1235_differentdateoftheyear".
% 
% **If you do not wish to use a lookup table, then press "cancel", and the software
% will allow you put in your own scale in UNITS/pixel.**
% 
% **This software will automatically adjust its window size to encompass up 
% 100 coordinates.**
% 
% However, if you wish to specify a sliding window size, input the size 
% (in the units you are going to in to the brackets of the variable "WINDOW_SIZE"
% on line ~92 of Coordinate_Mosaic_Metrics_MAP.m.
% 
% **Window size inclusion is governed by the following rules:**
% 
% 1) If the tif is present and windowsize is not specified, the analysis will 
% be done on everything within the dimensions of the image.
% 2) If the tif is present and windowsize is specified, the assumed center of 
% the image is calculated according to the borders of the tif. **In either case,
% it doesn�t �care� how many (or even if there are any) cells in the image.**
% 3) If the tif is not present and windowsize is not specified, the analysis will
% be done on everything within the min and max coordinates in both x and y directions.
% So if you have an image in which there is an absence of cells on one side, 
% for example, you might end up with a clipped area that is not a square.
% 4) If the tif is not present and windowsize is specified, the assumed center 
% of the image is calculated according to the min and max coordinates in both 
% x and y directions. So if you have an image in which there is an absence of 
% cells on one side, the center will shift towards the other side of the image.
%
% This script creates a map of metrics from a selected folder.
% Modified by Joe Carroll on 2/22/22 to add in some figure outputs
% Updated 339 & 340 to be pixelwindowsize(c), was previously
% pixelwindowsize
% Modified 10-11-23 to correct map scaling issue with Rob (was using max of
% the interped map, should use max of clims

clear;
close all force;

WINDOW_SIZE = [];

%% Crop the coordinates/image to this size in [scale], and calculate the area from it.
% If left empty, it uses the size of the image.


basePath = which('Coordinate_Mosaic_Metrics_MAP.m');

[basePath ] = fileparts(basePath);
path(path,fullfile(basePath,'lib')); % Add our support library to the path.

[basepath] = uigetdir(pwd);

[fnamelist, isdir ] = read_folder_contents(basepath,'csv');
[fnamelisttxt, isdirtxt ] = read_folder_contents(basepath,'txt');

fnamelist = [fnamelist; fnamelisttxt];
isdir = [isdir;isdirtxt];

liststr = {'microns (mm density)','degrees','arcmin'};
[selectedunit, oked] = listdlg('PromptString','Select output units:',...
                              'SelectionMode','single',...
                              'ListString',liststr);
if oked == 0
    error('Cancelled by user.');
end

selectedunit = liststr{selectedunit};                          

[scalingfname, scalingpath] = uigetfile(fullfile(basepath,'*.csv'),'Select scaling LUT, OR cancel if you want to input the scale directly.');

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
    [~, lutData] = load_scaling_file(fullfile(scalingpath,scalingfname));
end


%%
first = true;

proghand = waitbar(0,'Processing...');

for i=1:size(fnamelist,1)

    try
        if ~isdir{i}

            
            if length(fnamelist{i})>42
                waitbar(i/size(fnamelist,1), proghand, strrep(fnamelist{i}(1:42),'_','\_') );
            else
                waitbar(i/size(fnamelist,1), proghand, strrep(fnamelist{i},'_','\_') );
            end

            if isnan(scaleinput)
                % Calculate the scale for this identifier.                                
                LUTindex=find( cellfun(@(s) ~isempty(strfind(fnamelist{i},s )), lutData{1} ) );

                % Use whichever scale is most similar to our filename.
                sim = 1000*ones(length(LUTindex),1);
                for l=1:length(LUTindex)
                    sim(l) = lev(fnamelist{i}, lutData{1}{LUTindex(l)});
                end
                [~,simind]=min(sim);
                LUTindex = LUTindex(simind);
                
                axiallength = lutData{2}(LUTindex);
                pixelsperdegree = lutData{3}(LUTindex);

                micronsperdegree = (291*axiallength)/24;
                
                switch selectedunit
                    case 'microns (mm density)'
                        scaleval = 1 / (pixelsperdegree / micronsperdegree);
                    case 'degrees'
                        scaleval = 1/pixelsperdegree;
                    case 'arcmin'
                        scaleval = 60/pixelsperdegree;
                end
            else
                scaleval = scaleinput;
            end


            %Read in coordinates - assumes x,y
            coords=dlmread(fullfile(basepath,fnamelist{i}));
            
            % It should ONLY be a coordinate list, that means x,y, and
            % nothing else.
            if size(coords,2) ~= 2
                warning('Coordinate list contains more than 2 columns! Skipping...');
                continue;
            end

            if exist(fullfile(basepath, [fnamelist{i}(1:end-length('_coords.csv')) '.tif']), 'file')

                    im = imread( fullfile(basepath, [fnamelist{i}(1:end-length('_coords.csv')) '.tif']));

                    width = size(im,2);
                    height = size(im,1);
                    maxrowval = height;
                    maxcolval = width;
                else

                    coords = coords-min(coords)+1;
                    width  = ceil(max(coords(:,1)));
                    height = ceil(max(coords(:,2)));
                    maxrowval = max(coords(:,2));
                    maxcolval = max(coords(:,1));
                end

                statistics = cell(size(coords,1),1);
            
            if ~isempty(WINDOW_SIZE)
                
                pixelwindowsize = repmat(WINDOW_SIZE/scaleval,size(coords,1),1);
                
            else
                
                upper_bound = 150; %this is the number of BOUND cells to include
                
                if upper_bound > size(coords,1)
                    upper_bound = size(coords,1);
                end
                
                % Determine the window size dynamically for each coordinate
                pixelwindowsize = zeros(size(coords,1),1);

                parfor c=1:size(coords,1)               

                    thiswindowsize=1;
                    clipped_coords=[];
                    numbound=0;
                     while numbound < upper_bound
                        thiswindowsize = thiswindowsize+1;
                        rowborders = ([coords(c,2)-(thiswindowsize/2) coords(c,2)+(thiswindowsize/2)]);
                        colborders = ([coords(c,1)-(thiswindowsize/2) coords(c,1)+(thiswindowsize/2)]);

                        rowborders(rowborders<1) =1;
                        colborders(colborders<1) =1;
                        rowborders(rowborders>maxrowval) =maxrowval;
                        colborders(colborders>maxcolval) =maxcolval;

                        clipped_coords =coordclip(coords,colborders,...
                                                         rowborders,'i');
                        
                        % Ensure we're working with bound cells only.
                         if size(clipped_coords,1) > 5
                             % Next, create voronoi diagrams from the cells we've clipped.                             
                             [V,C] = voronoin(clipped_coords,{'QJ'}); % Returns the vertices of the Voronoi edges in VX and VY so that plot(VX,VY,'-',X,Y,'.')
 
                             bound = zeros(length(C),1);
                             for vc=1:length(C)
 
                                 vertices=V(C{vc},:);
 
                                 if (all(C{vc}~=1)  && all(vertices(:,1)<colborders(2)) && all(vertices(:,2)<rowborders(2)) ... % [xmin xmax ymin ymax] 
                                                  && all(vertices(:,1)>colborders(1)) && all(vertices(:,2)>rowborders(1))) 
                                     bound(vc) = 1;
                          
                                 end
                             end
 
                             numbound = sum(bound);
                         end

                     end
%                     axis([colborders rowborders])
                    pixelwindowsize(c) = thiswindowsize;
                end
            end
            disp('Determined window size.')
            
            %% Actually calculate the statistics
            for c=1:size(coords,1)
                
                rowborders = round([coords(c,2)-(pixelwindowsize(c)/2) coords(c,2)+(pixelwindowsize(c)/2)]); 
                colborders = round([coords(c,1)-(pixelwindowsize(c)/2) coords(c,1)+(pixelwindowsize(c)/2)]);

                rowborders(rowborders<1) =1;
                colborders(colborders<1) =1;
                rowborders(rowborders>maxrowval) =maxrowval;
                colborders(colborders>maxcolval) =maxcolval;
                
                clipped_coords =coordclip(coords,colborders,...
                                                 rowborders,'i');
                                             
                ccc=length(clipped_coords);
                % disp(ccc) if you want the number of clipped coordinates
                % displayed
                % [xmin xmax ymin ymax] 
                clip_start_end = [colborders rowborders];
                
                statistics{c} = determine_mosaic_stats( clipped_coords, scaleval, selectedunit, clip_start_end ,[colborders(2)-colborders(1) rowborders(2)-rowborders(1)], 4 );
                statistics{c}.Window_Size = pixelwindowsize(c)*scaleval;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Determine FFT Power Spectra %%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 if (exist('fit_fourier_spacing.m','file') == 2) && exist(fullfile(basepath, [fnamelist{i}(1:end-length('_coords.csv')) '.tif']), 'file')==2
%                     [pixel_spac, interped_map] = fit_fourier_spacing(im);
%                     statistics{c}.DFT_Spacing = pixel_spac*scaleval;                
%                 end


                warning off;
                [ success ] = mkdir(basepath,'Results');
                warning on;
            
            end
           
            
            %% Map output
%             metriclist = fieldnames(statistics{1});
%              [selectedmetric, oked] = listdlg('PromptString','Select map metric:',...
%                                            'SelectionMode','single',...
%                                            'ListString',metriclist);
%              
%              if oked == 0
%                  error('Cancelled by user.');
%             end

             %Hard code selection for bound density - added by JC 2/19/22
             metriclist = fieldnames(statistics{1});
             selectedmetric = 5; %5 = bound density, 7 = bound ICD      

            interped_map=zeros([height width]);
            sum_map=zeros([height width]);

            
            for c=1:size(coords,1)

                    thisval = statistics{c}.(metriclist{selectedmetric}); 

                    rowrange = round(coords(c,2)-(pixelwindowsize(c)/2):coords(c,2)+(pixelwindowsize(c)/2)); %fixedwith added c
                    colrange = round(coords(c,1)-(pixelwindowsize(c)/2):coords(c,1)+(pixelwindowsize(c)/2));

                    rowrange(rowrange<1) =[];
                    colrange(colrange<1) =[];
                    rowrange(rowrange>maxrowval) =[];
                    colrange(colrange>maxcolval) =[];
                  
                    interped_map(rowrange,colrange) = interped_map(rowrange,colrange) + thisval;
                    sum_map(rowrange, colrange) = sum_map(rowrange, colrange) + 1;

            end
                %

            interped_map = interped_map./sum_map;

            interped_map(isnan(interped_map)) =0;
            
            vmap=viridis; %calls viridis colormap function, added by Joe 2/19/22
            
            clims = [50000 225000]; % added to set limits of color scale, so all images use the same scale by Joe 2/19/22
            
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
            
            subjectID = lutData{1};% extract subject ID; added by Katie Litts in 2019
            disp([subjectID{LUTindex} ' Maximum value: ' num2str(round(maxval)) '(' num2str(maxcol) ',' num2str(maxrow) ')' ]) % display added by Katie Litts in 2019
                       
            title(['Minimum value: ' num2str(minval) '(' num2str(mincol) ',' num2str(minrow) ') Maximum value: ' num2str(maxval) '(' num2str(maxcol) ',' num2str(maxrow) ')'])
            
            result_fname = [fnamelist{i}(1:end-4) '_bound_map_' date '_' num2str(WINDOW_SIZE) metriclist{selectedmetric}];
            
            saveas(gcf,fullfile(basepath,'Results', [result_fname '_fig.png']));
            %updated to scale to the max of clims 10/11/23     
            scaled_map = interped_map-min(clims);
            scaled_map(scaled_map <0) =0; %in case there are min values below this
            scaled_map = uint8(255*scaled_map./(max(clims)-min(clims)));
            scaled_map(scaled_map  >255) = 255; %in case there are values above this
            imwrite(scaled_map, vmap, fullfile(basepath,'Results',[result_fname '_raw5.tif'])); %added by Joe Carroll 
            
            %Adding an output image with the marked location of peak density, added by Joe Carroll 2/19/22, updated to scale to the
            %max of clims 10/11/23
            scaled_map_mark = uint8(255*interped_map./max(clims));
            MARK = insertShape(scaled_map_mark,'circle',[maxcol maxrow 2], 'LineWidth' ,3, 'Color' , 'red');
            imwrite(MARK, vmap, fullfile(basepath,'Results',[result_fname '_marked.tif']));

            %save matrix of density values, added by Jenna Cava
            filename = fullfile(basepath,'Results',[subjectID{LUTindex} '_bounddensity_matrix_' date '.csv']);
            writematrix(interped_map, filename);
             
            %%
        end
    catch ex
        warning(['Unable to analyze ' fnamelist{i} ':']);
        warning([ex.message ', In file: ' ex.stack(1).file '  Line: ' num2str(ex.stack(1).line)]);
    end
end
close(proghand);

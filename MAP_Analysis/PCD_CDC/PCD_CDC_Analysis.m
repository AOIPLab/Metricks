%Creates coordinate files for density contours in AOSLO images and returns PCD & CDC data 
%NOTE, this requires a bounddensitymatrix file, NOT a cone coordinate file!

%Based on code created by JAC 9 May 2019
%Edited by J Grieshop March-12-2021 to automate file loading/saving
%Edited by Joe Carroll February-26-2022 to remove hard coded file folder
%save location and change output file structure.

clear all
close all
clc

basepath = which('PCD_CDC_Analysis.m');
[basepath] = fileparts(basepath);
path(path,fullfile(basepath,'lib2')); % Add our support library to the path.
%[basepath] = uigetdir(pwd); %select folder
[fnamelist] = read_folder_contents(basepath,'csv');
[scalingfname, scalingpath] = uigetfile(fullfile(basepath,'*.csv'),'Select scaling LUT.');

% Threshold percentile selection by the user
list = {'95', '90', '85', '80', '75', '70', '65', '60', '55', '50', '45', '40', '35', '30', '25', '20', '15', '10', '5'};
[indx, tf] = listdlg('PromptString', 'Select the desired threshold percentile.', 'SelectionMode', 'single', 'ListString', list);
if tf == 1
    threshold_percentile = str2num(list{indx})/100;
    thresh_str = list{indx};
else
    % canceled dialog box - end the program
    return
end

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

count = 1;

for i=1:size(fnamelist,1)
    
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
                
                scaleval = 1 / (pixelsperdegree / micronsperdegree);
            else
                scaleval = scaleinput;
    end

    
    densitymap = csvread(fnamelist{i});
    peak = max(densitymap(:));
    
    [maxval, maxind] = max(densitymap(:));

    % check that the max value is unique
    max_indices = find(densitymap == maxval);           
    [max_y_coords, max_x_coords] = find(densitymap == maxval);
    all_max_coords = [max_x_coords, max_y_coords];
    for j = 1:length(max_x_coords)
        all_maxes{count,1} = {fnamelist{i}, all_max_coords(j,1), all_max_coords(j,2)};
        count = count +1;
    end

    % finding the weighted (mean) centoid if multiple max locations found
    centroid_x = mean(max_x_coords);
    centroid_y = mean(max_y_coords);
   
    threshold = (densitymap >= (threshold_percentile*peak));
    contour = edge(threshold);
    
    %added for area
    pxareaAboveThresh = sum(sum(threshold == 1)); %Area in total pixels above % threshold using matrix
    umareaAboveThresh = (pxareaAboveThresh*(scaleval^2)); %Area in um2 above % threshold using matrix
    
    % added for ellipse
    [y_thresh, x_thresh] = find(contour);  % x and y are column vectors.
    ellipsefit_thresh = fit_ellipse(x_thresh,y_thresh);
    coord_thresh = [ellipsefit_thresh.X0_in, ellipsefit_thresh.Y0_in];
    contour = double(contour);
    CMap = [0,0,0; 0,1,0];
    contour  = ind2rgb(contour + 1, CMap);  
    
    % rotation matrix to rotate the axes with respect to an angle phi
    cos_phi = cos( ellipsefit_thresh.phi );
    sin_phi = sin( ellipsefit_thresh.phi );
    R = [ cos_phi sin_phi; -sin_phi cos_phi ];

    % the ellipse
    theta_r         = linspace(0,2*pi);
    ellipse_x_r     = ellipsefit_thresh.X0 + ellipsefit_thresh.a*cos( theta_r );
    ellipse_y_r     = ellipsefit_thresh.Y0 + ellipsefit_thresh.b*sin( theta_r );
    rotated_ellipse = R * [ellipse_x_r;ellipse_y_r];
    %unrotated_ellipse = [ellipse_x_r;ellipse_y_r];

    figure;
    imshow(contour);
    hold on;
    plot( rotated_ellipse(1,:),rotated_ellipse(2,:),'r' );
    %plot( unrotated_ellipse(1,:),unrotated_ellipse(2,:),'r' );
    plot(ellipsefit_thresh.X0_in, ellipsefit_thresh.Y0_in, '*b');
    plot(centroid_x, centroid_y, '*r');
    
    hold off;
    axis off;
    
    result_fname = [fnamelist{i} '_bestFitEllipse_'];
    f=getframe;
    imwrite(f.cdata, fullfile(basepath,[result_fname thresh_str '.tif']));
    
       
    % writing only the contour
    imwrite(contour, fullfile(basepath,[result_fname thresh_str '_only.tif']));
          
    %Joe's modification
    [y, x] = find(contour(:,:,2) == 1); %This seems to be correct
    coords = [x,y];
    writematrix(coords, fullfile(basepath,[result_fname 'contour_' thresh_str '.csv'])); %save coordinates of the percent contour
    
    %code to find densty at CDC
    ellipsefit_thresh.X0_rnd =  round(ellipsefit_thresh.X0_in);
    ellipsefit_thresh.Y0_rnd =  round(ellipsefit_thresh.Y0_in);
    densityatCDC = densitymap(ellipsefit_thresh.Y0_rnd, ellipsefit_thresh.X0_rnd);
    
    if (i == 1)
        data = [peak, centroid_x, centroid_y, pxareaAboveThresh, umareaAboveThresh, densityatCDC, ellipsefit_thresh.X0_rnd, ellipsefit_thresh.Y0_rnd, scaleval];
    else
        data = [data; peak, centroid_x, centroid_y,pxareaAboveThresh, umareaAboveThresh, densityatCDC, ellipsefit_thresh.X0_rnd, ellipsefit_thresh.Y0_rnd, scaleval];
    end
    
% Adding an output image with the marked location of peak density, added by Joe Carroll 2/19/22
vmap=viridis; %calls viridis colormap function, added by Joe 2/19/22
density_map_mark = densitymap-min(densitymap(:));
density_map_mark = uint8(255*density_map_mark./max(density_map_mark(:)));

MARK = insertShape(density_map_mark,'circle',[centroid_x centroid_y 2], 'LineWidth' ,3, 'Color' , 'red');
MARK = insertShape(MARK,'circle',[ellipsefit_thresh.X0_in ellipsefit_thresh.Y0_in 2], 'LineWidth' ,3, 'Color' , 'blue');
imwrite(MARK, vmap, fullfile(basepath,[result_fname 'marked.tif']));
end


% Write summary data to file
data = num2cell(data);
header = {'File Name', 'Peak', 'Centroid(max)_x', 'Centroid(max)_y',['PixelAreaAbove_0.' thresh_str], ['um2AreaAbove_0.' thresh_str], 'Density at CDC', 'EliCenter_0.8_x', 'EliCenter_0.8_y', 'um_per_pixel'};
EllipseCenterCoords = cat(2,fnamelist, data);
EllipseCenterCoords = cat(1, header, EllipseCenterCoords);
writecell(EllipseCenterCoords, fullfile(scalingpath, ['PCD_CDC_Analysis_Summary_', datestr(now, 'dd-mmm-yyyy'), '.csv']));

% Write al max value locations to file
unpacked_maxes = vertcat(all_maxes{:});
writecell(unpacked_maxes, fullfile(scalingpath, ['All_max_coords_' thresh_str '_percentile_' date '.csv']));



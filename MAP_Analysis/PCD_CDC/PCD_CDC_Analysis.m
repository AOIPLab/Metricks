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
    [maxrow,maxcol]=ind2sub(size(densitymap),maxind);       
    max_x = maxcol;
    max_y = maxrow;

    % check that the max value is unique
    max_indices = find(densitymap == maxval);           
    [max_y_coords, max_x_coords] = find(densitymap == maxval);
    all_max_coords = [max_x_coords, max_y_coords];
    for j = 1:length(max_x_coords)
        all_maxes{count,1} = {fnamelist{i}, all_max_coords(j,1), all_max_coords(j,2)};
        count = count +1;
    end

    % two ways of finding the centroid of the max values
    mean_x = mean(max_x_coords);
    mean_y = mean(max_y_coords);

    mid_x = (min(max_x_coords) + max(max_x_coords))/2;
    mid_y = (min(max_y_coords) + max(max_y_coords))/2;
    % unsure yet which one to go with
   
    threshold80 = (densitymap >= (0.8*peak));
    contour80 = edge(threshold80);
    
    %added for area of 80
    pxareaAboveThresh80 = sum(sum(threshold80 == 1)); %Area in total pixels above 80% threshold using matrix
    umareaAboveThresh80 = (pxareaAboveThresh80*(scaleval^2)); %Area in um2 above 80% threshold using matrix
    
    % added for ellipse
    [y80, x80] = find(contour80);  % x and y are column vectors.
    ellipsefit80 = fit_ellipse(x80,y80);
    coord80 = [ellipsefit80.X0_in, ellipsefit80.Y0_in];
    contour80 = double(contour80);
    CMap = [0,0,0; 0,1,0];
    contour80  = ind2rgb(contour80 + 1, CMap);  
    
    % rotation matrix to rotate the axes with respect to an angle phi
    cos_phi = cos( ellipsefit80.phi );
    sin_phi = sin( ellipsefit80.phi );
    R = [ cos_phi sin_phi; -sin_phi cos_phi ];

    % the ellipse
    theta_r         = linspace(0,2*pi);
    ellipse_x_r     = ellipsefit80.X0 + ellipsefit80.a*cos( theta_r );
    ellipse_y_r     = ellipsefit80.Y0 + ellipsefit80.b*sin( theta_r );
    rotated_ellipse = R * [ellipse_x_r;ellipse_y_r];
    %unrotated_ellipse = [ellipse_x_r;ellipse_y_r];

    figure;
    imshow(contour80);
    hold on;
    plot( rotated_ellipse(1,:),rotated_ellipse(2,:),'r' );
    %plot( unrotated_ellipse(1,:),unrotated_ellipse(2,:),'r' );
    plot(ellipsefit80.X0_in, ellipsefit80.Y0_in, '*b');
    plot(max_x, max_y, '*r');
    
    hold off;
    axis off;
    
    result_fname = [fnamelist{i} '_bestFitEllipse_'];
    f=getframe;
    imwrite(f.cdata, fullfile(basepath,[result_fname '80.tif']));
    
       
    % writing only the contour
    imwrite(contour80, fullfile(basepath,[result_fname '80_only.tif']));
          
    %Joe's modification
    [y, x] = find(contour80(:,:,2) == 1); %This seems to be correct
    coords = [x,y];
    writematrix(coords, fullfile(basepath,[result_fname 'contour_80.csv'])); %save coordinates of the 80 percent contour
    
    %code to find densty at CDC
    ellipsefit80.X0_rnd =  round(ellipsefit80.X0_in);
    ellipsefit80.Y0_rnd =  round(ellipsefit80.Y0_in);
    densityatCDC = densitymap(ellipsefit80.Y0_rnd, ellipsefit80.X0_rnd);
    
    if (i == 1)
        data = [peak, max_x, max_y, pxareaAboveThresh80, umareaAboveThresh80, densityatCDC, ellipsefit80.X0_rnd, ellipsefit80.Y0_rnd, scaleval];
    else
        data = [data; peak, max_x, max_y,pxareaAboveThresh80, umareaAboveThresh80, densityatCDC, ellipsefit80.X0_rnd, ellipsefit80.Y0_rnd, scaleval];
    end
    
% Adding an output image with the marked location of peak density, added by Joe Carroll 2/19/22
vmap=viridis; %calls viridis colormap function, added by Joe 2/19/22
density_map_mark = densitymap-min(densitymap(:));
density_map_mark = uint8(255*density_map_mark./max(density_map_mark(:)));

MARK = insertShape(density_map_mark,'circle',[max_x max_y 2], 'LineWidth' ,3, 'Color' , 'red');
MARK = insertShape(MARK,'circle',[ellipsefit80.X0_in ellipsefit80.Y0_in 2], 'LineWidth' ,3, 'Color' , 'blue');
imwrite(MARK, vmap, fullfile(basepath,[result_fname 'marked.tif']));
end


% Write summary data to file
data = num2cell(data);
header = {'File Name', 'Peak', 'Max_x', 'Max_y','PixelAreaAbove_0.80', 'um2AreaAbove_0.80', 'Density at CDC', 'EliCenter_0.8_x', 'EliCenter_0.8_y', 'um_per_pixel'};
EllipseCenterCoords = cat(2,fnamelist, data);
EllipseCenterCoords = cat(1, header, EllipseCenterCoords);
writecell(EllipseCenterCoords, fullfile(scalingpath, ['PCD_CDC_Analysis_Summary_', datestr(now, 'dd-mmm-yyyy'), '.csv']));

% Write al max value locations to file
unpacked_maxes = vertcat(all_maxes{:});
writecell(unpacked_maxes, fullfile(basepath, 'Results', ['all_max_coords_' date '.csv']));



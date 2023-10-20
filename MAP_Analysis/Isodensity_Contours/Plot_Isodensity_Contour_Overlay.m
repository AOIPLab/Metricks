%Creates coordinate files for density contours in AOSLO images. Created by
%JAC 9 May 2019
%Edited by N Wynne March-08-2021 to display multiple contours on the same
%Edited by J Grieshop March-12-2021 to automate file loading/saving and add colors to contours
%figure
%Edited by Joe Carroll to remove hard coded file folder save location
%NOTE, this requires a bounddensitymatrix file, NOT a cone coordinate file!

clear all
close all
clc

basePath = which('Plot_Isodensity_Contours_Overlay_v2.m');
[basePath] = fileparts(basePath);
path(path,fullfile(basePath,'lib3')); % Add our support library to the path.
[basepath] = uigetdir(pwd);
[fnamelist] = read_folder_contents(basepath,'csv');
[scalingfname, scalingpath] = uigetfile(fullfile(basepath,'*.csv'),'Select scaling LUT.');

for i = 1:length(fnamelist)
    
    densitymap = csvread(fnamelist{i});
    peak = max(densitymap(:));
    
    [maxval, maxind] = max(densitymap(:));  
    [maxrow,maxcol]=ind2sub(size(densitymap),maxind);       
    max_x = maxcol;
    max_y = maxrow;
    
    threshold75 = (densitymap >= (0.75*peak));
    contour75 = edge(threshold75); 
    % added for ellipse
    [y75, x75] = find(contour75);  % x and y are column vectors.
    ellipsefit75 = fit_ellipse(x75,y75);
    coord75 = [ellipsefit75.X0_in, ellipsefit75.Y0_in];
    %Uncomment to graph ellipse
    %graphEllipse('75', fnamelist{i},contour75, ellipsefit75);
    %comment out below lines to just do b&w
    contour75 = double(contour75);
    %red
    %CMap = [0,0,0; 0,0,1];
    %contour75  = ind2rgb(contour75 + 1, CMap);
    
    threshold80 = (densitymap >= (0.8*peak));
    contour80 = edge(threshold80);

    % added for ellipse
    [y80, x80] = find(contour80);  % x and y are column vectors.
    ellipsefit80 = fit_ellipse(x80,y80);
    coord80 = [ellipsefit80.X0_in, ellipsefit80.Y0_in];
    %Uncomment to graph ellipse
    %graphEllipse('80', fnamelist{i}, contour80, ellipsefit80)
    %comment out below lines to just do b&w
    contour80 = double(contour80);
    %orange
    %CMap = [0,0,0; 0,1,0];
    %contour80  = ind2rgb(contour80 + 1, CMap);
    
    
    threshold85 = (densitymap >= (0.85*peak));
    contour85 = edge(threshold85);
    % added for ellipse
    [y85, x85] = find(contour85);  % x and y are column vectors.
    ellipsefit85 = fit_ellipse(x85,y85);
    coord85 = [ellipsefit85.X0_in, ellipsefit85.Y0_in];
    %Uncomment to graph ellipse
    %graphEllipse('85', fnamelist{i}, contour85, ellipsefit85)
    %comment out below lines to just do b&w
    contour85 = double(contour85);
    %yellow
    %CMap = [0,0,0; 1,1,0];
    %contour85  = ind2rgb(contour85 + 1, CMap);
    
    
    threshold90 = (densitymap >= (0.9*peak));
    contour90 = edge(threshold90);
    % added for ellipse
    [y90, x90] = find(contour90);  % x and y are column vectors.
    ellipsefit90 = fit_ellipse(x90,y90);
    coord90 = [ellipsefit90.X0_in, ellipsefit90.Y0_in];
    %Uncomment to graph ellipse
    %graphEllipse('90', fnamelist{i}, contour90, ellipsefit90)
    %comment out below lines to just do b&w
    contour90 = double(contour90);
    %green
    %CMap = [0,0,0; 1,0.5,0];
    %contour90  = ind2rgb(contour90 + 1, CMap);
    
    
    threshold95 = (densitymap >= (0.95*peak));
    contour95 = edge(threshold95);
    % added for ellipse
    [y95, x95] = find(contour95);  % x and y are column vectors.
    ellipsefit95 = fit_ellipse(x95,y95);
    coord95 = [ellipsefit95.X0_in, ellipsefit95.Y0_in];
    %Uncomment to graph ellipse
    %graphEllipse('95', fnamelist{i}, contour95, ellipsefit95)
    %comment out below lines to just do b&w
    contour95 = double(contour95);
    %blue
    %CMap = [0,0,0; 1,0,0];
    %contour95  = ind2rgb(contour95 + 1, CMap);
    
    %switch which line is commented out to change to b&w
    %combined_contours = contour75 + contour80 + contour85 + contour90 + contour95;
    combined_contours = contour75 | contour80 | contour85 | contour90 | contour95;
    
    figure(1);
    imshow(combined_contours);
    saveas(gcf, fullfile(basepath, [fnamelist{i}(1:end-length('_bounddensity_matrix_xx-xxx-xxxx.csv')) '_combined_contours.tif']));
    
%   fixed x and y coordinates to be correct - JG
%   [y, x] = find(combined_contours == 1);
%   coords = [x,y];
    
    %result_fname = [fnamelist{i} '_combined_contours_'];

    %[y, x] = find(contour80(:,:,2) == 1); %This seems to be correct for RGB
    [y, x] = find(contour80 == 1); %This seems to be correct for B&W
    coords = [x,y];
    writematrix(coords, fullfile(basepath,[fnamelist{i}(1:end-length('_bounddensity_matrix_xx-xxx-xxxx.csv')) '_combined_contours.csv'])); %save combined contours
    
    if (i == 1)
        data = [peak, max_x, max_y, ellipsefit75.X0_in, ellipsefit75.Y0_in, ellipsefit80.X0_in, ellipsefit80.Y0_in, ellipsefit85.X0_in, ellipsefit85.Y0_in, ellipsefit90.X0_in, ellipsefit90.Y0_in, ellipsefit95.X0_in, ellipsefit95.Y0_in];
    else
        data = [data; peak, max_x, max_y, ellipsefit75.X0_in, ellipsefit75.Y0_in, ellipsefit80.X0_in, ellipsefit80.Y0_in, ellipsefit85.X0_in, ellipsefit85.Y0_in, ellipsefit90.X0_in, ellipsefit90.Y0_in, ellipsefit95.X0_in, ellipsefit95.Y0_in];
    end
end

data = num2cell(data);
header = {'File Name', 'Peak', 'Max_x', 'Max_y','EliCenter_0.75_x', 'EliCenter_0.75_y', 'EliCenter_0.8_x', 'EliCenter_0.8_y', 'EliCenter_0.85_x', 'EliCenter_0.85_y', 'EliCenter_0.9_x', 'EliCenter_0.9_y', 'EliCenter_0.95_x', 'EliCenter_0.95_y'};
EllipseCenterCoords = cat(2,fnamelist, data);
EllipseCenterCoords = cat(1, header, EllipseCenterCoords);
writecell(EllipseCenterCoords, fullfile(scalingpath, ['EllipseCenterCoords_', datestr(now, 'dd-mmm-yyyy'), '.csv']));



% Note this function has not been validated
function graphEllipse(percentile, filename, contour, ellipsefit)

 % rotation matrix to rotate the axes with respect to an angle phi
    cos_phi = cos( ellipsefit.phi );
    sin_phi = sin( ellipsefit.phi );
    R = [ cos_phi sin_phi; -sin_phi cos_phi ];

    % the ellipse
    theta_r         = linspace(0,2*pi);
    ellipse_x_r     = ellipsefit.X0 + ellipsefit.a*cos( theta_r );
    ellipse_y_r     = ellipsefit.Y0 + ellipsefit.b*sin( theta_r );
    rotated_ellipse = R * [ellipse_x_r;ellipse_y_r];
    %unrotated_ellipse = [ellipse_x_r;ellipse_y_r];

    figure;
    imshow(contour);
    hold on;
    plot( rotated_ellipse(1,:),rotated_ellipse(2,:),'r' );
    %plot( unrotated_ellipse(1,:),unrotated_ellipse(2,:),'r' );
    plot(ellipsefit.X0_in, ellipsefit.Y0_in, '*r');
    hold off;
    
    axis off;
    
    [filepath,name,ext] = fileparts(filename); %Doesn't work, need to update filenames and paths
    baseFileName = sprintf('%s.tiff', strcat(name, '_bestFitEllipse_', percentile));
    % Specify the specific folder:
    fullFileName = fullfile(pathnames, baseFileName);  
    f=getframe;
    imwrite(f.cdata,fullFileName)
    
    close all;
   
end
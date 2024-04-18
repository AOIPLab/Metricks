% Creates coordinate files for density contours in AOSLO images. Created by
% JAC 9 May 2019
% NOTE, this requires a bounddensitymatrix file, NOT a cone coordinate file!
%
% Input: bound density matrix .csv files from Metrick MAP output, same LUT
% table from Metricks.
%
% Output: overlay image for each subject, csv with contour coordinates for
% each subject, csv with summary data for all subjects.


clear all
close all
clc

% Add our support library to the path.
basePath = which('Plot_Isodensity_Contours_Overlay.m');
[basePath] = fileparts(basePath);
path(path,fullfile(basePath,'lib')); 

% User selects folder with data
dataPath = uigetdir('.','Select directory containing analyses');

% Read in csv names and then have user select the LUT
[fnameList] = read_folder_contents(dataPath,'csv');
[scalingFname, scalingPath] = uigetfile(fullfile(dataPath,'*.csv'),'Select scaling LUT.');

% Remove LUT file from fnameList
fnameList(ismember(fnameList,scalingFname))=[];

for i = 1:length(fnameList)
    
    densityMap = csvread(fullfile(dataPath, fnameList{i}));
    peak = max(densityMap(:));
    
    [maxVal, maxInd] = max(densityMap(:));  
    [maxRow,maxCol]=ind2sub(size(densityMap),maxInd);       
    maxX = maxCol;
    maxY = maxRow;
    
    threshold75 = (densityMap >= (0.75*peak));
    contour75 = edge(threshold75); 
    % added for ellipse
    [y75, x75] = find(contour75);  % x and y are column vectors.
    ellipseFit75 = fit_ellipse(x75,y75);
    coord75 = [ellipseFit75.X0_in, ellipseFit75.Y0_in];
    contour75 = double(contour75);
    
    threshold80 = (densityMap >= (0.8*peak));
    contour80 = edge(threshold80);
    % added for ellipse
    [y80, x80] = find(contour80);  % x and y are column vectors.
    ellipseFit80 = fit_ellipse(x80,y80);
    coord80 = [ellipseFit80.X0_in, ellipseFit80.Y0_in];
    contour80 = double(contour80);
        
    threshold85 = (densityMap >= (0.85*peak));
    contour85 = edge(threshold85);
    % added for ellipse
    [y85, x85] = find(contour85);  % x and y are column vectors.
    ellipseFit85 = fit_ellipse(x85,y85);
    coord85 = [ellipseFit85.X0_in, ellipseFit85.Y0_in];
    contour85 = double(contour85);
    
    threshold90 = (densityMap >= (0.9*peak));
    contour90 = edge(threshold90);
    % added for ellipse
    [y90, x90] = find(contour90);  % x and y are column vectors.
    ellipseFit90 = fit_ellipse(x90,y90);
    coord90 = [ellipseFit90.X0_in, ellipseFit90.Y0_in];
    contour90 = double(contour90);
    
    threshold95 = (densityMap >= (0.95*peak));
    contour95 = edge(threshold95);
    % added for ellipse
    [y95, x95] = find(contour95);  % x and y are column vectors.
    ellipseFit95 = fit_ellipse(x95,y95);
    coord95 = [ellipseFit95.X0_in, ellipseFit95.Y0_in];
    contour95 = double(contour95);
    
    % combine the contours
    combinedContours = contour75 | contour80 | contour85 | contour90 | contour95;
    
    % create and save the figure
    figure(1);
    imshow(combinedContours);
    saveas(gcf, fullfile(dataPath, [fnameList{i}(1:end-length('_bounddensity_matrix_xx-xxx-xxxx.csv')) '_combined_contours.tif']));

    % save matrix
    [y, x] = find(contour80 == 1);
    coords = [x,y];
    writematrix(coords, fullfile(dataPath,[fnameList{i}(1:end-length('_bounddensity_matrix_xx-xxx-xxxx.csv')) '_combined_contours.csv'])); %save combined contours
    
    % compile the data for each subject
    if (i == 1)
        data = [peak, maxX, maxY, ellipseFit75.X0_in, ellipseFit75.Y0_in, ellipseFit80.X0_in, ellipseFit80.Y0_in, ellipseFit85.X0_in, ellipseFit85.Y0_in, ellipseFit90.X0_in, ellipseFit90.Y0_in, ellipseFit95.X0_in, ellipseFit95.Y0_in];
    else
        data = [data; peak, maxX, maxY, ellipseFit75.X0_in, ellipseFit75.Y0_in, ellipseFit80.X0_in, ellipseFit80.Y0_in, ellipseFit85.X0_in, ellipseFit85.Y0_in, ellipseFit90.X0_in, ellipseFit90.Y0_in, ellipseFit95.X0_in, ellipseFit95.Y0_in];
    end
end

% combine the compiled data with the headers and save matrix
data = num2cell(data);
header = {'File Name', 'Peak', 'Max_x', 'Max_y','EliCenter_0.75_x', 'EliCenter_0.75_y', 'EliCenter_0.8_x', 'EliCenter_0.8_y', 'EliCenter_0.85_x', 'EliCenter_0.85_y', 'EliCenter_0.9_x', 'EliCenter_0.9_y', 'EliCenter_0.95_x', 'EliCenter_0.95_y'};
EllipseCenterCoords = cat(2,fnameList, data);
EllipseCenterCoords = cat(1, header, EllipseCenterCoords);
writecell(EllipseCenterCoords, fullfile(dataPath, ['EllipseCenterCoords_', datestr(now, datestr(now, 'dd_mmm_yyyy')), '.csv']));



% Note this function has not been validated - was used to trouble shoot
% during creation of this script. Not used now
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
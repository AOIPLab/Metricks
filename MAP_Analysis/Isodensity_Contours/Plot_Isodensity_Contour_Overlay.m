% Creates Isodensity Contour Overlay combined percentiles (75-95%) plots
% and coordinates
%
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

% Configure output header
header = {'File Name', 'Peak', 'Max_x', 'Max_y'};

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

% Threshold percentile selection by the user
list = {'75', '80', '85', '90', '95'};
[indx, tf] = listdlg('PromptString', 'Select the desired threshold percentile.', 'SelectionMode', 'multiple', 'ListString', list);

if tf ~= 1
    % Canceled dialog box - end the program
    return
end

for i = 1:length(fnameList)

    densityMap = csvread(fullfile(dataPath, fnameList{i}));
    peak = max(densityMap(:));
    
    [maxVal, maxInd] = max(densityMap(:));  
    [maxRow,maxCol]=ind2sub(size(densityMap),maxInd);       
    maxX = maxCol;
    maxY = maxRow;

    for j=1:length(indx)
        
        % Get the percentile
        thresholdPercentile = str2num(list{indx(j)})/100;
        threshStr{j} = list{indx(j)};


        threshold = (densityMap >= (thresholdPercentile*peak));
        contour = edge(threshold); 
        % Added for ellipse
        [y{j}, x{j}] = find(contour);  % x and y are column vectors.
        ellipseFit{j} = fit_ellipse(x{j},y{j});
        % coord = [ellipseFit{j}.X0_in, ellipseFit{j}.Y0_in];
        contour_cell{j} = double(contour);


        % Combine the contours
        if j == 1
            combinedContours = contour_cell{j};
            elipseCombined = [ellipseFit{j}.X0_in, ellipseFit{j}.Y0_in];
        else
            combinedContours = or(combinedContours, contour_cell{j});
            elipseCombined = [elipseCombined, ellipseFit{j}.X0_in, ellipseFit{j}.Y0_in];
        end
    end
      
    % Create and save the figure
    figure(1);
    imshow(combinedContours);
    saveas(gcf, fullfile(dataPath, [fnameList{i}(1:end-length('_bounddensity_matrix_xx-xxx-xxxx.csv')) '_combined_contours.tif']));

    % save coordinates - this is from a previous version; unsure about why
    % it was needed so commented out
    % [y, x] = find(contour80 == 1);
    % coords = [x,y];
    % writematrix(coords, fullfile(dataPath,[fnameList{i}(1:end-length('_bounddensity_matrix_xx-xxx-xxxx.csv')) '_combined_contours.csv'])); %save combined contours
    
    % Compile the data for each subject
    if (i == 1)
        data = [peak, maxX, maxY, elipseCombined];
    else
        data = [data; peak, maxX, maxY, elipseCombined];
    end
end

% Compile the rest of the header information
count = 4;
for j=1:length(indx)
    header{1,j+count} = sprintf('EliCenter_%s_x', threshStr{j});
    header{1,j+count+1} = sprintf('EliCenter_%s_y', threshStr{j});
    count = count+1;
end

% Combine the compiled data with the headers and save matrix
data = num2cell(data);
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
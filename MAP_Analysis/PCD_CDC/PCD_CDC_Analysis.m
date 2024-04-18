% Creates coordinate files for density contours in AOSLO images and returns PCD & CDC data 
% NOTE, this requires a bounddensitymatrix file, NOT a cone coordinate file!
%
% Input: bound density matrices (.csv) files from metrics MAP output. Also
% need the same LUT from the metrics MAP to use for this script.
% 
% Script will ask user for the folder where the data is stored (bound density matrices) and also to
% select the LUT file (LUT file can be in the same folder as the data or in
% a different folder). The script will ask wheat isodensity contour to
% report. 80% is what we typically use.
% 
% Outputs: individual contours and contours with markings images and csv, analysis
% summary csv, all max coords csv.
%
%Based on code created by JAC 9 May 2019

clear all
close all
clc

% Add our support library to the path.
basePath = which('PCD_CDC_Analysis.m');
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
list = {'95', '90', '85', '80', '75', '70', '65', '60', '55', '50', '45', '40', '35', '30', '25', '20', '15', '10', '5'};
[indx, tf] = listdlg('PromptString', 'Select the desired threshold percentile.', 'SelectionMode', 'single', 'ListString', list);
if tf == 1
    thresholdPercentile = str2num(list{indx})/100;
    threshStr = list{indx};
else
    % Canceled dialog box - end the program
    return
end

% Finds scale from the LUT or user enters it manually
scaleInput = NaN;
if scalingFname == 0        
    
    while isnan(scaleInput)                
        
        scaleInput = inputdlg('Input the scale in UNITS/PIXEL:','Input the scale in UNITS/PIXEL:');
        
        scaleInput = str2double(scaleInput);
        
        if isempty(scaleInput)
            error('Cancelled by user.');
        end
    end
else
    [~, lutData] = load_scaling_file(fullfile(scalingPath,scalingFname));
end


count = 1;
for i=1:size(fnameList,1) % Go through all files in list 
    
     if isnan(scaleInput)
                % Calculate the scale for this identifier.                                
                LUTindex=find( cellfun(@(s) ~isempty(strfind(fnameList{i},s )), lutData{1} ) );

                % Use whichever scale is most similar to our filename.
                sim = 1000*ones(length(LUTindex),1);
                for l=1:length(LUTindex)
                    sim(l) = lev(fnameList{i}, lutData{1}{LUTindex(l)});
                end
                [~,simInd]=min(sim);
                LUTindex = LUTindex(simInd);
                
                axialLength = lutData{2}(LUTindex);
                pixelsPerDegree = lutData{3}(LUTindex);

                micronsPerDegree = (291*axialLength)/24; % These numbers are from a book in Joe's office, Joe has confirmed it is correct
                
                scaleVal = 1 / (pixelsPerDegree / micronsPerDegree);
            else
                scaleVal = scaleInput;
    end

    
    % Load in the density map and find the max value (peak)
    densityMap = csvread(fullfile(dataPath,fnameList{i}));
    peak = max(densityMap(:));


    % Check that the max value/peak is unique         
    [maxYCoords, maxXCoords] = find(densityMap == peak);
    allMaxCoords = [maxXCoords, maxYCoords];
    for j = 1:length(maxXCoords)
        allMaxes{count,1} = {fnameList{i}, allMaxCoords(j,1), allMaxCoords(j,2)}; % store all maxes
        count = count +1;
    end

    % Finding the weighted (mean) centoid if multiple max locations found
    centroidX = mean(maxXCoords);
    centroidY = mean(maxYCoords);
   
    % Get threshold and contour
    threshold = (densityMap >= (thresholdPercentile*peak));
    contour = edge(threshold);
    
    % Added for area
    pxareaAboveThresh = sum(sum(threshold == 1)); %Area in total pixels above % threshold using matrix
    umareaAboveThresh = (pxareaAboveThresh*(scaleVal^2)); %Area in um2 above % threshold using matrix
    
    % Added for ellipse
    [yThresh, xThresh] = find(contour);  % x and y are column vectors.
    ellipsefitThresh = fit_ellipse(xThresh,yThresh);
    contour = double(contour);
    CMap = [0,0,0; 0,1,0];
    contour  = ind2rgb(contour + 1, CMap);  
    
    % Rotation matrix to rotate the axes with respect to an angle phi
    cosPhi = cos( ellipsefitThresh.phi );
    sinPhi = sin( ellipsefitThresh.phi );
    R = [ cosPhi sinPhi; -sinPhi cosPhi ];

    % The ellipse
    thetaR = linspace(0,2*pi);
    ellipseXR = ellipsefitThresh.X0 + ellipsefitThresh.a*cos( thetaR );
    ellipseYR = ellipsefitThresh.Y0 + ellipsefitThresh.b*sin( thetaR );
    rotatedEllipse = R * [ellipseXR;ellipseYR];

    % Create the plots
    figure;
    imshow(contour);
    hold on;
    plot( rotatedEllipse(1,:),rotatedEllipse(2,:),'r' );
    plot(ellipsefitThresh.X0_in, ellipsefitThresh.Y0_in, '*b');
    plot(centroidX, centroidY, '*r');
    
    hold off;
    axis off;
    
    % Writing the image
    resultFname = [fnameList{i} '_bestFitEllipse_'];
    f=getframe;
    imwrite(f.cdata, fullfile(dataPath,[resultFname threshStr '.tif']));
       
    % Writing only the contour
    imwrite(contour, fullfile(dataPath,[resultFname threshStr '_only.tif']));
          
    % Joe's modification
    [y, x] = find(contour(:,:,2) == 1); %This seems to be correct
    coords = [x,y];
    writematrix(coords, fullfile(dataPath,[resultFname 'contour_' threshStr '.csv'])); %save coordinates of the percent contour
    
    % Code to find densty at CDC
    ellipsefitThresh.X0_rnd =  round(ellipsefitThresh.X0_in);
    ellipsefitThresh.Y0_rnd =  round(ellipsefitThresh.Y0_in);
    densityatCDC = densityMap(ellipsefitThresh.Y0_rnd, ellipsefitThresh.X0_rnd);
    
    if (i == 1)
        data = [peak, centroidX, centroidY, pxareaAboveThresh, umareaAboveThresh, densityatCDC, ellipsefitThresh.X0_rnd, ellipsefitThresh.Y0_rnd, scaleVal];
    else
        data = [data; peak, centroidX, centroidY,pxareaAboveThresh, umareaAboveThresh, densityatCDC, ellipsefitThresh.X0_rnd, ellipsefitThresh.Y0_rnd, scaleVal];
    end
    
% Adding an output image with the marked location of peak density, added by Joe Carroll 2/19/22
vmap=viridis; %calls viridis colormap function, added by Joe 2/19/22
densityMapMark = densityMap-min(densityMap(:));
densityMapMark = uint8(255*densityMapMark./max(densityMapMark(:)));

MARK = insertShape(densityMapMark,'circle',[centroidX centroidY 2], 'LineWidth' ,3, 'Color' , 'red');
MARK = insertShape(MARK,'circle',[ellipsefitThresh.X0_in ellipsefitThresh.Y0_in 2], 'LineWidth' ,3, 'Color' , 'blue');
imwrite(MARK, vmap, fullfile(dataPath,[resultFname 'marked.tif']));

end


% Write summary data to file
data = num2cell(data);
header = {'File Name', 'Peak', 'Centroid(max)_x', 'Centroid(max)_y',['PixelAreaAbove_0.' threshStr], ['um2AreaAbove_0.' threshStr], 'Density at CDC', 'EliCenter_0.8_x', 'EliCenter_0.8_y', 'um_per_pixel'};
EllipseCenterCoords = cat(2,fnameList, data);
EllipseCenterCoords = cat(1, header, EllipseCenterCoords);
writecell(EllipseCenterCoords, fullfile(dataPath, ['PCD_CDC_Analysis_Summary_', datestr(now, 'dd_mmm_yyyy'), '.csv']));

% Write all max value locations to file
unpackedMaxes = vertcat(allMaxes{:});
writecell(unpackedMaxes, fullfile(dataPath, ['All_max_coords_' threshStr '_percentile_' datestr(now, 'dd_mmm_yyyy') '.csv']));



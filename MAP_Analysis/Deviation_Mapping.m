% Deviation Mapping
% AOIP
% Created by: Jenna Grieshop
% Date created: 7/18/24
%
% Description: This script compares metrics maps (e.g. density maps) to a
% normative averaged map to output a deviation map of the inputted maps for
% comparrison. Note: All maps including normative average & standard
% deviation maps must be the same scale. 
%
% Input: A folder containing: A LUT file with CDC coordinates for normative maps, Average
% normative map, standard deviation normative map (outputs from
% Density_Matrix_Averaging.m). Another folder containing Comparrison metric
% maps and LUT file with the CDC coordinates.
% The code will ask you for things in this order:
% 1. The normative averaged map
% 2. The directory of comparision maps
% 3. CDC LUT file for normative averaged map
% 4. CDC LUT file for comparison maps
% 5. Normative standard deviation map
%
% Outputs: Outputs are saved in the comparison metric maps folder. Raw
% data deviation map csv (for each), deviation map image and svg (for each),
% deviation percentage analysis csv (summary for all).
%


clear all;
close all;
clc;


basepath = which('Density_Matrix_Averaging.m');
[basepath] = fileparts(basepath);
path(path,fullfile(basepath,'lib')); % Add our support library to the path.


% User selects data
[normFname, normDataPath] = uigetfile('*.csv', 'Select Normative Average Map');
devDataPath = uigetdir('.','Select directory containing maps to analyze');

% Read in csv names and then have user select the LUT
[devFnameList] = read_folder_contents(devDataPath,'csv');

% Select and load in filename of the LUT with CDC
[normLUTfilename, normLUTpathname] = uigetfile('*.csv', 'Select file with Normative Map CDC coords', normDataPath);
[devLUTfilename, devLUTpathname] = uigetfile('*.csv', 'Select file with CDC coords from maps to analyze', devDataPath);

[stdfilename, stdpathname] = uigetfile('*.csv', 'Select Normative Stdev map', normDataPath);

% Remove LUT file from fnameList
devFnameList(ismember(devFnameList,devLUTfilename))=[];

% Load in the LUT file
[~, normLutData] = load_LUT_file(fullfile(normLUTpathname,normLUTfilename));
[~, devLutData] = load_LUT_file(fullfile(devLUTpathname,devLUTfilename));


devNumFiles = size(devFnameList,1);
 
%% Normative Data

% Find the information from the LUT file for the data                                
normLUTindex=find( cellfun(@(s) contains(normFname,s ), normLutData{1} ) );

normCDC_x = normLutData{2}(normLUTindex);
normCDC_y = normLutData{3}(normLUTindex);

% Load in the standard deviation and average normative maps
normAvg = readmatrix(fullfile(normDataPath,normFname));
normStdev = readmatrix(fullfile(stdpathname,stdfilename));

% Get the orignal coordinates for the corners of the matrices
l = size(normAvg);
ntl = [2,2];
ntr = [l(2)+1,2];
nbl = [2,l(1)+1];
nbr = [l(2)+1,l(1)+1];

% Get the adjusted coordinates for the corners of the matrix. Offset
% by the CDC coords
ntla = ntl-[normCDC_x,normCDC_y];
ntra = ntr-[normCDC_x,normCDC_y];
nbla = nbl-[normCDC_x,normCDC_y];
nbra = nbr-[normCDC_x,normCDC_y];

% Subtract from cdc coords too
normCDC_x = normCDC_x - normCDC_x+1;
normCDC_y = normCDC_y - normCDC_y+1;

narray{1} = ntla(1);
narray{2} = ntla(2);
narray{3} = ntra(1);
narray{4} = ntra(2);
narray{5} = nbla(1);
narray{6} = nbla(2);
narray{7} = nbra(1);
narray{8} = nbra(2);

narray = cell2mat(narray);

%% Compare maps
for j=1:devNumFiles % Go through all files in list 

    % Find the information from the LUT file for the data                                
    devLUTindex=find( cellfun(@(s) ~isempty(strfind(devFnameList{j},s )), devLutData{1} ) );
    
    devCDC_x{j} = devLutData{2}(devLUTindex);
    devCDC_y{j} = devLutData{3}(devLUTindex);
    devdata{j} = readmatrix(fullfile(devDataPath,devFnameList{j}));

    % Flip OS maps to match the OD maps - Average Maps are in OD
    % orientation
    if contains(devFnameList{j}, 'OS')
        devdata{j} = fliplr(devdata{j});
        devCDC_x{j} = length(devdata{j})-(devCDC_x{j}-1);
    end

    % Get the orignal coordinates for the corners of the matrices
    l = size(devdata{j});
    dtl{j} = [2,2];
    dtr{j} = [l(2)+1,2];
    dbl{j} = [2,l(1)+1];
    dbr{j} = [l(2)+1,l(1)+1];

    % Get the adjusted coordinates for the corners of the matrix. Offset
    % by the CDC coords
    dtla{j} = dtl{j}-[devCDC_x{j},devCDC_y{j}];
    dtra{j} = dtr{j}-[devCDC_x{j},devCDC_y{j}];
    dbla{j} = dbl{j}-[devCDC_x{j},devCDC_y{j}];
    dbra{j} = dbr{j}-[devCDC_x{j},devCDC_y{j}];

    % Subtract from cdc coords too
    devCDC_x{j} = devCDC_x{j} - devCDC_x{j}+1;
    devCDC_y{j} = devCDC_y{j} - devCDC_y{j}+1;

    darray{j,1} = dtla{j}(1);
    darray{j,2} = dtla{j}(2);
    darray{j,3} = dtra{j}(1);
    darray{j,4} = dtra{j}(2);
    darray{j,5} = dbla{j}(1);
    darray{j,6} = dbla{j}(2);
    darray{j,7} = dbra{j}(1);
    darray{j,8} = dbra{j}(2);


end

darray = cell2mat(darray);

count = 1;
for n = 1:devNumFiles

    % Find the minimum of all the coordinates
    minimumx = min(narray(1), darray(n,1));
    minimumy = min(narray(2), darray(n,2));

    if narray(1) < darray(n,1)
        minIx = 1;
    else
        minIx = 2;
    end

    if narray(2) < darray(n,2)
        minIy = 1;
    else
        minIy = 2;
    end
    
    maximumx = max(narray(7), darray(n,7));
    maximumy = max(narray(8), darray(n,8));

    if narray(7) > darray(n,7)
        maxIx = 1;
    else
        maxIx = 2;
    end

    if narray(8) > darray(n,8)
        maxIy = 1;
    else
        maxIy = 2;
    end
    
    
    offset = abs([minimumx-1, minimumy-1]);

    % Adjust coordinates by the minimum by adding the offset
    tlao{1} = ntla + offset;
    trao{1} = ntra + offset;
    blao{1} = nbla + offset;
    brao{1} = nbra + offset;
    
    nCDC_x = normCDC_x + offset(1);
    nCDC_y = normCDC_y + offset(2);

    % Adjust coordinates by the minimum by adding the offset
    tlao{2} = dtla{n} + offset;
    trao{2} = dtra{n} + offset;
    blao{2} = dbla{n} + offset;
    brao{2} = dbra{n} + offset;     
    
    dCDC_x = devCDC_x{n} + offset(1);
    dCDC_y = devCDC_y{n} + offset(2);

    

    % Figure out how much padding is needed for normative
    nleft = abs(tlao{minIx}(1)-tlao{1}(1));
    ntop = abs(tlao{minIy}(2)-tlao{1}(2));
    nright = abs(brao{maxIx}(1)-brao{1}(1));
    nbottom = abs(brao{maxIy}(2)-brao{1}(2));

    npaddedData = normAvg;
    nstdpaddedData = normStdev;

    nleftPad = zeros(size(npaddedData,1), nleft);
    npaddedData = horzcat(nleftPad, npaddedData);
    nstdpaddedData = horzcat(nleftPad, nstdpaddedData);

    ntopPad = zeros(ntop, size(npaddedData,2));
    npaddedData = vertcat(ntopPad, npaddedData);
    nstdpaddedData = vertcat(ntopPad, nstdpaddedData);

    nrightPad = zeros(size(npaddedData,1), nright(1));
    npaddedData = horzcat(npaddedData, nrightPad);
    nstdpaddedData = horzcat(nstdpaddedData, nrightPad);

    nbottomPad = zeros(nbottom, size(npaddedData, 2));
    npaddedData = vertcat(npaddedData, nbottomPad);
    nstdpaddedData = vertcat(nstdpaddedData, nbottomPad);
  
    % Change the padded portion to be Nan
    npaddedData(npaddedData==0) = NaN;
    nstdpaddedData(nstdpaddedData==0) = NaN;

    A(:,:,1) = npaddedData;


    % Figure out how much padding is needed for data

    dleft = abs(tlao{minIx}(1)-tlao{2}(1));
    dtop = abs(tlao{minIy}(2)-tlao{2}(2));
    dright = abs(brao{maxIx}(1)-brao{2}(1));
    dbottom = abs(brao{maxIy}(2)-brao{2}(2));

    dpaddedData = devdata{n};

    dleftPad = zeros(size(dpaddedData,1), dleft);
    dpaddedData = horzcat(dleftPad, dpaddedData);

    dtopPad = zeros(dtop, size(dpaddedData,2));
    dpaddedData = vertcat(dtopPad, dpaddedData);

    drightPad = zeros(size(dpaddedData,1), dright(1));
    dpaddedData = horzcat(dpaddedData, drightPad);

    dbottomPad = zeros(dbottom, size(dpaddedData, 2));
    dpaddedData = vertcat(dpaddedData, dbottomPad);
  
    % Change the padded portion to be Nan
    dpaddedData(dpaddedData==0) = NaN;

    A(:,:,2) = dpaddedData;

    CDC = [nCDC_x, nCDC_y];

    % Get rid of columns that include any NaNs

    for j=1:size(A,1)
        for g = 1:size(A,2)
            NanCount = sum(isnan(A(j,g,:)));
            if NanCount > 0
                % If less than N get rid of the spot for the N map one
                A(j,g,:) = NaN;
                nstdpaddedData(j,g) = NaN;
            end
        end
    end


    % Figure out how much will need to adjust CDC coords for new maps once full
    % NaN rows and columns deleted

    for s=1:size(A,2)
        if sum(A(:,s,1), "omitnan") > 0
            break;
        else
            CDC(1) = CDC(1)- 1;
        end
    end

    for t=1:size(A,1)
        if sum(A(t,:,1), "omitnan") > 0
            break;
        else
            CDC(2) = CDC(2)- 1;
        end
    end

    normativeAverage = A(:,:,1);
    comparisonData =  A(:,:,2);

    normativeAverage = normativeAverage(:,~all(isnan(normativeAverage))); 
    normativeAverage = normativeAverage(~all(isnan(normativeAverage),2), :); 

    comparisonData = comparisonData(:,~all(isnan(comparisonData))); 
    comparisonData = comparisonData(~all(isnan(comparisonData),2), :); 


    normativeStdev = nstdpaddedData(:,~all(isnan(nstdpaddedData))); 
    normativeStdev = normativeStdev(~all(isnan(normativeStdev),2), :); 




% TODO might need to force the matrices to be a square. Center the CDC and
% then get rid of the excess

    nstdev{1} = normativeStdev;
    pstd{1} = normativeAverage + nstdev{1};
    mstd{1} = normativeAverage - nstdev{1};

    nstdev{2} = normativeStdev.*2;
    pstd{2} = normativeAverage + nstdev{2};
    mstd{2} = normativeAverage - nstdev{2};

    % Additional standard deviation options
    % nstdev{3} = normativeStdev.*3;
    % pstd{3} = normativeAverage + nstdev{3};
    % mstd{3} = normativeAverage - nstdev{3};
    % 
    % nstdev{4} = normativeStdev.*4;
    % pstd{4} = normativeAverage + nstdev{4};
    % mstd{4} = normativeAverage - nstdev{4};


    for q=1:size(normativeStdev, 1)
        for u=1:size(normativeStdev,2)
            % Additional standard deviation options
            % if comparisonData(q,u) > pstd{4}(q,u) || comparisonData(q,u) < mstd{4}(q,u)
            %     resultMap(q,u) = 4;
            % elseif comparisonData(q,u) > pstd{3}(q,u) || comparisonData(q,u) < mstd{3}(q,u)
            %     resultMap(q,u) = 3;
            if comparisonData(q,u) <= pstd{1}(q,u) && comparisonData(q,u) >= mstd{1}(q,u)
                resultMap(q,u) = 1; % within +-1 stdev
            elseif comparisonData(q,u) <= pstd{2}(q,u) && comparisonData(q,u) >= mstd{2}(q,u)
                resultMap(q,u) = 2; % within +-2 stdev
            else
                resultMap(q,u) = 3; % more than +-2 stdev
            end
        end
    end

    % clims for the devaition map
    clims = [1 3];
    % Color scheme for the deviaiton map
    vmap = [0 1 0 
        1 1 0 
        1 0 0];

    % Determine what regions of the map are within the standard deviations
    std1 = sum(resultMap(:) == 1)/numel(resultMap) * 100;
    std2 = sum(resultMap(:) == 2)/numel(resultMap) * 100;
    std2p = sum(resultMap(:) == 3)/numel(resultMap) * 100;

    if (n == 1)
        data = [std1, std2, std2p];
    else
        data = [data; std1, std2, std2p];
    end

    % Plot the map
    vmap=viridis;
    dispfig=figure(count); 
    imagesc(resultMap, clims); % Added to use limits of color scale
    axis image;
    colormap(vmap); 
    colorbar;

    count = count + 1;


    % Save deviation map
    saveas(gcf,fullfile(devLUTpathname, ['Deviation_' devFnameList{n} '_' datestr(now, 'yyyymmdd') '_fig.png']));
    saveas(gcf,fullfile(devLUTpathname, ['Deviation_' devFnameList{n} '_' datestr(now, 'yyyymmdd') '_fig.svg']));
    writematrix(resultMap, fullfile(devLUTpathname,['Deviation_' devFnameList{n} '_' datestr(now, 'yyyymmdd') '_raw.csv']));


    % Clear variables to not have interference in next iteration
    clear A;
    clear normativeAverage;
    clear comparisonData;
    clear resultMap;
    clear nstdev;
    clear pstd;
    clear mstd;



end


% Write standard deviation percentages to file
data = num2cell(data);
header = {'File Name', '% within +-1 Stdev', '% within +-2 Stdev', '% greater than +-2 Stdev'};
StdevResults = cat(2,devFnameList, data);
StdevResults = cat(1, header, StdevResults);
writecell(StdevResults, fullfile(devDataPath, ['Deviation_Percentage_Analysis_Summary_', datestr(now, 'dd_mmm_yyyy'), '.csv']));










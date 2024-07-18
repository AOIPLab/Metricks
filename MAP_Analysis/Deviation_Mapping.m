% Deviation Mapping
% 
% Created by: Jenna Grieshop
% Date created: 7/18/24
%
%

clear all;
close all;
clc;


basepath = which('Density_Matrix_Averaging.m');
[basepath] = fileparts(basepath);
path(path,fullfile(basepath,'lib')); % Add our support library to the path.



% User selects folder with data
normDataPath = uigetdir('.','Select directory containing normative maps');
devDataPath = uigetdir('.','Select directory containing maps to analyze');


% Read in normative csv names and then have user select the LUT
[normFnameList] = read_folder_contents(normDataPath,'csv');

% Read in csv names and then have user select the LUT
[devFnameList] = read_folder_contents(devDataPath,'csv');

% select and load in filename of the LUT with CDC
[normLUTfilename, normLUTpathname] = uigetfile('*.csv', 'Select file with Normative Map CDC coords');
[devLUTfilename, devLUTpathname] = uigetfile('*.csv', 'Select file with CDC coords from maps to analyze');

[stdfilename, stdpathname] = uigetfile('*.csv', 'Select Normative Stdev map');

% Remove LUT file from fnameList
normFnameList(ismember(normFnameList,normLUTfilename))=[];
devFnameList(ismember(devFnameList,devLUTfilename))=[];

% load in the LUT file
[~, normLutData] = load_LUT_file(fullfile(normLUTpathname,normLUTfilename));
[~, devLutData] = load_LUT_file(fullfile(devLUTpathname,devLUTfilename));


normNumFiles = size(normFnameList,1);
normNumFilesDim = size(normFnameList);

devNumFiles = size(devFnameList,1);
devNumFilesDim = size(devFnameList);


% Get normative data sorted and ready
for i=1:normNumFiles % Go through all files in list 

    % Find the information from the LUT file for the data                                
    normLUTindex=find( cellfun(@(s) ~isempty(strfind(normFnameList{i},s )), normLutData{1} ) );
    
    normCDC_x{i} = normLutData{2}(normLUTindex);
    normCDC_y{i} = normLutData{3}(normLUTindex);
    normdata{i} = readmatrix(fullfile(normDataPath,normFnameList{i}));

    % figure out the sizes of the matricies
    normSz{i} = size(normdata{i},1);

    % get the orignal coordinates for the corners of the matrices
    l = size(normdata{i});
    ntl{i} = [2,2];
    ntr{i} = [l(2)+1,2];
    nbl{i} = [2,l(1)+1];
    nbr{i} = [l(2)+1,l(1)+1];

    % get the adjusted coordinates for the corners of the matrix. Offset
    % by the CDC coords
    ntla{i} = ntl{i}-[normCDC_x{i},normCDC_y{i}];
    ntra{i} = ntr{i}-[normCDC_x{i},normCDC_y{i}];
    nbla{i} = nbl{i}-[normCDC_x{i},normCDC_y{i}];
    nbra{i} = nbr{i}-[normCDC_x{i},normCDC_y{i}];

    %subtract from cdc coords too
    normCDC_x{i} = normCDC_x{i} - normCDC_x{i}+1;
    normCDC_y{i} = normCDC_y{i} - normCDC_y{i}+1;

    narray{i,1} = ntla{i}(1);
    narray{i,2} = ntla{i}(2);
    narray{i,3} = ntra{i}(1);
    narray{i,4} = ntra{i}(2);
    narray{i,5} = nbla{i}(1);
    narray{i,6} = nbla{i}(2);
    narray{i,7} = nbra{i}(1);
    narray{i,8} = nbra{i}(2);


end

narray = cell2mat(narray);

% Get data sorted and ready
for j=1:devNumFiles % Go through all files in list 

    % Find the information from the LUT file for the data                                
    devLUTindex=find( cellfun(@(s) ~isempty(strfind(devFnameList{j},s )), devLutData{1} ) );
    
    devCDC_x{j} = devLutData{2}(devLUTindex);
    devCDC_y{j} = devLutData{3}(devLUTindex);
    devdata{j} = readmatrix(fullfile(devDataPath,devFnameList{j}));

    % figure out the sizes of the matricies
    devSz{j} = size(devdata{j},1);

    % get the orignal coordinates for the corners of the matrices
    l = size(devdata{j});
    dtl{j} = [2,2];
    dtr{j} = [l(2)+1,2];
    dbl{j} = [2,l(1)+1];
    dbr{j} = [l(2)+1,l(1)+1];

    % get the adjusted coordinates for the corners of the matrix. Offset
    % by the CDC coords
    dtla{j} = dtl{j}-[devCDC_x{j},devCDC_y{j}];
    dtra{j} = dtr{j}-[devCDC_x{j},devCDC_y{j}];
    dbla{j} = dbl{j}-[devCDC_x{j},devCDC_y{j}];
    dbra{j} = dbr{j}-[devCDC_x{j},devCDC_y{j}];

    %subtract from cdc coords too
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
for m = 1:normNumFiles
    for n = 1:devNumFiles

        % find the minimum of all the coordinates
        minimumx = min(narray(m,1), darray(n,1));
        minimumy = min(narray(m,2), darray(n,2));

        if narray(m,1) < darray(n,1)
            minIx = 1;
        else
            minIx = 2;
        end

        if narray(m,2) < darray(n,2)
            minIy = 1;
        else
            minIy = 2;
        end
        
        maximumx = max(narray(m,7), darray(n,7));
        maximumy = max(narray(m,8), darray(n,8));

        if narray(m,7) > darray(n,7)
            maxIx = 1;
        else
            maxIx = 2;
        end

        if narray(m,8) > darray(n,8)
            maxIy = 1;
        else
            maxIy = 2;
        end
        
        
        offset = abs([minimumx-1, minimumy-1]);

        % adjust coordinates by the minimum by adding the offset
        tlao{1} = ntla{m} + offset;
        trao{1} = ntra{m} + offset;
        blao{1} = nbla{m} + offset;
        brao{1} = nbra{m} + offset;
        
        nCDC_x = normCDC_x{m} + offset(1);
        nCDC_y = normCDC_y{m} + offset(2);

        % adjust coordinates by the minimum by adding the offset
        tlao{2} = dtla{n} + offset;
        trao{2} = dtra{n} + offset;
        blao{2} = dbla{n} + offset;
        brao{2} = dbra{n} + offset;     
        
        dCDC_x = devCDC_x{n} + offset(1);
        dCDC_y = devCDC_y{n} + offset(2);

        % load in the standard deviation normative map
        nstdev = readmatrix(fullfile(stdpathname,stdfilename));

        % figure out how much padding is needed for normative
        nleft = abs(tlao{minIx}(1)-tlao{1}(1));
        ntop = abs(tlao{minIy}(2)-tlao{1}(2));
        nright = abs(brao{maxIx}(1)-brao{1}(1));
        nbottom = abs(brao{maxIy}(2)-brao{1}(2));
    
        npaddedData = normdata{m};
        nstdpaddedData = nstdev;
    
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
      
        % change the padded portion to be Nan
        npaddedData(npaddedData==0) = NaN;
        nstdpaddedData(nstdpaddedData==0) = NaN;

        A(:,:,1) = npaddedData;


        % figure out how much padding is needed for data

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
      
        % change the padded portion to be Nan
        dpaddedData(dpaddedData==0) = NaN;

        A(:,:,2) = dpaddedData;

        CDC = [nCDC_x, nCDC_y];

        %get rid of columns that include any NaNs

        for j=1:size(A,1)
            for g = 1:size(A,2)
                NanCount = sum(isnan(A(j,g,:)));
                if NanCount > 0
                    % if less than N get rid of the spot for the N map one
                    A(j,g,:) = NaN;
                    nstdpaddedData(j,g) = NaN;
                end
            end
        end


        % % figure out how much will need to adjust CDC coords for new maps once full
        % % NaN rows and columns deleted

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

        B(:,:,1) = A(:,~all(isnan(A(:,:,1)))); 
        C(:,:,1) = B(~all(isnan(B(:,:,1)),2), :); 

        D(:,:,2) = A(:,~all(isnan(A(:,:,2))),2); 
        E(:,:,2) = D(~all(isnan(D(:,:,2)),2), :,2); 

        clear A;
        A(:,:,1) = C(:,:,1);
        A(:,:,2) =  E(:,:,2);

        nstdev = nstdpaddedData(:,~all(isnan(nstdpaddedData))); 
        nstdev = nstdev(~all(isnan(nstdev),2), :); 

        % JG TODO check that cdc values are still in the right spot

        nstdev_2 = nstdev.*2;
        p2std = A(:,:,1) + nstdev_2;
        m2std = A(:,:,1) - nstdev_2;

        p1std = A(:,:,1) + nstdev;
        m1std = A(:,:,1) - nstdev;


        for q=1:size(nstdev_2, 1)
            for u=1:size(nstdev_2,1)
                if A(q,u,2) > p2std(q,u) || A(q,u,2) < m2std(q,u)
                    A(q,u,3) = 2;
                elseif A(q,u,2) > p1std(q,u) || A(q,u,2) < m1std(q,u)
                    A(q,u,3) = 1;
                else
                     A(q,u,3) = 0;
                end
            end
        end


        vmap=viridis;
        dispfig=figure(count); 
        imagesc(A(:,:,3)); % added to use limits of color scale
        axis image;
        colormap(vmap); 
        colorbar;

        count = count + 1;

        saveas(gcf,fullfile(devLUTpathname, ['Stdev_' devFnameList{n} '_' datestr(now, 'yyyymmdd') '_fig.png']));

        writematrix(A(:,:,3), fullfile(devLUTpathname,['stdev_' devFnameList{n} '_' datestr(now, 'yyyymmdd') '_raw.csv']));
        


        clear A;
        clear B;
        clear C;
        clear D;

    end
end











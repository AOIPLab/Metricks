% Author: Jenna Grieshop
% Date of creation: 2/22/22
%
% Description: Script that performs subtraction between two density
% matricies. They do not need to be the same size but they do need to be
% scaled to one another. 
%
% Input: Smaller matrix first, then larger matrix, then LUT file
% (PCD_CDC_Alalysis_Summary with only the information of the matrices being
% used & in the same order as selected).
%
% Output: csv, SVG, and tif of the resulting subtraction matrix. Note that
% svgs may only be readable in the original location they were saved to.
%
% Method: CDCs between teh two matrices are aligned and the smaller matrix
% is padded appropriately to be the same size as the larger one. The
% subtraction is performed and the excess padding is removed. The resulting
% matrix is the size of the smaller matrix.

clear all;
close all;
clc;

basepath = which('Density_Matrix_Subtraction.m');
[basepath] = fileparts(basepath);
path(path,fullfile(basepath,'lib')); % Add our support library to the path.


% select and load in first matrix name
[filename1, pathname1] = uigetfile('*.csv', 'MultiSelect', 'off');
% select and load in second matrix name
[filename2, pathname2] = uigetfile('*.csv', 'MultiSelect', 'off');

% load in the data
data{1} = readmatrix(fullfile(pathname1,filename1));
data{2} = readmatrix(fullfile(pathname2,filename2));


% select and load in filename of the LUT with CDC
[LUTfilename, LUTpathname] = uigetfile('*.csv', 'Select file with CDC coords');

% load in the LUT file
LUT = readtable(fullfile(LUTpathname, LUTfilename));

% figure out the sizes of the matricies
sz{1} = size(data{1},1);
sz{2} = size(data{2},1);

% get CDC coords out of table
x{1} = LUT{1,8}; 
x{2} = LUT{2,8};
y{1} = LUT{1,9};
y{2} = LUT{2,9};


% get the orignal coordinates for the corners of the second matrix
l2 = size(data{2});
tl2 = [2,2];
tr2 = [l2(2)+1,2];
bl2 = [2,l2(1)+1];
br2 = [l2(2)+1,l2(1)+1];

% get the adjusted coordinates for the corners of the second matrix. Offset
% by the CDC coords
tl2a = tl2-[x{2},y{2}];
tr2a = tr2-[x{2},y{2}];
bl2a = bl2-[x{2},y{2}];
br2a = br2-[x{2},y{2}];
%subtract from cdc coords too
x{2} = x{2} - x{2}+1;
y{2} = y{2} - y{2}+1;


% get the orignal coordinates for the corners of the first matrix
l = size(data{1});
tl1 = [2,2];
tr1 = [l(2)+1,2];
bl1 = [2,l(1)+1];
br1 = [l(2)+1,l(1)+1];

% get the adjusted coordinates for the corners of the first matrix. Offset
% by the CDC coords
tl1a = tl1-[x{1},y{1}];
tr1a = tr1-[x{1},y{1}];
bl1a = bl1-[x{1},y{1}];
br1a = br1-[x{1},y{1}];
%subtract from cdc coords too
x{1} = x{1} - x{1}+1;
y{1} = y{1} - y{1}+1;

% put the coordinates in an array
array = [tl2a; tr2a; bl2a; br2a; tl1a; tr1a; bl1a; br1a];

% find the minimum of all the coordinates
[minimumx, indexx] = min(array(:,1));
[minimumy, indexy] = min(array(:,2));

if indexx ~= indexy
    print("yo they don't match");
else
    offset = abs([minimumx-1, minimumy-1]);
end

% adjust all coordinates by the minimum by adding the offset
tl2a2 = tl2a + offset;
tr2a2 = tr2a + offset;
bl2a2 = bl2a + offset;
br2a2 = br2a + offset;

tl1a2 = tl1a + offset;
tr1a2 = tr1a + offset;
bl1a2 = bl1a + offset;
br1a2 = br1a + offset;

x{2} = x{2} + offset(1);
y{2} = y{2} + offset(2);

x{1} = x{1} + offset(1);
y{1} = y{1} + offset(2);


% determine which is the larger matrix
if sz{2} > sz{1}
    larger = 2;
    smaller = 1;
else
    larger = 1;
    smaller = 2;
end

% figure out how much padding is needed

front_top = abs(tl2a2-tl1a2);
back_bottom = abs(br2a2-br1a2);

paddedData = data{smaller};

% add padding to the top
top = zeros(front_top(2), size(paddedData,2));
paddedData = vertcat(top, paddedData);

% add padding to the front
front = zeros(size(paddedData,1), front_top(1));
paddedData = horzcat(front, paddedData);

% add padding to the bottom
bottom = zeros(back_bottom(2), size(paddedData, 2));
paddedData = vertcat(paddedData, bottom);

% add padding to the back
back = zeros(size(paddedData,1), back_bottom(1));
paddedData = horzcat(paddedData, back);

% change the padded portion to be Nan
paddedData(paddedData==0) = NaN;

data{smaller} = paddedData;

% subtract, but check equaland provide error message to user
if size(data{1}) == size(data{2})

    % perform subtraction
    resultMatrix = data{larger} - data{smaller};
    resultMatrix = resultMatrix(:,~all(isnan(resultMatrix))); % get rid of columns with only nans
    resultMatrix = resultMatrix(~all(isnan(resultMatrix),2), :); % get rid of row with only nans
    
    % remove extension from name
    [folder1, baseName1, extension1] = fileparts(filename1);
    [folder2, baseName2, extension2] = fileparts(filename2);
    
    % save averaged matrix
    n1 = split(baseName1, '_0p');
    n2 = split(baseName2, '_0p');
    filename = [n2{1} '_MINUS_' n1{1} '_' date];
    csvwrite(fullfile(LUTpathname, [filename '.csv']),resultMatrix);
    
    % save difference image
    imwrite(resultMatrix, parula(256),fullfile(LUTpathname, [filename '.tif']));
    
    % display difference map as a 3D plot
    vis = mesh(resultMatrix);
    view(0,90);
    set(gca, 'Visible', 'on')
    f = gcf;
    exportgraphics(f,fullfile(LUTpathname, [filename '_2.tif']),'Resolution',300)
    
    % save as svg for figure making
    matrix = resultMatrix;
    f = figure('visible', 'off');
    colormap(parula);
    image(matrix, 'CDataMapping', 'scaled');
    print(f, '-dsvg', fullfile(LUTpathname, [filename '.svg']));


else
    disp 'Matrices not same size'
end

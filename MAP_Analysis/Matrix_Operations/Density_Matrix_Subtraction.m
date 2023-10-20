% Author: Jenna Grieshop
% Date of creation: 2/22/22, v0

% Description: Script that performs element-wise subtraction on two density
% matrices. Matrices MUST be same size. First selected matrix-second selected matrix. 
% The resulting matrix is saved as an output.

% Description of edits (include date & author):
% Added difference image save and 3D display (2/22/22, Joe Carroll), v1

clear all;
close all;
clc;

% select and load in first matrix name
[filename1, pathname1] = uigetfile('*.csv', 'MultiSelect', 'off');
% select and load in second matrix name
[filename2, pathname2] = uigetfile('*.csv', 'MultiSelect', 'off');

% load in the data
data1 = readmatrix(filename1);
data2 = readmatrix(filename2);

% subtract, but check equaland provide error message to user
if size(data1) == size(data2)

% perform subtraction
resultMatrix = data1 - data2;

% remove extension from name
[folder1, baseName1, extension1] = fileparts(filename1);
[folder2, baseName2, extension2] = fileparts(filename2);

% save averaged matrix
filename = [baseName1 '_MINUS_' baseName2 '_' date '.csv'];
csvwrite(filename,resultMatrix);

% save difference image
imwrite(resultMatrix, parula(256),[filename '.tif']);

%display difference map asa 3D plot
mesh(resultMatrix)
colormap(parula(256))

else
    disp 'Matrices not same size'
end

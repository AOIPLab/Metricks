% Author: Jenna Grieshop
% Date of creation: 7/3/2024
%
% Description: Script that scales Foveal ROI to a set scale 0.25mpp (image
% size of 2400x2400 for 600um, 2000x2000 for 500um, 1200x1200 for 300um)
%
% Input: Folder containing foveal ROI images, image name must be in either
% of these formats:
% JC_XXXXX_date_OD/S_XpXXXXmpp_XXXpXXXppd_XXXum_date_jc.tif
% or XXXXX_date_OD/S_XpXXXXmpp_XXXpXXXppd_XXXum_date_jc.tif
%
% Output: scaled .tifs with new scales embedded in the file names with
% _scaled appended at end of name, csv that is the LUT file to use for the 
% coordinate_scaling.m script


clear all
close all
clc


basePath = which('Image_Scaling.m');
[basePath] = fileparts(basePath);
path(path,fullfile(basePath,'lib')); 

% User selects folder with data
dataPath = uigetdir('.','Select directory containing original images');

%User selects results folder
resultPath = uigetdir('.','Select directory to save scaled images in');

% Read in tif names
[fnameList] = read_folder_contents(dataPath,'tif');


for i=1:size(fnameList,1)
    spl = split(fnameList{i}, '_');
    if length(spl{1}) == 2
        og_mpp = split(spl{5},"mpp");
        og_mpp = str2double(strrep(og_mpp{1}, "p", "."));
        og_ppd = split(spl{6},"ppd");
        og_ppd = str2double(strrep(og_ppd{1}, "p", "."));
        roi_size = split(spl{7}, "um");
        roi_size = str2double(roi_size{1});
        file_desc = [spl{1} "_" spl{2} "_" spl{3} "_" spl{4}];
        file_desc = strrep(strjoin(file_desc), ' ' , '');
        spl{9} = strrep(spl{9}, '.tif', '');
        % else to account for other naming conventions other than AOIP
        % Beware of AOIP data that starts with initials other than JC, will
        % need to update if above
    else 
        og_mpp = split(spl{4},"mpp");
        og_mpp = str2double(strrep(og_mpp{1}, "p", "."));
        og_ppd = split(spl{5},"ppd");
        og_ppd = str2double(strrep(og_ppd{1}, "p", "."));
        roi_size = split(spl{6}, "um");
        roi_size = str2double(roi_size{1});
        file_desc = [spl{1} "_" spl{2} "_" spl{3}];
        file_desc = strrep(strjoin(file_desc), ' ' , '');
        spl{8} = strrep(spl{8}, '.tif', '');
    end
    

    if roi_size == 300
        new_dim = 1200;
        new_center = 600;
    elseif roi_size == 500
        new_dim = 2000;
        new_center = 1000; 
    elseif roi_size == 600
        new_dim = 2400;
        new_center = 1200;
    end

    % Load in image
    og_image = imread(fullfile(dataPath,fnameList{i}));
    og_dim = size(og_image, 1);
    og_center = og_dim/2;

    % calculate scale factor and image size in microns
    scale_factor = new_dim/og_dim;
    im_size_um = og_mpp * og_dim;

    % resize image and get new scales
    new_image = imresize(og_image, [new_dim new_dim]);
    new_ppd = og_ppd * scale_factor;
    new_mpp = im_size_um/new_dim;

    % save scaled image
    % format new scales to be correct for file name
    str_mpp = num2str(new_mpp);
    str_mpp = strrep(str_mpp, ".", "p");

    str_ppd = num2str(new_ppd);
    str_ppd = strrep(str_ppd, ".", "p");
    
    if length(spl{1}) == 2
        new_name = [spl{1}, "_" spl{2} "_" spl{3} "_" spl{4} "_" str_mpp "mpp_" str_ppd "ppd_" spl{7} "_" spl{8} "_" spl{9} "_scaled.tif"];
    else
        new_name = [spl{1}, "_" spl{2} "_" spl{3} "_" str_mpp "mpp_" str_ppd "ppd_" spl{6} "_" spl{7} "_" spl{8} "_scaled.tif"];
    end
    new_name = strjoin(new_name);
    new_name = strrep(new_name,' ','');

    imwrite(new_image, fullfile(resultPath, new_name));
    

    % compile data
    if i ==1 
        data = [file_desc, og_center, scale_factor, new_center, new_mpp, new_ppd];
    else
        data = [data; file_desc, og_center, scale_factor, new_center, new_mpp, new_ppd];
    end


end

data = num2cell(data);
writecell(data, fullfile(resultPath, "LUT_for_coordinate_scaling.csv"));







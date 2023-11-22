% Finds density at a specific distance away from the PCD in a density
% matrix. User enters the x and y distance in um (microns) from the PCD. 
% Negative distances are to the left and up.

% Inputs: Denisty matrices in the folder the code is running from, a LUT file
% (file name, axial length, and ppd) in a separate folder, user input about
% distance in x and y from the PCD

% Outputs: .csv saved to folder containing the LUT file.

% Created by Jenna Grieshop 11/22/23 - most code was recycled from
% PCD_CDC_Analysis.m code.


clear all
close all
clc

basepath = which('Density_at_Distance_from_PCD.m');
[basepath] = fileparts(basepath);
path(path,fullfile(basepath,'lib')); % Add our support library to the path.

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

x_input = NaN;
while isnan(x_input)

    x_input = inputdlg('Input the X distance from the PCD in um:', 'Input the X distance from the PCD in um:');
    
    x_input = str2double(x_input);

    if isempty(x_input)
            error('X input Cancelled by user.');
    end
end

y_input = NaN;
while isnan(y_input)
    
    y_input = inputdlg('Input the Y distance from the PCD in um:', 'Input the Y distance from the PCD in um:');
    
    y_input = str2double(y_input);
    if isempty(y_input)
            error('Y input cancelled by user.');
    end
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
                
                scaleval = 1 / (pixelsperdegree / micronsperdegree); %pixels/um
            else
                scaleval = scaleinput;
    end

    
    densitymap = csvread(fnamelist{i});
    peak = max(densitymap(:));
    
    [maxval, maxind] = max(densitymap(:));

    % check that the max value is unique
    max_indices = find(densitymap == maxval);           
    [max_y_coords, max_x_coords] = find(densitymap == maxval);
    all_max_coords = [max_x_coords, max_y_coords];
    for j = 1:length(max_x_coords)
        all_maxes{count,1} = {fnamelist{i}, all_max_coords(j,1), all_max_coords(j,2)};
        count = count +1;
    end

    % finding the weighted (mean) centoid if multiple max locations found
    centroid_x = mean(max_x_coords);
    centroid_y = mean(max_y_coords);

    x_dist = round(x_input * scaleval);
    y_dist = round(y_input * scaleval);

    new_x = centroid_x + x_dist;
    new_y = centroid_y + y_dist;

    d_at_point = densitymap(new_y, new_x) % x and y appear flipped bc of row and column rules in matlab
    
    if (i ==1)
        data = [d_at_point];
    else
        data = [data; d_at_point];
    end


end

data = num2cell(data);
header = {'File Name', 'Density at Point'};
compiled = cat(2, fnamelist, data);
compiled2 = cat(1, header, compiled);
writecell(compiled2, fullfile(scalingpath, ['Density_at_point_', num2str(x_input), 'um_', num2str(y_input), 'um_from_PCD_', datestr(now, 'dd-mmm-yyyy'), '.csv']))


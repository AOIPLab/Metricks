function Radial_Density_Batch(csv_file, lut_file, output_folder)
    % Radial Density - run through Python for batch processing
    % Inputs from Python:
    %   csv_file: path to the CSV to process
    %   lut_file: path to the LUT file
    %   output_folder: folder to save results
    %
    % This function allows batch processing.
    
    scaleFactor = 0.25; %umpp
    spacing = 1; %5 microns
    
    basepath = which('Radial_Density_Batch.m');
    [basepath] = fileparts(basepath);
    path(path,fullfile(basepath,'lib')); % Add support library
    
    disp(['Processing file: ', csv_file]);
    disp(['Using LUT file: ', lut_file]);
    
    % Load matrix
    data = readmatrix(csv_file);
    
    % Load LUT file
    [~, lutData] = load_LUT_file(lut_file);
    
    % Extract CDC info from LUT
    [~, filename_only, ~] = fileparts(csv_file);
    
    % Extract identifier: first two parts of filename
    underscore_idx = strfind(filename_only, '_');
    if numel(underscore_idx) >= 2
        identifier = filename_only(1:underscore_idx(2)-1);
    else
        identifier = filename_only;
    end
    disp(['Identifier: ', identifier]);
    
    LUTindex = find(cellfun(@(s) contains(s, identifier), lutData{1}));
    
    if isempty(LUTindex)
        warning('No match found in LUT for identifier: %s', identifier);
        return; % skip this file safely
    end
    
    CDC_x = lutData{2}(LUTindex);
    CDC_y = lutData{3}(LUTindex);
    CDC_density = data(CDC_y, CDC_x);  
    CDC = [CDC_y, CDC_x];

    
    if CDC_x ~= CDC_y
        warning('CDC is not centered in matrix');
    end
    
    px_num = (1/scaleFactor)*spacing;
    dimension = size(data, 1);
    points = floor((dimension/px_num)/2);
    
    for i = 1:points
        [X,Y] = meshgrid(1:size(data,1), 1:size(data,2));
        disk_locations = sqrt((X-CDC(1)).^2 + (Y-CDC(2)).^2) <= (px_num*i);
        outerboundary = imdilate(disk_locations, strel('disk',1)) & ~disk_locations;
        outerboundary = outerboundary';
        
        density_ring = data(outerboundary == 1);
        Avg(i) = mean(density_ring);
        eccentricity(i) = i;
    end
    
    Avg = Avg';
    densityveccentricity = [num2cell(eccentricity)', num2cell(Avg)];
    
    % Save results
    result_fname = fullfile(output_folder, [filename_only, '_', num2str(spacing), '_um_spacing_radial_density.csv']);
    header = {"Eccentricity", filename_only};
    header2 = {0, CDC_density};
    final_results = [header; header2; densityveccentricity];
    
    writecell(final_results, result_fname);
    
    disp(['Saved results to: ', result_fname]);
end
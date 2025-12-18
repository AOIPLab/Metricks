%% Batch_Run_Horizontal_Vertical_Density
% Extract horizontal and vertical density profiles relative to CDC location
clear; clc;

maltese = 0;
straight = 0;

pixel_size = 0.25; % assume 0.25 umpp

basePath = which('Batch_Run_Horizontal_Vertical_Density.m');

[basePath ] = fileparts(basePath);
path(path,fullfile(basePath,'lib')); % Add our support library to the path.

% --- User selects data directory, LUT file, output location ---
data_folder   = uigetdir('', 'Select folder containing density CSV files');
if data_folder == 0, error('No data folder selected'); end

[lut_file, lut_path] = uigetfile('*.csv', 'Select LUT CSV file with CDC coordinates');
if lut_file == 0, error('No LUT file selected'); end
lut_file = fullfile(lut_path, lut_file);

output_base = uigetdir('', 'Select folder to save results');
if output_base == 0, error('No output folder selected'); else output_root = fullfile(output_base, '\Results'); end %JG added else statement

%% Make sure results directory exists
if ~exist(output_root, 'dir')
    mkdir(output_root);
end

% Have user select what kind of row/column analysis they'd like to do
list = {'Straight Cross', 'Maltese Cross'};
[indx, tf] = listdlg('PromptString', 'Select the following row/column analyses to perform.', 'SelectionMode', 'multiple', 'ListString', list);
if tf == 1
    selections = size(indx,2);
    if selections == 2
        maltese = 1;
        straight = 1;
    else
        indx_max = max(indx);
        if indx_max == 1
            straight = 1;
        else
            maltese = 1;
        end
    end
    
else
    % Canceled dialog box - end the program
    return
end


%% ---- Load LUT ----
[hdr, lutData] = load_LUT_file_HV_Density(lut_file);

% lutData = {SubjectID cell, CDCx double, CDCy double}
subjectIDs = lutData{1};
CDCx_all   = lutData{2};
CDCy_all   = lutData{3};
eye_all    = lutData{4};

%% ---- Get all CSV files ----
files = dir(fullfile(data_folder, '*_bounddensity_matrix_*.csv'));

fprintf('Found %d files.\n', numel(files));

% Have the user select the thickness of the row/column extraction
if straight == 1
    m1 = 'Please enter the desired straight extraction thickness (# of rows/columns) MUST be an odd number:';
    thickness = inputdlg(m1);
    thickness = str2double(thickness{1});
    while mod(thickness, 2) ~= 1 
        warndlg('WARNING: Must enter odd number for thickness.')
        thickness = inputdlg(m1);
        thickness = str2double(thickness{1});
    end
end

if maltese ==1
    % Have the user select the thickness of the row/column extraction
    m2 = 'Please enter the desired maltese cross angle (degrees):';
    angle = inputdlg(m2);
    angle = str2double(angle{1});
end

%% ---- Loop through files ----
for i = 1:numel(files)
    fname = files(i).name;
    fprintf('\nProcessing %d/%d: %s\n', i, numel(files), fname);

    % Extract subject ID (before first two underscores)
    % e.g., AB_12636_bounddensity_matrix_* → AB_12636
    tokens = regexp(fname, '^[A-Z]{2}_[0-9]{4,5}', 'match');
    if isempty(tokens)
        warning('Skipping file — subject ID not found: %s', fname);
        continue;
    end
    subj = tokens{1};
    disp(['Looking for subject: ', subj]);
    
    % find match
    LUTindex = find(strcmp(subjectIDs, subj));
    
    if isempty(LUTindex)
        % display all LUT IDs for debugging
        disp('Available LUT IDs in LUT:');
        for k = 1:length(subjectIDs)
            disp(['"' subjectIDs{k} '"']);
        end
        warning('Subject %s not found in LUT. Skipping...', subj);
        continue;
    end

    % Pull CDC coords
    CDC_x = CDCx_all(LUTindex);
    CDC_y = CDCy_all(LUTindex);

    % Pull Eye
    eye = eye_all{LUTindex};

    % Load CSV image
    filepath = fullfile(files(i).folder, fname);
    M = readmatrix(filepath);

    if straight == 1

        straight_output = straight_cross(M, CDC_x, CDC_y, thickness, pixel_size, eye);

        %% ---- Save output ----
        % unit is microns per pixel
        % create the output file name
        outname = strcat(subj,'_', eye, '_Straight_Density_', string(datetime('now','TimeZone','local','Format','yyyyMMdd_hhmmss')), '.csv');
        outfile = fullfile(output_root, outname);
        writecell(straight_output, outfile);

    end
    if maltese ==1

        [distance_density_mean, degree_center_points_full] = Extract_radial_densities_one(M, [CDC_x, CDC_y], angle);
        % now need to save the density curves from 0, 90, 180, 270

        right = distance_density_mean{degree_center_points_full == 0};
        top = distance_density_mean{degree_center_points_full == 90};
        left = distance_density_mean{degree_center_points_full == 180};
        bottom = distance_density_mean{degree_center_points_full == 270};

        % Flip left and right if OS so that matches orientation of OD
         % ======= Orientation correction =======
        if strcmpi(eye, 'OS') 
            temp = left;
            left = right;
            right = temp;
        end

        ecc_l = (0:size(left,2)-1) * pixel_size;
        ecc_r = (0:size(right,2)-1) * pixel_size;
        ecc_t = (0:size(top,2)-1) * pixel_size;
        ecc_b = (0:size(bottom,2)-1) * pixel_size;

         % find minimum eccentricity for the average
        if size(ecc_r, 2) > size(ecc_l,2)
            h_ecc_avg = ecc_l;
        else
            h_ecc_avg = ecc_r;
        end

         % find minimum eccentricity for the average
        if size(ecc_b, 2) > size(ecc_t,2)
            v_ecc_avg = ecc_t;
        else
            v_ecc_avg = ecc_b;
        end

         v_dens{1} = top;
         v_dens{2} = bottom;
         h_dens{1} = left;
         h_dens{2} = right;

         % prepare for vertical averaging
        maxNumCol_v = max(cellfun(@(c) size(c,2), v_dens));  % max number of columns
        combined_v = cell2mat(cellfun(@(c){padarray(c,[0,maxNumCol_v-size(c,2)],NaN,'Post')}, v_dens)');

        % prepare for horizontal averaging
        maxNumCol_h = max(cellfun(@(c) size(c,2), h_dens));  % max number of columns
        combined_h = cell2mat(cellfun(@(c){padarray(c,[0,maxNumCol_h-size(c,2)],NaN,'Post')}, h_dens)');

        % agerage the vertical and horizontal densities
        avg_v_dens = rmmissing(mean(combined_v, 1));
        avg_h_dens = rmmissing(mean(combined_h, 1));
        
        results = {ecc_l, h_dens{1}', ecc_r, h_dens{2}', ecc_t, v_dens{1}', ecc_b, v_dens{2}', h_ecc_avg, avg_h_dens', v_ecc_avg, avg_v_dens'};
        headers = {'Ecc_Horz_L_um', 'Density_Horz_L', 'Ecc_Horz_R_um', 'Density_Horz_R', 'Ecc_Vert_T_um,', 'Density_Vert_T', 'Ecc_Vert_B_um', 'Density_Vert_B', 'Ecc_Avg_Horz_um', 'Density_Avg_Horz', 'Ecc_Avg_Vert_um', 'Density_Avg_Vert'};

        D = {};
        for k = 1:numel(results)
            D(1:numel(results{k}),k)=num2cell(results{k});
        end
        maltese_output = [headers; D];

        outname = strcat(subj,'_', eye, '_Maltese_Density_', string(datetime('now','TimeZone','local','Format','yyyyMMdd_hhmmss')), '.csv');
        outfile = fullfile(output_root, outname);
        writecell(maltese_output, outfile);


        
    end

 
end

fprintf('\nBatch complete.\n');
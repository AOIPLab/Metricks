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
        % Have the user select the thickness of the row/column extraction
        m1 = 'Please enter the desired row/column extraction thickness (# of rows/columns) MUST be an odd number:';
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

    %% Extract meridian profiles
    % ======= Extract vertical profile =======
    dens_vert = M(:, CDC_x);

    % Pixel offsets (top = negative, bottom = positive)
    pix_v = ( (1:length(dens_vert))' - CDC_y );
    ecc_vert = pix_v * pixel_size;   % microns

    % ======= Extract horizontal profile =======
    dens_horz = M(CDC_y, :);

    % Pixel offsets (left = negative, right = positive)
    pix_h = ( (1:length(dens_horz)) - CDC_x );
    ecc_horz = pix_h(:) * pixel_size;   % microns

    % ======= Orientation correction =======
    if strcmpi(eye, 'OS')  
        dens_horz = fliplr(dens_horz); % flip horizontally for OS to match OD
        ecc_horz = -flipud(ecc_horz);
    end



    %% split data at the cdc point

    % find the 0 point in the eccentricity
    v_0_index = find(ecc_vert == 0); 
    h_0_index = find(ecc_horz == 0); 

    % split the vertical data into top and bottom
    v_ecc_t = -flipud(ecc_vert(1:v_0_index));
    v_dens{1} = flipud(dens_vert(1:v_0_index))'; % top
    v_ecc_b = ecc_vert(v_0_index:end);
    v_dens{2} = dens_vert(v_0_index:end)'; % bottom
    % find minimum eccentricity for the average
    if size(v_ecc_b, 1) > size(v_ecc_t,1)
        v_ecc_avg = v_ecc_t;
    else
        v_ecc_avg = v_ecc_b;
    end

    % split the horizontal data into left and right
    h_ecc_l = -flipud(ecc_horz(1:h_0_index));
    h_dens{1} = fliplr(dens_horz(1:h_0_index)); % left
    h_ecc_r = ecc_horz(h_0_index:end);
    h_dens{2} = dens_horz(h_0_index:end); % right
    % find minimum eccentricity for the average
    if size(h_ecc_r, 1) > size(h_ecc_l,1)
        h_ecc_avg = h_ecc_l;
    else
        h_ecc_avg = h_ecc_r;
    end


    % prepare for vertical averaging
    maxNumCol_v = max(cellfun(@(c) size(c,2), v_dens));  % max number of columns
    combined_v = cell2mat(cellfun(@(c){padarray(c,[0,maxNumCol_v-size(c,2)],NaN,'Post')}, v_dens)');

    % prepare for horizontal averaging
    maxNumCol_h = max(cellfun(@(c) size(c,2), h_dens));  % max number of columns
    combined_h = cell2mat(cellfun(@(c){padarray(c,[0,maxNumCol_h-size(c,2)],NaN,'Post')}, h_dens)');

    % agerage the vertical and horizontal densities
    avg_v_dens = rmmissing(mean(combined_v, 1));
    avg_h_dens = rmmissing(mean(combined_h, 1));


    % compile results and create headers for output sheet
    results = {h_ecc_l, h_dens{1}', h_ecc_r, h_dens{2}', v_ecc_t, v_dens{1}', v_ecc_b, v_dens{2}', h_ecc_avg, avg_h_dens', v_ecc_avg, avg_v_dens'};
    headers = {'Ecc_Horz_L_um', 'Density_Horz_L', 'Ecc_Horz_R_um', 'Density_Horz_R', 'Ecc_Vert_T_um,', 'Density_Vert_T', 'Ecc_Vert_B_um', 'Density_Vert_B', 'Ecc_Avg_Horz_um', 'Density_Avg_Horz', 'Ecc_Avg_Vert_um', 'Density_Avg_Vert'};

    %% Plotting for sanity check
    plot(v_ecc_t, v_dens{1});
    hold on
    plot(v_ecc_b, v_dens{2});
    plot(h_ecc_l, h_dens{1});
    plot(h_ecc_r, h_dens{2});
    % 
    % plot(v_ecc_avg, avg_v_dens);
    % plot(h_ecc_avg, avg_h_dens);

    %% ---- Save output ----

    % create the output file name
    outname = strcat(subj,'_', eye, '_HorizVertDensity_', string(datetime('now','TimeZone','local','Format','yyyyMMdd_hhmmss')), '.csv');
    outfile = fullfile(output_root, outname);

    % organize all the data into a format that can be saved to the csv
    D = {};
    for k = 1:numel(results)
      D(1:numel(results{k}),k)=num2cell(results{k});
    end
    combined_output = [headers; D];
    writecell(combined_output, outfile);
 
end

fprintf('\nBatch complete.\n');
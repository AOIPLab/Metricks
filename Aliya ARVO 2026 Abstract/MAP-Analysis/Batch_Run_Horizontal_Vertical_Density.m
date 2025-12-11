%% Batch_Run_Horizontal_Vertical_Density
% Extract horizontal and vertical density profiles relative to CDC location
clear; clc;

% % --- User selects folders and LUT file ---
% data_folder   = uigetdir('', 'Select folder containing density CSV files');
% if data_folder == 0, error('No data folder selected'); end
% 
% [lut_file, lut_path] = uigetfile('*.csv', 'Select LUT CSV file with CDC coordinates');
% if lut_file == 0, error('No LUT file selected'); end
% lut_file = fullfile(lut_path, lut_file);
% 
% output_base = uigetdir('', 'Select folder to save results');
% if output_base == 0, error('No output folder selected'); end


% %% ---- ALIYA INPUTS ----
data_folder = 'C:\Users\ne61322\Box\PRO-23898\WIP_Siddiqui_FovealCone\1-Data and Analysis\5-Analysis';
lut_file    = 'C:\Users\ne61322\Box\PRO-23898\WIP_Siddiqui_FovealCone\1-Data and Analysis\5-Analysis\CDC_20251125.csv';
output_root = 'C:\Users\ne61322\Box\PRO-23898\WIP_Siddiqui_FovealCone\1-Data and Analysis\5-Analysis\Results';

%% Make sure results directory exists
if ~exist(output_root, 'dir')
    mkdir(output_root);
end

%% ---- Load LUT ----
[hdr, lutData] = load_LUT_file(lut_file);

% lutData = {SubjectID cell, CDCx double, CDCy double}
subjectIDs = lutData{1};
CDCx_all   = lutData{2};
CDCy_all   = lutData{3};

%% ---- Get all CSV files ----
files = dir(fullfile(data_folder, '*_bound_density_map_*.csv'));

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

    % Load CSV image
    filepath = fullfile(files(i).folder, fname);
    M = readmatrix(filepath);

    % ---- Call density extraction function ----
    [ dens_horiz, dens_vert, ecc_horiz, ecc_vert ] = ...
        extract_meridian_profiles(M, CDC_x, CDC_y, 0.25, 'OD');
    % default eye is OD; assume our scale of 0.25 µm per pixel

    %% ---- Save output ----
    outname = sprintf('%s_HorizVertDensity.csv', subj);
    outfile = fullfile(output_root, outname);
    
    orientation = [repmat({'Horizontal'}, numel(ecc_horiz), 1); ...
                   repmat({'Vertical'},   numel(ecc_vert), 1)];
    
    ecc_all  = [ecc_horiz(:); ecc_vert(:)];
    dens_all = [dens_horiz(:); dens_vert(:)];
    
    T = table(orientation, ecc_all, dens_all, ...
        'VariableNames', {'Orientation','Eccentricity_um','Density'});
    
    writetable(T, outfile);
    fprintf('Saved: %s\n', outfile);
end
fprintf('\nBatch complete.\n');
%% Radial Density Extraction for Multiple Subjects
%
% THIS CODE IS A WORK IN PROGRESS. DO NOT USE.
%
% Inputs:
%   - Directory containing raw CSV files of images
%   - Lookup table (CSV or XLSX) with Subject IDs and CDC coordinates
%
% Outputs:
%   - Per-subject radial density CSV
%   - All-subject radial density CSV
%

clear all; clc;

% Select directory containing raw CSV analyses
root_path_raw = uigetdir('.', 'Select directory containing density analyses');

% Select lookup table with CDC coordinates
[filename_lut, pathname_lut] = uigetfile({'*.csv;*.xlsx'}, 'Select LUT with CDC coordinates');
[~, ~, ext] = fileparts(filename_lut);
if strcmp(ext, '.csv')
    lut_table = readtable(fullfile(pathname_lut, filename_lut));
elseif strcmp(ext, '.xlsx')
    lut_table = readtable(fullfile(pathname_lut, filename_lut));
else
    error('Unsupported LUT format');
end

% User sets spacing and window in microns
spacing = str2double(inputdlg('Enter desired SPACING (um)'));
window  = str2double(inputdlg('Enter desired WINDOW (um)'));

% Get all CSV file paths in the directory
csv_files = dir(fullfile(root_path_raw, '*.csv'));
num_subjects = length(csv_files);

% Initialize storage for all-subject radial densities
all_radial_data = cell(num_subjects,1);

%% Loop through subjects
for i = 1:num_subjects
    % Load CSV
    data_csv = readmatrix(fullfile(csv_files(i).folder, csv_files(i).name));
    
    % Get subject ID from filename (truncate after second underscore)
    [~, fname, ~] = fileparts(csv_files(i).name);  % remove folder and extension
    parts = strsplit(fname, '_');                   % split filename at underscores
    subject_id = strjoin(parts(1:2), '_');         % keep only first two parts

    % print to verify
    fprintf('Processing file: %s → subject_id: %s\n', csv_files(i).name, subject_id);

    
    % Find CDC coordinates in LUT
    idx = find(strcmp(lut_table.SubjectID, subject_id));
    if isempty(idx)
        warning('Subject %s not found in LUT. Skipping.', subject_id);
        continue;
    end
    
    cdc_x = lut_table.CDC_X(idx);
    cdc_y = lut_table.CDC_Y(idx);
    
    % Convert to radial coordinates
    [X, Y] = meshgrid(1:size(data_csv,2), 1:size(data_csv,1));
    R = sqrt((X - cdc_x).^2 + (Y - cdc_y).^2);
    
    % Flatten for binning
    R_flat = R(:);
    values_flat = data_csv(:);
    
    % Define bins
    r_max = max(R_flat);
    bin_centers = 0:spacing:r_max;
    radial_avg = zeros(length(bin_centers),1);
    
    % Bin the values
    for b = 1:length(bin_centers)
        in_bin = values_flat(R_flat >= (bin_centers(b) - window/2) & R_flat < (bin_centers(b) + window/2));
        if ~isempty(in_bin)
            radial_avg(b) = mean(in_bin);
        else
            radial_avg(b) = NaN; % handle empty bins
        end
    end
    
    % Save per-subject CSV
    T = table(bin_centers', radial_avg, 'VariableNames', {'Radius_um', 'Density'});
    writetable(T, fullfile(root_path_raw, strcat(subject_id, '_RadialDensity.csv')));
    
    % Store for all-subject aggregation
    all_radial_data{i} = radial_avg;
end

%% Combine all subjects
min_length = min(cellfun(@length, all_radial_data));
radial_matrix = NaN(num_subjects, min_length);

for i = 1:num_subjects
    radial_matrix(i,:) = all_radial_data{i}(1:min_length)';
end

% Compute average across subjects
avg_radial = mean(radial_matrix,1, 'omitnan');

% Save combined CSV
T_all = table(bin_centers(1:min_length)', avg_radial', 'VariableNames', {'Radius_um', 'AverageDensity'});
writetable(T_all, fullfile(root_path_raw, 'AllSubjects_RadialDensity.csv'));


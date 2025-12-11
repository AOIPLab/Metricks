%% Batch Run Radial Density with GUI file selection

% Pick CSV folder
csv_folder = uigetdir(pwd, 'Select folder containing CSV files');
if csv_folder == 0
    error('No folder selected. Exiting.');
end

% Pick LUT file
[lut_filename, lut_pathname] = uigetfile({'*.csv','CSV files (*.csv)'}, 'Select LUT file', pwd);
if lut_filename == 0
    error('No LUT file selected. Exiting.');
end
lut_file = fullfile(lut_pathname, lut_filename);

% Output folder
output_folder = fullfile(csv_folder, 'RadialDensity_Results');
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Get list of CSV files
csv_files = dir(fullfile(csv_folder, '*.csv'));

% Loop through each CSV and call your function
for k = 1:length(csv_files)
    csv_path = fullfile(csv_folder, csv_files(k).name);
    disp(['Processing: ', csv_files(k).name]);
    
    try
        Radial_Density_Batch(csv_path, lut_file, output_folder);
    catch ME
        warning('Error processing file %s: %s', csv_files(k).name, ME.message);
    end
end

disp('Batch processing complete.');

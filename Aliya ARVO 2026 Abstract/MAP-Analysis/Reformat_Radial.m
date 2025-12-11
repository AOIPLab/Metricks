%Reformat outputs from Batch Radial Script

%load CSVs
folderPath = uigetdir(pwd, 'Select folder with CSV files');
csvFiles = dir(fullfile(folderPath, '*.csv'));

masterTable = table();

for k = 1:length(csvFiles)

    fileName = csvFiles(k).name;

    underscores = strfind(fileName, '_');
    if length(underscores) >= 2
        subjectID = fileName(1:underscores(2)-1);
    else
        subjectID = fileName(1:end-4); % fallback
    end
    
    T = readtable(fullfile(folderPath, fileName));
    
    T.Properties.VariableNames{2} = subjectID;
    
    if isempty(masterTable)
 
        masterTable = T;
    else
        % Merge by Eccentricity
        masterTable = outerjoin(masterTable, T, 'Keys', 'Eccentricity', 'MergeKeys', true);
    end
    
    fprintf('Processed %d/%d: %s\n', k, length(csvFiles), fileName);
end

% Save tables
writetable(masterTable, fullfile(folderPath, 'Density_Extractions_Radial.csv'));
fprintf('radial density table saved.\n');
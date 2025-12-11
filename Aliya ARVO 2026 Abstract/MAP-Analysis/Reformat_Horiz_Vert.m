%Reformat outputs from Batch Horizontal Vertical Density Script
%load CSVs
folderPath = uigetdir(pwd, 'Select folder with CSV files');
csvFiles = dir(fullfile(folderPath, '*.csv'));
masterHoriz = table();
masterVert = table();

for k = 1:length(csvFiles)
    fileName = csvFiles(k).name;
    
    % Extract subject ID
    underscores = strfind(fileName, '_');
    if length(underscores) >= 2
        subjectID = fileName(1:underscores(2)-1);
    else
        subjectID = fileName(1:end-4);%debug from stack search
    end

    T = readtable(fullfile(folderPath, fileName));
    
    % Horizontal table
    T_h = table();
    T_h.Subject = repmat({subjectID}, height(T), 1);
    T_h.Eccentricity_um = T.Eccentricity_Horiz_um;
    T_h.Density = T.Density_Horiz;
    masterHoriz = [masterHoriz; T_h];
    
    % Vertical table
    T_v = table();
    T_v.Subject = repmat({subjectID}, height(T), 1);
    T_v.Eccentricity_um = T.Eccentricity_Vert_um;
    T_v.Density = T.Density_Vert;
    masterVert = [masterVert; T_v];
    
    fprintf('Processed %d/%d: %s\n', k, length(csvFiles), fileName);
end

% Save tables
writetable(masterHoriz, fullfile(folderPath, 'Density_Extractions_Horizontal.csv'));
writetable(masterVert, fullfile(folderPath, 'Density_Extractions_Vertical.csv'));
fprintf('Horizontal and Vertical master tables saved.\n');
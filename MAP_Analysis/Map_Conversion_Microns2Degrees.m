% Script to convert density matrix in microns to degrees
% 1/25/2024
% Jenna Grieshop


root_path = uigetdir('.','Select directory containing analyses');
root_dir = dir(root_path);
root_dir = struct2cell(root_dir)';  

% looks for all the MATFILE matrices
MATFILE_dir = root_dir(...
    ~cellfun(@isempty, strfind(root_dir(:,1), 'MATFILE')),:);

[scalingfname, scalingpath] = uigetfile(fullfile(root_path,'*.csv'),'Select scaling LUT, OR cancel if you want to input the scale directly.');

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



for i=1:size(MATFILE_dir,1)
    

    if isnan(scaleinput)
    % load in the data
    data = load(fullfile(MATFILE_dir{i,2}, MATFILE_dir{i,1}));

    converted_interped_map = zeros([size(data.interped_map,1), size(data.interped_map,2)]);
    
    % Calculate the scale for this identifier.                                
    LUTindex=find( cellfun(@(s) ~isempty(strfind(MATFILE_dir{i},s )), lutData{1} ) );
        for x=1:size(LUTindex, 1)
            if x == size(LUTindex, 1) % if it is the last/only item in the LUT - if only matches with the eye and not subID will have axial length as NAN (would happen if LUT doesn't have info needed for this dataset)
                LUTindex = LUTindex(x);
                break
            end
            val = LUTindex(x+1) - LUTindex(x); % checking if there are two eyes from the same subject in LUT
            if val == 1
                LUTindex = LUTindex(x);
                break
            end
        end
    end
    
    axiallength = lutData{2}(LUTindex);
    pixelsperdegree = lutData{3}(LUTindex);
    micronsperdegree = (291*axiallength)/24;
    subjectID = lutData{1};
    
    mic_scaleval = 1 / (pixelsperdegree / micronsperdegree);
    deg_scaleval = 1/pixelsperdegree;
    
    % convert from microns to degrees
    for j=1:size(data.interped_map,1)
        for k=1:size(data.interped_map,2)
            area = (size(data.interped_map,1)*1000^2)/data.interped_map(j,k);
            px_area = area/mic_scaleval^2;
            deg_area = px_area*(deg_scaleval^2);
            converted_density = size(data.interped_map,1)/deg_area;
            converted_interped_map(j,k) = converted_density;
        end
    end

    filename = fullfile(root_path,'Results',[subjectID{LUTindex} '_bounddensity_matrix_converted2DEG_' date '.csv']);
    writematrix(converted_interped_map, filename);




end



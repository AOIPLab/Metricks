function [ LUT_row LUT_col ] = load_LUT_file_HV_Density( fileloc )
% Robert Cooper 06-18-2012
%   This function loads the needed scaling information for our metrics.

fid = fopen(fileloc,'r');

LUT_col = textscan(fid , '%s %f %f %s', 'delimiter', ',');

i=1;
% Reform into cells that contain each row
for i=1:size(LUT_col{1},1)
   
    
    ID = LUT_col{1}(i);
    CDC_x{1} = LUT_col{2}(i);
    CDC_y{1} = LUT_col{3}(i);  
    Eye{1} = LUT_col{4}(i);
    
    LUT_row{i} = [ID CDC_x CDC_y Eye];
    
end

fclose(fid);
end


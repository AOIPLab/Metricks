function [ scaling_row scaling_col ] = load_scaling_file_ddd( fileloc )
% Robert Cooper 06-18-2012
%   This function loads the needed scaling information for our metrics.
% modified for Density at a Distance Degrees

fid = fopen(fileloc,'r');

scaling_col = textscan(fid , '%s %f %f %f %f', 'delimiter', ',');

i=1;
% Reform into cells that contain each row
for i=1:size(scaling_col{1},1)
   
    
    ID = scaling_col{1}(i);
    axial{1} = scaling_col{2}(i);
    pix_per_deg{1} = scaling_col{3}(i);  
    PCD_x{1} = scaling_col{4}(i);
    PCD_y{1} = scaling_col{5}(i);
    
    scaling_row{i} = [ID axial pix_per_deg];
    
end

fclose(fid);
end


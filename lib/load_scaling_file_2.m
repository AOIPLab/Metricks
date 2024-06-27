function [ scaling_row scaling_col ] = load_scaling_file2( fileloc )
% Robert Cooper 06-18-2012
%   This function loads the needed scaling information for our metrics.

fid = fopen(fileloc,'r');

scaling_col = textscan(fid , '%s %f %f %f %f %f', 'delimiter', ',');

i=1;
% Reform into cells that contain each row
for i=1:size(scaling_col{1},1)
   
    
    ID = scaling_col{1}(i);
    og_center{1} = scaling_col{2}(i);
    scaling_factor{1} = scaling_col{3}(i);
    new_center{1} = scaling_col{4}(i);
    mpp{1} = scaling_col{5}(i);
    ppd{1} = scaling_col{6}(i);
    
    scaling_row{i} = [ID og_center scaling_factor, new_center, mpp, ppd];
    
end

fclose(fid);
end


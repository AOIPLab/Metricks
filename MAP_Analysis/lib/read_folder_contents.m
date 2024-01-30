function [ file_list, isdir, numfiles ] = read_folder_contents( root_dir, extension, wildcard )
% function [ file_list, isdir, numfiles ] = READ_FOLDER_CONTENTS( root_dir, extension, wildcard )
% Robert Cooper 12 - 05 - 09
% This function extracts all filenames from a folder into a cell array.
%
% [ file_list, isdir, numfiles ] = READ_FOLDER_CONTENTS( root_dir )
%       Returns an Nx1 cell array of character vectors (file_list);
%       each cell contains the filename of a file in the directory
%       specified (root_dir).
%
%       'isdir' is a corresponding Nx1 vector that contains boolean values
%       if a given path is a directory (true) or a file (false).
%
%       'numfiles' is the total number of files that are detected in the
%       directory.
%       
%
% [ ... ] = READ_FOLDER_CONTENTS( root_dir, extension )
%       Including an character vector with JUST the extension of the files of
%       interest will only return the above for files that match that
%       extension.
%
%       For example: READ_FOLDER_CONTENTS( 'C:\Windows', 'dll'); will
%       return a cell array of all dll files in the windows directory.
%
% [ ... ] = READ_FOLDER_CONTENTS( root_dir, extension, wildcard )
%       Including an character vector with the extension and wildcard
%       of the files of interest will only return the above for files that
%       match that wildcard and extension.
%
%       For example: READ_FOLDER_CONTENTS( 'C:\Windows', 'dll', 'py'); will
%       return a cell array of all dll files in the Windows directory that
%       contain the string 'py'.
%



x=1;
if ~exist('wildcard','var') || strcmp(wildcard, '*')
    wildcard = '*';
else
    wildcard = ['*' wildcard '*'];
end

if ~exist('extension','var')
    file_list=dir(fullfile(root_dir, wildcard) );
    file_list=file_list(3:end);
else
    file_list=dir(fullfile(root_dir,[wildcard extension]));
end
numoffiles=length(file_list);


% This should NOT be used with large numoffiles sizes (large memory
% footprint)
filenames=cell(numoffiles,1);
is_dir=false(numoffiles,1);

for x=1:1:numoffiles

    temp=file_list(x,1).name;
    temp_isdir=file_list(x,1).isdir;
    filenames{x}=temp; % -2 offset to correct for . and .. in the beginning
    is_dir(x)=temp_isdir;
    x=x+1;

end

numfiles=numoffiles;
file_list=filenames(1:numoffiles); 
isdir=is_dir(1:numoffiles);

end


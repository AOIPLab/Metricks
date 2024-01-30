function [ file_paths ] = read_folder_contents_rec( basedir, extension, wildcard )
% function [ file_list, isdir, numfiles ] = READ_FOLDER_CONTENTS_REC( root_dir, extension, wildcard )
% Robert Cooper 17 - 06
% This function recursively extracts all filenames from a folder into a cell array.
%
% [ file_list, isdir, numfiles ] = READ_FOLDER_CONTENTS_REC( root_dir )
%       Returns an Nx1 cell array of character vectors (file_list);
%       each cell contains the filename of a file in the subtree of
%       directory specified specified (root_dir).
%       
%
% [ ... ] = READ_FOLDER_CONTENTS_REC( root_dir, extension )
%       Including an character vector with JUST the extension of the files of
%       interest will only return the above for files that match that
%       extension.
%
%       For example: READ_FOLDER_CONTENTS( 'C:\Windows', 'dll'); will
%       return a cell array of all dll files in the Windows subdirectory.
%
% [ ... ] = READ_FOLDER_CONTENTS_REC( root_dir, extension, wildcard )
%       Including an character vector with the extension and wildcard
%       of the files of interest will only return the above for files that
%       match that wildcard and extension.
%
%       For example: READ_FOLDER_CONTENTS( 'C:\Windows', 'dll', 'py'); will
%       return a cell array of all dll files in the Windows subdirectory that
%       contain the string 'py'.
%

if ~exist('wildcard','var') || strcmp(wildcard, '*')
    wildcard = '*';
end

if ~exist('extension','var')
    file_list=dir(fullfile(basedir, wildcard) );
    file_list=file_list(3:end);
else
    file_list=dir(fullfile(basedir,[wildcard extension]));

    % If we have defined an extension, then make sure to look for only
    % directories too, as they won't trigger here.
    [dir_list, onlydir, numdirs] = read_folder_contents(basedir);
    
    dir_list = dir_list(onlydir);
    numdirs = numdirs - sum(~onlydir);
    onlydir = onlydir(onlydir);
end

file_paths={};

[file_list, isdir, numfiles] = read_folder_contents(basedir, extension, wildcard);

if exist('extension','var')
    file_list = [file_list; dir_list];
    isdir = [isdir; onlydir];
    numfiles = numfiles+numdirs;
end

for i=1:numfiles   
    if isdir(i)
        file_paths = [file_paths; read_folder_contents_rec(fullfile(basedir,file_list{i}), extension, wildcard)];
    else
        file_paths = [file_paths; {fullfile(basedir,file_list{i})}];
    end    
end

end


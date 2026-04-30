function write_metrics_results(basepath, result_fname, fname, statistics, first)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    % Grab the names of the fields we're working with
    datafields = fieldnames(statistics);

    if first

        warning off;
        [ success ] = mkdir(basepath,'Results');
        warning on;
        if ~success
            error('Failed to make results folder! Exiting...');
        end
        
        fid= fopen(fullfile(basepath,'Results', result_fname),'w');
    
        % If it is the first time writing the file, then write the
        % header
        fprintf(fid,'Filename'); 
  
        numfields = size(datafields,1);                
    
        k=1;
    
        while k <= numfields
    
            val = statistics.(datafields{k});
    
            % If it is a multi-dimensional field, remove it
            % from our csv, and write it separately.
            if size(val,1) ~= 1 || size(val,2) ~= 1   
                disp([datafields{k} ' removed!']);
                datafields = datafields([1:k-1 k+1:end]);                        
                numfields = numfields-1;                        
            else
    %           disp([fields{k} ' added!']);
                fprintf(fid,',%s',datafields{k});
                k = k+1;
            end 
    
    
        end  
        fprintf(fid,'\n');
    
    else % If it isn't the first entry, then append.
        fid= fopen(fullfile(basepath,'Results',result_fname ),'a');
    end
    
    % Write the file we've worked on as the first column
    fprintf(fid,'%s', fname);
    
    for k=1:size(datafields,1)
        val = statistics.(datafields{k});
        if size(val,1) == 1 || size(val,2) == 1
            val = statistics.(datafields{k});
    
            fprintf(fid,',%1.2f',val);
        end
    end
    
    fprintf(fid,'\n');
    fclose(fid);



end
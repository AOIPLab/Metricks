%% Robert F Cooper
% 1/13/21
%
% This script's sole job is to run multiple montages in series; to save time, 
% instead of operating on individual files though, it operates on a series
% of folders.


thisfolder = pwd;

[foveafile, foveafold] = uigetfile(fullfile(pwd,'*.xlsx'), 'Select the folders containing the foveas you wish to assess.');

[~,filelist,~] = xlsread(fullfile(foveafold,foveafile));

[lutfname, lutfolder] = uigetfile(fullfile(pwd,'*.csv'),'Select scaling LUT, OR cancel if you want to input the scale directly.');

%%
restartf=1;
endf=length(filelist);

%%
for f=restartf:endf
    restartf=f;

    something=false;
    if ~isempty(filelist{f,1})
        
        splitlist=split(filelist{f,1},'_');
        filedir=fullfile(lutfolder, [splitlist{2} '_' splitlist{4}],'confocal');
        
        [ scalinginfo, ~, lut ]=determine_scaling(filedir, filelist(f,1), fullfile(lutfolder,lutfname) ,'degrees');

        foveaim = imread(fullfile(filedir,filelist{f,1}));
        
        [spac, spacmap, confmap, summap, imbox] = fit_fourier_spacing(foveaim,96,true);
        
        densim = zeros(size(foveaim));
        
        weightedmap = spacmap./confmap;
        
        density_map = sqrt(3)./ (2*(weightedmap*scalinginfo).^2);
        density_map(isinf(density_map))=0;
        
        densim(imbox(2):imbox(2)+imbox(4), imbox(1):imbox(1)+imbox(3)) = density_map;
        
        result_path = fullfile(filedir,'Results_fovea_map');
        mkdir(result_path)

       
        lower01 = quantile(density_map(~isnan(density_map)),0.01);
        
        scaled_density = density_map-lower01;
        upper99 = quantile(scaled_density(~isnan(scaled_density)),0.99);
        
        scaled_density = 255.*(scaled_density./ upper99 );
        scaled_density(scaled_density>255) = 255;
        
        imwrite(uint8(scaled_density),parula(256), fullfile(result_path, [splitlist{2} '_' splitlist{4} '_Fovea.tif']));
        save(fullfile(result_path, [splitlist{2} '_' splitlist{4} '_Fovea.mat']),...
            'weightedmap','density_map','densim', 'spacmap', 'confmap', 'summap', 'imbox','-v7.3')
        
    end
end
function [avg_pixel_spac, interped_spac_map, interped_err_map, sum_map, imbox ] = fit_fourier_spacing(test_image, roi_size)


if ~exist('test_image','var') || isempty(test_image)
    [filename, pathname] = uigetfile('*.tif', 'Pick an image to segment');

    test_image = imread( fullfile(pathname, filename) );
    if size(test_image,3) >1
        test_image = test_image(:,:,1);
    end
end
% tic;

im_size = size(test_image);
if ~exist('roi_size','var') 
    roi_size = im_size;
end

roi_step = floor(roi_size/4);
interped_spac_map=[];

imcomps = bwconncomp( imclose(test_image>0,ones(5)) );
imbox = regionprops(imcomps, 'BoundingBox');


boxsizes = zeros(size(imbox,1),1);
for i=1:size(imbox,1)
    boxsizes(i)= imbox(i).BoundingBox(3)*imbox(i).BoundingBox(4);
end   
[~, maxsizeind]=max(boxsizes);
imbox = floor(imbox(maxsizeind).BoundingBox);

imbox(imbox<=0) = 1;
width_diff = im_size(2)-(imbox(1)+imbox(3));
if width_diff  < 0 
    imbox(3) = imbox(3)+width_diff;
end
height_diff = im_size(1)-(imbox(2)+imbox(4));
if height_diff  < 0 
    imbox(4) = imbox(4)+height_diff;
end

if any( im_size <= roi_size)        
    % Our roi size should always be divisible by 2 (for simplicity).
    if rem(min(roi_size),2) ~= 0
        roi_size = min(roi_size)-1;
    end
    roi = {test_image(1:end-1,1:end-1)};
else
    % Our roi size should always be divisible by 2 (for simplicity).
    if rem(roi_size,2) ~= 0
        roi_size = roi_size-1;
    end
    roi = cell(round((size(test_image)-roi_size)/roi_step));

    for i=imbox(2):roi_step:imbox(2)+imbox(4)-roi_size
        for j=imbox(1):roi_step:imbox(1)+imbox(3)-roi_size

            numzeros = sum(sum(test_image(i:i+roi_size-1, j:j+roi_size-1)<=10));
            
            if numzeros < (roi_size*roi_size)*0.05
                roi{round(i/roi_step)+1,round(j/roi_step)+1} = test_image(i:i+roi_size-1, j:j+roi_size-1);
            else
                roi{round(i/roi_step)+1,round(j/roi_step)+1} =[];
            end
        end
    end
end

numind = size(roi,1)*size(roi,2);
pixel_spac = nan(size(roi));
err = nan(size(roi));

  %         power_spect_export = power_spect-min(power_spect(:));
%         power_spect_export = power_spect_export./max(power_spect_export(:));
%         power_spect_export = power_spect_export.*255;
% %         
%         imwrite(uint8(power_spect_export),['pwr_spect ' num2str(r) '.tif']);

for r=1:length(pixel_spac(:))
    if ~isempty(roi{r})        
        
%         if roi_size <= 256 % We don't want this run on massive images (RAM sink)
% 
%             padsize = 2^(nextpow2(roi_size)+1);
%             padsize = (padsize-roi_size)/2;
% 
%             power_spect = fftshift(fft2( padarray(roi{r}, [padsize padsize]) ));
%         else
            power_spect = fftshift(fft2( roi{r} ));
%         end
        power_spect = log10(abs(power_spect).^2);
       
%         figure(100); imagesc(power_spect); axis image;

        rhosampling = .5;
        thetasampling = 1;

        [polarroi, power_spect_radius] = imcart2pseudopolar(power_spect,rhosampling,thetasampling,[],'linear');
        polarroi = circshift(polarroi,-90/thetasampling,1);
%         figure(101); imagesc(polarroi); axis image;
        
        upper_n_lower = [thetasampling:45 136:225 316:360]/thetasampling;
        left_n_right = [46:135 226:315]/thetasampling;
        upper_n_lower_fourierProfile = mean(polarroi(upper_n_lower,:));
        left_n_right_fourierProfile = mean(polarroi(left_n_right,:));
        fullfourierProfile = mean(polarroi);
%         figure(101); plot(upper_n_lower_fourierProfile); axis image;

        if ~all(isinf(upper_n_lower_fourierProfile)) && ~all(isnan(upper_n_lower_fourierProfile))

            [pixel_spac(r), ~, err(r)] = fourierFit(upper_n_lower_fourierProfile,[], false);
            pixel_spac(r) = 1/ (pixel_spac(r) / ((power_spect_radius*2)/rhosampling));
%             if err(r) < 0.2
%                 pixel_spac(r)
%             end
        else
            pixel_spac(r) = NaN;
        end
    end
end

avg_pixel_spac = mean(pixel_spac(~isnan(pixel_spac)) );
std_pixel_spac = std(pixel_spac(~isnan(pixel_spac)));
interped_spac_map = avg_pixel_spac;
interped_err_map = err;


%% If we've sampled over the region, then create the heat map
if length(roi) > 1
    interped_spac_map=zeros(im_size);
    interped_err_map=zeros(im_size);
    interped_corrected_err_map=zeros(im_size);
    sum_map=zeros(im_size);
    
    for i=imbox(2):roi_step:imbox(2)+imbox(4)-roi_size
        for j=imbox(1):roi_step:imbox(1)+imbox(3)-roi_size

            if ~isnan( pixel_spac(round(i/roi_step)+1,round(j/roi_step)+1) )
                
                thiserr = err(round(i/roi_step)+1,round(j/roi_step)+1)^2;
%                 if thiserr > .44
                    interped_err_map(i:i+roi_size-1, j:j+roi_size-1) = interped_err_map(i:i+roi_size-1, j:j+roi_size-1) + thiserr;                
                    thisspac = pixel_spac(round(i/roi_step)+1,round(j/roi_step)+1);
                

                    interped_spac_map(i:i+roi_size-1, j:j+roi_size-1) = interped_spac_map(i:i+roi_size-1, j:j+roi_size-1) + thiserr*thisspac;

                    sum_map(i:i+roi_size-1, j:j+roi_size-1) = sum_map(i:i+roi_size-1, j:j+roi_size-1) + 1;

%                 end
            else

            end
        end
    end
    
    
    
%     [X,Y]=meshgrid( 1:roi_step:(size(test_image,2)-roi_size-1), 1:roi_step:(size(test_image,1)-roi_size-1));
%     [Xq,Yq]=meshgrid( 1:(size(test_image,2)-roi_size-1), 1:(size(test_image,1)-roi_size-1));
%     interped_spac_map = interp2( X,Y, pixel_spac, Xq, Yq);
    
    interped_spac_map = interped_spac_map( imbox(2):imbox(2)+imbox(4), imbox(1):imbox(1)+imbox(3) );
    interped_err_map = interped_err_map( imbox(2):imbox(2)+imbox(4), imbox(1):imbox(1)+imbox(3) );
    sum_map = sum_map( imbox(2):imbox(2)+imbox(4), imbox(1):imbox(1)+imbox(3) );
    
    figure(1);clf; imagesc((2/sqrt(3)).*interped_spac_map./interped_err_map); axis image;
    figure(2);clf; imagesc(interped_err_map./sum_map); axis image; colormap(flipud(jet(256)));
    caxis([0 1])
    figure(3);clf; imagesc(sum_map); axis image; colormap gray;
end


% pause;
        
end
function [combined_output] = straight_cross(M, CDC_x, CDC_y, thickness, pixel_size, eye)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    thick = floor(thickness/2);


    %% Extract meridian profiles
    % ======= Extract vertical profile =======
    dens_vert = M(:, CDC_x-thick:CDC_x+thick); %add thickness here
    dens_vert = mean(dens_vert,2);

    % Pixel offsets (top = negative, bottom = positive)
    pix_v = ( (1:length(dens_vert))' - CDC_y );
    ecc_vert = pix_v * pixel_size;   % microns

    % ======= Extract horizontal profile =======
    dens_horz = M(CDC_y-thick:CDC_y+thick, :); %add thickness here
    dens_horz = mean(dens_horz,1);

    % Pixel offsets (left = negative, right = positive)
    pix_h = ( (1:length(dens_horz)) - CDC_x );
    ecc_horz = pix_h(:) * pixel_size;   % microns

    % ======= Orientation correction =======
    if strcmpi(eye, 'OS')  
        dens_horz = fliplr(dens_horz); % flip horizontally for OS to match OD
        ecc_horz = -flipud(ecc_horz);
    end

    %% split data at the cdc point

    % find the 0 point in the eccentricity
    v_0_index = find(ecc_vert == 0); 
    h_0_index = find(ecc_horz == 0); 

    % split the vertical data into top and bottom
    v_ecc_t = -flipud(ecc_vert(1:v_0_index));
    v_dens{1} = flipud(dens_vert(1:v_0_index))'; % top
    v_ecc_b = ecc_vert(v_0_index:end);
    v_dens{2} = dens_vert(v_0_index:end)'; % bottom
    % find minimum eccentricity for the average
    if size(v_ecc_b, 1) > size(v_ecc_t,1)
        v_ecc_avg = v_ecc_t;
    else
        v_ecc_avg = v_ecc_b;
    end

    % split the horizontal data into left and right
    h_ecc_l = -flipud(ecc_horz(1:h_0_index));
    h_dens{1} = fliplr(dens_horz(1:h_0_index)); % left
    h_ecc_r = ecc_horz(h_0_index:end);
    h_dens{2} = dens_horz(h_0_index:end); % right
    % find minimum eccentricity for the average
    if size(h_ecc_r, 1) > size(h_ecc_l,1)
        h_ecc_avg = h_ecc_l;
    else
        h_ecc_avg = h_ecc_r;
    end


    % prepare for vertical averaging
    maxNumCol_v = max(cellfun(@(c) size(c,2), v_dens));  % max number of columns
    combined_v = cell2mat(cellfun(@(c){padarray(c,[0,maxNumCol_v-size(c,2)],NaN,'Post')}, v_dens)');

    % prepare for horizontal averaging
    maxNumCol_h = max(cellfun(@(c) size(c,2), h_dens));  % max number of columns
    combined_h = cell2mat(cellfun(@(c){padarray(c,[0,maxNumCol_h-size(c,2)],NaN,'Post')}, h_dens)');

    % agerage the vertical and horizontal densities
    avg_v_dens = rmmissing(mean(combined_v, 1));
    avg_h_dens = rmmissing(mean(combined_h, 1));


    % compile results and create headers for output sheet
    results = {h_ecc_l, h_dens{1}', h_ecc_r, h_dens{2}', v_ecc_t, v_dens{1}', v_ecc_b, v_dens{2}', h_ecc_avg, avg_h_dens', v_ecc_avg, avg_v_dens'};
    headers = {'Ecc_Horz_L_um', 'Density_Horz_L', 'Ecc_Horz_R_um', 'Density_Horz_R', 'Ecc_Vert_T_um,', 'Density_Vert_T', 'Ecc_Vert_B_um', 'Density_Vert_B', 'Ecc_Avg_Horz_um', 'Density_Avg_Horz', 'Ecc_Avg_Vert_um', 'Density_Avg_Vert'};

    %% Plotting for sanity check
    % plot(v_ecc_t, v_dens{1});
    % hold on
    % plot(v_ecc_b, v_dens{2});
    % plot(h_ecc_l, h_dens{1});
    % plot(h_ecc_r, h_dens{2});
    % 
    % plot(v_ecc_avg, avg_v_dens);
    % plot(h_ecc_avg, avg_h_dens);

    

    % organize all the data into a format that can be saved to the csv
    D = {};
    for k = 1:numel(results)
      D(1:numel(results{k}),k)=num2cell(results{k});
    end
    combined_output = [headers; D];
    
end
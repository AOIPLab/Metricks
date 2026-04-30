function [ mosaic_stats_um, mosaic_stats_deg, mosaic_stats_arcmin  ] = determine_mosaic_stats( coords, pixelsperdegree, micronsperdegree, bounds, ignore_idx, reliability )
% Robert Cooper 09-24-14
% This function takes in a list of coordinates in a m-2 matrix, and
% calculates the mean nearest neighbor, cell area created by the
% coordinates, and calculates the density of the coordinates

%% Computing different scaling variables - added by MG 4/3/2026

% microns or cones/mm^2 for density
scaleval_um = 1 / (pixelsperdegree / micronsperdegree);

% degrees
scaleval_deg = 1/pixelsperdegree;

% arcmin
scaleval_arcmin = 60/pixelsperdegree;


%% Coords are in X,Y!

clipped_row_col = [bounds(2)-bounds(1) bounds(4)-bounds(3)];

clipped_coords = coordclip(coords,bounds(1:2),bounds(3:4),'i');

%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine Mean N-N %%
%%%%%%%%%%%%%%%%%%%%%%%%

dist_between_pts=pdist2(clipped_coords,clipped_coords); % Measure the distance from each set of points to the other
max_ident=eye(length(dist_between_pts)).*max(dist_between_pts(:)); % Make diagonal not the minimum for any observation

[minval minind]=min(dist_between_pts+max_ident); % Find the minimum distance from one set of obs to another

mean_nn_dist_um = mean(minval.*scaleval_um); % Distance in units
mean_nn_dist_deg = mean(minval.*scaleval_deg); % Distance in units
mean_nn_dist_arcmin = mean(minval.*scaleval_arcmin); % Distance in units

% std(minval.*um_per_pix)
regularity_nn_index_um = mean_nn_dist_um/std(minval.*scaleval_um); % Units cancel out but we'll still save everything separately 
regularity_nn_index_deg = mean_nn_dist_deg/std(minval.*scaleval_deg);
regularity_nn_index_arcmin = mean_nn_dist_arcmin/std(minval.*scaleval_arcmin);

% wb = waitbar(.2,'Determining Voronoi');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine Voronoi Cell Area %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sixsided=0;
bound = false(size(coords,1),1);
cellarea = zeros(size(coords,1),1);
numedges = zeros(size(coords,1),1);

if size(coords,1) > 2

    [V,C] = voronoin(coords,{'QJ'}); % Returns the vertices of the Voronoi edges in VX and VY so that plot(VX,VY,'-',X,Y,'.')
    fastbound = (V(:,1)<bounds(2) & V(:,1)>bounds(1) & V(:,2)<bounds(4) & V(:,2)>bounds(3));
    
    % figure(10); hold on;
    for i=1:length(C)

        vertices=V(C{i},:);
      
        if all(fastbound(C{i})) && all(i ~= ignore_idx) && all(C{i}~=1)  

            cellarea(i) = polyarea(vertices(:,1),vertices(:,2));

            % Code to display number of sides for each voronoi domain
            numedges(i)=size(V(C{i},1),1);
            switch(numedges(i))
    %             case 4
    %                 color = 'm';
    %             case 5
    %                 color = 'c';
                case 6
    %                 color = 'g';
                  sixsided = sixsided+1;
    %             case 7
    %                 color = 'y';
    %             case 8
    %                 color = 'r';
    %             case 9
    %                 color = 'b';
            end
            % figure(10);
            % patch(V(C{i},1),V(C{i},2),ones(size(V(C{i},1))),'FaceColor',color);
            % hold on;

            bound(i) = true;
        end

    end

end
% hold off;
% toc
% figure(2);
% voronoi(coords(:,1),coords(:,2));
if sum(bound) ~= 0
    coords_bound= coords(bound,:); % Clip out the unbounded cells

    cellarea_deg = cellarea((cellarea~=0)).*(scaleval_deg.^2);
    cellarea_arcmin = cellarea((cellarea~=0)).*(scaleval_arcmin.^2);
    cellarea_um= cellarea((cellarea~=0)).*(scaleval_um.^2); % Clip out unbounded cells, convert to square microns
    
    numedges = numedges(numedges~=0);
    
    mean_cellarea_deg=mean(cellarea_deg);
    mean_cellarea_arcmin=mean(cellarea_arcmin);
    mean_cellarea_um=mean(cellarea_um);


    regularity_voro_index_deg = mean_cellarea_deg/std(cellarea_deg);
    regularity_voro_index_arcmin = mean_cellarea_arcmin/std(cellarea_arcmin);
    regularity_voro_index_um = mean_cellarea_um/std(cellarea_um);



    regularity_voro_sides = mean(numedges)/std(numedges);
    
%     disp([ 'Mean: ' num2str(mean(numedges))  ' Std deviation: ' num2str(std(numedges)) ] );
    percent_six_sided = 100*sixsided/size(coords_bound,1);
else
    cellarea_um = 0;
    cellarea_arcmin = 0;
    cellarea_deg = 0;

    mean_cellarea_um = 0;
    mean_cellarea_arcmin = 0;
    mean_cellarea_deg = 0;

    regularity_voro_index_deg=0;
    regularity_voro_index_arcmin=0;
    regularity_voro_index_um=0;


    regularity_voro_sides=0;
    percent_six_sided=0;
    coords_bound = [];
end
% waitbar(0.4,wb,'Determining Density');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine Number of Cells, Density Direct Count (D_dc) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numcells=length(clipped_coords); % Total number of cells
total_cell_area_um = sum(cellarea_um); % Total cell area in um
total_cell_area_arcmin = sum(cellarea_arcmin);
total_cell_area_deg = sum(cellarea_deg); 

% microns (mm density)
    total_coord_area_um=(((clipped_row_col(1)*clipped_row_col(2))*((scaleval_um^2)/(1000^2))))*1000^2; 
% degrees
    total_coord_area_deg=((clipped_row_col(1)*clipped_row_col(2))*((scaleval_deg^2)));    
% arcmin 
    total_coord_area_arcmin=((clipped_row_col(1)*clipped_row_col(2))*((scaleval_arcmin^2))); 


pixel_density = numcells/(clipped_row_col(1)*clipped_row_col(2));

density_dc_um=(1000^2)*numcells/total_coord_area_um; % cells/um^2
density_dc_deg=numcells/total_coord_area_deg; 
density_dc_arcmin=numcells/total_coord_area_arcmin; 

if ~isempty(coords_bound)
    
        density_bound_um = (1000^2)*size(coords_bound,1)./total_cell_area_um;
        density_bound_deg = size(coords_bound,1)./total_cell_area_deg;
        density_bound_arcmin = size(coords_bound,1)./total_cell_area_arcmin;
else
    density_bound_um = 0;
    density_bound_deg = 0;
    density_bound_arcmin = 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine Inter-Cell Distance %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


inter_cell_dist_um = zeros(size(clipped_coords,1),1);
inter_cell_dist_deg = zeros(size(clipped_coords,1),1);
inter_cell_dist_arcmin = zeros(size(clipped_coords,1),1);

max_cell_dist_um = zeros(size(clipped_coords,1),1);
max_cell_dist_deg = zeros(size(clipped_coords,1),1);
max_cell_dist_arcmin = zeros(size(clipped_coords,1),1);

correct_inter_cell_dist_um = zeros(sum(bound),1);
correct_inter_cell_dist_deg = zeros(sum(bound),1);
correct_inter_cell_dist_arcmin = zeros(sum(bound),1);

correct_max_cell_dist_um = zeros(sum(bound),1);
correct_max_cell_dist_deg = zeros(sum(bound),1);
correct_max_cell_dist_arcmin = zeros(sum(bound),1);

correct_nn_cell_dist_um = zeros(sum(bound),1);
correct_nn_cell_dist_deg = zeros(sum(bound),1);
correct_nn_cell_dist_arcmin = zeros(sum(bound),1);


if size(coords,1) > 2

    dt = DelaunayTri(coords);    

    % Find all instances of each bound cell.
    boundinds = find(bound);
    for k=1:numel(boundinds)

        % If its bound, then we've flagged it as such, and can use it in the triangulation
        % Only take the first row because that is the cell of interest's
        % relative distance to its neighboring cells
        ind = boundinds(k);

        [i, j] =find(dt.Triangulation == ind);

        conn_ind = dt.Triangulation(i,:);

        coord_row = unique(conn_ind( conn_ind ~= ind)); % Find all of the unique coordinate points that isn't the "center" coordinate

        if(size(i,1)~=1)
            coord_row = [ind; coord_row]; % Add the "center" to the top, so we know the order for the distances
        else
            coord_row = [ind; coord_row']; 
        end

        cell_dist = squareform(pdist([coords(coord_row,1) coords(coord_row,2)]));
            
        correct_inter_cell_dist_um(k) = scaleval_um*(sum(cell_dist(1,:)) / (length(cell_dist(1,:))-1));
        correct_inter_cell_dist_deg(k) = scaleval_deg*(sum(cell_dist(1,:)) / (length(cell_dist(1,:))-1));
        correct_inter_cell_dist_arcmin(k) = scaleval_arcmin*(sum(cell_dist(1,:)) / (length(cell_dist(1,:))-1));

        correct_max_cell_dist_um(k)   = scaleval_um*max(cell_dist(1,:));
        correct_max_cell_dist_deg(k)   = scaleval_deg*max(cell_dist(1,:));
        correct_max_cell_dist_arcmin(k)   = scaleval_arcmin*max(cell_dist(1,:));

        correct_nn_cell_dist_um(k)    = scaleval_um*min(cell_dist(1,2:end));
        correct_nn_cell_dist_deg(k)    = scaleval_deg*min(cell_dist(1,2:end));
        correct_nn_cell_dist_arcmin(k)    = scaleval_arcmin*min(cell_dist(1,2:end));

    end

    % Repeat the above, but with all cells in the unbound region..
    dt = DelaunayTri(clipped_coords);
    for k=1:size(clipped_coords,1)

        [i, j] =find(dt.Triangulation == k);

        conn_ind = dt.Triangulation(i,:);

        coord_row = unique(conn_ind( conn_ind ~= k)); % Find all of the unique coordinate points that isn't the "center" coordinate

        if(size(i,1)~=1)
            coord_row = [k; coord_row]; % Add the "center" to the top, so we know the order for the distances
        else
            coord_row = [k; coord_row']; 
        end

        cell_dist = squareform(pdist([coords(coord_row,1) coords(coord_row,2)]));

        inter_cell_dist_um(k) = scaleval_um*(sum(cell_dist(1,:)) / (length(cell_dist(1,:))-1));
        inter_cell_dist_deg(k) = scaleval_deg*(sum(cell_dist(1,:)) / (length(cell_dist(1,:))-1));
        inter_cell_dist_arcmin(k) = scaleval_arcmin*(sum(cell_dist(1,:)) / (length(cell_dist(1,:))-1));

        max_cell_dist_um(k)   = scaleval_um*max(cell_dist(1,:));
        max_cell_dist_deg(k)   = scaleval_deg*max(cell_dist(1,:));
        max_cell_dist_arcmin(k)   = scaleval_arcmin*max(cell_dist(1,:));
    end

    mean_inter_cell_dist_um = mean(inter_cell_dist_um);
    mean_inter_cell_dist_deg = mean(inter_cell_dist_deg);
    mean_inter_cell_dist_arcmin = mean(inter_cell_dist_arcmin);

    mean_max_cell_dist_um   = mean( max_cell_dist_um);
    mean_max_cell_dist_deg   = mean( max_cell_dist_deg);
    mean_max_cell_dist_arcmin   = mean( max_cell_dist_arcmin);

else
    mean_inter_cell_dist_um = scaleval_um*pdist(coords);
    mean_inter_cell_dist_deg = scaleval_deg*pdist(coords);
    mean_inter_cell_dist_arcmin = scaleval_arcmin*pdist(coords);

    mean_max_cell_dist_um = mean_inter_cell_dist_um;
    mean_max_cell_dist_deg = mean_inter_cell_dist_deg;
    mean_max_cell_dist_arcmin = mean_inter_cell_dist_arcmin;

end
    
if ~isempty(coords_bound)
    mean_correct_nn_dist_um = mean( correct_nn_cell_dist_um );
    mean_correct_nn_dist_deg = mean( correct_nn_cell_dist_deg );
    mean_correct_nn_dist_arcmin = mean( correct_nn_cell_dist_arcmin );

    mean_correct_inter_cell_dist_um = mean(correct_inter_cell_dist_um);
    mean_correct_inter_cell_dist_deg = mean(correct_inter_cell_dist_deg);
    mean_correct_inter_cell_dist_arcmin = mean(correct_inter_cell_dist_arcmin);

    regularity_ic_index_um = mean(correct_inter_cell_dist_um)./std(correct_inter_cell_dist_um);
    regularity_ic_index_deg = mean(correct_inter_cell_dist_deg)./std(correct_inter_cell_dist_deg);
    regularity_ic_index_arcmin = mean(correct_inter_cell_dist_arcmin)./std(correct_inter_cell_dist_arcmin);


    mean_correct_max_cell_dist_um   = mean( correct_max_cell_dist_um );
    mean_correct_max_cell_dist_deg   = mean( correct_max_cell_dist_deg );
    mean_correct_max_cell_dist_arcmin   = mean( correct_max_cell_dist_arcmin );
    
else
    regularity_ic_index_um = 0;
    regularity_ic_index_deg = 0;
    regularity_ic_index_arcmin = 0;

    mean_correct_nn_dist_um = 0;
    mean_correct_nn_dist_deg = 0;
    mean_correct_nn_dist_arcmin = 0;

    mean_correct_inter_cell_dist_um = 0;
    mean_correct_inter_cell_dist_deg = 0;
    mean_correct_inter_cell_dist_arcmin = 0;

    mean_correct_max_cell_dist_um = 0;
    mean_correct_max_cell_dist_deg = 0;
    mean_correct_max_cell_dist_arcmin = 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine Density Recovery Profile %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% [ density_per_rad, um_drp_sizes, drp_spac]=calculate_DRP(coords, [bounds(1:2); bounds(3:4)], scale, pixel_density, reliability );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output List Formatting %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make the returned structs


mosaic_stats_um = struct('Number_Unbound_Cells', numcells,'Number_Bound_Cells', sum(bound), 'Total_Area', total_coord_area_um, 'Total_Bound_Area',total_cell_area_um,...                      
                      'Bound_Density',density_bound_um, 'Bound_NN_Distance',mean_correct_nn_dist_um,'Bound_IC_Distance',mean_correct_inter_cell_dist_um,'Bound_Furthest_Distance',mean_correct_max_cell_dist_um,...
                      'Bound_Mean_Voronoi_Area', mean_cellarea_um,'Bound_Percent_Six_Sided_Voronoi',percent_six_sided,'Unbound_DRP_Distance', 0,...
                      'Bound_Voronoi_Area_RI',regularity_voro_index_um,'Bound_Voronoi_Sides_RI',regularity_voro_sides, 'Bound_NN_RI', regularity_nn_index_um, 'Bound_IC_RI', regularity_ic_index_um,...
                      'Unbound_Density', density_dc_um ,'Unbound_NN_Distance', mean_nn_dist_um, 'Unbound_IC_Distance',mean_inter_cell_dist_um, 'Unbound_Furthest_Distance',mean_max_cell_dist_um);

mosaic_stats_deg = struct('Number_Unbound_Cells', numcells,'Number_Bound_Cells', sum(bound), 'Total_Area', total_coord_area_deg, 'Total_Bound_Area',total_cell_area_deg,...                      
                      'Bound_Density',density_bound_deg, 'Bound_NN_Distance',mean_correct_nn_dist_deg,'Bound_IC_Distance',mean_correct_inter_cell_dist_deg,'Bound_Furthest_Distance',mean_correct_max_cell_dist_deg,...
                      'Bound_Mean_Voronoi_Area', mean_cellarea_deg,'Bound_Percent_Six_Sided_Voronoi',percent_six_sided,'Unbound_DRP_Distance', 0,...
                      'Bound_Voronoi_Area_RI',regularity_voro_index_deg,'Bound_Voronoi_Sides_RI',regularity_voro_sides, 'Bound_NN_RI', regularity_nn_index_deg, 'Bound_IC_RI', regularity_ic_index_deg,...
                      'Unbound_Density', density_dc_deg ,'Unbound_NN_Distance', mean_nn_dist_deg, 'Unbound_IC_Distance',mean_inter_cell_dist_deg, 'Unbound_Furthest_Distance',mean_max_cell_dist_deg);

mosaic_stats_arcmin = struct('Number_Unbound_Cells', numcells,'Number_Bound_Cells', sum(bound), 'Total_Area', total_coord_area_arcmin, 'Total_Bound_Area',total_cell_area_arcmin,...                      
                      'Bound_Density',density_bound_arcmin, 'Bound_NN_Distance',mean_correct_nn_dist_arcmin,'Bound_IC_Distance',mean_correct_inter_cell_dist_arcmin,'Bound_Furthest_Distance',mean_correct_max_cell_dist_arcmin,...
                      'Bound_Mean_Voronoi_Area', mean_cellarea_arcmin,'Bound_Percent_Six_Sided_Voronoi',percent_six_sided,'Unbound_DRP_Distance', 0,...
                      'Bound_Voronoi_Area_RI',regularity_voro_index_arcmin,'Bound_Voronoi_Sides_RI',regularity_voro_sides, 'Bound_NN_RI', regularity_nn_index_arcmin, 'Bound_IC_RI', regularity_ic_index_arcmin,...
                      'Unbound_Density', density_dc_arcmin ,'Unbound_NN_Distance', mean_nn_dist_arcmin, 'Unbound_IC_Distance',mean_inter_cell_dist_arcmin, 'Unbound_Furthest_Distance',mean_max_cell_dist_arcmin);




end


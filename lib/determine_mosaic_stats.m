function [ mosaic_stats ] = determine_mosaic_stats( coords, scale, scaledeg, unit, bounds , ignore_idx, reliability )
% Robert Cooper 09-24-14
% This function takes in a list of coordinates in a m-2 matrix, and
% calculates the mean nearest neighbor, cell area created by the
% coordinates, and calculates the density of the coordinates

%% Coords are in X,Y!


clipped_row_col = [bounds(2)-bounds(1) bounds(4)-bounds(3)];

clipped_coords = coordclip(coords,bounds(1:2),bounds(3:4),'i');

%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine Mean N-N %%
%%%%%%%%%%%%%%%%%%%%%%%%

dist_between_pts=pdist2(clipped_coords,clipped_coords); % Measure the distance from each set of points to the other
max_ident=eye(length(dist_between_pts)).*max(dist_between_pts(:)); % Make diagonal not the minimum for any observation

[minval minind]=min(dist_between_pts+max_ident); % Find the minimum distance from one set of obs to another

mean_nn_dist=mean(minval.*scale); % Distance in units

% std(minval.*um_per_pix)
regularity_nn_index = mean_nn_dist/std(minval.*scale);

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
    cellarea_deg = cellarea((cellarea~=0)).*(scaledeg.^2);
    cellarea= cellarea((cellarea~=0)).*(scale.^2); % Clip out unbounded cells, convert to square microns
    numedges = numedges(numedges~=0);
    
    mean_cellarea=mean(cellarea);
    regularity_voro_index = mean_cellarea/std(cellarea);
    regularity_voro_sides = mean(numedges)/std(numedges);
    
%     disp([ 'Mean: ' num2str(mean(numedges))  ' Std deviation: ' num2str(std(numedges)) ] );
    percent_six_sided = 100*sixsided/size(coords_bound,1);
else
    cellarea=0;
    cellarea_deg = 0;
    mean_cellarea=0;
    regularity_voro_index=0;
    regularity_voro_sides=0;
    percent_six_sided=0;
    coords_bound = [];
end
% waitbar(0.4,wb,'Determining Density');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine Number of Cells, Density Direct Count (D_dc) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numcells=length(clipped_coords); % Total number of cells
total_cell_area=sum(cellarea); % Total cell area in units
total_cell_area_deg = sum(cellarea_deg);

if strcmp(unit,'microns (mm density)')
    total_coord_area=((clipped_row_col(1)*clipped_row_col(2))*((scale^2)/(1000^2))); 
else
    total_coord_area=((clipped_row_col(1)*clipped_row_col(2))*((scale^2)));    
end

pixel_density = numcells/(clipped_row_col(1)*clipped_row_col(2));
density_dc=numcells/total_coord_area; % cells/mm^2

if ~isempty(coords_bound)
    if strcmp(unit,'microns (mm density)')
        density_bound = (1000^2)*size(coords_bound,1)./total_cell_area;
        density_bound_deg = size(coords_bound,1)./total_cell_area_deg;
    else
        density_bound = size(coords_bound,1)./total_cell_area;
        density_bound_deg = size(coords_bound,1)./total_cell_area;
    end
else
    density_bound = 0;
    density_bound_deg = 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine Inter-Cell Distance %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


inter_cell_dist = zeros(size(clipped_coords,1),1);
max_cell_dist = zeros(size(clipped_coords,1),1);

correct_inter_cell_dist = zeros(sum(bound),1);
correct_max_cell_dist = zeros(sum(bound),1);
correct_nn_cell_dist = zeros(sum(bound),1);
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
            
        correct_inter_cell_dist(k) = scale*(sum(cell_dist(1,:)) / (length(cell_dist(1,:))-1));
        correct_max_cell_dist(k)   = scale*max(cell_dist(1,:));
        correct_nn_cell_dist(k)    = scale*min(cell_dist(1,2:end));        
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

        inter_cell_dist(k) = scale*(sum(cell_dist(1,:)) / (length(cell_dist(1,:))-1));
        max_cell_dist(k)   = scale*max(cell_dist(1,:));
    end


    mean_inter_cell_dist = mean(inter_cell_dist);    
    mean_max_cell_dist   = mean( max_cell_dist );
else
    mean_inter_cell_dist = scale*pdist(coords);
    mean_max_cell_dist = mean_inter_cell_dist;
end
    
if ~isempty(coords_bound)
    mean_correct_nn_dist = mean( correct_nn_cell_dist );
    mean_correct_inter_cell_dist = mean(correct_inter_cell_dist);
    regularity_ic_index = mean(correct_inter_cell_dist)./std(correct_inter_cell_dist);
    mean_correct_max_cell_dist   = mean( correct_max_cell_dist );    
    
else
    regularity_ic_index = 0;
    mean_correct_nn_dist=0;
    mean_correct_inter_cell_dist=0;
    mean_correct_max_cell_dist=0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine Density Recovery Profile %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% [ density_per_rad, um_drp_sizes, drp_spac]=calculate_DRP(coords, [bounds(1:2); bounds(3:4)], scale, pixel_density, reliability );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output List Formatting %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make the returned struct
if strcmp(unit,'microns (mm density)')
    total_coord_area=total_coord_area*1000^2;      
end

mosaic_stats = struct('Number_Unbound_Cells', numcells,'Number_Bound_Cells', sum(bound), 'Total_Area', total_coord_area, 'Total_Bound_Area',total_cell_area,...                      
                      'Bound_Density',density_bound, 'Bound_NN_Distance',mean_correct_nn_dist,'Bound_IC_Distance',mean_correct_inter_cell_dist,'Bound_Furthest_Distance',mean_correct_max_cell_dist,...
                      'Bound_Mean_Voronoi_Area', mean_cellarea,'Bound_Percent_Six_Sided_Voronoi',percent_six_sided,'Unbound_DRP_Distance', 0,...
                      'Bound_Voronoi_Area_RI',regularity_voro_index,'Bound_Voronoi_Sides_RI',regularity_voro_sides, 'Bound_NN_RI', regularity_nn_index, 'Bound_IC_RI', regularity_ic_index,...
                      'Unbound_Density', density_dc ,'Unbound_NN_Distance', mean_nn_dist, 'Unbound_IC_Distance',mean_inter_cell_dist, 'Unbound_Furthest_Distance',mean_max_cell_dist, 'Bound_Density_DEG',density_bound_deg);

end


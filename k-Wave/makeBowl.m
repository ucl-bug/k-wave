function bowl = makeBowl(grid_size, bowl_pos, radius, diameter, focus_pos, varargin)
%MAKEBOWL Create a binary map of a bowl within a 3D grid.
%
% DESCRIPTION:
%     makeBowl creates a binary map of a bowl from a symmetric spherical
%     shell within a three-dimensional grid. The bowl position is denoted
%     by 1's in the matrix with 0's elsewhere. The bowl surface within the
%     Cartesian grid is found using a bi-directional line search along each
%     row and column to find the grid point with radius closest to the
%     desired radius.
%
%     The position of the bowl is set by bowl_pos, which corresponds to the
%     center of the rear bowl surface. The orientation of the bowl is set
%     by focus_pos, which corresponds to any point on the axis of the bowl
%     (note, this must not be equal to bowl_pos). It is assumed that the
%     solid angle of the bowl is equal to or less than 2*pi steradians. If
%     the radius is set to inf, a disc is generated.
%
%     In some cases, the generated bowl will not be simply connected, and
%     there will be a small number of overlapping grid points. To remove
%     these points, set the optional input 'RemoveOverlap' to true.
%
% USAGE:
%     bowl = makeBowl(grid_size, bowl_pos, radius, diameter, focus_pos)
%     bowl = makeBowl(grid_size, bowl_pos, radius, diameter, focus_pos, ...)
%
% INPUTS:
%     grid_size       - size of the 3D grid given as a three element
%                       vector [Nx, Ny, Nz] [grid points] 
%     bowl_pos        - centre of the rear surface of the bowl given as a
%                       three element vector [bx, by, bz] [grid points] 
%     radius          - radius of curvature of the bowl [grid points]
%     diameter        - aperture diameter of the bowl [grid points]
%     focus_pos       - any point on the beam axis of the bowl given as a
%                       three element vector [fx, fy, fz] [grid points] 
%
% OPTIONAL INPUTS:
%     Optional 'string', value pairs that may be used to modify the default
%     computational settings. 
%
%     'Binary'        - Boolean controlling whether the bowl map is
%                       returned as a double precision matrix (false) or
%                       a logical matrix (true) (default = false).
%     'Plot'          - Boolean controlling whether the bowl is plotted
%                       using voxelPlot (default = false).
%     'RemoveOverlap' - Boolean controlling whether overlapped grid points
%                       are removed (default = false).
%
% OUTPUTS:
%     bowl            - 3D binary map of a bowl
%
% ABOUT:
%     author          - Bradley Treeby and Yan To Ling
%     date            - 17th November 2014
%     last update     - 7th April 2017
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2014-2017 Bradley Treeby and Yan To Ling
%
% See also makeArc, makeBall, makeCartBowl, makeMultiBowl, makeSphere,
% makeSphericalSection

% This file is part of k-Wave. k-Wave is free software: you can
% redistribute it and/or modify it under the terms of the GNU Lesser
% General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% 
% k-Wave is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
% more details. 
% 
% You should have received a copy of the GNU Lesser General Public License
% along with k-Wave. If not, see <http://www.gnu.org/licenses/>.

% =========================================================================
% DEFINE LITERALS
% =========================================================================

% threshold used to find the closest point to the radius
THRESHOLD               = 0.5;

% number of grid points to expand the bounding box compared to
% sqrt(2)*diameter 
BOUNDING_BOX_EXP        = 2;

% optional input defaults
NUM_REQUIRED_INPUTS     = 5;
PLOT_BOWL_DEF           = false;
BINARY_DEF              = false;
REMOVE_OVERLAP_DEF      = false;

% =========================================================================
% INPUT CHECKING
% =========================================================================

% assign optional inputs with defaults
plot_bowl               = PLOT_BOWL_DEF;
binary                  = BINARY_DEF;
remove_overlap          = REMOVE_OVERLAP_DEF;

% replace with user defined values if provided
if nargin < NUM_REQUIRED_INPUTS
    error('Incorrect number of inputs.');
elseif ~isempty(varargin)
    for input_index = 1:2:length(varargin)
        switch varargin{input_index}
            case 'Binary'
                
                % assign input
                binary = varargin{input_index + 1};
                
                % check value
                if ~islogical(binary)
                    error('Optional input ''Binary'' must be Boolean.');
                end
                
            case 'Plot'
                
                % assign input
                plot_bowl = varargin{input_index + 1};
                
                % check value
                if ~islogical(plot_bowl)
                    error('Optional input ''Plot'' must be Boolean.');
                end
                
            case 'RemoveOverlap'
                
                % assign input
                remove_overlap = varargin{input_index + 1};
                
                % check value
                if ~islogical(remove_overlap)
                    error('Optional input ''RemoveOverlap'' must be Boolean.');
                end
                
            otherwise
                error('Unknown optional input.');
        end
    end
end

% assign appropriate flipfunc function
if exist('flip', 'builtin')
    flipfunc = @(x,dim) flip(x, dim);
else
    flipfunc = @(x,dim) flipdim(x, dim); %#ok<DFLIPDIM>
end

% force integer input values
grid_size = round(grid_size);
bowl_pos  = round(bowl_pos);
focus_pos = round(focus_pos);
diameter  = round(diameter);
radius    = round(radius);

% check the input ranges
if any(grid_size < 1)
    error('The grid size must be positive.');
end
if any(bowl_pos < 1) || any(bowl_pos > grid_size)
    error('The centre of the bowl must be within the grid.');
end
if radius <= 0
    error('The radius must be positive.');
end
if diameter <= 0
    error('The diameter must be positive.');
end
if diameter > 2*radius
    error('The diameter of the bowl must be less than twice the radius of curvature.');
end
if ~rem(diameter, 2)
    error('The diameter must be an odd number of grid points.');
end
if all(bowl_pos == focus_pos)
    error('The focus_pos must be different to the bowl_pos.');
end

% =========================================================================
% BOUND THE GRID TO SPEED UP CALCULATION
% =========================================================================

% create bounding box slightly larger than bowl diameter * sqrt(2)
Nx = round(sqrt(2) * diameter) + BOUNDING_BOX_EXP;
Ny = Nx;
Nz = Nx;
grid_size_sm = [Nx, Ny, Nz];

% set the bowl position to be the centre of the bounding box
bx = ceil(Nx/2);
by = ceil(Ny/2);
bz = ceil(Nz/2);
bowl_pos_sm = [bx, by, bz];

% set the focus position to be in the direction specified by the user
fx = bx + (focus_pos(1) - bowl_pos(1));
fy = by + (focus_pos(2) - bowl_pos(2));
fz = bz + (focus_pos(3) - bowl_pos(3));
focus_pos_sm = [fx, fy, fz];

% preallocate storage variable
if binary
    bowl_sm = false(Nx, Ny, Nz);
else
    bowl_sm = zeros(Nx, Ny, Nz);
end

% =========================================================================
% CREATE DISTANCE MATRIX
% =========================================================================

if ~isinf(radius)

    % find half the arc angle
    half_arc_angle = asin(diameter ./ (2 .* radius));
    
    % find centre of sphere on which the bowl lies
    distance_cf = sqrt( (bx - fx).^2 + (by - fy).^2 + (bz - fz).^2 );
    cx = round(radius ./ distance_cf .* (fx - bx) + bx);
    cy = round(radius ./ distance_cf .* (fy - by) + by);
    cz = round(radius ./ distance_cf .* (fz - bz) + bz);
    c = [cx, cy, cz];    
    
    % generate matrix with distance from the centre
    pixel_map = makePixelMapPoint(grid_size_sm, c);
    
    % set search radius to bowl radius
    search_radius = radius;    
    
else
    
    % generate matrix with distance from the centre
    pixel_map = makePixelMapPlane(grid_size_sm, bowl_pos_sm - focus_pos_sm, bowl_pos_sm);
    
    % set search radius to 0 (the disc is flat)
    search_radius = 0;
    
end

% calculate distance from search radius
pixel_map = abs(pixel_map - search_radius);

% =========================================================================
% DIMENSION 1
% =========================================================================

% find the grid point that corresponds to the outside of the bowl in the
% first dimension in both directions (the index gives the distance along
% this dimension) 
[value_forw, index_forw] = min(pixel_map, [], 1);
[value_back, index_back] = min(flipfunc(pixel_map, 1), [], 1);

% extract the linear index in the y-z plane of the values that lie on the
% bowl surface 
yz_ind_forw = find(value_forw < THRESHOLD);
yz_ind_back = find(value_back < THRESHOLD);

% use these subscripts to extract the x-index of the grid points that lie
% on the bowl surface 
x_ind_forw = index_forw(yz_ind_forw);
x_ind_back = index_back(yz_ind_back);        

% convert the linear index to equivalent subscript values
[y_ind_forw, z_ind_forw] = ind2sub([Ny, Nz], yz_ind_forw);
[y_ind_back, z_ind_back] = ind2sub([Ny, Nz], yz_ind_back);

% combine x-y-z indices into a linear index
linear_index_forw = sub2ind([Nx, Ny, Nz], x_ind_forw, y_ind_forw, z_ind_forw);
linear_index_back = sub2ind([Nx, Ny, Nz], Nx + 1 - x_ind_back, y_ind_back, z_ind_back);

% assign these values to the bowl
bowl_sm(linear_index_forw) = 1;
bowl_sm(linear_index_back) = 1;

% set existing bowl values to a distance of zero in the pixel map (this
% avoids problems with overlapping pixels)
pixel_map(bowl_sm == 1) = 0;

% =========================================================================
% DIMENSION 2
% =========================================================================

% find the grid point that corresponds to the outside of the bowl in the
% second dimension in both directions (the pixel map is first re-ordered to
% [X, Y, Z] -> [Y, Z, X]) 
[value_forw, index_forw] = min(permute(pixel_map, [2, 3, 1]), [], 1);
[value_back, index_back] = min(flipfunc(permute(pixel_map, [2, 3, 1]), 1), [], 1);

% extract the linear index in the y-z plane of the values that lie on the
% bowl surface 
zx_ind_forw = find(value_forw < THRESHOLD);
zx_ind_back = find(value_back < THRESHOLD);

% use these subscripts to extract the y-index of the grid points that lie
% on the bowl surface 
y_ind_forw = index_forw(zx_ind_forw);
y_ind_back = index_back(zx_ind_back);        

% convert the linear index to equivalent subscript values
[z_ind_forw, x_ind_forw] = ind2sub([Nz, Nx], zx_ind_forw);
[z_ind_back, x_ind_back] = ind2sub([Nz, Nx], zx_ind_back);

% combine x-y-z indices into a linear index
linear_index_forw = sub2ind([Nx, Ny, Nz], x_ind_forw, y_ind_forw, z_ind_forw);
linear_index_back = sub2ind([Nx, Ny, Nz], x_ind_back, Ny + 1 - y_ind_back, z_ind_back);

% assign these values to the bowl
bowl_sm(linear_index_forw) = 1;
bowl_sm(linear_index_back) = 1;        

% set existing bowl values to a distance of zero in the pixel map (this
% avoids problems with overlapping pixels)
pixel_map(bowl_sm == 1) = 0;

% =========================================================================
% DIMENSION 3
% =========================================================================

% find the grid point that corresponds to the outside of the bowl in the
% third dimension in both directions (the pixel map is first re-ordered to
% [X, Y, Z] -> [Z, X, Y]) 
[value_forw, index_forw] = min(permute(pixel_map, [3, 1, 2]), [], 1);
[value_back, index_back] = min(flipfunc(permute(pixel_map, [3, 1, 2]), 1), [], 1);

% extract the linear index in the y-z plane of the values that lie on the
% bowl surface 
xy_ind_forw = find(value_forw < THRESHOLD);
xy_ind_back = find(value_back < THRESHOLD);

% use these subscripts to extract the z-index of the grid points that lie
% on the bowl surface 
z_ind_forw = index_forw(xy_ind_forw);
z_ind_back = index_back(xy_ind_back);        

% convert the linear index to equivalent subscript values
[x_ind_forw, y_ind_forw] = ind2sub([Nx, Ny], xy_ind_forw);
[x_ind_back, y_ind_back] = ind2sub([Nx, Ny], xy_ind_back);

% combine x-y-z indices into a linear index
linear_index_forw = sub2ind([Nx, Ny, Nz], x_ind_forw, y_ind_forw, z_ind_forw);
linear_index_back = sub2ind([Nx, Ny, Nz], x_ind_back, y_ind_back, Nz + 1 - z_ind_back);

% assign these values to the bowl
bowl_sm(linear_index_forw) = 1;
bowl_sm(linear_index_back) = 1;

% =========================================================================
% RESTRICT SPHERE TO BOWL
% =========================================================================

% remove grid points within the sphere that do not form part of the bowl
if ~isinf(radius)
        
    % form vector from the geometric bowl centre to the back of the bowl 
    v1 = bowl_pos_sm - c;
    
    % calculate length of vector
    l1 = sqrt(sum((bowl_pos_sm - c).^2));

    % loop through the non-zero elements in the bowl matrix
    bowl_ind = find(bowl_sm == 1);
    for ind = 1:length(bowl_ind)
        
        % extract the indices of the current point
        [x_ind, y_ind, z_ind] = ind2sub([Nx, Ny, Nz], bowl_ind(ind));
        p = [x_ind, y_ind, z_ind];

        % form vector from the geometric bowl centre to the current point
        % on the bowl 
        v2 = p - c;
        
        % calculate length of vector
        l2 = sqrt(sum((p - c).^2));
        
        % find the angle between the two vectors using the dot product,
        % normalised using the vector lengths
        theta = acos(sum( v1 .* v2 ./ (l1 .* l2) ));
        
%         % alternative calculation normalised using radius of curvature
%         theta2 = acos(sum( v1 .* v2 ./ radius.^2 ));
        
        % if the angle is greater than the half angle of the bowl, remove
        % it from the bowl
        if theta > half_arc_angle
            bowl_sm(bowl_ind(ind)) = 0;
        end

    end    
    
else
    
    % form a distance map from the centre of the disc
    pixelMapPoint = makePixelMapPoint(grid_size_sm, bowl_pos_sm);
    
    % set all points in the disc greater than the diameter to zero
    bowl_sm(pixelMapPoint > (diameter ./ 2)) = 0;  
    
end

% =========================================================================
% REMOVE OVERLAPPED POINTS
% =========================================================================

if remove_overlap

    % define the shapes that capture the overlapped points, along with the
    % corresponding mask of which point to delete
    overlap_shapes{1}               = zeros(3, 3, 3);
    overlap_shapes{1}(1, 1, :)      = 1;
    overlap_shapes{1}(2, 2, :)      = 1;
    overlap_shapes{1}(3, 3, :)      = 1;
    overlap_shapes{1}(1, 2, 2)      = 1;
    overlap_shapes{1}(2, 3, 2)      = 1;

    overlap_delete{1}               = zeros(3, 3, 3);
    overlap_delete{1}(1, 2, 2)      = 1;
    overlap_delete{1}(2, 3, 2)      = 1;

    overlap_shapes{2}               = zeros(3, 3, 3);
    overlap_shapes{2}(1, 1, :)      = 1;
    overlap_shapes{2}(2, 2, :)      = 1;
    overlap_shapes{2}(3, 3, :)      = 1;
    overlap_shapes{2}(1, 2, 2)      = 1;

    overlap_delete{2}               = zeros(3, 3, 3);
    overlap_delete{2}(1, 2, 2)      = 1;

    overlap_shapes{3}               = zeros(3, 3, 3);
    overlap_shapes{3}(1:2, 1, :)    = 1;
    overlap_shapes{3}(3, 2, :)      = 1;
    overlap_shapes{3}(2, 2, 2)      = 1;

    overlap_delete{3}               = zeros(3, 3, 3);
    overlap_delete{3}(2, 2, 2)      = 1;

    overlap_shapes{4}               = zeros(3, 3, 3);
    overlap_shapes{4}(1, 1, :)      = 1;
    overlap_shapes{4}(2, 2, :)      = 1;
    overlap_shapes{4}(3, 3, :)      = 1;
    overlap_shapes{4}(1, 2, 1)      = 1;

    overlap_delete{4}               = zeros(3, 3, 3);
    overlap_delete{4}(1, 2, 1)      = 1;

    overlap_shapes{5}               = zeros(3, 3, 3);
    overlap_shapes{5}(1:2, 2, :)    = 1;
    overlap_shapes{5}(3, 3, :)      = 1;
    overlap_shapes{5}(3, 2, 1)      = 1;

    overlap_delete{5}               = zeros(3, 3, 3);
    overlap_delete{5}(3, 2, 1)      = 1;

    overlap_shapes{6}               = zeros(3, 3, 3);
    overlap_shapes{6}(1, :, 3)      = 1;
    overlap_shapes{6}(2, :, 2)      = 1;
    overlap_shapes{6}(2, :, 1)      = 1;
    overlap_shapes{6}(3, 2, 1)      = 1;

    overlap_delete{6}               = zeros(3, 3, 3);
    overlap_delete{6}(3, 2, 1)      = 1;

    overlap_shapes{7}               = zeros(3, 3, 3);
    overlap_shapes{7}(1, 3, :)      = 1;
    overlap_shapes{7}(2, 1:2, :)      = 1;
    overlap_shapes{7}(3, 1, 1)      = 1;

    overlap_delete{7}               = zeros(3, 3, 3);
    overlap_delete{7}(3, 1, 1)      = 1;

    overlap_shapes{8}               = zeros(3, 3, 3);
    overlap_shapes{8}(:, :, 2)      = 1;
    overlap_shapes{8}(1, 1, 1)      = 1;

    overlap_delete{8}               = zeros(3, 3, 3);
    overlap_delete{8}(1, 1, 1)      = 1;

    overlap_shapes{9}               = zeros(3, 3, 3);
    overlap_shapes{9}(1, :, 1)      = 1;
    overlap_shapes{9}(2, :, 2)      = 1;
    overlap_shapes{9}(2, :, 3)      = 1;
    overlap_shapes{9}(2, 2, 1)      = 1;

    overlap_delete{9}               = zeros(3, 3, 3);
    overlap_delete{9}(2, 2, 1)      = 1;

    overlap_shapes{10}               = zeros(3, 3, 3);
    overlap_shapes{10}(2:3, 3, 1)    = 1;
    overlap_shapes{10}(1, 3, 2:3)    = 1;
    overlap_shapes{10}(1, 2, 3)      = 1;
    overlap_shapes{10}(2, 2, 2)      = 1;
    overlap_shapes{10}(3, 2, 1)      = 1;
    overlap_shapes{10}(2:3, 1, 2)    = 1;
    overlap_shapes{10}(2, 1, 3)      = 1;


    overlap_delete{10}               = zeros(3, 3, 3);
    overlap_delete{10}(2, 1, 2)      = 1;

    % set loop flag
    points_remaining = true;
    
    % initialise deleted point counter
    deleted_points = 0;

    % set list of possible permutations
    perm_list = [...
        1 2 3; ...
        1 3 2; ...
        2 3 1; ...
        2 1 3; ...
        3 1 2; ...
        3 2 1];

    while points_remaining

        % get linear index of non-zero bowl elements
        index_mat = find(bowl_sm > 0);

        % set Boolean delete variable
        delete_point = false;    

        % loop through all points on the bowl, and find the all the points with
        % more than 8 neighbours
        for index = 1:length(index_mat)

            % extract subscripts for current point
            [cx, cy, cz] = ind2sub([Nx, Ny, Nz], index_mat(index));    

            % ignore edge points
            if (cx > 1) && (cx < Nx) && (cy > 1) && (cy < Ny) && (cz > 1) && (cz < Nz)

                % extract local region around current point
                region = bowl_sm(cx-1:cx+1, cy-1:cy+1, cz-1:cz+1); 

                % if there's more than 8 neighbours, check the point for
                % deletion
                if (sum(region(:)) - 1) > 8

                    % loop through the different shapes
                    for shape_index = 1:length(overlap_shapes)

                        % check every permutation of the shape, and apply the
                        % deletion mask if the pattern matches

                        % loop through possible shape permutations
                        for ind1 = 1:size(perm_list, 1)

                            % get shape and delete mask
                            overlap_s = overlap_shapes{shape_index};
                            overlap_d = overlap_delete{shape_index};

                            % permute
                            overlap_s = permute(overlap_s, perm_list(ind1, :));
                            overlap_d = permute(overlap_d, perm_list(ind1, :));

                            % loop through possible shape reflections
                            for ind2 = 1:7

                                % flipfunc the shape
                                switch ind2
                                    case 2
                                        overlap_s = flipfunc(overlap_s, 1);
                                        overlap_d = flipfunc(overlap_d, 1);
                                    case 3
                                        overlap_s = flipfunc(overlap_s, 2);
                                        overlap_d = flipfunc(overlap_d, 2);
                                    case 4
                                        overlap_s = flipfunc(overlap_s, 3);
                                        overlap_d = flipfunc(overlap_d, 3);
                                    case 5
                                        overlap_s = flipfunc(flipfunc(overlap_s, 1), 2);
                                        overlap_d = flipfunc(flipfunc(overlap_d, 1), 2);
                                    case 6
                                        overlap_s = flipfunc(flipfunc(overlap_s, 1), 3);
                                        overlap_d = flipfunc(flipfunc(overlap_d, 1), 3);
                                    case 7
                                        overlap_s = flipfunc(flipfunc(overlap_s, 2), 3);
                                        overlap_d = flipfunc(flipfunc(overlap_d, 2), 3);
                                end

                                % rotate the shape 4 x 90 degrees
                                for ind3 = 1:4

                                    % check if the shape matches
                                    if all(overlap_s(:) == region(:))
                                        delete_point = true;
                                    end

                                    % break from loop if a match is found
                                    if delete_point
                                        break;
                                    end

                                    % rotate shape
                                    overlap_s = rot90(overlap_s);
                                    overlap_d = rot90(overlap_d);

                                end

                                % break from loop if a match is found
                                if delete_point
                                    break;
                                end

                            end

                            % break from loop if a match is found
                            if delete_point
                                break;
                            end      

                        end                

                        % remove point from bowl if required, and update
                        % counter
                        if delete_point
                            bowl_sm(cx-1:cx+1, cy-1:cy+1, cz-1:cz+1) = bowl_sm(cx-1:cx+1, cy-1:cy+1, cz-1:cz+1) .* double(~overlap_d);
                            deleted_points = deleted_points + 1;
                            break
                        end

                    end

                end

            end

            % break from loop if a match is found
            if delete_point
                break;
            end

        end    

        % break from while loop if the outer for loop has completed
        % without deleting a point
        if index == length(index_mat)
            points_remaining = false;
        end

    end
    
    % display status
    if deleted_points
        disp([num2str(deleted_points) ' overlapped points removed from bowl']);
    end
    
end
    
% =========================================================================
% PLACE BOWL WITHIN LARGER GRID
% =========================================================================

% preallocate storage variable
if binary
    bowl = false(grid_size);
else
    bowl = zeros(grid_size);
end

% calculate position of bounding box within larger grid
x1 = bowl_pos(1) - bx + 1;
x2 = x1 + Nx - 1;
y1 = bowl_pos(2) - by + 1;
y2 = y1 + Ny - 1;
z1 = bowl_pos(3) - bz + 1;
z2 = z1 + Nz - 1;

% truncate bounding box if it falls outside the grid
if x1 < 1
    bowl_sm(1:abs(x1) + 1, :, :) = [];
    x1 = 1;
end
if y1 < 1
    bowl_sm(:, 1:abs(y1) + 1, :) = [];
    y1 = 1;
end
if z1 < 1
    bowl_sm(:, :, 1:abs(z1) + 1) = [];
    z1 = 1;
end
if x2 > grid_size(1)
    bowl_sm(end - (x2 - grid_size(1)) + 1:end, :, :) = [];
    x2 = grid_size(1);
end
if y2 > grid_size(2)
    bowl_sm(:, end - (y2 - grid_size(2)) + 1:end, :) = [];
    y2 = grid_size(2);
end
if z2 > grid_size(3)
    bowl_sm(:, :, end - (z2 - grid_size(3)) + 1:end) = [];
    z2 = grid_size(3);
end

% place bowl into grid
bowl(x1:x2, y1:y2, z1:z2) = bowl_sm;

% =========================================================================
% PLOT
% =========================================================================

% plot results using voxelPlot
if plot_bowl
    voxelPlot(double(bowl));
end
function [grid_data, order_index, reorder_index] = cart2grid(kgrid, cart_data, axisymmetric)
%CART2GRID Interpolate a set of Cartesian points onto a binary grid.
%
% DESCRIPTION:
%     cart2grid interpolates the set of Cartesian points defined by
%     cart_data onto a binary matrix defined by the kWaveGrid object
%     kgrid using nearest neighbour interpolation. An error is returned if
%     the Cartesian points are outside the computational domain defined by
%     kgrid.
%
% USAGE:
%     [grid_data, order_index, reorder_index] = cart2grid(kgrid, cart_data)
%
% INPUTS:
%     kgrid         - k-Wave grid object returned by kWaveGrid
%     cart_data     - 1 x N, 2 x N, or 3 x N (for 1, 2, and 3 dimensions)
%                     array of Cartesian sensor points 
%
% OPTIONAL INPUTS:
%     axisymmetric  - Boolean defining whether the grid in 2D is for an
%                     axisymmetric coordinate system (default = false)
%
% OUTPUTS:
%     grid_data     - binary matrix with the same dimensions as kgrid
%     order_index   - the order that the Cartesian points appear in
%                     grid_data according to MATLAB's standard column-wise
%                     linear matrix index ordering  
%     reorder_index - the order that the binary points in grid_data
%                     (according to MATLAB's standard column-wise linear
%                     matrix index ordering) appear in the original
%                     Cartesian data   
%
% ABOUT:
%     author        - Bradley Treeby
%     date          - 5th June 2009
%     last update   - 16th April 2018
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2018 Bradley Treeby
%
% See also interpCartData, kWaveGrid

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

% check for axisymmetric input
if nargin < 3 || isempty(axisymmetric)
    axisymmetric = false;
end
if axisymmetric && kgrid.dim ~= 2
    error('Axisymmetric flag only supported in 2D.');
end

% detect whether the inputs are for one, two, or three dimensions
switch kgrid.dim
    case 1
        
        % one-dimensional
        data_x = cart_data(1, :);     
    
        % scale position values to grid centered pixel coordinates using
        % nearest neighbour interpolation
        data_x = round(data_x./kgrid.dx);

        % shift pixel coordinates to coincide with matrix indexing
        data_x = data_x + floor(kgrid.Nx/2) + 1;
        
        % check if the points all lie within the grid
        if max(data_x(:)) > kgrid.Nx || min(data_x(:)) < 1
            error('Cartesian points must lie within the grid defined by kgrid.');
        end
        
        % create empty grid
        grid_data = zeros(kgrid.Nx, 1);

        % create index variable
        point_index = 1:length(data_x);

        % map values
        for data_index = 1:length(data_x)
            grid_data(data_x(data_index)) = point_index(data_index);
        end

        % extract reordering index
        reorder_index = reshape(grid_data(grid_data ~= 0), [], 1);        
        
    case 2
        
        % two-dimensional
        data_x = cart_data(1, :);
        data_y = cart_data(2, :);
        
        % scale position values to grid centered pixel coordinates using
        % nearest neighbour interpolation
        data_x = round(data_x./kgrid.dx);
        data_y = round(data_y./kgrid.dy);

        % shift pixel coordinates to coincide with matrix indexing (leave
        % y-direction = radial-direction if axisymmetric)
        data_x = data_x + floor(kgrid.Nx/2) + 1;
        if ~axisymmetric
            data_y = data_y + floor(kgrid.Ny/2) + 1;
        else
            data_y = data_y + 1;
        end
        % check if the points all lie within the grid
        if max(data_x(:)) > kgrid.Nx || max(data_y(:)) > kgrid.Ny || min([data_x(:); data_y(:)]) < 1
            error('Cartesian points must lie within the grid defined by kgrid.');
        end
        
        % create empty grid
        grid_data = zeros(kgrid.Nx, kgrid.Ny);

        % create index variable
        point_index = 1:length(data_x);

        % map values
        for data_index = 1:length(data_x)
            grid_data(data_x(data_index), data_y(data_index)) = point_index(data_index);
        end
        
        % extract reordering index
        reorder_index = reshape(grid_data(grid_data ~= 0), [], 1);
       
    case 3
        
        % three dimensional
        data_x = cart_data(1, :);
        data_y = cart_data(2, :);
        data_z = cart_data(3, :);
        
        % scale position values to grid centered pixel coordinates using
        % nearest neighbour interpolation
        data_x = round(data_x./kgrid.dx);
        data_y = round(data_y./kgrid.dy);
        data_z = round(data_z./kgrid.dz);

        % shift pixel coordinates to coincide with matrix indexing
        data_x = data_x + floor(kgrid.Nx/2) + 1;
        data_y = data_y + floor(kgrid.Ny/2) + 1;
        data_z = data_z + floor(kgrid.Nz/2) + 1;

        % check if the points all lie within the grid
        if max(data_x(:)) > kgrid.Nx || max(data_y(:)) > kgrid.Ny || max(data_z(:)) > kgrid.Nz || min([data_x(:); data_y(:); data_z(:)]) < 1
            error('Cartesian points must lie within the grid defined by kgrid.');
        end        
        
        % create empty grid
        grid_data = zeros(kgrid.Nx, kgrid.Ny, kgrid.Nz);

        % create index variable
        point_index = 1:length(data_x);

        % map values
        for data_index = 1:length(data_x)
            grid_data(data_x(data_index), data_y(data_index), data_z(data_index)) = point_index(data_index);
        end

        % extract reordering index
        reorder_index = reshape(grid_data(grid_data ~= 0), [], 1, 1);
      
    otherwise
        error('Input cart_data must be a 1, 2, or 3 dimensional matrix.');
end

% compute the reverse ordering index (i.e., what is the index of each point
% in the reordering vector)
order_index = ones(length(reorder_index), 2);
order_index(:, 1) = reorder_index;
order_index(:, 2) = 1:length(reorder_index);
order_index = sortrows(order_index, 1);
order_index = order_index(:, 2);

% reset binary grid values
grid_data(grid_data ~= 0) = 1;

% check if any Cartesian points have been mapped to the same grid point,
% thereby reducing the total number of points
num_discarded_points = size(cart_data, 2) - sum(grid_data(:));
if num_discarded_points ~= 0
    disp(['  cart2grid: ' num2str(num_discarded_points) ' Cartesian points mapped to overlapping grid points']);
end  
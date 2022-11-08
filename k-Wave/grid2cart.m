function [cart_data, order_index] = grid2cart(kgrid, grid_data)
%GRID2CART Return the Cartesian coordinates of the non-zero points of a binary grid.
%
% DESCRIPTION:
%     grid2cart returns the set of Cartesian coordinates corresponding to
%     the non-zero elements in the binary matrix grid_data, in the 
%     coordinate framework defined in kgrid.
%
% USAGE:
%     [cart_data, order_index] = grid2cart(kgrid, grid_data)
%
% INPUTS:
%     kgrid       - k-Wave grid object returned by kWaveGrid
%     grid_data   - binary grid with the same dimensions as the k-Wave
%                   grid kgrid 
%
% OUTPUTS:
%     cart_data   - 1 x N, 2 x N, or 3 x N (for 1, 2, and 3
%                   dimensions) array of Cartesian sensor points
%     order_index - the order that the Cartesian points appear in
%                   grid_data according to MATLAB's standard column-wise
%                   linear matrix index ordering  
%
% ABOUT:
%     author      - Ben Cox and Bradley Treeby
%     date        - 18th August 2009
%     last update - 4th June 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2017 Ben Cox and Bradley Treeby
%
% See also cart2grid, interpCartData, kWaveGrid

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

% set all non-zero elements in grid_data to 1
grid_data = (grid_data ~= 0);

% detect whether the inputs are for one, two, or three dimensions
% then return the Cartesian coordinate of grid_data's non-zero points
switch kgrid.dim
    case 1 
        
        % one-dimensional
        cart_data(1, :) = kgrid.x(grid_data);
        
    case 2 
        
        % two-dimensional
        cart_data(1, :) = kgrid.x(grid_data);
        cart_data(2, :) = kgrid.y(grid_data);
        
    case 3 
        
        % three dimensional
        cart_data(1, :) = kgrid.x(grid_data);
        cart_data(2, :) = kgrid.y(grid_data);
        cart_data(3, :) = kgrid.z(grid_data);
        
    otherwise
        error('Unknown input parameters.');
end

% find the index numbers of grid_data's non-zero points
order_index = find(grid_data ~= 0);
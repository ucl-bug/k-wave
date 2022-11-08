function pixel_map = makePixelMapPoint(grid_size, centre_pos)
%MAKEPIXELMAPPOINT Create matrix of grid point distances from a point.
%
% DESCRIPTION:
%     makePixelMapPoint generates a matrix populated with values of how far
%     each pixel in a grid is from the centre (given in pixel coordinates).
%     The centre position is controlled by centre_pos.
%
% USAGE:
%     pixel_map = makePixelMapPoint(grid_size, centre_pos)
%
% INPUTS:
%     grid_size       - size of Cartesian grid given as a two [Nx, Ny] or
%                       three [Nx, Ny, Nz] element vector [grid points] 
%     centre_pos      - coordinates of the centre given as a two [cx, cy]
%                       or three [cx, cy, cz] element vector [grid
%                       points]  
%
% OUTPUTS:
%     pixel_map       - pixel map
%
% ABOUT:
%     author          - Yan To Ling and Bradley Treeby
%     date            - 25th January 2015
%     last update     - 9th April 2017
%       
% See also makePixelMap, makePixelMapPlane
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2015-2017 Yan To Ling and Bradley Treeby

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

% check for number of dimensions
num_dim = length(grid_size);

% check that centre_pos has the same dimensions
if length(grid_size) ~= length(centre_pos)
    error('The inputs centre_pos and grid_size must have the same number of dimensions.');
end

switch num_dim
    case 2

        % assign inputs and force to be integers
        Nx = round(grid_size(1));
        Ny = round(grid_size(2));
        cx = round(centre_pos(1));
        cy = round(centre_pos(2));

        % generate index vectors in each dimension
        nx = reshape(0:Nx - 1, [Nx, 1]) - cx + 1;
        ny = reshape(0:Ny - 1, [1, Ny]) - cy + 1;

        % combine index matrices
        pixel_map = zeros(Nx, Ny);
        pixel_map = bsxfun(@plus, nx.^2, pixel_map);
        pixel_map = bsxfun(@plus, ny.^2, pixel_map);
        pixel_map = sqrt(pixel_map);

    case 3

        % assign inputs and force to be integers
        Nx = round(grid_size(1));
        Ny = round(grid_size(2));
        Nz = round(grid_size(3));
        cx = round(centre_pos(1));
        cy = round(centre_pos(2));
        cz = round(centre_pos(3));

        % generate index vectors in each dimension
        nx = reshape(0:Nx - 1, [Nx, 1, 1]) - cx + 1;
        ny = reshape(0:Ny - 1, [1, Ny, 1]) - cy + 1;
        nz = reshape(0:Nz - 1, [1, 1, Nz]) - cz + 1;

        % combine index matrices
        pixel_map = zeros(Nx, Ny, Nz);
        pixel_map = bsxfun(@plus, nx.^2, pixel_map);
        pixel_map = bsxfun(@plus, ny.^2, pixel_map);
        pixel_map = bsxfun(@plus, nz.^2, pixel_map);
        pixel_map = sqrt(pixel_map);
        
    otherwise
        
        % throw error
        error('Grid size must be 2 or 3D.');

end

function pixel_map = makePixelMapPlane(grid_size, normal, point)
%MAKEPIXELMAPPLANE Create matrix of grid point distances from a plane.
%
% DESCRIPTION:
%     Function to generate a matrix populated with values of how far each
%     pixel in a grid is from a plane.
%
% USAGE:
%     pixel_map = makePixelMapPlane(grid_size, normal, point)
%
% INPUTS:
%     grid_size     - size of Cartesian grid given as a two [Nx, Ny] or
%                     three [Nx, Ny, Nz] element vector [grid points]
%     normal        - a vector normal to the plane given as a three element
%                     vector [au] 
%     point         - a point on the plane given as a three element vector
%                     [px, py, pz] [grid points]  
%
% OUTPUTS:
%     pixel_map     - pixel map
%
% ABOUT:
%     author        - Yan To Ling and Bradley Treeby
%     date          - 9th December 2014
%     last update   - 9th April 2017
%       
% See also makePixelMap, makePixelMapPlane
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2014-2017 Yan To Ling and Bradley Treeby

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

% error checking
if all(normal == 0)
    error('Normal vector should not be zero.');
end

% check for number of dimensions
num_dim = length(grid_size);

switch num_dim
    case 2

        % assign inputs and force to be integers
        Nx = round(grid_size(1));
        Ny = round(grid_size(2));

        % create coordinate meshes
        [px, py]         = ndgrid(1:Nx, 1:Ny);
        [pointx, pointy] = ndgrid(ones(1, Nx) * point(1),  ones(1, Ny) * point(2));
        [nx, ny]         = ndgrid(ones(1, Nx) * normal(1), ones(1, Ny) * normal(2));
        
        % calculate distance according to Eq. (6) at
        % http://mathworld.wolfram.com/Point-PlaneDistance.html
        pixel_map = abs((px - pointx) .* nx + (py - pointy) .* ny) ./ sqrt(sum(normal.^2));
        
    case 3
        
        % assign inputs and force to be integers
        Nx = round(grid_size(1));
        Ny = round(grid_size(2));
        Nz = round(grid_size(3));

        % create coordinate meshes
        [px, py, pz]             = ndgrid(1:Nx, 1:Ny, 1:Nz);
        [pointx, pointy, pointz] = ndgrid(ones(1, Nx) * point(1),  ones(1, Ny) * point(2),  ones(1, Nz) * point(3));
        [nx, ny, nz]             = ndgrid(ones(1, Nx) * normal(1), ones(1, Ny) * normal(2), ones(1, Nz) * normal(3));
        
        % calculate distance according to Eq. (6) at
        % http://mathworld.wolfram.com/Point-PlaneDistance.html
        pixel_map = abs((px - pointx) .* nx + (py - pointy) .* ny + (pz - pointz) .* nz) ./ sqrt(sum(normal.^2));

    otherwise
        
        % throw error
        error('Grid size must be 2 or 3D.');
        
end
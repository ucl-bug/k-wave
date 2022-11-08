function sphere = makeSphere(Nx, Ny, Nz, radius, plot_sphere, binary)
%MAKESPHERE Create a binary map of a sphere within a 3D grid.
%
% DESCRIPTION:
%     makeSphere creates a binary map of a spherical shell (using an
%     extension of the midpoint circle algorithm) within a
%     three-dimensional grid. The sphere position is denoted by 1's in the
%     matrix with 0's elsewhere. If the Boolean input parameter 'binary' is
%     set to false (the default), the sphere map is returned as a double
%     precision matrix. If it is set to true, the map is returned as a
%     logical matrix. 
%
% USAGE:
%     sphere = makeSphere(Nx, Ny, Nz, radius)
%     sphere = makeSphere(Nx, Ny, Nz, radius, plot_sphere)
%     sphere = makeSphere(Nx, Ny, Nz, radius, plot_sphere, binary)
%     sphere = makeSphere(Nx, Ny, Nz, radius, [], binary)
%
% INPUTS:
%     Nx, Ny, Nz      - size of the 3D grid [grid points]
%     radius          - sphere radius [grid points]
%
% OPTIONAL INPUTS:
%     plot_sphere     - Boolean controlling whether the sphere is
%                       plotted using voxelPlot (default = false)
%     binary          - Boolean controlling whether the sphere map is
%                       returned as a double precision matrix (false) or
%                       a logical matrix (true) (default = false)
%
% OUTPUTS:
%     sphere          - 3D binary map of a sphere
%
% ABOUT:
%     author          - Bradley Treeby
%     date            - 16th June 2009
%     last update     - 8th June 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2017 Bradley Treeby
%
% See also makeBall, makeCartSphere, makeCircle, makeSphericalSection

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

% check for plot_sphere input
if nargin < 5 || isempty(plot_sphere)
    plot_sphere = false;
end

% check for binary input
if nargin < 6 || isempty(binary)
    binary = false;
end

% enforce a centered sphere
cx = floor(Nx / 2) + 1;
cy = floor(Ny / 2) + 1;
cz = floor(Nz / 2) + 1;

% preallocate the storage variable
if binary
    sphere = false(Nx, Ny, Nz);
else
    sphere = zeros(Nx, Ny, Nz);
end

% create a guide circle from which the individal radii can be extracted
guide_circle = makeCircle(Ny, Nx, cy, cx, radius);

% step through the guide circle points and create partially filled discs
centerpoints = (cx - radius):cx;
reflection_offset = length(centerpoints):-1:2;
for centerpoint_index = 1:length(centerpoints)
   
    % extract the current row from the guide circle
    row_data = guide_circle(:, centerpoints(centerpoint_index));

    % add an index to the grid points in the current row
    row_index = row_data .* (1:length(row_data)).';
    
    % calculate the radius 
    swept_radius = (max(row_index) - min(row_index(row_index ~= 0))) / 2;
    
    % create a circle to add to the sphere
    circle = makeCircle(Ny, Nz, cy, cz, swept_radius);

    % make an empty fill matrix
    if binary
        circle_fill = false(Ny, Nz);
    else
        circle_fill = zeros(Ny, Nz);
    end
    
    % fill in the circle line by line
    fill_centerpoints = (cz - swept_radius):(cz + swept_radius);
    for fill_centerpoint_index = 1:length(fill_centerpoints)
       
        % extract the first row
        row_data = circle(:, fill_centerpoints(fill_centerpoint_index));
        
        % add an index to the grid points in the current row
        row_index = row_data .* (1:length(row_data)).';
        
        % calculate the diameter
        start_index = min(row_index(row_index ~= 0));
        stop_index = max(row_index);
        
        % count how many points on the line
        num_points = sum(row_data);
        
        % fill in the line
        if start_index ~= stop_index && (stop_index - start_index) >= num_points
            circle_fill((start_index + num_points / 2):(stop_index - num_points / 2), fill_centerpoints(fill_centerpoint_index)) = 1;
        end
    end
        
    % remove points from the filled circle that existed in the previous
    % layer
    if centerpoint_index == 1
        sphere(centerpoints(centerpoint_index), :, :) = circle + circle_fill;
        prev_circle = circle + circle_fill;
    else
        prev_circle_alt = circle + circle_fill;
        circle_fill = circle_fill - prev_circle;
        circle_fill(circle_fill < 0) = 0;
        sphere(centerpoints(centerpoint_index), :, :) = circle + circle_fill;
        prev_circle = prev_circle_alt;
    end
    
    % create the other half of the sphere at the same time
    if centerpoint_index ~= length(centerpoints)
        sphere(cx + reflection_offset(centerpoint_index) - 1, :, :) = sphere(centerpoints(centerpoint_index), :, :);
    end
end

% plot results
if plot_sphere
    voxelPlot(double(sphere));
end
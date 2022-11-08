function [ss, dist_map] = makeSphericalSection(radius, height, width, plot_section, binary)
%MAKESPHERICALSECTION Create a binary map of a sphere segment within a 3D grid.
%
% DESCRIPTION:
%     makeSphericalSection creates a binary map of a section of a spherical
%     surface within a three-dimensional matrix. The sphere is created
%     using an extension of the midpoint circle algorithm. A single grid
%     point is taken as the sphere center so the total diameter will always
%     be an odd number of grid points. The sphere is then truncated based
%     on the values for height and width. The face of the spherical section
%     faces in the positive x-direction and the optional width parameter
%     truncates the size in the y-direction. 
%
%     If the optional input parameter "binary" is set to false (the
%     default), the section map is returned as a double precision matrix.
%     If it is set to true, the map is returned as a logical matrix. The
%     average distance between each grid point in the spherical section and
%     its contiguous neighbours can also be returned. This is given as a
%     ratio compared to the average neighbour distance for a flat surface.  
%
% SYNTAX:
%     ss = makeSphericalSection(radius, height)
%     ss = makeSphericalSection(radius, height, width)
%     ss = makeSphericalSection(radius, height, width, plot_section)
%     ss = makeSphericalSection(radius, height, [], plot_section)
%     ss = makeSphericalSection(radius, height, width, plot_section, binary)
%     ss = makeSphericalSection(radius, height, [], [] binary)
%
%     [ss, dist_map] = makeSphericalSection(radius, height)
%     [ss, dist_map] = makeSphericalSection(radius, height, width)
%     [ss, dist_map] = makeSphericalSection(radius, height, [], plot_section)
%     [ss, dist_map] = makeSphericalSection(radius, height, width, plot_section)
%     [ss, dist_map] = makeSphericalSection(radius, height, width, plot_section, binary)
%     [ss, dist_map] = makeSphericalSection(radius, height, [], [] binary)
%
% INPUTS:
%     radius          - radius of curvature [grid points]
%     height          - transducer height [grid points]
%
% OPTIONAL INPUTS:
%     width           - section width (must be specified as an odd
%                       number) [grid points] 
%     plot_section    - Boolean controlling whether the spherical section
%                       is plotted using voxelPlot (default = false) 
%     binary          - Boolean controlling whether the spherical section
%                       is returned as a double precision matrix (false)
%                       or a logical matrix (true) (default = false)
%
% OUTPUTS:
%     ss              - binary matrix containing spherical section
%     dist_map        - ratio of average neighbour distance for each grid
%                       point within the spherical section compared to a
%                       flat surface 
%
% ABOUT:
%     author          - Bradley Treeby
%     date            - 3rd February 2012
%     last update     - 9th April 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2017 Bradley Treeby
%
% See also makeBall, makeBowl, makeSphere

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

%#ok<*UNRCH>

% option to use the spherical sections or square sections of the sphere
% returned by makeSphere (true by default)
use_spherical_sections = true;

% force inputs to be integers
radius = round(radius);
height = round(height);

% check for width input
if nargin < 3 || (nargin == 4 && isempty(width))
    
    % set width truncation flag to false
    use_width = false;
    
else
    
    % set width truncation flag to true
    use_width = true;
    
    % force input to an integer
    width = round(width);
    
    % check that it's an odd number
    if ~rem(width, 2)
        error('Input width must be an odd number.');
    end
    
end

% check for plot input
if nargin < 4 || isempty(plot_section)
    plot_section = false;
end  

% check for data type input
if nargin < 5 || isempty(binary)
    binary = false;
end  

% calculate minimum grid dimensions to fit entire sphere
Nx = 2*radius + 1;

% create sphere
ss = makeSphere(Nx, Nx, Nx, radius, [], binary);

% truncate to given height
if use_spherical_sections
    ss = ss(1:height, :, :);
else
    ss = permute(ss(:, 1:height, :), [2, 3, 1]);
end

% flatten transducer and store the maximum and indices
mx = squeeze(max(ss, [], 1));

% calculate the total length/width of the transducer
length = sum(mx(ceil(end/2), :), 2);

% truncate transducer grid based on length (removes empty rows and columns)
offset = (Nx - length)/2;
ss = ss(:, 1 + offset:end - offset, 1 + offset:end - offset);

% also truncate to given width if defined by user
if use_width
    
    % check the value is appropriate
    if width > length
        error('Input for width must be less than or equal to transducer length.');
    end
    
    % calculate offset
    offset = (length - width)/2;
    
    % truncate transducer grid
    ss = ss(:, 1 + offset:end - offset, :); 
    
end
        
% compute average distance between each grid point and its contiguous
% neighbours if dist_map output is defined
if nargout == 2    

    % calculate x-index of each grid point in the spherical section, create
    % mask and remove singleton dimensions 
    [mx, mx_ind] = max(ss, [], 1);
    mask = squeeze(mx ~= 0);
    mx_ind = squeeze(mx_ind).*mask;     
    
    % double check there there is only one value of spherical section in
    % each matrix column
    if sum(mx(:)) ~= sum(ss(:))
        error('mean neighbour distance cannot be calculated uniquely due to overlapping points in the x-direction');
    end  
    
    % calculate average distance to grid point neighbours in the flat case
    x_dist = repmat([1 0 1], 3, 1);
    y_dist = x_dist.';
    flat_dist = sqrt(x_dist.^2 + y_dist.^2);
    flat_dist = mean(flat_dist(:));
    
    % compute distance map 
    dist_map = zeros(size(mx_ind));
    sz = size(mx_ind);
    for m = 1:sz(1)
        for n = 1:sz(2)

            % clear map
            local_heights = zeros(3, 3);

            % extract the height (x-distance) of the 8 neighbouring grid
            % points
            if m == 1 && n == 1
                local_heights(2:3, 2:3) = mx_ind(m:m + 1, n:n + 1);
            elseif m == sz(1) && n == sz(2)
                local_heights(1:2, 1:2) = mx_ind(m - 1:m, n - 1:n);
            elseif m == 1 && n == sz(2)
                local_heights(2:3, 1:2) = mx_ind(m:m + 1, n - 1:n);
            elseif m == sz(1) && n == 1
                local_heights(1:2, 2:3) = mx_ind(m - 1:m, n:n + 1);
            elseif m == 1
                local_heights(2:3, :) = mx_ind(m:m + 1, n - 1:n + 1);
            elseif m == sz(1)
                local_heights(1:2, :) = mx_ind(m - 1:m, n - 1:n + 1);
            elseif n == 1
                local_heights(:, 2:3) = mx_ind(m - 1:m + 1, n:n + 1);
            elseif n == sz(2)
                local_heights(:, 1:2) = mx_ind(m - 1:m + 1, n - 1:n);
            else
                local_heights = mx_ind(m - 1:m + 1, n - 1:n + 1);
            end

            % compute average variation from center
            local_heights_var = abs(local_heights - local_heights(2, 2));      

            % threshold no neighbours
            local_heights_var(local_heights == 0) = 0;

            % calculate total distance from centre
            dist = sqrt(x_dist.^2 + y_dist.^2 + local_heights_var.^2);

            % average and store as a ratio
            dist_map(m, n) = 1 + (mean(dist(:)) - flat_dist)./flat_dist;

        end 
    end

    % threshold out the non-transducer grid points
    dist_map(mask ~= 1) = 0;
    
    % flattened plot and distance values
    if plot_section 
        figure;
        subplot(2, 1, 1);
        imagesc(mx_ind);
        axis image;
        colorbar;
        title('Height of Transducer Face');
        subplot(2, 1, 2);
        dist_map_plot = dist_map;
        dist_map_plot(dist_map_plot ~= 0) = dist_map_plot(dist_map_plot ~= 0) -1;
        imagesc(100*dist_map_plot);
        axis image;
        colorbar;
        title('Percentage Increase in Average Neighbour Distance Compared to a Flat Surface');
    end

end

% plot if required 
if plot_section 

    voxelPlot(double(ss));
    view(150, 20);
    
%     % create surface plot
%     figure;
%     sz = size(ss);
%     [x, y] = meshgrid(1:sz(3),1:sz(2));
%     surface_ind = find(mx ~= 0);   
%     plot3(x(surface_ind), y(surface_ind), mx_ind(mx ~= 0), 'k.');
%     axis image;

end
function binary_sensor_data = interpCartData(kgrid, cart_sensor_data, cart_sensor_mask, binary_sensor_mask, interp)
%INTERPCARTDATA Interpolate data from a Cartesian to a binary sensor mask.
%
% DESCRIPTION:
%     interpCartData takes a matrix of time-series data recorded over a set
%     of Cartesian sensor points given by cart_sensor_mask and computes the
%     equivalent time-series at each sensor position on the binary sensor
%     mask binary_sensor_mask using interpolation. The properties of
%     binary_sensor_mask are defined by the k-Wave grid object kgrid.
%     Two and three dimensional data are supported.
%
% USAGE:
%     binary_sensor_data = interpCartData(kgrid, cart_sensor_data, cart_sensor_mask, binary_sensor_mask)
%     binary_sensor_data = interpCartData(kgrid, cart_sensor_data, cart_sensor_mask, binary_sensor_mask, interp)
%
% INPUTS:
%     kgrid               - k-Wave grid object returned by kWaveGrid
%     cart_sensor_data    - original sensor data measured over 
%                           cart_sensor_mask indexed as
%                           cart_sensor_data(sensor position, time)  
%     cart_sensor_mask    - Cartesian sensor mask over which
%                           cart_sensor_data is measured 
%     binary_sensor_mask  - binary sensor mask at which equivalent
%                           time-series are computed via interpolation 
%
% OPTIONAL INPUTS:
%     interp              - interpolation mode used to compute the
%                           time-series, both 'nearest' and 'linear'
%                           (two-point) modes are supported 
%                           (default = 'nearest') 
%
% OUTPUTS:
%     binary_sensor_data  - array of time-series corresponding to the
%                           sensor positions given by binary_sensor_mask  
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 14th July 2009
%     last update - 4th June 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2017 Bradley Treeby
%
% See also cart2grid, kspaceFirstOrder2D, kspaceFirstOrder3D, kWaveGrid

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

% start the clock
tic;

% assign defaults to optional inputs
if nargin < 5
    interp = 'nearest';
end

% extract the number of data points
[num_cart_data_points, num_time_points] = size(cart_sensor_data);

% update command line status
disp('Interpolating Cartesian sensor data...');
disp(['  interpolation mode: ' interp]);
disp(['  number of Cartesian sensor points: ' num2str(num_cart_data_points)]);
disp(['  number of binary sensor points: ' num2str(sum(binary_sensor_mask(:)))]);

switch numDim(kgrid.k)
    case 3
        
        % extract the cartesian coordinates of the binary sensor mask
        cart_bsm(1, :) = reshape(kgrid.x(binary_sensor_mask ~= 0), [], 1, 1);
        cart_bsm(2, :) = reshape(kgrid.y(binary_sensor_mask ~= 0), [], 1, 1);
        cart_bsm(3, :) = reshape(kgrid.z(binary_sensor_mask ~= 0), [], 1, 1);
        
        % preallocate storage variable
        num_binary_sensor_points = length(cart_bsm(1, :));
        binary_sensor_data = zeros(num_binary_sensor_points, num_time_points);
        
        % nearest neighbour interpolation of the data points
        if strcmp(interp, 'nearest')
            for point_index = 1:num_binary_sensor_points

                % find the measured data point that is closest
                dist = sqrt( (cart_bsm(1, point_index) - cart_sensor_mask(1, :)).^2 + (cart_bsm(2, point_index) - cart_sensor_mask(2, :)).^2 + (cart_bsm(3, point_index) - cart_sensor_mask(3, :)).^2 );
                [~, dist_min_index] = min(dist);

                % assign value
                binary_sensor_data(point_index, :) = cart_sensor_data(dist_min_index, :);

            end
        elseif strcmp(interp, 'linear')
            for point_index = 1:num_binary_sensor_points

                % find the distances to the measured data points
                dist = sqrt( (cart_bsm(1, point_index) - cart_sensor_mask(1, :)).^2 + (cart_bsm(2, point_index) - cart_sensor_mask(2, :)).^2 + (cart_bsm(3, point_index) - cart_sensor_mask(3, :)).^2 );

                % append the distance information onto the data set
                new_col_pos = length(cart_sensor_data(1, :)) + 1;
                cart_sensor_data_ro = cart_sensor_data;
                cart_sensor_data_ro(:, new_col_pos) = dist;

                % reorder the data set based on distance information
                cart_sensor_data_ro = sortrows(cart_sensor_data_ro, new_col_pos);

                % linearly interpolate between the two closest points
                perc = cart_sensor_data_ro(2, new_col_pos) / ( cart_sensor_data_ro(1, new_col_pos) + cart_sensor_data_ro(2, new_col_pos) );
                binary_sensor_data(point_index, :) = perc * cart_sensor_data_ro(1, 1:new_col_pos - 1) + (1 - perc)*cart_sensor_data_ro(2, 1:new_col_pos - 1);

            end
        else
            error('Unknown interpolation option.');
        end        
        
    case 2

        % extract the cartesian coordinates of the binary sensor mask
        cart_bsm(1, :) = reshape(kgrid.x(binary_sensor_mask ~= 0), [], 1);
        cart_bsm(2, :) = reshape(kgrid.y(binary_sensor_mask ~= 0), [], 1);

        % preallocate storage variable
        num_binary_sensor_points = length(cart_bsm(1, :));
        binary_sensor_data = zeros(num_binary_sensor_points, num_time_points);

        % nearest neighbour interpolation of the data points
        if strcmp(interp, 'nearest')
            for point_index = 1:num_binary_sensor_points

                % find the measured data point that is closest
                dist = sqrt( (cart_bsm(1, point_index) - cart_sensor_mask(1, :)).^2 + (cart_bsm(2, point_index) - cart_sensor_mask(2, :)).^2 );
                [~, dist_min_index] = min(dist);

                % assign value
                binary_sensor_data(point_index, :) = cart_sensor_data(dist_min_index, :);

            end
        elseif strcmp(interp, 'linear')
            for point_index = 1:num_binary_sensor_points

                % find the distances to the measured data points
                dist = sqrt( (cart_bsm(1, point_index) - cart_sensor_mask(1, :)).^2 + (cart_bsm(2, point_index) - cart_sensor_mask(2, :)).^2 );

                % append the distance information onto the data set
                new_col_pos = length(cart_sensor_data(1, :)) + 1;
                cart_sensor_data_ro = cart_sensor_data;
                cart_sensor_data_ro(:, new_col_pos) = dist;

                % reorder the data set based on distance information
                cart_sensor_data_ro = sortrows(cart_sensor_data_ro, new_col_pos);

                % linearly interpolate between the two closest points
                perc = cart_sensor_data_ro(2, new_col_pos) / ( cart_sensor_data_ro(1, new_col_pos) + cart_sensor_data_ro(2, new_col_pos) );
                binary_sensor_data(point_index, :) = perc * cart_sensor_data_ro(1, 1:new_col_pos - 1) + (1 - perc) * cart_sensor_data_ro(2, 1:new_col_pos - 1);

            end
        else
            error('Unknown interpolation option.');
        end
        
    otherwise
        error('Data must be two- or three-dimensional.');
end

% update command line status
disp(['  computation completed in ' scaleTime(toc)]);
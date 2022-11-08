function data = reorderBinarySensorData(data, reorder_index)
%REORDERBINARYSENSORDATA Reorder data from a binary sensor mask.
%
% DESCRIPTION:
%     reorderBinarySensorData reorders sensor data returned by the k-Wave
%     simulation functions (e.g., kspaceFirstOrder2D and
%     kspaceFirstOrder3D) according to a given reordering index. This is
%     useful when using a binary sensor mask in place of a Cartesian sensor
%     mask with nearest neighbour interpolation (Cartesian sensor masks are
%     not supported by the C++ codes).
%
%     Example syntax (assuming sensor.record is not defined)
%
%         % convert Cartesian sensor mask to binary sensor mask
%         [sensor.mask, ~, reorder_index] = cart2grid(kgrid, sensor.mask);
%
%         % run simulation
%         sensor_data = kspaceFirstOrder3DC(kgrid, medium, source, sensor)
%
%         % reorder output to order used in Cartesian sensor mask
%         sensor_data = reorderBinarySensorData(sensor_data, reorder_index);
%
%     If recording multiple outputs, then pass the individual outputs
%     (e.g., sensor_data.ux) to reorderBinarySensorData.
%
% USAGE:
%     data = reorderBinarySensorData(data, reorder_index)
%
% INPUTS:
%     data          - sensor data returned by the k-Wave simulation
%                     functions indexed as (sensor_point_index, time_index)
%     reorder_index - order of the Cartesian sensor points in the binary
%                     mask (this is returned by cart2grid)
%
% OUTPUTS:
%     data          - reordered data according to the order of the sensor
%                     points defined in the original Cartesian sensor mask
%
% ABOUT:
%     author        - Bradley E. Treeby
%     date          - 28th February 2019
%     last update   - 21st March 2019
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2019 Bradley Treeby

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

% calculate the position for reordering data
new_col_pos = length(data(1, :)) + 1;

% append the reordering data
data(:, new_col_pos) = reorder_index;

% reorder based on the order_index
data = sortrows(data, new_col_pos);

% remove the reordering data
data = data(:, 1:new_col_pos - 1);
function [reordered_sensor_data, indices_new] = reorderSensorData(kgrid, sensor, sensor_data)
%REORDERSENSORDATA Reorder sensor data from kspaceFirstOrder2D based on angle. 
%
% DESCRIPTION:
%     reorderSensorData reorders the time series from kspaceFirstOrder2D
%     based on the angle that each sensor point makes with the centre of
%     the grid. The sensor mask must be a binary mask, and the angles are
%     defined from the upper left quadrant or negative y-axis in the same
%     way as within makeCircle and makeCartCircle.
%
% USAGE:
%     [reordered_sensor_data, indices_new] = reorderSensorData(kgrid, sensor, sensor_data)
%
% INPUTS:
%     kgrid         - k-Wave grid object returned by kWaveGrid
%     sensor        - k-Wave sensor structure where sensor.mask is
%                     defined as binary grid
%     sensor_data   - sensor data returned by kspaceFirstOrder2D ordered
%                     using MATLAB's standard column-wise linear matrix
%                     indexing
%
% OUTPUTS:
%     reordered_sensor_data 
%                   - time varying sensor data reordered by the angle that
%                     each sensor point makes with the centre of the grid
%     indices_new   - the indices of the reordered time series, where
%                     reordered_sensor_data = sensor_data(indices_new, :) 
%
% ABOUT:
%     author        - Ben Cox
%     date          - 30th May 2012
%     last update   - 21st March 2019
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2012-2019 Ben Cox
%
% See also kspaceFirstOrder2D, reorderBinarySensorData, unmaskSensorData

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

% check simulation is 2D
if kgrid.dim ~= 2
    error('The simulation must be 2D.')
end

% check sensor.mask is a binary mask
if (size(kgrid(:), 2) ~= size(sensor.mask(:), 2))
    error('The sensor must be defined as a binary mask.')
end

% find the coordinates of the sensor points
x_sensor = kgrid.x(sensor.mask == 1);
y_sensor = kgrid.y(sensor.mask == 1);

% find the angle of each sensor point (from the centre)
angle = atan2(-x_sensor, -y_sensor);
angle(angle < 0) = 2 * pi + angle(angle < 0);

% sort the sensor points in order of increasing angle
[angle_sorted, indices_new] = sort(angle, 'ascend'); %#ok<ASGLU>

% reorder the measure time series so that adjacent time series correspond
% to adjacent sensor points.
reordered_sensor_data = sensor_data(indices_new, :);
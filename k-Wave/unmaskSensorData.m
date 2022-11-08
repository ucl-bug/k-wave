function unmasked_sensor_data = unmaskSensorData(kgrid, sensor, sensor_data)
%UNMASKSENSORDATA Reorder data recorded using a binary sensor mask.
%
% DESCRIPTION:
%     unmaskSensorData restores the grid position of the data recorded at a
%     single time-step within the time-series data returned by
%     kspaceFirstOrder1D, kspaceFirstOrder2D, or kspaceFirstOrder3D when
%     using a binary sensor mask.
%
% USAGE:
%     unmasked_sensor_data = unmaskSensorData(kgrid, sensor, sensor_data)
%
% INPUTS:
%     kgrid       - k-Wave grid object returned by kWaveGrid
%     sensor      - k-Wave sensor structure where sensor.mask is defined
%                   as a binary grid 
%     sensor_data - sensor data (returned by the first order simulation
%                   functions) at a single time-step ordered using
%                   MATLAB's standard column-wise linear matrix indexing 
%
% OUTPUTS:
%     unmasked_sensor_data 
%                 - matrix with the sensor data reordered to its original
%                   position on the sensor mask 
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 6th April 2009
%     last update - 7th June 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2017 Bradley Treeby
%
% See also kWaveGrid, kspaceFirstOrder1D, kspaceFirstOrder2D,
% kspaceFirstOrder3D

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

% check for reversed inputs (B.0.5 and earlier) in the form 
% (sensor_data, sensor_mask, kgrid), and call the function recursively
if ~isa(kgrid, 'kWaveGrid')
    disp('WARNING: Functions inputs for unmaskSensorData have changed. Please see documentation.');
    sens.mask = sensor;
    unmasked_sensor_data = unmaskSensorData(sensor_data, sens, kgrid);
    return
end

% create an empty matrix
switch numDim(kgrid.k)
    case 1
        unmasked_sensor_data = zeros(kgrid.Nx, 1);
    case 2
        unmasked_sensor_data = zeros(kgrid.Nx, kgrid.Ny);
    case 3
        unmasked_sensor_data = zeros(kgrid.Nx, kgrid.Ny, kgrid.Nz);
end

% reorder input data
unmasked_sensor_data(sensor.mask ~= 0) = sensor_data;
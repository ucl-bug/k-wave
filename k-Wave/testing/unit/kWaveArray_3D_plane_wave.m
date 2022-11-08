function test_pass = kWaveArray_3D_plane_wave(plot_comparisons, plot_simulations)
% DESCRIPTION:
%       Unit test to compare an off-grid line, disc, and rectangle
%       detectors created using the kWaveArray class with an on-grid point
%       detector for an incident plane wave in 3D.
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 15th July 2019
%       last update - 15th July 2019
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

% check for plot inputs, and set to true if nargin is zero (to allow the
% test to be run independent of runUnitTests)
if nargin == 0
    plot_comparisons = true;
    plot_simulations = true;
end

% set pass variable
test_pass = true;

% set comparison threshold
comparison_thresh = 2e-12;

% =========================================================================
% DEFINE GRID PROPERTIES
% =========================================================================

% grid properties
Nx = 96;
Ny = 48;
Nz = 48;
dx = 0.5e-3;
dy = dx;
dz = dx;
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% medium properties
medium.sound_speed = 1500;
medium.density = 1000;

% time array
kgrid.makeTime(medium.sound_speed);

% =========================================================================
% SIMULATION WITH GRID-BASED DETECTOR
% =========================================================================

% define sensor mask
sensor.mask = zeros(Nx, Ny, Nz);
x_pos = 20;
sensor.mask(x_pos, Ny/2, Nz/2) = 1;

% set source as a plane wave
source.u_mask = zeros(Nx, Ny, Nz);
source.u_mask(80, :, :) = 1;
source.ux = toneBurst(1/kgrid.dt, 300e3, 3) / 1.5e6;

% run k-Wave simulation
sensor_data_on_grid = kspaceFirstOrder3D(kgrid, medium, source, sensor, 'PlotSim', plot_simulations);

% =========================================================================
% SIMULATION WITH OFF-GRID LINE SENSOR
% =========================================================================

% create empty array
karray = kWaveArray('BLIType', 'exact', 'UpsamplingRate', 100);

% add line element
lat_pos1 = 12;
lat_pos2 = 36;
start_point = [kgrid.x_vec(x_pos), kgrid.y_vec(lat_pos1), kgrid.z_vec(lat_pos1)];
end_point   = [kgrid.x_vec(x_pos), kgrid.y_vec(lat_pos2), kgrid.z_vec(lat_pos2)];
karray.addLineElement(start_point, end_point);

% assign binary mask from karray to the sensor mask
sensor.mask = karray.getArrayBinaryMask(kgrid);

% run k-Wave simulation
sensor_data_off_grid = kspaceFirstOrder3D(kgrid, medium, source, sensor, ...
    'PlotSim', plot_simulations, ...
    'DisplayMask', 'off');

% combine data to give one trace per physical array element
sensor_data_line = karray.combineSensorData(kgrid, sensor_data_off_grid);

% =========================================================================
% SIMULATION WITH OFF-GRID RECT SENSOR
% =========================================================================

% create empty array
karray = kWaveArray('BLIType', 'exact', 'UpsamplingRate', 100);

% add rect element
lat_pos1 = 12;
lat_pos2 = 36;
position = [kgrid.x_vec(x_pos), kgrid.y_vec(Ny/2), kgrid.z_vec(Nz/2)];
Lx = (lat_pos1 - lat_pos2) * kgrid.dx;
Ly = Lx;
theta = [0, 90, 0];
karray.addRectElement(position, Lx, Ly, theta)

% assign binary mask from karray to the sensor mask
sensor.mask = karray.getArrayBinaryMask(kgrid);

% run k-Wave simulation
sensor_data_off_grid = kspaceFirstOrder3D(kgrid, medium, source, sensor, ...
    'PlotSim', plot_simulations, ...
    'DisplayMask', 'off');

% combine data to give one trace per physical array element
sensor_data_rect = karray.combineSensorData(kgrid, sensor_data_off_grid);

% =========================================================================
% SIMULATION WITH OFF-GRID DISC SENSOR
% =========================================================================

% create empty array
karray = kWaveArray('BLIType', 'exact', 'UpsamplingRate', 100);

% add disc element
position = [kgrid.x_vec(x_pos), kgrid.y_vec(Ny/2), kgrid.z_vec(Nz/2)];
diameter = 12 * kgrid.dx;
focus_pos = [0, kgrid.y_vec(Ny/2), kgrid.z_vec(Nz/2)];
karray.addDiscElement(position, diameter, focus_pos)

% assign binary mask from karray to the sensor mask
sensor.mask = karray.getArrayBinaryMask(kgrid);

% run k-Wave simulation
sensor_data_off_grid = kspaceFirstOrder3D(kgrid, medium, source, sensor, ...
    'PlotSim', plot_simulations, ...
    'DisplayMask', 'off');

% combine data to give one trace per physical array element
sensor_data_disc = karray.combineSensorData(kgrid, sensor_data_off_grid);

% =========================================================================
% CHECK ERROR
% =========================================================================

% compute error
err_line = 100 * max(abs(sensor_data_on_grid - sensor_data_line)) / max(sensor_data_on_grid(:));
err_rect = 100 * max(abs(sensor_data_on_grid - sensor_data_rect)) / max(sensor_data_on_grid(:));
err_disc = 100 * max(abs(sensor_data_on_grid - sensor_data_disc)) / max(sensor_data_on_grid(:));

% check values are the same
if (err_line > comparison_thresh) || (err_rect > comparison_thresh) || (err_disc > comparison_thresh)
    test_pass = false;
end

% =========================================================================
% VISUALISATION
% =========================================================================

% plot sensor traces
if plot_comparisons
    
    figure;
    subplot(2, 1, 1);
    plot(kgrid.t_array, sensor_data_on_grid, 'k-');
    hold on;
    plot(kgrid.t_array, sensor_data_line, 'r--');
    legend('on grid', 'off grid');
    title('Line Sensor');

    subplot(2, 1, 2);
    plot(kgrid.t_array, 100 * abs(sensor_data_on_grid - sensor_data_line) ./ max(sensor_data_on_grid(:)), 'k-');
    ylabel('Error [%]');
    
    figure;
    subplot(2, 1, 1);
    plot(kgrid.t_array, sensor_data_on_grid, 'k-');
    hold on;
    plot(kgrid.t_array, sensor_data_rect, 'r--');
    legend('on grid', 'off grid');
    title('Rectangle Sensor');

    subplot(2, 1, 2);
    plot(kgrid.t_array, 100 * abs(sensor_data_on_grid - sensor_data_rect) ./ max(sensor_data_on_grid(:)), 'k-');
    ylabel('Error [%]');
    
    figure;
    subplot(2, 1, 1);
    plot(kgrid.t_array, sensor_data_on_grid, 'k-');
    hold on;
    plot(kgrid.t_array, sensor_data_disc, 'r--');
    legend('on grid', 'off grid');
    title('Disc Sensor');

    subplot(2, 1, 2);
    plot(kgrid.t_array, 100 * abs(sensor_data_on_grid - sensor_data_disc) ./ max(sensor_data_on_grid(:)), 'k-');
    ylabel('Error [%]');
    
end
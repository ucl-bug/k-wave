function test_pass = kWaveArray_on_vs_off_grid_bowl_detector(plot_comparisons, plot_simulations)
% DESCRIPTION:
%       Unit test to compare an off-grid bowl detector created using the
%       kWaveArray class with a point detector.
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 5th July 2019
%       last update - 6th December 2019
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
comparison_thresh = 5;

% =========================================================================
% DEFINE GRID PROPERTIES
% =========================================================================

% grid properties
Nx = 64;
dx = 1e-3;
Ny = Nx;
dy = dx;
Nz = Nx;
dz = dx;
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% medium properties
medium.sound_speed = 1500;

% time array
kgrid.makeTime(medium.sound_speed, [], 40e-6);

% =========================================================================
% SIMULATION WITH GRID-BASED SOURCE
% =========================================================================

% define positions
src_x_pos = Nx/2;
src_y_pos = Ny/2;
src_z_pos = Nz/2;
sens_y_offset = Ny/4;

% set source as a ball shaped initial pressure
source.p0 = makeBall(Nx, Ny, Nz, src_x_pos, src_y_pos, src_z_pos, 3);

% define sensor mask
sensor.mask = zeros(Nx, Ny, Nz);
sensor.mask(src_x_pos, src_y_pos + sens_y_offset, src_z_pos) = 1;

% run k-Wave simulation
sensor_data_on_grid = kspaceFirstOrder3D(kgrid, medium, source, sensor, ...
    'PlotSim', plot_simulations, ...
    'PlotPML', false);

% =========================================================================
% SIMULATION WITH OFF-GRID SOURCE
% =========================================================================

% create empty array
karray = kWaveArray('BLITolerance', 0.01);

% add bowl element with focal position coincident with the source
position = [kgrid.x_vec(src_x_pos), kgrid.y_vec(src_y_pos + sens_y_offset), kgrid.z_vec(src_z_pos)];
radius = sens_y_offset * kgrid.dx;
diameter = radius / 2;
focus_pos = [kgrid.x_vec(src_x_pos), kgrid.y_vec(src_y_pos), kgrid.z_vec(src_z_pos)];
karray.addBowlElement(position, radius, diameter, focus_pos)

% assign binary mask from karray to the sensor mask
sensor.mask = karray.getArrayBinaryMask(kgrid);

% run k-Wave simulation
sensor_data_off_grid = kspaceFirstOrder3D(kgrid, medium, source, sensor, ...
    'PlotSim', plot_simulations, ...
    'PlotPML', false, ...
    'PMLInside', false, ...
    'PMLSize', 'auto');

% combine data to give one trace per physical array element
sensor_data_off_grid = karray.combineSensorData(kgrid, sensor_data_off_grid);

% compute error
err = 100 * max(abs(sensor_data_on_grid - sensor_data_off_grid)) / max(sensor_data_on_grid(:));

% check values are the same
if err > comparison_thresh
    test_pass = false;
end

% =========================================================================
% VISUALISATION
% =========================================================================

if plot_comparisons
    
    % plot sensor traces
    figure;
    subplot(2, 1, 1);
    plot(kgrid.t_array, sensor_data_on_grid, 'k-');
    hold on;
    plot(kgrid.t_array, sensor_data_off_grid, 'r--');
    legend('on grid', 'off grid');

    subplot(2, 1, 2);
    plot(kgrid.t_array, 100 * abs(sensor_data_on_grid - sensor_data_off_grid) ./ max(sensor_data_on_grid(:)), 'k-');
    ylabel('Error [%]');
    
    % plot the offgrid source and on-grid source
    figure;
    src_mask = karray.getArrayGridWeights(kgrid);
    src_mask(src_x_pos, src_y_pos + sens_y_offset, src_z_pos) = 2;
    plotSlices(src_mask);
    
end
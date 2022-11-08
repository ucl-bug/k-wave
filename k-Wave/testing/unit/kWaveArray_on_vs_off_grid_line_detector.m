function test_pass = kWaveArray_on_vs_off_grid_line_detector(plot_comparisons, plot_simulations)
% DESCRIPTION:
%       Unit test to compare an off-grid line detector created using the
%       kWaveArray class with a conventional on-grid line detector.
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 4th July 2019
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
comparison_thresh = 0.2;

% =========================================================================
% DEFINE GRID PROPERTIES
% =========================================================================

% grid properties
Nx = 128;
dx = 0.5e-3;
Ny = 128;
dy = 0.5e-3;
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% medium properties
medium.sound_speed = 1500;

% time array
kgrid.makeTime(medium.sound_speed);

% =========================================================================
% SIMULATION WITH GRID-BASED DETECTOR
% =========================================================================

% define sensor mask
sensor.mask = zeros(Nx, Ny);
x_pos = 30;
y1_pos = 50;
y2_pos = 78;
sensor.mask(x_pos, y1_pos:y2_pos) = 1;

% set source as a disc shaped initial pressure
source.p0 = makeDisc(Nx, Ny, 3*Nx/4, Ny/2, 10);

% run k-Wave simulation
sensor_data_on_grid = kspaceFirstOrder2D(kgrid, medium, source, sensor, 'PlotSim', plot_simulations);

% sum up the data, and divide by the number of elements
sensor_data_on_grid = sum(sensor_data_on_grid, 1) ./ sum(sensor.mask(:));

% =========================================================================
% SIMULATION WITH OFF-GRID DETECTOR
% =========================================================================

% create empty array
karray = kWaveArray;

% add line element (length offset by half a grid point)
start_point = [kgrid.x_vec(x_pos), kgrid.y_vec(y1_pos) - kgrid.dx/2];
end_point = [kgrid.x_vec(x_pos), kgrid.y_vec(y2_pos) + kgrid.dx/2];
karray.addLineElement(start_point, end_point);

% assign binary mask from karray to the sensor mask
sensor.mask = karray.getArrayBinaryMask(kgrid);

% run k-Wave simulation
sensor_data_off_grid = kspaceFirstOrder2D(kgrid, medium, source, sensor, 'PlotSim', plot_simulations);

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

% plot sensor traces
if plot_comparisons
    
    figure;
    subplot(2, 1, 1);
    plot(kgrid.t_array, sensor_data_on_grid, 'k-');
    hold on;
    plot(kgrid.t_array, sensor_data_off_grid, 'r--');
    legend('on grid', 'off grid');

    subplot(2, 1, 2);
    plot(kgrid.t_array, 100 * abs(sensor_data_on_grid - sensor_data_off_grid) ./ max(sensor_data_on_grid(:)), 'k-');
    ylabel('Error [%]');
    
end
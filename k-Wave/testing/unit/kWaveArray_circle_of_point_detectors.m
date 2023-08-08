function test_pass = kWaveArray_circle_of_point_detectors(plot_comparisons, plot_simulations)
% DESCRIPTION:
%       Unit test to compare a circle of off-grid point detectors.
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 31st May 2023
%       last update - 31st May 2023
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2023 Bradley Treeby

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
comparison_thresh = 3;

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
kgrid.makeTime(medium.sound_speed, [], 40e-6);

% set source
source.p_mask = zeros(Nx, Ny);
source.p_mask(Nx/2+1, Ny/2+1) = 1;
source.p = createCWSignals(kgrid.t_array, 500e3, 1, 0, 5);

% =========================================================================
% SIMULATION WITH OFF-GRID SOURCE
% =========================================================================

% create empty array
karray = kWaveArray('BLITolerance', 0.01);

% add elements
positions = makeCartCircle(15e-3, 100, [0, 0]);
measure = 1;
element_dim = 1;
for ind = 1:length(positions)
    karray.addCustomElement(positions(:, ind), measure, element_dim, 'Point');
end

% assign binary mask from karray to the sensor mask
sensor.mask = karray.getArrayBinaryMask(kgrid);

% run k-Wave simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, ...
    'PlotSim', plot_simulations, ...
    'PMLInside', false, ...
    'PlotPML', false);

% combine data to give one trace per physical array element
sensor_data = karray.combineSensorData(kgrid, sensor_data);

% check all the traces are the same
maxErr = 0;
for ind = 2:length(positions)
    err = 100 * max(abs(sensor_data(ind, :) - sensor_data(1, :))) / max(sensor_data(1, :));
    if err > maxErr
        maxErr = err;
        maxErrInd = ind;
    end
end

% check values are the same
if maxErr > comparison_thresh
    test_pass = false;
end

% =========================================================================
% VISUALISATION
% =========================================================================

% plot sensor traces
if plot_comparisons
    
    figure;
    imagesc(sensor.mask);
    title('Sensor mask');
    axis image;

    figure;
    imagesc(sensor_data);
    title('Sensor data');
    colorbar;

    figure;
    subplot(2, 1, 1);
    plot(kgrid.t_array, sensor_data(1, :), 'k-');
    hold on;
    plot(kgrid.t_array, sensor_data(maxErrInd, :), 'r--');
    legend('Point 1', 'Point 2');

    subplot(2, 1, 2);
    plot(kgrid.t_array, 100 * abs(sensor_data(1, :) - sensor_data(maxErrInd, :)) ./ max(sensor_data(1, :)), 'k-');
    ylabel('Error [%]');
    
end
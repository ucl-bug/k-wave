function test_pass = kspaceFirstOrderAS_compare_binary_and_cartesian_sensor_mask(plot_comparisons, plot_simulations)
% DESCRIPTION:
%     Unit test to compare cartesian and binary sensor masks for the
%     axisymmetric code.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 7th December 2017
%     last update - 21st March 2019
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2017-2019 Bradley Treeby

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

%#ok<*NOPRT>

% check for plot inputs, and set to true if nargin is zero (to allow the
% test to be run independent of runUnitTests)
if nargin == 0
    plot_comparisons = true;
    plot_simulations = true;
end

% set additional literals to give further permutations of the test
NN_COMPARISON_THRESH  = 1e-14;
LIN_COMPARISON_THRESH = 0.3;

% set pass variable
test_pass = true;

% create the computational grid
Nx = 128;           
Ny = 128;          
dx = 25e-3/Nx;    
dy = dx;          
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
medium.sound_speed = 1500 * ones(Nx, Ny);  % [m/s]
medium.sound_speed(Nx/2:end, :) = 2000;
medium.sound_speed_ref = 1500;

medium.density = 1000 * ones(Nx, Ny);
medium.density(Nx/2:end, :) = 1200;

% define source mask
source.p0 = zeros(Nx, Ny);
source.p0(22, 1:Ny/4) = 1;

% define Cartesian sensor points
num_points = 50;
radius = 7e-3;
sensor.mask = makeCartCircle(radius, num_points, [0, 0]);
sensor.mask(:, sensor.mask(2, :) < 0) = [];
sensor.mask(:, sensor.mask(1, :) < 0) = [];
sensor.mask = fliplr(sensor.mask);

% run the simulation as normal
sensor_data_c_lin = kspaceFirstOrderAS(kgrid, medium, source, sensor, ...
    'RadialSymmetry', 'WSWA-FFT', 'PlotSim', plot_simulations);

% run the simulation using nearest-neighbour interpolation
sensor_data_c_nn = kspaceFirstOrderAS(kgrid, medium, source, sensor, ...
    'RadialSymmetry', 'WSWA-FFT', 'CartInterp', 'nearest', 'PlotSim', plot_simulations);

% convert sensor mask
sensor.mask = cart2grid(kgrid, sensor.mask, true);

% run the simulation again
sensor_data_b = kspaceFirstOrderAS(kgrid, medium, source, sensor, ...
    'RadialSymmetry', 'WSWA-FFT', 'CartInterp', 'nearest', 'PlotSim', plot_simulations);

% compute errors
diff_nn = max(abs(sensor_data_c_nn(:) - sensor_data_b(:))) / max(abs(sensor_data_c_nn(:))) 
diff_lin = max(abs(sensor_data_c_lin(:) - sensor_data_b(:))) / max(abs(sensor_data_c_lin(:))) 

% check for test pass
if (diff_nn > NN_COMPARISON_THRESH) || (diff_lin > LIN_COMPARISON_THRESH)
    test_pass = false;
end

% plot
if plot_comparisons

    figure;
    subplot(5, 1, 1);
    imagesc(sensor_data_c_lin);
    colorbar;
    title('Cartesian - Linear');

    subplot(5, 1, 2);
    imagesc(sensor_data_c_nn);
    colorbar;
    title('Cartesian - Nearest Neighbour');

    subplot(5, 1, 3);
    imagesc(sensor_data_b);
    colorbar;
    title('Binary');

    subplot(5, 1, 4);
    imagesc(abs(sensor_data_b - sensor_data_c_lin));
    colorbar;
    title('Diff (Linear - Binary)');
    
    subplot(5, 1, 5);
    imagesc(abs(sensor_data_b - sensor_data_c_nn));
    colorbar;
    title('Diff (Nearest Neighbour - Binary)');
    
end
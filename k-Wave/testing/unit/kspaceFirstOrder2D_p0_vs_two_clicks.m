function test_pass = kspaceFirstOrder2D_p0_vs_two_clicks(plot_comparisons, plot_simulations)
% DESCRIPTION:
%     Compare source.p0 and source.p with two clicks scaled by CFL. These
%     shuold be identical if there is no absorption, the PML is turned off,
%     and the k-space source correction is turned off.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 15th May 2018
%     last update - 15th May 2018
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2018 Bradley Treeby

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

% set pass threshold
comparison_threshold = 1e-14;

% set pass variable
test_pass = true;

% =========================================================================
% 2D SIMULATION
% =========================================================================

% create the computational grid
Nx = 64;            % number of grid points in the x (row) direction
Ny = 64;            % number of grid points in the y (column) direction
dx = 0.1e-3;        % grid point spacing in the x direction [m]
dy = 0.1e-3;        % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
medium.sound_speed = 1500;  % [m/s]

% create time array
CFL = 0.234;
kgrid.makeTime(medium.sound_speed, CFL);

% create initial pressure distribution using makeDisc
disc_x_pos = 25;    % [grid points]
disc_y_pos = 25;    % [grid points]
disc_radius = 4;    % [grid points]
disc_1 = makeDisc(Nx, Ny, disc_x_pos, disc_y_pos, disc_radius);

disc_x_pos = 40;    % [grid points]
disc_y_pos = 30;    % [grid points]
disc_radius = 2;    % [grid points]
disc_2 = makeDisc(Nx, Ny, disc_x_pos, disc_y_pos, disc_radius);

p0 = disc_1 + disc_2;
source.p0 = p0;

% define a centered circular sensor
sensor.mask = makeCircle(Nx, Ny, 0, 0, Nx/4);

% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, ...
    'Smooth', false, ...
    'PMLAlpha', 0, ...
    'PlotSim', plot_simulations);

% redefine using two clicks, turning off the k-space source correction
clear source
source.p_mask = zeros(Nx, Ny);
source.p_mask(p0 ~= 0) = 1;
source.p_mode = 'additive-no-correction';
source.p = [0.5, 0.5] ./ (2 * medium.sound_speed * kgrid.dt / kgrid.dx);

% run the simulation
sensor_data2 = kspaceFirstOrder2D(kgrid, medium, source, sensor, ...
    'PMLAlpha', 0, ...
    'PlotSim', plot_simulations);

% check error
err = max(sensor_data(:) - sensor_data2(:));

% check for test pass
if err > comparison_threshold
    test_pass = false;
end

% =========================================================================
% PLOT
% =========================================================================

% plot the simulated sensor data
if plot_comparisons
    
    figure;

    subplot(3, 1, 1);
    imagesc(sensor_data, [-1, 1]);
    title('p0');
    colorbar;

    subplot(3, 1, 2);
    imagesc(sensor_data2, [-1, 1]);
    title('two clicks');
    colorbar;

    subplot(3, 1, 3);
    imagesc(abs(sensor_data - sensor_data2));
    title('difference');
    colorbar;
   
    colormap(getColorMap);
    
end
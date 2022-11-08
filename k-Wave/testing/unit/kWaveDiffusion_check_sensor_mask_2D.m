function test_pass = kWaveDiffusion_check_sensor_mask_2D(plot_comparisons, plot_simulations) 
% DESCRIPTION:
%       Unit test to compare the output from the k-Wave diffusion class
%       when using a sensor mask with the output manually extracted at each
%       time step.
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 16th August 2015
%       last update - 11th December 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2015-2017 Bradley Treeby

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

% set pass variable
test_pass = true;

% set comparison threshold
COMPARISON_THRESH = 1e-15;

% check for plot inputs, and set to true if nargin is zero (to allow the
% test to be run independent of runUnitTests)
if nargin == 0
    plot_comparisons = true;
    plot_simulations = true;
end

% =========================================================================
% SIMULATION
% =========================================================================

% create grid
Nx = 128;
Ny = 128;
dx = 1e-3;
dy = 1e-3;
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define medium properties
medium.density              = 1079;     % [kg/m^3]
medium.thermal_conductivity = 0.52;     % [W/(m.K)]
medium.specific_heat        = 3540;     % [J/(kg.K)]

% set Gaussian initial temperature distribution
offset = 0;
width = 4*dx;
source.T0 = 37 + exp(-( ((kgrid.x-offset)/width).^2 + ((kgrid.y-offset)/width).^2 ));

% define a circular sensor mask
sensor.mask = makeCircle(Nx, Ny, Nx/2, Ny/4, 10);

% create kWaveDiffusion object without sensor
kdiff = kWaveDiffusion(kgrid, medium, source, [], 'PlotSim', false);

% pre-allocate storage variable
dt = 0.5;
Nt = 100;
sensor_data_manual = zeros(sum(sensor.mask(:)), Nt);

% take time steps one by one
for ind = 1:Nt
    kdiff.takeTimeStep(1, dt);
    sensor_data_manual(:, ind) = kdiff.T(sensor.mask == 1);
end

% create kWaveDiffusion object with sensor
kdiff = kWaveDiffusion(kgrid, medium, source, sensor, 'PlotSim', plot_simulations);

% take time steps
kdiff.takeTimeStep(Nt, dt);

% compute error
error = abs(kdiff.sensor_data - sensor_data_manual);

% compute the maximum error
L_inf = max(error);

% compute pass
if (L_inf > COMPARISON_THRESH)
    test_pass = false;
end

% =========================================================================
% PLOT COMPARISONS
% =========================================================================

if plot_comparisons

    figure;
    
    subplot(3, 1, 1);
    imagesc(kdiff.sensor_data);
    colorbar;
    title('Sensor mask');
    
    subplot(3, 1, 2);
    imagesc(sensor_data_manual);
    colorbar;
    title('Manual');
    
    subplot(3, 1, 3);
    imagesc(error);
    colorbar;
    title('Error');
    
end
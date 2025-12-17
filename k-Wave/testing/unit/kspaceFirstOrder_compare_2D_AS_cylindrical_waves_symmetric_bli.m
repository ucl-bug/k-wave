function test_pass = kspaceFirstOrder_compare_2D_AS_cylindrical_waves_symmetric_bli(plot_comparisons, plot_simulations)
% DESCRIPTION:
%     Unit test to compare single frequency cylindrical waves in 2D and AS
%     when the BLI in the 2D code is forced to be symmetric.
%
% ABOUT:
%     author        - Bradley Treeby
%     date          - 10th October 2018
%     last update   - 18th December 2018
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2017- Bradley Treeby

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

% set additional literals to give further permutations of the test
USE_KSPACE  = true;
USE_PML     = true;
CFL         = 0.1;

% set pass threshold
comparison_threshold = 1e-2;

% set pass variable
test_pass = true;

% =========================================================================
% AXISYMMETRIC SIMULATIONS
% =========================================================================

% setup grid
sc = 1;
pml_size = 10;
Nx = 128 * sc;
Ny = 128 * sc;
source_sensor_distance = 32;
dx = 0.1e-3;
dy = dx;
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% medium
medium.sound_speed = 1500;

% time
t_end = 7e-6;
kgrid.makeTime(medium.sound_speed, CFL, t_end);

% source
source.p_mask = zeros(Nx, Ny);
source.p_mask(:, 1) = 1;

% source signal
freq = 2e6;
source_sig = createCWSignals(kgrid.t_array, freq, 1, 0);
source.p = source_sig;

% sensor
sensor.mask = zeros(Nx, Ny);
sensor.mask(Nx/2, source_sensor_distance) = 1;

% pml
pml_alpha = 2;
if ~USE_PML
    pml_alpha = 0;
end

% with source correction
sensor_data_as = kspaceFirstOrderAS(kgrid, medium, source, sensor, ...
    'PMLSize', [2, pml_size], ...
    'PMLAlpha', [0, pml_alpha], ...
    'RadialSymmetry', 'WSWA-FFT', ...
    'UsekSpace', USE_KSPACE, ...
    'PlotSim', plot_simulations);

% =========================================================================
% 2D k-Wave
% =========================================================================

% create a 2D BLI centered on the grid
r = sqrt(kgrid.x.^2 + kgrid.y.^2);
mask = sin(pi * r ./ dx) ./ (pi * r ./ dx);
mask(r == 0) = 1;

figure;
imagesc(mask);
axis image;
title('Symmetric BLI');

% assign source
source.p_mask = ones(Nx, Ny);
source.p = bsxfun( @times, reshape(mask(:), [], 1), reshape(source_sig, 1, []) );

% sensor
sensor.mask = zeros(Nx, Ny);
sensor.mask(Nx/2 + 1, Ny/2 + source_sensor_distance) = 1;

% record final pressure
sensor.record = {'p', 'p_final'};

% simulations
sensor_data_2D = kspaceFirstOrder2D(kgrid, medium, source, sensor, ...
    'PMLSize', pml_size, ...
    'PMLAlpha', pml_alpha, ...
    'UsekSpace', USE_KSPACE);

% =========================================================================
% ANALYSIS
% =========================================================================

% compute error
err = max(abs(sensor_data_as(:) - sensor_data_2D.p(:))) / max(abs(sensor_data_2D.p(:)));
    
% check for test pass
if (err > comparison_threshold)
    test_pass = false;
end   

% =========================================================================
% PLOTTING
% =========================================================================

if plot_comparisons

    % plot final pressure
    figure;
    subplot(4, 1, 1);
    imagesc(sensor_data_2D.p_final);
    axis image;
    title('Final Pressure');
    
    % plot k-space
    subplot(4, 1, 2);
    spect_2D = fftshift(abs(fft2(sensor_data_2D.p_final)));
    imagesc(spect_2D);
    axis image;
    title('k-space');

    % plot signal at the end of the domain
    subplot(4, 1, 3);
    plot(sensor_data_2D.p);
    hold on;
    plot(sensor_data_as);
    legend('2D', 'AS');
    title('Time series');
    
    % plot error
    subplot(4, 1, 4);
    plot(100 * abs(sensor_data_2D.p - sensor_data_as) ./ max(sensor_data_2D.p(:)));
    hold on;
    title('Error in time series [%]');
    
end
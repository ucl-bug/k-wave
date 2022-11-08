function test_pass = kspaceFirstOrderAS_compare_FFT_implementations(plot_comparisons, plot_simulations)
% DESCRIPTION:
%     Unit test to compare the WSWA and WSWS implementations with
%     WSWA-FFT and WSWS-FFT.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 7th December 2017
%     last update - 7th December 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2017 Bradley Treeby

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
USE_KSPACE          = true;
COMPARISON_THRESH   = 1e-14;

% set pass variable
test_pass = true;

% check for dtt functions
if ~exist('dtt1D', 'file') 
    disp('WARNING: kspaceFirstOrderAS_compare_FFT_implementations not tested, as dtt library is not present on path');
    return
end

% create the computational grid
Nx = 128;           
Ny = 64;          
dx = 25e-3/Nx;    
dy = dx;          
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
medium.sound_speed = 1500 * ones(Nx, Ny);  % [m/s]
medium.sound_speed(Nx/2:end, :) = 2000;
medium.sound_speed_ref = 1500;

medium.density = 1000 * ones(Nx, Ny);
medium.density(Nx/2:end) = 1200;

% define source mask
source.p0 = zeros(Nx, Ny);
source.p0(22, 1:Ny/4) = 1;

% define a single sensor point
sensor.mask = zeros(Nx, Ny);
sensor.mask(Nx/4, 10) = 1;

% run the simulations
sensor_data_wswa_fft = kspaceFirstOrderAS(kgrid, medium, source, sensor, ...
    'RadialSymmetry', 'WSWA-FFT', 'UsekSpace', USE_KSPACE, 'PlotSim', plot_simulations);
sensor_data_wswa     = kspaceFirstOrderAS(kgrid, medium, source, sensor, ...
    'RadialSymmetry', 'WSWA'    , 'UsekSpace', USE_KSPACE, 'PlotSim', plot_simulations);

sensor_data_wsws_fft = kspaceFirstOrderAS(kgrid, medium, source, sensor, ...
    'RadialSymmetry', 'WSWS-FFT', 'UsekSpace', USE_KSPACE, 'PlotSim', plot_simulations);
sensor_data_wsws     = kspaceFirstOrderAS(kgrid, medium, source, sensor, ...
    'RadialSymmetry', 'WSWS'    , 'UsekSpace', USE_KSPACE, 'PlotSim', plot_simulations);

% compute errors
diff_wswa = max(abs(sensor_data_wswa_fft(:) - sensor_data_wswa(:))) / max(abs(sensor_data_wswa_fft(:)))
diff_wsws = max(abs(sensor_data_wsws_fft(:) - sensor_data_wsws(:))) / max(abs(sensor_data_wsws_fft(:)))

% check for test pass
if (diff_wswa > COMPARISON_THRESH) || ...
   (diff_wsws > COMPARISON_THRESH)
    test_pass = false;
end

% plot
if plot_comparisons
    
    figure;
    subplot(2, 1, 1);
    plot(sensor_data_wswa_fft, 'k-');
    hold on;
    plot(sensor_data_wswa, 'r:');
    title('WSWA');
    subplot(2, 1, 2);
    plot(abs(sensor_data_wswa_fft - sensor_data_wswa), 'k-');

    figure;
    subplot(2, 1, 1);
    plot(sensor_data_wsws_fft, 'k-');
    hold on;
    plot(sensor_data_wsws, 'r:');
    title('WSWS');
    subplot(2, 1, 2);
    plot(abs(sensor_data_wsws_fft - sensor_data_wsws), 'k-');
    
end
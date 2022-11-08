function test_pass = kspaceFirstOrderAS_compare_symmetries(plot_comparisons, plot_simulations)
% DESCRIPTION:
%     Unit test to compare the axisymmetric code for a cicular piston with
%     different symmetries.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 7th December 2017
%     last update - 26th April 2018
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2017-2018 Bradley Treeby

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
COMPARISON_THRESH   = 1e-6;

% set pass variable
test_pass = true;

% create the computational grid
Nx = 128;            % number of grid points in the axial dimension
Ny = 64;           % number of grid points in the radial dimension
dx = 25e-3/Nx;    	% grid point spacing [m]
dy = dx;            % grid point spacing [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
medium.sound_speed = 1500;  % [m/s]
medium.alpha_coeff = 0.1;
medium.density = 1000;

% create the time array
CFL = 0.3;
t_end = 15e-6;
kgrid.makeTime(medium.sound_speed, CFL, t_end);

% define source mask
source.p_mask = zeros(Nx, Ny);
source.p_mask(25, 1:Ny/4) = 1;

% define a gaussian time varying source
source.p = exp(-((kgrid.t_array - 2e-6)/2.5e-7).^2);

% define a single sensor point
sensor.mask = zeros(Nx, Ny);
sensor.mask(Nx/2, Ny/4) = 1;

% input arguements
input_args = {'UsekSpace', true, 'PMLAlpha', 2, 'Smooth', false};

% run the simulation
sensor_data1 = kspaceFirstOrderAS(kgrid, medium, source, sensor, input_args{:}, ...
    'RadialSymmetry', 'WSWA-FFT', 'PlotSim', plot_simulations);
sensor_data2 = kspaceFirstOrderAS(kgrid, medium, source, sensor, input_args{:}, ...
    'RadialSymmetry', 'WSWS-FFT', 'PlotSim', plot_simulations);

% compute errors
diff_12 = max(abs(sensor_data1(:) - sensor_data2(:))) / max(abs(sensor_data1(:)));

% check for test pass
if (diff_12 > COMPARISON_THRESH)
    test_pass = false;
end

% plot
if plot_comparisons
    
    figure;
    subplot(2, 1, 1)
    plot(sensor_data1);
    hold on;
    plot(sensor_data2, '--');
    legend('WSWA', 'WSWS');
    title('Recorded signals');

    subplot(2, 1, 2)
    plot(abs(sensor_data1 - sensor_data2));
    title('Difference: WSWA to WSWS');
    
end
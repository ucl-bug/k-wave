function test_pass = kspaceFirstOrderAS_time_reversal(plot_comparisons, plot_simulations)
% DESCRIPTION:
%     Unit test to compare time reversal image reconstruction using the
%     axisymmetric code.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 28th January 2019
%     last update - 25th July 2019
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2019- Bradley Treeby

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

% set comparison threshold
COMPARISON_THRESH = 3;

% set pass variable
test_pass = true;

% create the computational grid
PML_size = 10;
Nx = 128;
Ny = 64;
dx = 25e-3/Nx;    
dy = dx;          
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define an expanded grid for creating the source mask
kgrid_exp = kWaveGrid(Nx, dx, Ny*2, dy);

% define the properties of the propagation medium
medium.sound_speed = 1500;

% define time
CFL = 0.3;
t_end = 20e-6;
kgrid.makeTime(medium.sound_speed, CFL, t_end);

% define a smoothed source mask
p0 = makeDisc(Nx, 2*Ny, Nx/2 + 5, Ny + 1, Nx/8);
p0 = smooth(p0, true, 'Blackman');
p0(:, 1:Ny) = [];
source.p0 = p0;

% define a circular sensor mask
sensor.mask = makeCircle(Nx, 2*Ny, Nx/2, Ny + 1, Nx/2 - 2 * PML_size);
sensor.mask(:, 1:Ny) = [];

% set input options
input_args = {...
    'PlotSim', plot_simulations, ...
    'DataCast', 'single', ...
    'RadialSymmetry', 'WSWA-FFT', ...
    'PMLSize', PML_size, ...
    'Smooth', false};

% run the simulation
sensor_data = kspaceFirstOrderAS(kgrid, medium, source, sensor, input_args{:});

% run time reversal simulation
source = rmfield(source, 'p0');
sensor.time_reversal_boundary_data = sensor_data;
p0_recon = kspaceFirstOrderAS(kgrid, medium, source, sensor, input_args{:});

% trim the data to avoid the reconstruction outside the sensor surface
trim = 20;
p0 = p0(1 + trim:end - trim, 1:end - trim);
p0_recon = p0_recon(1 + trim:end - trim, 1:end - trim);

% compute relative L_inf error [%]
error = 100 * max(abs(p0(:) - p0_recon(:))) ./ max(abs(p0(:)));

% check for test pass
if error > COMPARISON_THRESH
    test_pass = false;
end

% plot
if plot_comparisons

    figure;
    subplot(1, 3, 1);
    imagesc(p0);
    colorbar;
    axis image;
    title('p0');

    subplot(1, 3, 2);
    imagesc(p0_recon);
    colorbar;
    axis image;
    title('Reconstruction');

    subplot(1, 3, 3);
    imagesc(100 * abs(p0_recon - p0));
    colorbar;
    axis image;
    title('Error [%]');
    
    scaleFig(1.5, 1);
    
end
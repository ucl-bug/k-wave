function test_pass = kspaceFirstOrderAS_check_PML_outside(plot_comparisons, plot_simulations)
% DESCRIPTION:
%     Unit test to compare the axisymmetric code with PML inside and
%     outside.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 14th December 2017
%     last update - 14th December 2017
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
COMPARISON_THRESH   = 1e-15;

% set pass variable
test_pass = true;

% =========================================================================
% PML INSIDE - DEFAULT SIZE
% =========================================================================

% set PML size
pml_x_size = 20;
pml_y_size = 20;

% create the computational grid
Nx = 128; 
Ny = 64;  
dx = 25e-3/Nx;
dy = dx;      
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
medium.sound_speed = 1500;  % [m/s]
medium.alpha_coeff = 0.1;
medium.density = 1000;

% create the time array
dt = 50e-9;
Nt = 300;
kgrid.setTime(Nt, dt);

% define source
src_rad = 16;
source.p0 = zeros(Nx, Ny);
source.p0(pml_x_size + 1, 1:src_rad) = 1;

% record the maximum pressure everywhere
sensor.mask = zeros(Nx, Ny);
sensor.record = {'p_max_all'};

% run the simulation
sensor_data_pml_inside = kspaceFirstOrderAS(kgrid, medium, source, sensor, ...
    'PlotSim', plot_simulations, 'Smooth', false);
sensor_data_pml_inside = sensor_data_pml_inside.p_max_all;

% trim the pml
sensor_data_pml_inside = sensor_data_pml_inside(1 + pml_x_size:end - pml_x_size, 1:end - pml_y_size);

% =========================================================================
% PML OUTSIDE - DEFAULT SIZE
% =========================================================================

% create the computational grid
Nx = 128 - 2 * pml_x_size;  
Ny = 64 - pml_y_size;      
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% create the time array
kgrid.setTime(Nt, dt);

% define source
source.p0 = zeros(Nx, Ny);
source.p0(1, 1:src_rad) = 1;

% record the maximum pressure everywhere
sensor.mask = zeros(Nx, Ny);
sensor.record = {'p_max_all'};

% run the simulation
sensor_data_pml_outside = kspaceFirstOrderAS(kgrid, medium, source, sensor, ...
    'PlotSim', plot_simulations, 'PMLInside', false, 'Smooth', false);
sensor_data_pml_outside = sensor_data_pml_outside.p_max_all;

% compute errors
diff = max(abs(sensor_data_pml_inside(:) - sensor_data_pml_outside(:))) / max(abs(sensor_data_pml_inside(:)));

% check for test pass
if (diff > COMPARISON_THRESH)
    test_pass = false;
end

% plot
if plot_comparisons
    
    figure;
    subplot(3, 1, 1)
    imagesc(sensor_data_pml_inside);
    axis image;
    colorbar;
    title('PML Inside');

    subplot(3, 1, 2)
    imagesc(sensor_data_pml_outside);
    axis image;
    colorbar;
    title('PML Outside');
    
    subplot(3, 1, 3)
    imagesc(abs(sensor_data_pml_outside - sensor_data_pml_inside));
    axis image;
    colorbar;
    title('Difference');
    
end

% =========================================================================
% PML INSIDE - CUSTOM SIZE, P SOURCE
% =========================================================================

% set PML size
pml_x_size = 16;
pml_y_size = 15;

% create the computational grid
Nx = 128; 
Ny = 64;    
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% create the time array
kgrid.setTime(Nt, dt);

% define source mask
src_rad = 16;
clear source;
source.p_mask = zeros(Nx, Ny);
source.p_mask(pml_x_size + 1, 1:src_rad) = 1;

% define source
freq = 1e6;
source.p = createCWSignals(kgrid.t_array, freq, 1, 0);

% record the maximum pressure everywhere
sensor.mask = zeros(Nx, Ny);
sensor.record = {'p_max_all'};

% run the simulation
sensor_data_pml_inside = kspaceFirstOrderAS(kgrid, medium, source, sensor, ...
    'PlotSim', plot_simulations, 'Smooth', false, 'PMLSize', [pml_x_size, pml_y_size]);
sensor_data_pml_inside = sensor_data_pml_inside.p_max_all;

% trim the pml
sensor_data_pml_inside = sensor_data_pml_inside(1 + pml_x_size:end - pml_x_size, 1:end - pml_y_size);

% =========================================================================
% PML OUTSIDE - CUSTOM SIZE, P SOURCE
% =========================================================================

% create the computational grid
Nx = 128 - 2 * pml_x_size;  
Ny = 64 - pml_y_size;      
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% create the time array
kgrid.setTime(Nt, dt);

% define source
source.p_mask = zeros(Nx, Ny);
source.p_mask(1, 1:src_rad) = 1;

% record the maximum pressure everywhere
sensor.mask = zeros(Nx, Ny);
sensor.record = {'p_max_all'};

% run the simulation
sensor_data_pml_outside = kspaceFirstOrderAS(kgrid, medium, source, sensor, ...
    'PlotSim', plot_simulations, 'PMLInside', false, 'Smooth', false, 'PMLSize', [pml_x_size, pml_y_size]);
sensor_data_pml_outside = sensor_data_pml_outside.p_max_all;

% compute errors
diff = max(abs(sensor_data_pml_inside(:) - sensor_data_pml_outside(:))) / max(abs(sensor_data_pml_inside(:)));

% check for test pass
if (diff > COMPARISON_THRESH)
    test_pass = false;
end

% plot
if plot_comparisons
    
    figure;
    subplot(3, 1, 1)
    imagesc(sensor_data_pml_inside);
    axis image;
    colorbar;
    title('PML Inside');

    subplot(3, 1, 2)
    imagesc(sensor_data_pml_outside);
    axis image;
    colorbar;
    title('PML Outside');
    
    subplot(3, 1, 3)
    imagesc(abs(sensor_data_pml_outside - sensor_data_pml_inside));
    axis image;
    colorbar;
    title('Difference');
    
end
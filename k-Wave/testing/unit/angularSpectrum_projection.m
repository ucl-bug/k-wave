function test_pass = angularSpectrum_projection(plot_comparisons, plot_simulations)
% DESCRIPTION:
%     Unit test to project between the data measured on two parallel
%     planes. The data in this case is generated using k-Wave.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 19th March 2018
%     last update - 19th February 2019
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2018-2019 Bradley Treeby

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

% set comparison threshold [%]
comparison_threshold_grid_exp = 2.5;
comparison_threshold_forward  = 1;
comparison_threshold_backward = 35;
comparison_threshold_back_alt = 2e-5;

% set pass variable
test_pass = true;

% =========================================================================
% GENERATE DATA
% =========================================================================

% grid expansion
grid_expansion = 50;

% simulation settings
pml_size = 10;
Nx = 64 - 2*pml_size;
Ny = 64 - 2*pml_size;
Nz = 48 - 2*pml_size;
dx = 1e-4;
kgrid = kWaveGrid(Nx, dx, Ny, dx, Nz, dx);

% define medium properties
c0 = 1500;
medium.sound_speed = c0;

% define sensor mask
pos1 = 10;
pos2 = Nz;
sensor.mask = [1, 1, pos1, Nx, Ny, pos1; 1, 1, pos2, Nx, Ny, pos2].';

% calculate distance between two planes
proj_dist = (pos2 - pos1) * dx;

% define time axis, modifying dt so that the z/c0 offset between the two
% planes is an integer number of time steps
kgrid.makeTime(c0);
time_offset = proj_dist / c0;
time_offset_steps = round(time_offset / kgrid.dt);
kgrid.setTime(kgrid.Nt, time_offset / time_offset_steps);

% define source mask as focused bowl
radius = 33;
diameter = 25;
source.p_mask = makeBowl([Nx, Ny, Nz], [Nx/2, Ny/2, 1], radius, diameter, [Nx/2, Ny/2, Nz]);

% define input as toneburst
f0 = 2e6;
cycles = 3;
source.p = toneBurst(1/kgrid.dt, f0, cycles);

% run simulation
sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, ...
    'PMLInside', false', ...
    'PMLSize', pml_size, ...
    'PlotSim', plot_simulations, ...
    'DataCast', 'single');

% =========================================================================
% ANGULAR SPECTRUM PROJECTION - FORWARD
% =========================================================================

% set the input plane
plane_1 = squeeze(sensor_data(1).p);
plane_2 = squeeze(sensor_data(2).p);

% run projection, saving the output on the last plane
[~, pressure_time_forw] = angularSpectrum(plane_1, dx, ...
    kgrid.dt, proj_dist, c0, ...
    'Plot', plot_simulations, ...
    'GridExpansion', grid_expansion);

% run projection, saving the output on the last plane without grid expansion
[~, pressure_time_forw_no_exp] = angularSpectrum(plane_1, dx, ...
    kgrid.dt, proj_dist, c0, ...
    'Plot', plot_simulations);

% =========================================================================
% ANGULAR SPECTRUM PROJECTION - BACKWARD
% =========================================================================

% run projection, saving the output on the last plane
[~, pressure_time_back] = angularSpectrum(plane_2, dx, ...
    kgrid.dt, proj_dist, c0, ...
    'Plot', plot_simulations, ...
    'GridExpansion', grid_expansion, ...
    'Reverse', true);

% run projection, saving the output on the last plane
[~, pressure_time_back_alt] = angularSpectrum(plane_2, dx, ...
    kgrid.dt, -proj_dist, c0, ...
    'Plot', plot_simulations, ...
    'GridExpansion', grid_expansion);

% =========================================================================
% COMPARISON
% =========================================================================

% compute relative L_inf error between two back projections
err_back_alt = 100 * max(pressure_time_back(:) - pressure_time_back_alt(:)) / max(pressure_time_back(:));

% compare relative L_inf error when using grid expansion
err_forw_exp = 100 * max(pressure_time_forw(:) - pressure_time_forw_no_exp(:)) / max(pressure_time_forw(:));

% trim the forward projection planes to the same length
plane_2_tr = plane_2(:, :, time_offset_steps + 1:end);
pressure_time_forw_tr = pressure_time_forw(:, :, 1:size(plane_2_tr, 3));

% trim the backward projection planes to the same length
pressure_time_back_tr = pressure_time_back(:, :, time_offset_steps + 1:end);
plane_1_tr = plane_1(:, :, 1:size(pressure_time_back_tr, 3));

% compute relative L_inf error (%)
err_forw = 100 * max (abs(squeeze(plane_2_tr(:)) - squeeze(pressure_time_forw_tr(:))) ) / max(abs(squeeze(plane_2_tr(:))));
err_back = 100 * max (abs(squeeze(plane_1_tr(:)) - squeeze(pressure_time_back_tr(:))) ) / max(abs(squeeze(plane_1_tr(:))));

% check for test pass
if (err_forw > comparison_threshold_forward) || ...
        (err_back > comparison_threshold_backward) || ...
        (err_back_alt > comparison_threshold_back_alt) || ...
        (err_forw_exp > comparison_threshold_grid_exp)
    test_pass = false;
end

% =========================================================================
% PLOTTING
% =========================================================================

if plot_comparisons

    % plot trace in the middle
    figure;
    subplot(2, 1, 1);
    plot(squeeze(plane_2_tr(Nx/2, Ny/2, :)));
    hold on;
    plot(squeeze(pressure_time_forw_tr(Nx/2, Ny/2, :)));
    xlabel('Time Index');
    ylabel('Pressure');
    legend('Reference', 'ASA', 'Location', 'Best');
    title('Forward Projection');
    
    subplot(2, 1, 2);
    plot(100 * abs(squeeze(plane_2_tr(Nx/2, Ny/2, :)) - squeeze(pressure_time_forw_tr(Nx/2, Ny/2, :))) ./ max(abs(squeeze(plane_2_tr(Nx/2, Ny/2, :)))));
    xlabel('Time Index');
    ylabel('Error [%]');
    
    % plot all of the data
    figure;
    subplot(3, 1, 1);
    imagesc(reshape(plane_2_tr, Nx*Ny, []));
    xlabel('Time Index');
    ylabel('Grid Index');
    cb = colorbar;
    title(cb, '[Pa]');
    title('Reference');
    
    subplot(3, 1, 2);
    imagesc(reshape(pressure_time_forw_tr, Nx*Ny, []));
    xlabel('Time Index');
    ylabel('Grid Index');
    cb = colorbar;
    title(cb, '[Pa]');
    title('Forward Projection');
    
    subplot(3, 1, 3);
    imagesc(100 * abs(reshape(plane_2_tr, Nx*Ny, []) - reshape(pressure_time_forw_tr, Nx*Ny, [])) ./ max(abs(plane_2_tr(:)))); 
    xlabel('Time Index');
    ylabel('Grid Index');
    cb = colorbar;
    title(cb, '[%]');
    title('Error');
    
    % plot trace in the middle
    figure;
    subplot(2, 1, 1);
    plot(squeeze(plane_1_tr(Nx/2, Ny/2, :)));
    hold on;
    plot(squeeze(pressure_time_back_tr(Nx/2, Ny/2, :)));
    xlabel('Time Index');
    ylabel('Pressure');
    legend('Reference', 'ASA', 'Location', 'Best');
    title('Backward Projection');
    
    subplot(2, 1, 2);
    plot(100 * abs(squeeze(plane_1_tr(Nx/2, Ny/2, :)) - squeeze(pressure_time_back_tr(Nx/2, Ny/2, :))) ./ max(abs(squeeze(plane_1_tr(Nx/2, Ny/2, :)))));
    xlabel('Time Index');
    ylabel('Error [%]');
    
    % plot all of the data
    figure;
    subplot(3, 1, 1);
    imagesc(reshape(plane_1_tr, Nx*Ny, []));
    xlabel('Time Index');
    ylabel('Grid Index');
    cb = colorbar;
    title(cb, '[Pa]');
    title('Reference');
    
    subplot(3, 1, 2);
    imagesc(reshape(pressure_time_back_tr, Nx*Ny, []));
    xlabel('Time Index');
    ylabel('Grid Index');
    cb = colorbar;
    title(cb, '[Pa]');
    title('Backward Projection');
    
    subplot(3, 1, 3);
    imagesc(100 * abs(reshape(plane_1_tr, Nx*Ny, []) - reshape(pressure_time_back_tr, Nx*Ny, [])) ./ max(abs(plane_1_tr(:)))); 
    xlabel('Time Index');
    ylabel('Grid Index');
    cb = colorbar;
    title(cb, '[%]');
    title('Error');
    
    % plot all of the data - two ways of going backwards
    figure;
    subplot(3, 1, 1);
    imagesc(reshape(pressure_time_back, Nx*Ny, []));
    xlabel('Time Index');
    ylabel('Grid Index');
    cb = colorbar;
    title(cb, '[Pa]');
    title('Backward Projection Reverse');
    
    subplot(3, 1, 2);
    imagesc(reshape(pressure_time_back_alt, Nx*Ny, []));
    xlabel('Time Index');
    ylabel('Grid Index');
    cb = colorbar;
    title(cb, '[Pa]');
    title('Backward Projection Negative Z');
    
    subplot(3, 1, 3);
    imagesc(100 * abs(reshape(pressure_time_back, Nx*Ny, []) - reshape(pressure_time_back_alt, Nx*Ny, [])) ./ max(abs(pressure_time_back(:)))); 
    xlabel('Time Index');
    ylabel('Grid Index');
    cb = colorbar;
    title(cb, '[%]');
    title('Error');
    
    % plot all of the data - forward with and without grid expansion
    figure;
    subplot(3, 1, 1);
    imagesc(reshape(pressure_time_forw, Nx*Ny, []));
    xlabel('Time Index');
    ylabel('Grid Index');
    cb = colorbar;
    title(cb, '[Pa]');
    title('Forward Projection - Expanded');
    
    subplot(3, 1, 2);
    imagesc(reshape(pressure_time_forw_no_exp, Nx*Ny, []));
    xlabel('Time Index');
    ylabel('Grid Index');
    cb = colorbar;
    title(cb, '[Pa]');
    title('Forward Projection - No Expansion');
    
    subplot(3, 1, 3);
    imagesc(100 * abs(reshape(pressure_time_forw, Nx*Ny, []) - reshape(pressure_time_forw_no_exp, Nx*Ny, [])) ./ max(abs(pressure_time_forw(:)))); 
    xlabel('Time Index');
    ylabel('Grid Index');
    cb = colorbar;
    title(cb, '[%]');
    title('Error');      
    
end
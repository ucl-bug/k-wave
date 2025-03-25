function test_pass = kWaveArray_multiple_hologram_elements_3D(plot_comparisons, plot_simulations)
% DESCRIPTION:
%       Unit test to verify the behaviour of multiple hologram elements
%       with different amplitudes and phases in the kWaveArray class.
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 4th March 2025
%       last update - 24th March 2025
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2025- Bradley Treeby

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

% check for image processing toolbox and fail gracefull
v = ver;
if ~any(strcmp('Image Processing Toolbox', {v.Name}))
    warning('Skipping test, image processing toolbox not installed.');
    return
end

% set comparison threshold
comparison_thresh = 5;  % 5% error threshold

% =========================================================================
% DEFINE GRID PROPERTIES
% =========================================================================

% grid properties
Nx = 63;
Ny = 63;
Nz = 63;
dx = 1e-3;
kgrid = kWaveGrid(Nx, dx, Ny, dx, Nz, dx);

% medium properties
c0 = 1500;  % [m/s]
medium.sound_speed = c0;

% time array
kgrid.makeTime(c0, 0.3, 40e-6);

% define source parameters
source_freq = 200e3;  % [Hz]

% define sensor to capture final pressure field
sensor.record = {'p_final'};

% =========================================================================
% SIMULATION WITH MULTIPLE HOLOGRAM ELEMENTS
% =========================================================================

% Create the kWaveArray
karray = kWaveArray();

% Define locations for two hologram elements
% First element - a plane offset in the x-direction
x_offset_1 = 15 * dx;
num_points_side = 15;
spacing = dx/2;

% Create a grid of points for the first hologram
[y_grid, z_grid] = meshgrid(-spacing * (num_points_side-1)/2:spacing:spacing * (num_points_side-1)/2);
x_grid = x_offset_1 * ones(size(y_grid));
integration_points_1 = [x_grid(:), y_grid(:), z_grid(:)]';

% Create amplitude and phase patterns for the first hologram
% Use a simple focused pattern (phase depends on distance from center)
center_point_1 = [0, 0, 0];
amp_1 = ones(1, size(integration_points_1, 2));
distances = sqrt(sum((integration_points_1 - repmat(center_point_1', 1, size(integration_points_1, 2))).^2, 1));
phase_1 = 2 * pi * distances / (c0 / source_freq);

% Second element - a plane offset in the negative x-direction
x_offset_2 = -15 * dx;
center_point_2 = center_point_1 + [x_offset_2, 0, 0];

% Create a grid of points for the second hologram
x_grid = x_offset_2 * ones(size(y_grid));
integration_points_2 = [x_grid(:), y_grid(:), z_grid(:)]';

% Create amplitude and phase patterns for the second hologram
% Use a different pattern (e.g., phase gradient in y direction)
amp_2 = ones(1, size(integration_points_2, 2)) * 0.7;  % 70% amplitude compared to first
phase_gradient = 2;  % phase gradient factor
phase_2 = phase_gradient * integration_points_2(2, :) / spacing;

% Add hologram elements to the array
hologram_area = (num_points_side * spacing)^2;
karray.addHologramElement(center_point_1, integration_points_1, source_freq, amp_1, phase_1, hologram_area);
karray.addHologramElement(center_point_2, integration_points_2, source_freq, amp_2, phase_2, hologram_area);

% Create source using binary mask from karray
source_holograms.p_mask = karray.getArrayBinaryMask(kgrid);

% Get distributed source signals
el_amp = [1.0, 3];  % Element-level amplitudes
el_phase = [0, pi/4]; % Element-level phases
source_holograms.p = karray.getDistributedSourceSignalCW(kgrid, el_amp, el_phase);

% Run hologram simulation
sensor_data_holograms = kspaceFirstOrder3D(kgrid, medium, source_holograms, sensor, 'PlotSim', plot_simulations);

% =========================================================================
% SIMULATION WITH SUPERPOSITION OF INDIVIDUAL HOLOGRAM ELEMENTS
% =========================================================================

% Run separate simulations for each hologram element and superimpose
sensor_data_individual = zeros(size(sensor_data_holograms.p_final));

for elem_idx = 1:2

    % Create array with just one element
    karray_single = kWaveArray();
    
    if elem_idx == 1
        karray_single.addHologramElement(center_point_1, integration_points_1, source_freq, amp_1, phase_1, hologram_area);
    else
        karray_single.addHologramElement(center_point_2, integration_points_2, source_freq, amp_2, phase_2, hologram_area);
    end
    
    % Create source
    source_single.p_mask = karray_single.getArrayBinaryMask(kgrid);
    source_single.p = karray_single.getDistributedSourceSignalCW(kgrid, el_amp(elem_idx), el_phase(elem_idx));
    
    % Run simulation
    sensor_data_single = kspaceFirstOrder3D(kgrid, medium, source_single, sensor, 'PlotSim', plot_simulations);
    
    % Add to superposition
    sensor_data_individual = sensor_data_individual + sensor_data_single.p_final;
    
end

% =========================================================================
% COMPARISON AND ERROR CALCULATION
% =========================================================================

% Compare the combined and superimposed fields
err = 100 * max(abs(sensor_data_holograms.p_final(:) - sensor_data_individual(:))) / max(abs(sensor_data_individual(:)));

% Check if error is within threshold
if err > comparison_thresh
    test_pass = false;
    disp(['Test failed: Error between combined and superimposed fields is ' num2str(err) '%, which exceeds threshold of ' num2str(comparison_thresh) '%']);
else
    disp(['Test passed: Error between combined and superimposed fields is ' num2str(err) '%']);
end

% plot pressure fields and error
if plot_comparisons   
    figure;
    
    % Plot central XY planes
    subplot(2, 3, 1);
    imagesc(squeeze(sensor_data_holograms.p_final(:, :, ceil(Nz/2))));
    title('Combined Holograms - XY plane');
    axis image;
    colorbar;
    
    subplot(2, 3, 2);
    imagesc(squeeze(sensor_data_individual(:, :, ceil(Nz/2))));
    title('Superimposed Holograms - XY plane');
    axis image;
    colorbar;
    
    subplot(2, 3, 3);
    imagesc(squeeze(abs(sensor_data_holograms.p_final(:, :, ceil(Nz/2)) - sensor_data_individual(:, :, ceil(Nz/2)))));
    title(['Difference - XY plane (Error: ' num2str(err, '%.2f') '%)']);
    axis image;
    colorbar;
    
    % Plot central XZ planes
    subplot(2, 3, 4);
    imagesc(squeeze(sensor_data_holograms.p_final(:, ceil(Ny/2), :)));
    title('Combined Holograms - XZ plane');
    axis image;
    colorbar;
    
    subplot(2, 3, 5);
    imagesc(squeeze(sensor_data_individual(:, ceil(Ny/2), :)));
    title('Superimposed Holograms - XZ plane');
    axis image;
    colorbar;
    
    subplot(2, 3, 6);
    imagesc(squeeze(abs(sensor_data_holograms.p_final(:, ceil(Ny/2), :) - sensor_data_individual(:, ceil(Ny/2), :))));
    title(['Difference - XZ plane (Error: ' num2str(err, '%.2f') '%)']);
    axis image;
    colorbar;
    
    sgtitle('Multiple Hologram Elements Test');
end

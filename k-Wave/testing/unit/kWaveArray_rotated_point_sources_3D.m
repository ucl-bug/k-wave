function test_pass = kWaveArray_rotated_point_sources_3D(plot_comparisons, plot_simulations)
% DESCRIPTION:
%       Unit test to compare the acoustic fields generated by on-grid point
%       sources with fields from a rotated array of off-grid point sources
%       created using the kWaveArray class.
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 26th February 2025
%       last update - 26th February 2025
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
kgrid.makeTime(c0, 0.3, 25e-6);

% define source parameters
source_grid_pos = 20:40;
num_source_points_side = length(source_grid_pos);
num_source_points_total = length(source_grid_pos).^2;
source_freq = 200e3;  % [Hz]
rotation_angle = 45;  % [degrees]

% define spatially varying source amplitudes and phases
amp = linspace(0.1, 2, num_source_points_total);
phase = linspace(0, pi, num_source_points_total);

% create source signals with varying amplitude and phase
source_signals = createCWSignals(kgrid.t_array, source_freq, amp, phase);

% define sensor to capture final pressure field
sensor.record = {'p_final'};

% =========================================================================
% SIMULATION WITH ON-GRID SOURCES
% =========================================================================

% define on-grid source mask with points along a line
source.p_mask = zeros(Nx, Ny, Nz);
source.p_mask(ceil(Nx/2), source_grid_pos, source_grid_pos) = 1;

% assign source signals to each source point
source.p = source_signals;

% run k-Wave simulation with on-grid sources
sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, 'PlotSim', plot_simulations);

% =========================================================================
% SIMULATION WITH ROTATED OFF-GRID SOURCES
% =========================================================================

% create empty kWaveArray
karray = kWaveArray();

% add custom point elements at same locations as on-grid sources
element_dim = 1;
measure = kgrid.dx.^element_dim;
for ind1 = 1:num_source_points_side
    for ind2 = 1:num_source_points_side
        z_pos = kgrid.z_vec(source_grid_pos(ind1));
        y_pos = kgrid.y_vec(source_grid_pos(ind2));
        karray.addCustomElement([0; y_pos; z_pos], measure, element_dim, ['ind 1: ' num2str(ind1) ' 2: ' num2str(ind2)]);
    end
end

% apply rotation to the array
translation = [0; 0; 0];  % no translation
rotation = [rotation_angle, 0, 0];
karray.setArrayPosition(translation, rotation);

% create source using binary mask from karray
source_off_grid.p_mask = karray.getArrayBinaryMask(kgrid);

% get distributed source signals
source_off_grid.p = karray.getDistributedSourceSignal(kgrid, source_signals);

% run k-Wave simulation with rotated off-grid sources
sensor_data_off_grid = kspaceFirstOrder3D(kgrid, medium, source_off_grid, sensor, 'PlotSim', plot_simulations);

% =========================================================================
% SIMULATION WITH HOLOGRAPHY ELEMENT
% =========================================================================

karray = kWaveArray();

% setup integration points
integration_points = [];
linear_ind = 1;
for ind1 = 1:num_source_points_side
    for ind2 = 1:num_source_points_side
        z_pos = kgrid.z_vec(source_grid_pos(ind1));
        y_pos = kgrid.y_vec(source_grid_pos(ind2));
        integration_points(:, linear_ind) = [0; y_pos; z_pos]; %#ok<AGROW>
        linear_ind = linear_ind + 1;
    end
end

karray.addHologramElement(integration_points, amp, phase);
karray.setArrayPosition(translation, rotation);

% create source using binary mask from karray
source_off_grid.p_mask = karray.getArrayBinaryMask(kgrid);

% get distributed source signals
el_amp = 1;
el_phase = 0;
source_off_grid.p = karray.getDistributedSourceSignalCW(kgrid, source_freq, el_amp, el_phase);

% run k-Wave simulation with rotated off-grid sources
sensor_data_holography = kspaceFirstOrder3D(kgrid, medium, source_off_grid, sensor, 'PlotSim', plot_simulations);

% =========================================================================
% COMPARISON AND ERROR CALCULATION
% =========================================================================

% final pressure fields for comparison
crop_range = 15:45;
ref = squeeze(sensor_data.p_final(ceil(Nx/2), :, :));
ref = ref(crop_range, crop_range);

for ind = 1:2
    if ind == 1
        rotated = squeeze(sensor_data_off_grid.p_final(ceil(Nx/2), :, :));
        str = 'Off grid';
    else
        rotated = squeeze(sensor_data_holography.p_final(ceil(Nx/2), :, :));
        str = 'Hologram';
    end

    % rotate the off-grid result back to match orientation of on-grid result
    rotated = imrotate(rotated, -rotation_angle, 'bicubic', 'crop');
    rotated = rotated(crop_range, crop_range);
    
    % compute error as percentage of maximum value
    err = 100 * max(abs(ref(:) - rotated(:))) / max(abs(ref(:)));
    
    % check if error is within threshold
    if err > comparison_thresh
        test_pass = false;
    end
    
    % plot pressure fields and error
    if plot_comparisons
    
        figure;
        
        subplot(1, 3, 1);
        imagesc(ref);
        title('On-grid Result');
        axis image;
        colorbar;
        
        subplot(1, 3, 2);
        imagesc(rotated);
        title(['Rotated Result (Rotated ' num2str(rotation_angle) '°)']);
        axis image;
        colorbar;
        
        subplot(1, 3, 3);
        imagesc(abs(ref - rotated));
        title(['Difference (Error: ' num2str(err, '%.2f') '%)']);
        axis image;
        colorbar;

        sgtitle(str);
        
    end

end

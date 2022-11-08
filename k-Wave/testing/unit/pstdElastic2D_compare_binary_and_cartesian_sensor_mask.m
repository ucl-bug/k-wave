function test_pass = pstdElastic2D_compare_binary_and_cartesian_sensor_mask(plot_comparisons, plot_simulations)
% DESCRIPTION:
%     Unit test to compare cartesian and binary sensor masks.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 21st March 2019
%     last update - 9th April 2019
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2019 Bradley Treeby

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

% set comparison threshold
comparison_thresh = 1e-14;

% set pass variable
test_pass = true;

% create the computational grid
Nx = 128;           
Ny = 128;          
dx = 25e-3/Nx;    
dy = dx;          
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
medium.sound_speed_compression = 1500 * ones(Nx, Ny);  % [m/s]
medium.sound_speed_compression(Nx/2:end, :) = 2000;
medium.sound_speed_shear = zeros(Nx, Ny);
medium.sound_speed_shear(Nx/2:end, :) = 1200;

medium.density = 1000 * ones(Nx, Ny);
medium.density(Nx/2:end, :) = 1200;

% define source mask
source.p0 = zeros(Nx, Ny);
source.p0(22, Ny/4:3*Ny/4) = 1;

% define Cartesian sensor points using points exactly on the grid
circ_mask = makeCircle(Nx, Ny, Nx/2, Ny/2, Nx/2 - 10);
x_points = kgrid.x(circ_mask == 1);
y_points = kgrid.y(circ_mask == 1);
sensor.mask = [x_points, y_points].';

% record all output variables
sensor.record = {...
    'p', ...
    'p_max', ...
    'p_min', ...
    'p_rms', ...
    'u', ...
    'u_max', ...
    'u_min', ...
    'u_rms', ...
    'u_non_staggered', ...
    'I', ...
    'I_avg'};
    
% run the simulation as normal
sensor_data_c_ln = pstdElastic2D(kgrid, medium, source, sensor, ...
    'PlotSim', plot_simulations);

% run the simulation using nearest-neighbour interpolation
sensor_data_c_nn = pstdElastic2D(kgrid, medium, source, sensor, ...
    'CartInterp', 'nearest', ...
    'PlotSim', plot_simulations);

% convert sensor mask
[sensor.mask, reorder_index] = cart2grid(kgrid, sensor.mask);

% run the simulation again
sensor_data_b = pstdElastic2D(kgrid, medium, source, sensor, ...
    'PlotSim', plot_simulations);

% reorder the binary sensor data
sensor_data_b.p                 = reorderBinarySensorData(sensor_data_b.p, reorder_index);
sensor_data_b.p_max             = reorderBinarySensorData(sensor_data_b.p_max, reorder_index);
sensor_data_b.p_min             = reorderBinarySensorData(sensor_data_b.p_min, reorder_index);
sensor_data_b.p_rms             = reorderBinarySensorData(sensor_data_b.p_rms, reorder_index);
sensor_data_b.ux                = reorderBinarySensorData(sensor_data_b.ux, reorder_index);
sensor_data_b.uy                = reorderBinarySensorData(sensor_data_b.uy, reorder_index);
sensor_data_b.ux_max            = reorderBinarySensorData(sensor_data_b.ux_max, reorder_index);
sensor_data_b.uy_max            = reorderBinarySensorData(sensor_data_b.uy_max, reorder_index);
sensor_data_b.ux_min            = reorderBinarySensorData(sensor_data_b.ux_min, reorder_index);
sensor_data_b.uy_min            = reorderBinarySensorData(sensor_data_b.uy_min, reorder_index);
sensor_data_b.ux_rms            = reorderBinarySensorData(sensor_data_b.ux_rms, reorder_index);
sensor_data_b.uy_rms            = reorderBinarySensorData(sensor_data_b.uy_rms, reorder_index);
sensor_data_b.ux_non_staggered	= reorderBinarySensorData(sensor_data_b.ux_non_staggered, reorder_index);
sensor_data_b.uy_non_staggered  = reorderBinarySensorData(sensor_data_b.uy_non_staggered, reorder_index);
sensor_data_b.Ix                = reorderBinarySensorData(sensor_data_b.Ix, reorder_index);
sensor_data_b.Iy                = reorderBinarySensorData(sensor_data_b.Iy, reorder_index);
sensor_data_b.Ix_avg            = reorderBinarySensorData(sensor_data_b.Ix_avg, reorder_index);
sensor_data_b.Iy_avg            = reorderBinarySensorData(sensor_data_b.Iy_avg, reorder_index);

% compute errors
err_p_nn = max(abs(sensor_data_c_nn.p(:) - sensor_data_b.p(:))) / max(abs(sensor_data_b.p(:)));
err_p_ln = max(abs(sensor_data_c_ln.p(:) - sensor_data_b.p(:))) / max(abs(sensor_data_b.p(:)));

err_p_max_nn = max(abs(sensor_data_c_nn.p_max(:) - sensor_data_b.p_max(:))) / max(abs(sensor_data_b.p_max(:)));
err_p_max_ln = max(abs(sensor_data_c_ln.p_max(:) - sensor_data_b.p_max(:))) / max(abs(sensor_data_b.p_max(:)));

err_p_min_nn = max(abs(sensor_data_c_nn.p_min(:) - sensor_data_b.p_min(:))) / max(abs(sensor_data_b.p_min(:)));
err_p_min_ln = max(abs(sensor_data_c_ln.p_min(:) - sensor_data_b.p_min(:))) / max(abs(sensor_data_b.p_min(:)));

err_p_rms_nn = max(abs(sensor_data_c_nn.p_rms(:) - sensor_data_b.p_rms(:))) / max(abs(sensor_data_b.p_rms(:)));
err_p_rms_ln = max(abs(sensor_data_c_ln.p_rms(:) - sensor_data_b.p_rms(:))) / max(abs(sensor_data_b.p_rms(:)));

err_ux_nn = max(abs(sensor_data_c_nn.ux(:) - sensor_data_b.ux(:))) / max(abs(sensor_data_b.ux(:)));
err_ux_ln = max(abs(sensor_data_c_ln.ux(:) - sensor_data_b.ux(:))) / max(abs(sensor_data_b.ux(:)));

err_uy_nn = max(abs(sensor_data_c_nn.uy(:) - sensor_data_b.uy(:))) / max(abs(sensor_data_b.uy(:)));
err_uy_ln = max(abs(sensor_data_c_ln.uy(:) - sensor_data_b.uy(:))) / max(abs(sensor_data_b.uy(:)));

err_ux_max_nn = max(abs(sensor_data_c_nn.ux_max(:) - sensor_data_b.ux_max(:))) / max(abs(sensor_data_b.ux_max(:)));
err_ux_max_ln = max(abs(sensor_data_c_ln.ux_max(:) - sensor_data_b.ux_max(:))) / max(abs(sensor_data_b.ux_max(:)));

err_uy_max_nn = max(abs(sensor_data_c_nn.uy_max(:) - sensor_data_b.uy_max(:))) / max(abs(sensor_data_b.uy_max(:)));
err_uy_max_ln = max(abs(sensor_data_c_ln.uy_max(:) - sensor_data_b.uy_max(:))) / max(abs(sensor_data_b.uy_max(:)));

err_ux_min_nn = max(abs(sensor_data_c_nn.ux_min(:) - sensor_data_b.ux_min(:))) / max(abs(sensor_data_b.ux_min(:)));
err_ux_min_ln = max(abs(sensor_data_c_ln.ux_min(:) - sensor_data_b.ux_min(:))) / max(abs(sensor_data_b.ux_min(:)));

err_uy_min_nn = max(abs(sensor_data_c_nn.uy_min(:) - sensor_data_b.uy_min(:))) / max(abs(sensor_data_b.uy_min(:)));
err_uy_min_ln = max(abs(sensor_data_c_ln.uy_min(:) - sensor_data_b.uy_min(:))) / max(abs(sensor_data_b.uy_min(:)));

err_ux_rms_nn = max(abs(sensor_data_c_nn.ux_rms(:) - sensor_data_b.ux_rms(:))) / max(abs(sensor_data_b.ux_rms(:)));
err_ux_rms_ln = max(abs(sensor_data_c_ln.ux_rms(:) - sensor_data_b.ux_rms(:))) / max(abs(sensor_data_b.ux_rms(:)));

err_uy_rms_nn = max(abs(sensor_data_c_nn.uy_rms(:) - sensor_data_b.uy_rms(:))) / max(abs(sensor_data_b.uy_rms(:)));
err_uy_rms_ln = max(abs(sensor_data_c_ln.uy_rms(:) - sensor_data_b.uy_rms(:))) / max(abs(sensor_data_b.uy_rms(:)));

err_ux_non_staggered_nn = max(abs(sensor_data_c_nn.ux_non_staggered(:) - sensor_data_b.ux_non_staggered(:))) / max(abs(sensor_data_b.ux_non_staggered(:)));
err_ux_non_staggered_ln = max(abs(sensor_data_c_ln.ux_non_staggered(:) - sensor_data_b.ux_non_staggered(:))) / max(abs(sensor_data_b.ux_non_staggered(:)));

err_uy_non_staggered_nn = max(abs(sensor_data_c_nn.uy_non_staggered(:) - sensor_data_b.uy_non_staggered(:))) / max(abs(sensor_data_b.uy_non_staggered(:)));
err_uy_non_staggered_ln = max(abs(sensor_data_c_ln.uy_non_staggered(:) - sensor_data_b.uy_non_staggered(:))) / max(abs(sensor_data_b.uy_non_staggered(:)));

err_Ix_nn = max(abs(sensor_data_c_nn.Ix(:) - sensor_data_b.Ix(:))) / max(abs(sensor_data_b.Ix(:)));
err_Ix_ln = max(abs(sensor_data_c_ln.Ix(:) - sensor_data_b.Ix(:))) / max(abs(sensor_data_b.Ix(:)));

err_Iy_nn = max(abs(sensor_data_c_nn.Iy(:) - sensor_data_b.Iy(:))) / max(abs(sensor_data_b.Iy(:)));
err_Iy_ln = max(abs(sensor_data_c_ln.Iy(:) - sensor_data_b.Iy(:))) / max(abs(sensor_data_b.Iy(:)));

err_Ix_avg_nn = max(abs(sensor_data_c_nn.Ix_avg(:) - sensor_data_b.Ix_avg(:))) / max(abs(sensor_data_b.Ix_avg(:)));
err_Ix_avg_ln = max(abs(sensor_data_c_ln.Ix_avg(:) - sensor_data_b.Ix_avg(:))) / max(abs(sensor_data_b.Ix_avg(:)));

err_Iy_avg_nn = max(abs(sensor_data_c_nn.Iy_avg(:) - sensor_data_b.Iy_avg(:))) / max(abs(sensor_data_b.Iy_avg(:)));
err_Iy_avg_ln = max(abs(sensor_data_c_ln.Iy_avg(:) - sensor_data_b.Iy_avg(:))) / max(abs(sensor_data_b.Iy_avg(:)));

% check for test pass
if (err_p_nn > comparison_thresh) || ...
        (err_p_ln > comparison_thresh) || ...
        (err_p_max_nn > comparison_thresh) || ...
        (err_p_max_ln > comparison_thresh) || ...
        (err_p_min_nn > comparison_thresh) || ...
        (err_p_min_ln > comparison_thresh) || ...
        (err_p_rms_nn > comparison_thresh) || ...
        (err_p_rms_ln > comparison_thresh) || ...
        (err_ux_nn > comparison_thresh) || ...
        (err_ux_ln > comparison_thresh) || ...
        (err_ux_ln > comparison_thresh) || ...        
        (err_ux_max_nn > comparison_thresh) || ...
        (err_ux_max_ln > comparison_thresh) || ...
        (err_ux_min_nn > comparison_thresh) || ...
        (err_ux_min_ln > comparison_thresh) || ...
        (err_ux_rms_nn > comparison_thresh) || ...
        (err_ux_rms_ln > comparison_thresh) || ...
        (err_ux_non_staggered_nn > comparison_thresh) || ...
        (err_ux_non_staggered_ln > comparison_thresh) || ...
        (err_uy_nn > comparison_thresh) || ...
        (err_uy_ln > comparison_thresh) || ...      
        (err_uy_max_nn > comparison_thresh) || ...
        (err_uy_max_ln > comparison_thresh) || ...
        (err_uy_min_nn > comparison_thresh) || ...
        (err_uy_min_ln > comparison_thresh) || ...
        (err_uy_rms_nn > comparison_thresh) || ...
        (err_uy_rms_ln > comparison_thresh) || ...
        (err_uy_non_staggered_nn > comparison_thresh) || ...
        (err_uy_non_staggered_ln > comparison_thresh) || ...
        (err_Ix_nn > comparison_thresh) || ...
        (err_Ix_ln > comparison_thresh) || ...
        (err_Ix_avg_nn > comparison_thresh) || ...
        (err_Ix_avg_ln > comparison_thresh) || ...        
        (err_Iy_nn > comparison_thresh) || ...
        (err_Iy_ln > comparison_thresh) || ...
        (err_Iy_avg_nn > comparison_thresh) || ...
        (err_Iy_avg_ln > comparison_thresh)
    test_pass = false;
end

% plot
if plot_comparisons

    figure;
    subplot(5, 1, 1);
    imagesc(sensor_data_c_ln.p);
    colorbar;
    title('Cartesian - Linear');

    subplot(5, 1, 2);
    imagesc(sensor_data_c_nn.p);
    colorbar;
    title('Cartesian - Nearest Neighbour');

    subplot(5, 1, 3);
    imagesc(sensor_data_b.p);
    colorbar;
    title('Binary');

    subplot(5, 1, 4);
    imagesc(abs(sensor_data_b.p - sensor_data_c_ln.p));
    colorbar;
    title('Diff (Linear - Binary)');
    
    subplot(5, 1, 5);
    imagesc(abs(sensor_data_b.p - sensor_data_c_nn.p));
    colorbar;
    title('Diff (Nearest Neighbour - Binary)');
    
end
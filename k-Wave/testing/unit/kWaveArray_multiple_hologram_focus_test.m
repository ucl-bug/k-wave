function test_pass = kWaveArray_multiple_hologram_focus_test(plot_comparisons, plot_simulations)
% DESCRIPTION:
%       Unit test to verify the focusing capabilities of multiple hologram 
%       elements arranged in a circle using the acoustic reciprocity 
%       principle.
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 25th March 2025
%       last update - 26th March 2025
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

% set comparison threshold (percentage error)
phase_comparison_thresh = 1e-3;

% =========================================================================
% DEFINE GRID PROPERTIES
% =========================================================================

% compute points per temporal period
source_freq = 250e3;  % [Hz]
ppw = 8;
c0 = 1500;  % [m/s]
cfl = 0.3;
dx = c0 / (ppw * source_freq);   % [m]
PPP = round(ppw / cfl);

% grid properties
Nx = 96;
Ny = 90;
Nz = 40;
kgrid = kWaveGrid(Nx, dx, Ny, dx, Nz, dx);

% medium properties
medium.sound_speed = c0;

% create the time array using an integer number of points per period
dt = 1 / (PPP * source_freq);
t_end = Nx * dx / c0;
Nt = round(t_end / dt);
kgrid.setTime(Nt, dt);

% =========================================================================
% DEFINE HOLOGRAM PROPERTIES
% =========================================================================

% Define parameters for the circle of hologram elements
num_holograms = 30;             % Number of hologram elements arranged in a circle
radius = Ny/3 * dx;             % Radius of the circle
hologram_size = 2.5 * dx;        % Size of each hologram element
num_points_per_hologram = 25;    % Number of integration points per hologram element

% Random seed for reproducibility
rng(42);

% =========================================================================
% FORWARD SIMULATION: POINT SOURCE AND HOLOGRAM SENSORS
% =========================================================================

% Create point source at the focus
source_pt.p_mask = zeros(Nx, Ny, Nz);
src_x = Nx/2 - 10;
src_y = Ny/2 + 5;
source_pt.p_mask(src_x, src_y, Nz/2 + 1) = 1;

% Create CW source
source_phase = 0;  % Can be any value
source_pt.p = createCWSignals(kgrid.t_array, source_freq, 1, source_phase);

% Create hologram array for sensing
karray = kWaveArray();

% Generate positions for hologram elements in a circle
theta = linspace(0, 2*pi*(num_holograms-1)/num_holograms, num_holograms);
hologram_positions = zeros(3, num_holograms);
hologram_points_list = cell(num_holograms, 1);

for hol_ind = 1:num_holograms

    % Calculate position of hologram element on the circle
    hologram_positions(:, hol_ind) = [radius * cos(theta(hol_ind)), radius * sin(theta(hol_ind)), 0];
    
    % Create tangential vectors in the plane perpendicular to normal_vector
    % For a circle in the XY plane, these are simple
    tangent1 = [-sin(theta(hol_ind)), cos(theta(hol_ind)), 0]';  % Tangent to the circle
    tangent2 = [0, 0, 1]';  % Vertical direction
    
    % Create a tighter grid of points
    point_spacing = hologram_size/(sqrt(num_points_per_hologram)+1);
    [u, v] = meshgrid(-hologram_size/2:point_spacing:hologram_size/2);
    
    % Initialize array for hologram points
    hologram_points = zeros(3, numel(u));
    
    % Place points on the plane perpendicular to the normal vector
    for j = 1:numel(u)
        % Calculate position relative to the hologram center
        offset = u(j) * tangent1 + v(j) * tangent2;
        hologram_points(:, j) = hologram_positions(:, hol_ind) + offset;
    end
    
    % Store the points for this hologram element
    hologram_points_list{hol_ind} = hologram_points;
    
    % Add hologram element to the array with random amplitudes and phases
    random_amp = 1 + rand(1, size(hologram_points, 2)) / 4;
    random_phase = 0 + rand(1, size(hologram_points, 2)) * pi / 2;
    hologram_area = hologram_size^2;
    
    karray.addHologramElement(hologram_positions(:, hol_ind), hologram_points, ...
        source_freq, random_amp, random_phase, hologram_area);
end

if plot_comparisons
    plotHologramPoints(hologram_points_list, hologram_positions);
end

% Set sensor mask using hologram array
sensor_hologram.mask = karray.getArrayBinaryMask(kgrid);

% Configure sensor to record time-varying pressure
sensor_hologram.record = {'p'};

% Run first simulation with point source and hologram sensors
sensor_data_pt_src = kspaceFirstOrder3D(kgrid, medium, source_pt, sensor_hologram, ...
    'PlotSim', plot_simulations, ...
    'PlotScale', 'auto');

% Extract amplitude and phase from the sensor data
[amp, phase] = extractAmpPhase(sensor_data_pt_src.p(:, end - PPP + 1:end), 1/kgrid.dt, source_freq, ...
    'Dim', 2, 'Window', 'Rectangular', 'FFTPadding', 1);

% Account for phase asymmetry between createCWSignals and extractAmpPhase
phase = phase - pi/2;

% Use combineSensorDataCW to get complex amplitude at each array element
sensor_data_complex = karray.combineSensorDataCW(kgrid, amp .* exp(1i * phase));

% =========================================================================
% REVERSE SIMULATION: HOLOGRAM SOURCE AND POINT SENSOR
% =========================================================================

% Set source mask using hologram array
source_hologram.p_mask = karray.getArrayBinaryMask(kgrid);

% Get distributed source signal for CW input
% Note: we use the complex conjugate of the recorded field to time-reverse
source_hologram.p = karray.getDistributedSourceSignalCW(kgrid, abs(sensor_data_complex), -angle(sensor_data_complex));

% Create sensor mask covering central plane
sensor_pt.mask = zeros(Nx, Ny, Nz);
sensor_pt.mask(:, :, Nz/2 + 1) = 1;
sensor_pt.record = {'p'};

% Run second simulation with hologram source
sensor_data_holo_src = kspaceFirstOrder3D(kgrid, medium, source_hologram, sensor_pt, ...
    'PlotSim', plot_simulations, ...
    'PlotScale', 'auto');

% Extract amplitude and phase at the point sensor, accounting for sin/cos
[amp, phase] = extractAmpPhase(sensor_data_holo_src.p(:, end - PPP + 1:end), 1/kgrid.dt, source_freq, ...
    'Dim', 2, 'Window', 'Rectangular', 'FFTPadding', 1);
phase = phase - pi/2;
amp = reshape(amp, Nx, Ny);
phase = reshape(phase, Nx, Ny);
phase_pt = squeeze(phase(src_x, src_y));

% Conjugate phase a final time to match the input phase
phase_pt = -phase_pt;

% =========================================================================
% FIND FOCUS AND CALCULATE ERROR
% =========================================================================

% Find the maximum amplitude in the pressure field
amp_masked = amp;
mask = makeDisc(Nx, Ny, 0, 0, Nx/3 - 5);
amp_masked(mask == 0) = 0;
[~, max_idx] = maxND(amp_masked);
focus_error = norm(max_idx - [src_x, src_y]);

% Calculate phase error
phase_error = abs(mod(phase_pt - source_phase, 2*pi)) * 180/pi;
phase_error = min(phase_error, 360 - phase_error);

% Output results
fprintf('Distance from expected focus: %.2f grid points\n', focus_error);
fprintf('Phase error: %.2f degrees\n', phase_error);

% Check if the focus is within threshold (should be within 2 grid points of the center)
focus_threshold = 2 * dx;
if (focus_error > focus_threshold)
    test_pass = false;
    fprintf('TEST FAILED: Focus error exceeds threshold\n');
else
    fprintf('TEST PASSED: Focus error within threshold\n');
end

% Check if the phase error is within threshold
if (phase_error > phase_comparison_thresh)
    test_pass = false;
    fprintf('TEST FAILED: Phase error exceeds threshold\n');
else
    fprintf('TEST PASSED: Phase error within threshold\n');
end

% =========================================================================
% PLOT COMPARISONS
% =========================================================================

if plot_comparisons

    figure;   
    imagesc(kgrid.y_vec*1000, kgrid.x_vec*1000, amp);
    axis image;
    colorbar;
    xlabel('y [mm]');
    ylabel('x [mm]');
    hold on;
    
    % Plot source and hologram positions
    plot(kgrid.y_vec(src_y)*1000, kgrid.x_vec(src_x)*1000, 'ro', 'MarkerSize', 20, 'LineWidth', 2);
    plot(hologram_positions(2, :)*1000, hologram_positions(1, :)*1000, 'ko', 'MarkerSize', 6);
    
    sgtitle(['Multiple Hologram Focus Test - ' num2str(num_holograms) ' elements']);
    
end

end

% Function to plot hologram_points_list
function plotHologramPoints(hologram_points_list, hologram_positions)

    figure;
    hold on;
    
    % Use different colors for each hologram element
    colors = jet(length(hologram_points_list));
    
    % Plot each hologram's integration points
    for hol_ind = 1:length(hologram_points_list)
        points = hologram_points_list{hol_ind};
        scatter3(points(2,:), points(1,:), points(3,:), 20, colors(hol_ind,:), 'filled');
        
        % Plot the hologram center position with a larger marker
        scatter3(hologram_positions(2,hol_ind), hologram_positions(1,hol_ind), hologram_positions(3,hol_ind), 100, colors(hol_ind,:), 'filled', 'MarkerEdgeColor', 'k');
    end
    
    % Plot focus point at origin
    scatter3(0, 0, 0, 200, 'k', 'filled', 'MarkerEdgeColor', 'w');
    text(0, 0, 0, '  Focus', 'FontSize', 12);
    
    % Connect each hologram center to the focus with a line
    for hol_ind = 1:length(hologram_points_list)
        line([hologram_positions(2,hol_ind), 0], [hologram_positions(1,hol_ind), 0], [hologram_positions(3,hol_ind), 0], 'Color', [0.5 0.5 0.5], 'LineStyle', '--');
    end
    
    % Set axis properties
    axis equal;
    grid on;
    view(45, 30);
    xlabel('y [m]');
    ylabel('x [m]');
    zlabel('z [m]');
    title('Hologram Integration Points');
    
    % Add legend
    legend_entries = cell(1, length(hologram_points_list));
    for hol_ind = 1:length(hologram_points_list)
        legend_entries{hol_ind} = ['Hologram ' num2str(hol_ind)];
    end
    legend(legend_entries, 'Location', 'eastoutside');

end

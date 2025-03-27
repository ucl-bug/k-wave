function test_pass = kWaveArray_complex_grid_weights(plot_comparisons, ~)
% DESCRIPTION:
%       Unit test to verify that getElementGridWeights and getArrayGridWeights
%       correctly return complex matrices for hologram elements that incorporate
%       amplitude and phase information.
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 25th March 2025
%       last update - 25th March 2025
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
end

% set pass variable
test_pass = true;

% set comparison threshold
comparison_thresh = 1e-12;  % Threshold for numerical precision differences

% =========================================================================
% DEFINE GRID PROPERTIES
% =========================================================================

% grid properties
Nx = 64;
Ny = 64;
Nz = 64;
dx = 1e-3;
kgrid = kWaveGrid(Nx, dx, Ny, dx, Nz, dx);

% define source frequency for hologram elements
source_freq = 200e3;  % [Hz]

% =========================================================================
% CREATE HOLOGRAM ELEMENTS
% =========================================================================

% Create the kWaveArray
karray = kWaveArray();

% Define locations for three hologram elements with different properties
% First element - a plane at x = 15*dx
x_offset_1 = 15 * dx;
num_points_side = 15;
spacing = dx/2;

% Create a grid of points for the first hologram
[y_grid, z_grid] = meshgrid(-spacing * (num_points_side-1)/2:spacing:spacing * (num_points_side-1)/2);
x_grid = x_offset_1 * ones(size(y_grid));
integration_points_1 = [x_grid(:), y_grid(:), z_grid(:)]';

% Create amplitude and phase patterns for the first hologram (focused)
center_point_1 = [15*dx, 0, 0];
amp_1 = ones(1, size(integration_points_1, 2));
distances = sqrt(sum((integration_points_1 - repmat(center_point_1', 1, size(integration_points_1, 2))).^2, 1));
phase_1 = 2 * pi * distances / (1500 / source_freq);

% Second element - a plane at x = -15*dx
x_offset_2 = -15 * dx;

% Create a grid of points for the second hologram
x_grid = x_offset_2 * ones(size(y_grid));
integration_points_2 = [x_grid(:), y_grid(:), z_grid(:)]';

% Create amplitude and phase patterns for the second hologram (phase gradient in y)
amp_2 = ones(1, size(integration_points_2, 2)) * 0.7;  % 70% amplitude
phase_gradient = 2;  % phase gradient factor
phase_2 = phase_gradient * integration_points_2(2, :) / spacing;

% Third element - a plane at z = 10*dx with varying amplitude pattern
z_offset_3 = 10 * dx;

% Create a grid of points for the third hologram
z_grid_3 = z_offset_3 * ones(size(y_grid));
integration_points_3 = [y_grid(:), x_grid(:), z_grid_3(:)]';

% Create amplitude and phase patterns for the third hologram (Gaussian amplitude)
center_y = 0;
center_x = 0;
sigma = 5 * spacing;
amp_3 = exp(-((integration_points_3(1, :) - center_y).^2 + (integration_points_3(2, :) - center_x).^2) / (2 * sigma^2));
phase_3 = zeros(1, size(integration_points_3, 2));  % constant phase

% Calculate hologram area
hologram_area = (num_points_side * spacing)^2;

% Add hologram elements to the array
karray.addHologramElement(center_point_1,     integration_points_1, source_freq, amp_1, phase_1, hologram_area);
karray.addHologramElement([x_offset_2, 0, 0], integration_points_2, source_freq, amp_2, phase_2, hologram_area);
karray.addHologramElement([0, 0, z_offset_3], integration_points_3, source_freq, amp_3, phase_3, hologram_area);

% =========================================================================
% TEST ELEMENT GRID WEIGHTS
% =========================================================================

% Get grid weights for each element
element_weights_1 = karray.getElementGridWeights(kgrid, 1);
element_weights_2 = karray.getElementGridWeights(kgrid, 2);
element_weights_3 = karray.getElementGridWeights(kgrid, 3);

% Check that the grid weights are complex for the first two elements, and
% real for the last element (0 phase)
is_complex_1 = ~isreal(element_weights_1);
is_complex_2 = ~isreal(element_weights_2);
is_real_3 = isreal(element_weights_3);

if ~is_complex_1 || ~is_complex_2 || ~is_real_3
    test_pass = false;
    disp('Test failed: Element grid weights should be complex, complex, and real matrices');
end

% =========================================================================
% TEST ARRAY GRID WEIGHTS
% =========================================================================

% Get array grid weights (sum of all element grid weights)
array_weights = karray.getArrayGridWeights(kgrid);

% Check that array weights are complex
is_array_complex = ~isreal(array_weights);
if ~is_array_complex
    test_pass = false;
    disp('Test failed: Array grid weights should be a complex matrix');
end

% =========================================================================
% TEST MANUAL SUPERPOSITION VS. getArrayGridWeights
% =========================================================================

% Manually superimpose element weights
manual_superposition = element_weights_1 + element_weights_2 + element_weights_3;

% Calculate difference between manual superposition and array weights
weights_diff = abs(manual_superposition - array_weights);
max_error = max(weights_diff(:));

if max_error > comparison_thresh
    test_pass = false;
    disp(['Test failed: Max error between manual superposition and array weights is ' num2str(max_error) ...
        ', which exceeds threshold of ' num2str(comparison_thresh)]);
end

% =========================================================================
% VISUALIZATION
% =========================================================================

if plot_comparisons
    figure;
    
    % Plot amplitude of element weights (middle slices)
    subplot(3, 3, 1);
    imagesc(squeeze(abs(element_weights_1(ceil(Nx/2), :, :))));
    title('Element 1 Amplitude (YZ slice)');
    axis image;
    colorbar;
    
    subplot(3, 3, 2);
    imagesc(squeeze(abs(element_weights_2(ceil(Nx/2), :, :))));
    title('Element 2 Amplitude (YZ slice)');
    axis image;
    colorbar;
    
    subplot(3, 3, 3);
    imagesc(squeeze(abs(element_weights_3(ceil(Nx/2), :, :))));
    title('Element 3 Amplitude (YZ slice)');
    axis image;
    colorbar;
    
    % Plot phase of element weights (middle slices)
    subplot(3, 3, 4);
    imagesc(squeeze(angle(element_weights_1(ceil(Nx/2), :, :))));
    title('Element 1 Phase (YZ slice)');
    axis image;
    colorbar;
    
    subplot(3, 3, 5);
    imagesc(squeeze(angle(element_weights_2(ceil(Nx/2), :, :))));
    title('Element 2 Phase (YZ slice)');
    axis image;
    colorbar;
    
    subplot(3, 3, 6);
    imagesc(squeeze(angle(element_weights_3(ceil(Nx/2), :, :))));
    title('Element 3 Phase (YZ slice)');
    axis image;
    colorbar;
    
    % Plot array weights and comparison
    subplot(3, 3, 7);
    imagesc(squeeze(abs(array_weights(ceil(Nx/2), :, :))));
    title('Array Weights Amplitude (YZ slice)');
    axis image;
    colorbar;
    
    subplot(3, 3, 8);
    imagesc(squeeze(angle(array_weights(ceil(Nx/2), :, :))));
    title('Array Weights Phase (YZ slice)');
    axis image;
    colorbar;
    
    subplot(3, 3, 9);
    imagesc(squeeze(weights_diff(ceil(Nx/2), :, :)));
    title(['Difference (Max Error: ' num2str(max_error, '%.2e') ')']);
    axis image;
    colorbar;
    
    sgtitle('Complex Grid Weights Test');
end

% =========================================================================
% TEST BINARY MASK CONSISTENCY
% =========================================================================

% Get binary masks for each element
mask_1 = karray.getElementBinaryMask(kgrid, 1);
mask_2 = karray.getElementBinaryMask(kgrid, 2);
mask_3 = karray.getElementBinaryMask(kgrid, 3);

% Get array binary mask
array_mask = karray.getArrayBinaryMask(kgrid);

% Manually combine element masks
manual_mask = mask_1 | mask_2 | mask_3;

% Check if array mask matches the manually combined mask
if ~isequal(array_mask, manual_mask)
    test_pass = false;
    disp('Test failed: Array binary mask does not match the combined element masks');
end

% Check that grid weights are non-zero only where mask is true
mask_check_1 = all(abs(element_weights_1(~mask_1)) < comparison_thresh);
mask_check_2 = all(abs(element_weights_2(~mask_2)) < comparison_thresh);
mask_check_3 = all(abs(element_weights_3(~mask_3)) < comparison_thresh);
array_mask_check = all(abs(array_weights(~array_mask)) < comparison_thresh);

if ~mask_check_1 || ~mask_check_2 || ~mask_check_3 || ~array_mask_check
    test_pass = false;
    disp('Test failed: Grid weights are non-zero outside the corresponding binary mask');
end

function test_pass = kWaveArray_combineSensorDataCW_sampling_test(~, ~)
% DESCRIPTION:
%       Unit test for the combineSensorDataCW.
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

% set pass variable
test_pass = true;

% set comparison threshold
comparison_thresh = 1e-12;

% =========================================================================
% DEFINE A SIMPLE GRID AND HOLOGRAM ARRAY
% =========================================================================

% Create a small grid
Nx = 20;
Ny = 20;
Nz = 20;
dx = 1e-3;
kgrid = kWaveGrid(Nx, dx, Ny, dx, Nz, dx);

% Set source frequency
source_freq = 250e3;

% Create a simple hologram array with two elements
karray = kWaveArray();

% Define integration points for first element
center1 = [kgrid.x_vec(5), kgrid.y_vec(10), kgrid.z_vec(10)];
[x1, y1, z1] = ndgrid(kgrid.x_vec(4:6), kgrid.y_vec(9:11), kgrid.z_vec(9:11));
integration_points1 = [x1(:)'; y1(:)'; z1(:)'];
amp1 = ones(1, size(integration_points1, 2));
phase1 = zeros(1, size(integration_points1, 2));
area1 = (3*dx)^2;

% Add first hologram element
karray.addHologramElement(center1, integration_points1, source_freq, amp1, phase1, area1);

% Define integration points for second element
center2 = [kgrid.x_vec(15), kgrid.y_vec(10), kgrid.z_vec(10)];
[x2, y2, z2] = ndgrid(kgrid.x_vec(14:16), kgrid.y_vec(9:11), kgrid.z_vec(9:11));
integration_points2 = [x2(:)'; y2(:)'; z2(:)'];
amp2 = 2 * ones(1, size(integration_points2, 2));  % Double amplitude
phase2 = pi/4 * ones(1, size(integration_points2, 2));  % Phase offset
area2 = (3*dx)^2;

% Add second hologram element
karray.addHologramElement(center2, integration_points2, source_freq, amp2, phase2, area2);

% =========================================================================
% TEST 1: DIRECT VERIFICATION WITH KNOWN INPUT-OUTPUT RELATIONSHIP
% =========================================================================

% Get the array binary mask
mask = karray.getArrayBinaryMask(kgrid);
mask_indices = find(mask);

% Get individual element grid weights to manually compute expected output
element1_weights = karray.getElementGridWeights(kgrid, 1);
element2_weights = karray.getElementGridWeights(kgrid, 2);

% Create synthetic complex sensor data
% Simple case: uniform field with known amplitude and phase
field_amplitude = 1.0;
field_phase = pi/3;
complex_sensor_data = field_amplitude * exp(1i * field_phase) * ones(length(mask_indices), 1);

% Calculate expected output manually for each element
% For element 1
element1_indices = find(element1_weights ~= 0);
element1_mask_ind = ismember(mask_indices, element1_indices);
weights_phase_rev = abs(element1_weights(element1_indices)).^2 ./ conj(element1_weights(element1_indices));
element1_expected = sum(complex_sensor_data(element1_mask_ind) .* weights_phase_rev);
element1_expected = element1_expected / (area1 / (dx^2));

% For element 2
element2_indices = find(element2_weights ~= 0);
element2_mask_ind = ismember(mask_indices, element2_indices);
weights_phase_rev = abs(element2_weights(element2_indices)).^2 ./ conj(element2_weights(element2_indices));
element2_expected = sum(complex_sensor_data(element2_mask_ind) .* weights_phase_rev);
element2_expected = element2_expected / (area2 / (dx^2));

% Expected combined output
expected_output = [element1_expected; element2_expected];

% Use combineSensorDataCW to get the actual output
actual_output = karray.combineSensorDataCW(kgrid, complex_sensor_data);

% Compare expected and actual outputs
error_abs = abs(expected_output - actual_output);
error_rel = error_abs ./ abs(expected_output);

fprintf('Element 1 - Expected: %.6f + %.6fi, Actual: %.6f + %.6fi\n', ...
    real(element1_expected), imag(element1_expected), real(actual_output(1)), imag(actual_output(1)));
fprintf('Element 2 - Expected: %.6f + %.6fi, Actual: %.6f + %.6fi\n', ...
    real(element2_expected), imag(element2_expected), real(actual_output(2)), imag(actual_output(2)));
fprintf('Relative error: [%.6e, %.6e]\n', error_rel(1), error_rel(2));

% Test passes if relative error is below threshold
if any(error_rel > comparison_thresh)
    test_pass = false;
    fprintf('TEST FAILED: Error exceeds threshold\n');
else
    fprintf('TEST PASSED: Error within threshold\n');
end

% =========================================================================
% TEST 2: ADDITIVE PROPERTY WITH MULTIPLE FIELD COMPONENTS
% =========================================================================

% Create a second synthetic field with different amplitude and phase
field2_amplitude = 0.5;
field2_phase = -pi/6;
complex_sensor_data2 = field2_amplitude * exp(1i * field2_phase) * ones(length(mask_indices), 1);

% Compute expected combined outputs separately
actual_output1 = karray.combineSensorDataCW(kgrid, complex_sensor_data);
actual_output2 = karray.combineSensorDataCW(kgrid, complex_sensor_data2);
actual_output_sum = actual_output1 + actual_output2;

% Compute actual combined output for the sum of fields
complex_sensor_data_sum = complex_sensor_data + complex_sensor_data2;
actual_output_combined = karray.combineSensorDataCW(kgrid, complex_sensor_data_sum);

% Compare results
error_sum = abs(actual_output_sum - actual_output_combined);
error_sum_rel = error_sum ./ abs(actual_output_sum);

fprintf('\nTesting additive property of combineSensorDataCW:\n');
fprintf('Max relative error for field addition: %.6e\n', max(error_sum_rel));

% Test passes if error is below threshold
if any(error_sum_rel > comparison_thresh)
    test_pass = false;
    fprintf('TEST FAILED: Additive property error exceeds threshold\n');
else
    fprintf('TEST PASSED: Additive property error within threshold\n');
end

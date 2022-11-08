function test_pass = acousticFieldPropagator_exp_vs_sine_propagator(plot_comparisons, ~)
% DESCRIPTION:
%     Unit test to compare methods for the extraction of amplitude and
%     phase using the exponential and sin propagators.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 9th January 2019
%     last update - 9th January 2019
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

% check for plot inputs, and set to true if nargin is zero (to allow the
% test to be run independent of runUnitTests)
if nargin == 0
    plot_comparisons = true;
end

% set comparison threshold
comparison_threshold = 0.01;

% set pass variable
test_pass = true;

% =========================================================================
% SETUP
% =========================================================================

% define grid size
Nx = 256;
Ny = 64;
Nz = 64;
dx = 0.1e-3;

% define source frequency and sound speed
f0 = 1e6;
c0 = 1500;

% input phase offset
phase_in = 0;

% angular frequency [rad/s]
w0 = 2 * pi * f0;

% compute time
t1 = 1.5 * dx * Nx / c0;

% =========================================================================
% TIME OFFSET 1D
% =========================================================================

% point source
amp_in = zeros(Nx, 1);
amp_in(1) = 1;

% compute field patterns using exponential propagator
[pressure_field_1, amp_exp_prop, ~] = acousticFieldPropagator(amp_in, phase_in, dx, f0, c0, 'Time', t1, 'UseRamp', false);
[pressure_field_2] = acousticFieldPropagator(amp_in, phase_in, dx, f0, c0, 'Time', t1 + pi / (2 * w0), 'UseRamp', false);

% take imaginary part of complex pressure field to give output for sin
% propagator
pressure_field_1 = imag(pressure_field_1);
pressure_field_2 = imag(pressure_field_2);

% compute amplitude and phase from pressure field at two times
amp_real_prop = sqrt(pressure_field_1.^2 + pressure_field_2.^2);

% compute error
err = max(abs(amp_real_prop(:) - amp_exp_prop(:))) / max(abs(amp_exp_prop(:)));

% check for test pass
if (err > comparison_threshold)
    test_pass = false;
end

% plot comparison
if plot_comparisons
    
    figure;

    subplot(1, 3, 1);
    plot(amp_exp_prop);
    title('Amplitude - Complex Propagator');
    colorbar

    subplot(1, 3, 2);
    plot(amp_real_prop);
    title('Amplitude - Real Propagator');
    colorbar

    subplot(1, 3, 3);
    plot(abs(amp_real_prop - amp_exp_prop));
    title('Amplitude - Diff');
    colorbar

    scaleFig(1.5, 1);
    
end

% =========================================================================
% TIME OFFSET 2D
% =========================================================================

% line with steered phase
aperture_width = 20;
amp_in = zeros(Nx, Ny);
amp_in(1, Ny/2 - aperture_width/2 + 1: Ny/2 + aperture_width/2) = 1;

% calculate phase offset
steering_angle = 11.32;
T = 1/f0;
phase_term = 2*pi/T * ( dx*(0:aperture_width - 1)*sind(steering_angle) / c0);
phase_in = zeros(Nx, Ny);
phase_in(1, Ny/2 - aperture_width/2 + 1: Ny/2 + aperture_width/2) = phase_term;

% compute field patterns using exponential propagator
[pressure_field_1, amp_exp_prop, ~] = acousticFieldPropagator(amp_in, phase_in, dx, f0, c0, 'Time', t1, 'UseRamp', false);
[pressure_field_2] = acousticFieldPropagator(amp_in, phase_in, dx, f0, c0, 'Time', t1 + pi / (2 * w0), 'UseRamp', false);

% take imaginary part of complex pressure field to give output for sin
% propagator
pressure_field_1 = imag(pressure_field_1);
pressure_field_2 = imag(pressure_field_2);

% compute amplitude and phase from pressure field at two times
amp_real_prop = sqrt(pressure_field_1.^2 + pressure_field_2.^2);

% compute error
err = max(abs(amp_real_prop(:) - amp_exp_prop(:))) / max(abs(amp_exp_prop(:)));

% check for test pass
if (err > comparison_threshold)
    test_pass = false;
end

% plot comparison
if plot_comparisons

    figure;

    subplot(1, 3, 1);
    imagesc(amp_exp_prop);
    title('Amplitude - Complex Propagator');
    colorbar

    subplot(1, 3, 2);
    imagesc(amp_real_prop);
    title('Amplitude - Real Propagator');
    colorbar

    subplot(1, 3, 3);
    imagesc(abs(amp_real_prop - amp_exp_prop));
    title('Amplitude - Diff');
    colorbar

    scaleFig(1.5, 1);
    
end

% =========================================================================
% PHASE OFFSET 2D
% =========================================================================

% line with steered phase
aperture_width = 20;
amp_in = zeros(Nx, Ny);
amp_in(1, Ny/2 - aperture_width/2 + 1: Ny/2 + aperture_width/2) = 1;

% calculate phase offset
steering_angle = 11.32;
T = 1/f0;
phase_term = 2*pi/T * ( dx*(0:aperture_width - 1)*sind(steering_angle) / c0);
phase_in = zeros(Nx, Ny);
phase_in(1, Ny/2 - aperture_width/2 + 1: Ny/2 + aperture_width/2) = phase_term;

% compute field patterns using exponential propagator
[pressure_field_1, amp_exp_prop, ~] = acousticFieldPropagator(amp_in, phase_in, dx, f0, c0, 'Time', t1, 'UseRamp', false);
[pressure_field_2] = acousticFieldPropagator(amp_in, phase_in + pi/2, dx, f0, c0, 'Time', t1, 'UseRamp', false);

% take imaginary part of complex pressure field to give output for sin
% propagator
pressure_field_1 = imag(pressure_field_1);
pressure_field_2 = imag(pressure_field_2);

% compute amplitude and phase from pressure field at two times
amp_real_prop = sqrt(pressure_field_1.^2 + pressure_field_2.^2);

% compute error
err = max(abs(amp_real_prop(:) - amp_exp_prop(:))) / max(abs(amp_exp_prop(:)));

% check for test pass
if (err > comparison_threshold)
    test_pass = false;
end

% plot comparison
if plot_comparisons

    figure;

    subplot(1, 3, 1);
    imagesc(amp_exp_prop);
    title('Amplitude - Complex Propagator');
    colorbar

    subplot(1, 3, 2);
    imagesc(amp_real_prop);
    title('Amplitude - Real Propagator');
    colorbar

    subplot(1, 3, 3);
    imagesc(abs(amp_real_prop - amp_exp_prop));
    title('Amplitude - Diff');
    colorbar

    scaleFig(1.5, 1);
    
end

% =========================================================================
% TIME OFFSET 3D
% =========================================================================

% line with steered phase
aperture_width = 20;
amp_in = zeros(Nx, Ny, Nz);
amp_in(1, Ny/2 - aperture_width/2 + 1: Ny/2 + aperture_width/2, Nz/2 - aperture_width/2 + 1: Nz/2 + aperture_width/2) = 1;

% calculate phase offset
phase_in = 0;

% compute field patterns using exponential propagator
[pressure_field_1, amp_exp_prop, ~] = acousticFieldPropagator(amp_in, phase_in, dx, f0, c0, 'Time', t1, 'UseRamp', false);
[pressure_field_2] = acousticFieldPropagator(amp_in, phase_in, dx, f0, c0, 'Time', t1 + pi / (2 * w0), 'UseRamp', false);

% take imaginary part of complex pressure field to give output for sin
% propagator
pressure_field_1 = imag(pressure_field_1);
pressure_field_2 = imag(pressure_field_2);

% compute amplitude and phase from pressure field at two times
amp_real_prop = sqrt(pressure_field_1.^2 + pressure_field_2.^2);

% compute error
err = max(abs(amp_real_prop(:) - amp_exp_prop(:))) / max(abs(amp_exp_prop(:)));

% check for test pass
if (err > comparison_threshold)
    test_pass = false;
end

% plot comparison
if plot_comparisons

    % take slice
    amp_exp_prop = squeeze(amp_exp_prop(:, :, Nz/2));
    amp_real_prop = squeeze(amp_real_prop(:, :, Nz/2));
    
    % plot
    figure;

    subplot(1, 3, 1);
    imagesc(amp_exp_prop);
    title('Amplitude - Complex Propagator');
    colorbar

    subplot(1, 3, 2);
    imagesc(amp_real_prop);
    title('Amplitude - Real Propagator');
    colorbar

    subplot(1, 3, 3);
    imagesc(abs(amp_real_prop - amp_exp_prop));
    title('Amplitude - Diff');
    colorbar

    scaleFig(1.5, 1);
    
end
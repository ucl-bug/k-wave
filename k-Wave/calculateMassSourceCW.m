function [source_estimate, output] = calculateMassSourceCW(measured_data, dx, freq, c0, source_offset, grid_expansion, varargin)
%CALCULATEMASSSOURCECW Compute k-Wave input plane from measured CW data.
%
% DESCRIPTION:
%     calculateMassSourceCW takes a measured 2D plane of complex pressure
%     values (e.g., measured using a hydrophone in a scanning tank under
%     steady state conditions) and calculates an equivalent CW additive
%     pressure source positioned a given distance away that recreates the
%     measured data when used as an input to k-Wave (i.e., by assigning the
%     equivalent source as an input to acousticFieldPropagator, or by
%     generating a time-varying input signal using createCWSignal and
%     assigning this to source.p with source.p_mode = 'additive'). The
%     measured data is assumed to be for a single-frequency continuous wave
%     source, specified by freq, and is given as a complex 2D matrix.
%
%     The equivalent source is calculated using an iterative optimisation
%     based on gradient descent, where functional gradients are calculated
%     using the adjoint. Both the forward and adjoint operators are
%     computed using acousticFieldPropagator. The calculation assumes the
%     propagation is linear and the medium is lossless. The algorithm and
%     approach are described in detail in [1].
%
%     If the source is larger than the measured input plane, for example,
%     if measuring a focused bowl transducer, a suitable value should be
%     specified for the value of grid_expansion. Note, the value of
%     source_offset does not need to match the position of the real source
%     in the experiment. 
%
%     An alternative approach to project measured data using k-Wave
%     is to directly use the measured data as a pressure source with a
%     Dirichlet boundary condition (i.e., by generating a time-varying
%     input signal from the measured data using createCWSignal and
%     assigning this to source.p and setting source.p_mode = 'dirichlet').
%     However, this approach leads to errors in the imposed spatial
%     gradient, which manifests as errors in the projected field. Thus, for
%     accurate holographic projections using k-Wave, it is recommended to
%     use calculateMassSourceCW to first calculate the input data (see [1]
%     for a comparison).
%
% USAGE:
%     source_estimate = calculateMassSourceCW(measured_data, dx, freq, c0, source_offset)
%     source_estimate = calculateMassSourceCW(measured_data, dx, freq, c0, source_offset, grid_expansion)
%     source_estimate = calculateMassSourceCW(measured_data, dx, freq, c0, source_offset, [], ...)
%     source_estimate = calculateMassSourceCW(measured_data, dx, freq, c0, source_offset, grid_expansion, ...)
%     [source_estimate, output] = calculateMassSourceCW(measured_data, dx, freq, c0, source_offset)
%     ...
%
% INPUTS:
%     measured_data       - 2D matrix of complex pressure values [Pa].
%     dx                  - Spatial step between grid points in the input
%                           plane [m]. 
%     freq                - Source frequency [Hz].
%     c0                  - Speed of sound in the medium [m/s].
%     source_offset       - Offset between the measured input plane and the
%                           source plane [grid points]. For example, if
%                           source_offset = 1, the input plane and source
%                           plane are on adjacent grid points.
%     grid_expansion      - Number of grid points used to expand the size
%                           of the estimated source plane in each lateral
%                           dimension relative to the measured input plane
%                           (set to 0 if not defined).
%
% OPTIONAL INPUTS:
%     Optional 'string', value pairs that may be used to modify the default
%     computational settings. 
%     
%     'Mask'              - Binary mask to select a subset of the measured
%                           data to compute the error functional (default =
%                           all the data, i.e., ones(size(measured_data))).
%     'NumSteps'          - Number of gradient descent steps (default =
%                           20).
%     'StepSize'          - Starting size of gradient descent step 
%                           (default = 0.5). 
%     'StepSizeIncrement' - Multiplicative factor used to increase the step
%                           size when the error is reduced (default = 1.1).
%     'StepSizeDecrement' - Multiplicative factor used to decrease the step
%                           size when the error is increased 
%                           (default = 0.5). 
%     'Plot'              - Boolean controlling whether the update steps
%                           are displayed (default = true).
%     'ReturnIterations'  - Boolean controlling whether the source estimate
%                           at each gradient descent step is returned. If
%                           set to true, the source_estimate output is
%                           given as a 3D matrix, each x-y plane
%                           corresponds to the source estimate at each step
%                           (default = false).
%     'UseCpp'            - Boolean controlling whether the simulations are
%                           run using the C++ implementation of 
%                           acousticFieldPropagator (default = false).
%
% OUTPUTS:
%     source_estimate     - If 'ReturnIterations' is false (the default),
%                           source_estimate is given as a 2D matrix of
%                           complex pressure values [Pa]. If
%                           'ReturnIterations' is true, source_estimate is
%                           given as a 3D matrix containing the source
%                           estimate after  each iteration, indexed as (x,
%                           y, iteration). 
%     output              - Structure containing details of the
%                           optimisation: 
%
%                           .linf_error
%                           .l2_error
%                           .step_size
%                           .number_steps
%                           .number_function_calls
%                           .modelled_data
%
% ABOUT:
%     author              - Bradley Treeby, Ben Cox, Felix Lucka
%     date                - 20th February 2018
%     last update         - 22nd February 2020
%
% REFERENCES:
%     [1] Treeby, B., Lucka, F., Martin, E., & Cox, B. T. (2018).
%     Equivalent-Source Acoustic Holography for Projecting Measured
%     Ultrasound Fields Through Complex Media. IEEE Transactions on
%     Ultrasonics, Ferroelectrics, and Frequency Control, 65(10),
%     1857-1864.
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2018-2020 Bradley Treeby, Ben Cox, Felix Lucka
%
% See also acousticFieldPropagator, calculateMassSource

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

% =========================================================================
% INPUT CHECKING
% =========================================================================

% define defaults
step_size               = 0.5;
step_size_incr          = 1.1;
step_size_decr          = 0.5;
num_update_steps        = 20;
plot_updates            = true;
use_cpp_code            = false;
return_source_updates   = false;
use_mask                = false;

% check for grid expansion
if nargin < 6 || isempty(grid_expansion)
    grid_expansion = 0;
end

% replace with user defined values if provided
if ~isempty(varargin)
    for input_index = 1:2:length(varargin)
        switch varargin{input_index}
            case 'Mask'
                use_mask = true;
                data_mask = varargin{input_index + 1};
            case 'NumSteps'
                num_update_steps = varargin{input_index + 1};
            case 'StepSize'
                step_size = varargin{input_index + 1};
            case 'StepSizeDecrement'
                step_size_decr = varargin{input_index + 1};
            case 'StepSizeIncrement'
                step_size_incr = varargin{input_index + 1};
            case 'Plot'
                plot_updates = logical(varargin{input_index + 1});
            case 'ReturnIterations'
                return_source_updates = logical(varargin{input_index + 1});
            case 'UseCpp'
                use_cpp_code = logical(varargin{input_index + 1});
            otherwise
                error('Unknown optional input.');
        end
    end
end

% grid size parameters
Nx = size(measured_data, 1) + 2 * grid_expansion;
Ny = size(measured_data, 2) + 2 * grid_expansion;
Nz = 1 + source_offset;

% set indices for grid expansion
x1 = 1 + grid_expansion;
x2 = Nx - grid_expansion;
y1 = 1 + grid_expansion;
y2 = Ny - grid_expansion;

% counter for the number of function evaluations
num_function_calls = 0;

% create mask to use all data if not specified
if ~use_mask
    data_mask = ones(size(measured_data));
end

% =========================================================================
% SOURCE OPTIMISATION
% =========================================================================

% create input plane for initial guess using AFP
amp_in                    = zeros(Nx, Ny, Nz);
amp_in(x1:x2, y1:y2, 1)   = abs(measured_data);
phase_in                  = zeros(Nx, Ny, Nz);
phase_in(x1:x2, y1:y2, 1) = angle(measured_data);

% get initial guess using AFP, conjugating the phase before and after
if use_cpp_code
    source_estimate = acousticFieldPropagatorC(amp_in, -phase_in, dx, freq, c0);
else
    source_estimate = acousticFieldPropagator(amp_in, -phase_in, dx, freq, c0);
end
source_estimate = source_estimate(:, :, end);
source_estimate = conj(source_estimate);

% update counter
num_function_calls = num_function_calls + 1;

% store source estimate
if return_source_updates
    
    % preallocate
    source_estimate_matrix = zeros(Nx, Ny, num_update_steps);
    
    % store
    source_estimate_matrix(:, :, 1) = source_estimate;
    
end

% compute first gradient
[current_error, gradient, modelled_data] = errorFunctional(...
    measured_data, source_estimate, data_mask, dx, freq, c0, ...
    source_offset, grid_expansion, use_cpp_code);

% update counter
num_function_calls = num_function_calls + 2;

% preallocate vectors
linf_err_conv     = zeros(num_update_steps, 1);
l2_err_conv       = zeros(num_update_steps, 1);
step_size_history = zeros(num_update_steps, 1);

% compute error in data
linf_err_conv(1) = 100 * max(abs(abs(measured_data(:)) - abs(modelled_data(:))) ) ./ max(abs(measured_data(:)));
l2_err_conv(1)   = 100 * sqrt( sum((abs(measured_data(:)) - abs(modelled_data(:))).^2 ) ./ sum( abs(measured_data(:)).^2 ) );

% open figure and update
if plot_updates
    f1 = figure;
    updatePlot(f1, measured_data, modelled_data, data_mask, linf_err_conv, l2_err_conv, step_size_history, 1);
end

% loop through update steps
for update_index = 2:num_update_steps

    % take steps
    step_error = current_error;
    while step_error >= current_error

        % take a step in the gradient direction
        source_estimate_updated = source_estimate - step_size .* gradient;

        % compute error
        [step_error, step_gradient, modelled_data] = errorFunctional(...
            measured_data, source_estimate_updated, data_mask, dx, ...
            freq, c0, source_offset, grid_expansion, use_cpp_code);
        
        % update counter
        num_function_calls = num_function_calls + 2;
        
        % if the error is increased, reduce the step size
        if step_error >= current_error
            step_size = step_size .* step_size_decr;
            disp(['Error increased, step size decreased to ' num2str(step_size)]);
        end

    end

    % keep the current value of source estimate and gradient
    source_estimate = source_estimate_updated;
    gradient = step_gradient;
    current_error = step_error;

    % store source estimate
    if return_source_updates
        source_estimate_matrix(:, :, update_index) = source_estimate;
    end
    
    % store step size
    step_size_history(update_index) = step_size;
    
    % increase step size
    step_size = step_size .* step_size_incr;
    disp(['Error reduced, step size increased to ' num2str(step_size)]);
    
    % compute relative error norms (not used for calculation)
    data_error = modelled_data - measured_data;
    data_error(data_mask == 0) = 0;
    linf_err_conv(update_index) = 100 * max(abs(data_error(:))) ./ max(abs(measured_data(:)));
    l2_err_conv(update_index)   = 100 * sqrt( sum(abs(data_error(:)).^ 2) ./ sum( abs(measured_data(:)).^2 ) );
    
    % update plot
    if plot_updates
        updatePlot(f1, measured_data, modelled_data, data_mask, linf_err_conv, l2_err_conv, step_size_history, update_index);
    end
    
end

% assign output
if nargout == 2
    output.linf_error            = linf_err_conv;
    output.l2_error              = l2_err_conv;
    output.step_size             = step_size_history;
    output.number_steps          = num_update_steps;
    output.number_function_calls = num_function_calls;
    output.modelled_data         = modelled_data;
end

% return all source updates if optional input is set
if return_source_updates
    source_estimate = source_estimate_matrix;
end

% =========================================================================

function [error, gradient, modelled_data] = errorFunctional(measured_data, source_estimate, data_mask, dx, freq, c0, source_offset, grid_expansion, use_cpp_code)
% author: Bradley Treeby, Ben Cox, Felix Lucka
% date: 19th February 2018
% last update: 9th September 2019

% grid size parameters
Nx = size(measured_data, 1) + 2 * grid_expansion;
Ny = size(measured_data, 2) + 2 * grid_expansion;
Nz = 1 + source_offset;

% set indices for grid expansion
x1 = 1 + grid_expansion;
x2 = Nx - grid_expansion;
y1 = 1 + grid_expansion;
y2 = Ny - grid_expansion;

% assign input plane
amp_in            = zeros(Nx, Ny, Nz);
amp_in(:, :, 1)   = abs(source_estimate);
phase_in          = zeros(Nx, Ny, Nz);
phase_in(:, :, 1) = angle(source_estimate);

% use source plane to compute hologram plane using AFP
if use_cpp_code
    projected_field = acousticFieldPropagatorC(amp_in, phase_in, dx, freq, c0);
else
    projected_field = acousticFieldPropagator(amp_in, phase_in, dx, freq, c0);
end
modelled_data = projected_field(1 + grid_expansion:end - grid_expansion, 1 + grid_expansion:end - grid_expansion, end);

% compute difference between complex fields
data_error = modelled_data - measured_data;

% discard error values outside mask
data_error(data_mask == 0) = 0;

% compute error
error = sum(abs(data_error(:)).^ 2);

% assign input plane
amp_in                    = zeros(Nx, Ny, Nz);
amp_in(x1:x2, y1:y2, 1)   = abs(data_error);
phase_in                  = zeros(Nx, Ny, Nz);
phase_in(x1:x2, y1:y2, 1) = angle(data_error);

% compute gradient, conjugating the phase before and after
if use_cpp_code
    field = acousticFieldPropagatorC(amp_in, -phase_in, dx, freq, c0);
else
    field = acousticFieldPropagator(amp_in, -phase_in, dx, freq, c0);
end
gradient = field(:, :, end);
gradient = conj(gradient);

% =========================================================================

function updatePlot(figure_handle, measured_data, modelled_data, data_mask, linf_err_conv, l2_err_conv, step_size_history, update_index)
% author: Bradley Treeby
% date: 26th February 2018
% last update: 9th September 2019

% compute amplitude error
amp_error = abs(abs(measured_data) - abs(modelled_data));
amp_error(data_mask == 0) = 0;

% compute phase error
phase_error = angle(measured_data) - angle(modelled_data);
phase_error(phase_error > pi) = phase_error(phase_error > pi) - 2*pi;
phase_error(phase_error < -pi) = phase_error(phase_error < -pi) + 2*pi;
phase_error(data_mask == 0) = 0;

% select figure
figure(figure_handle);

% plot
subplot(3, 3, 1);
[~, scale, prefix] = scaleSI(max(abs(measured_data(:))));
imagesc(scale * abs(measured_data));
axis image;
cb = colorbar;
ylabel(cb, ['[' prefix 'Pa]']);
title('Measured Data Amplitude');
set(gca, 'XTick', [], 'YTick', []);

subplot(3, 3, 2);
imagesc(angle(measured_data));
axis image;
cb = colorbar;
ylabel(cb, '[rad]');
title('Measured Data Phase');
set(gca, 'XTick', [], 'YTick', []);

subplot(3, 3, 4);
[~, scale, prefix] = scaleSI(max(abs(modelled_data(:))));
imagesc(scale * abs(modelled_data));
axis image;
cb = colorbar;
ylabel(cb, ['[' prefix 'Pa]']);
title('Modelled Data Amplitude');
set(gca, 'XTick', [], 'YTick', []);

subplot(3, 3, 5);
imagesc(angle(modelled_data));
axis image;
colorbar;
cb = colorbar;
ylabel(cb, '[rad]');
title('Modelled Data Phase');
set(gca, 'XTick', [], 'YTick', []);

subplot(3, 3, 7);
[~, scale, prefix] = scaleSI(max(amp_error(:)));
imagesc(scale * amp_error);
axis image;
cb = colorbar;
ylabel(cb, ['[' prefix 'Pa]']);
title('Data Amplitude Error');
set(gca, 'XTick', [], 'YTick', []);

subplot(3, 3, 8);
imagesc(phase_error);
axis image;
cb = colorbar;
ylabel(cb, '[rad]');
title('Data Phase Error');
set(gca, 'XTick', [], 'YTick', []);

subplot(3, 3, 3);
semilogy(linf_err_conv(1:update_index));
xlabel('Iteration');
title('L\infty Error [%]');

subplot(3, 3, 6);
semilogy(l2_err_conv(1:update_index));
xlabel('Iteration');
title('L2 Error [%]');

subplot(3, 3, 9);
plot(step_size_history(1:update_index));
xlabel('Iteration');
title('Step Size');

drawnow;
function [source_estimate, output] = calculateMassSource(measured_data, dx, dt, c0, source_offset, grid_expansion, varargin)
%CALCULATEMASSSOURCE Compute k-Wave input plane from measured time-varying data.
%
% DESCRIPTION:
%     calculateMassSource takes a measured 2D plane of time-varying
%     pressure signals (e.g., measured using a hydrophone in a scanning
%     tank) and calculates an equivalent time-varying additive pressure
%     source positioned a given distance away that recreates the measured
%     data when used as an input to k-Wave (i.e., by assigning the
%     equivalent source to source.p with source.p_mode = 'additive').
%
%     The equivalent source is calculated using an iterative optimisation
%     based on gradient descent, where functional gradients are calculated
%     using the adjoint. Both the forward and adjoint operators are
%     computed using kspaceFirstOrder3D. The calculation assumes the
%     propagation is linear and the medium is lossless. The algorithm and
%     approach are described in detail in [1].
%
%     If the source is larger than the measured input plane, for example,
%     if measuring a focused bowl transducer, a suitable value should be
%     specified for the value of grid_expansion. Note, the value of
%     source_offset does not need to match the position of the real source
%     in the experiment. 
%
%     An alternative approach to project measured data using k-Wave is to
%     directly use the measured data as a pressure source with a Dirichlet
%     boundary condition (i.e., by assigning the measured data to source.p
%     and setting source.p_mode = 'dirichlet'). However, this approach
%     leads to errors in the imposed spatial gradient, which manifests as
%     errors in the projected field. Thus, for accurate holographic
%     projections using k-Wave, it is recommended to use
%     calculateMassSource to first calculate the input data (see [1] for a
%     comparison). 
%
% USAGE:
%     source_estimate = calculateMassSource(measured_data, dx, dt, c0, source_offset)
%     source_estimate = calculateMassSource(measured_data, dx, dt, c0, source_offset, grid_expansion)
%     [source_estimate, output] = calculateMassSource(measured_data, dx, dt, c0, source_offset, grid_expansion)
%     [source_estimate, output] = calculateMassSource(measured_data, dx, dt, c0, source_offset, grid_expansion, ...)
%     ...
%
% INPUTS:
%     measured_data       - 3D matrix containing the time varying pressure
%                           over a 2D input plane indexed as (x, y, t)
%                           [Pa]. 
%     dx                  - Spatial step between grid points in the input
%                           plane [m]. 
%     dt                  - Temporal step between time points in the input
%                           plane [s]. 
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
%                           given as a 4D matrix, each x-y plane
%                           corresponds to the source estimate at each step
%                           (default = false).
%     'UseCpp'            - Integer controlling whether the simulations are
%                           run using the C++ or CUDA implementations of
%                           kspaceFirstOrder3D, where
%                               0: Matlab code
%                               1: C++ CPU code
%                               2: C++ GPU code
%                           (default = 0).
%
% OUTPUTS:
%     source_estimate     - If 'ReturnIterations' is false (the default),
%                           source_estimate is given as a 3D matrix of
%                           pressure [Pa] indexed as (x, y, t). If
%                           'ReturnIterations' is true, source_estimate is
%                           given as a 4D matrix containing the source
%                           estimate after each iteration, indexed as 
%                           (x, y, t, iteration).
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
%     author              - Bradley Treeby
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
% Copyright (C) 2018-2020 Bradley Treeby
%
% See also calculateMassSourceCW, kspaceFirstOrder3D

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
use_cpp_code            = 0;
return_source_updates   = false;

% check for grid expansion
if nargin < 6 || isempty(grid_expansion)
    grid_expansion = 0;
end

% replace with user defined values if provided
if ~isempty(varargin)
    for input_index = 1:2:length(varargin)
        switch varargin{input_index}
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
                use_cpp_code = varargin{input_index + 1};
            otherwise
                error('Unknown optional input.');
        end
    end
end

% grid size parameters
Nx = size(measured_data, 1) + 2 * grid_expansion;
Ny = size(measured_data, 2) + 2 * grid_expansion;
Nt = size(measured_data, 3);

% set indices for grid expansion
x1 = 1 + grid_expansion;
x2 = Nx - grid_expansion;
y1 = 1 + grid_expansion;
y2 = Ny - grid_expansion;

% counter for the number of function evaluations
num_function_calls = 0;

% =========================================================================
% SOURCE OPTIMISATION
% =========================================================================

% create input plane for initial guess
input_plane = zeros(Nx, Ny, Nt);
input_plane(x1:x2, y1:y2, :) = measured_data;

% get initial guess, using adjoint
adjoint_flag = true;
source_estimate = kWaveProjection(input_plane, dx, dt, c0, ...
    source_offset, adjoint_flag, use_cpp_code);

% update counter
num_function_calls = num_function_calls + 1;

% store source estimate
if return_source_updates
    
    % preallocate
    source_estimate_matrix = zeros(Nx, Ny, Nt, num_update_steps);
    
    % store
    source_estimate_matrix(:, :, :, 1) = source_estimate;
    
end

% compute first gradient
[current_error, gradient, modelled_data] = errorFunctional(measured_data, ...
    source_estimate, dx, dt, c0, source_offset, grid_expansion, use_cpp_code);

% update counter
num_function_calls = num_function_calls + 2;

% preallocate vectors
linf_err_conv        = zeros(num_update_steps, 1);
l2_err_conv          = zeros(num_update_steps, 1);
step_size_history    = zeros(num_update_steps, 1);

% compute error in data
linf_err_conv(1) = 100 * max(abs( measured_data(:) - modelled_data(:) )) ./ max(abs(measured_data(:)));
l2_err_conv(1)   = 100 * sqrt( sum( (measured_data(:) - modelled_data(:)).^2 ) ./ sum( measured_data(:).^2 ));

% open figure and update
if plot_updates
    f1 = figure;
    updatePlot(f1, measured_data, modelled_data, linf_err_conv, l2_err_conv, step_size_history, 1, source_estimate);
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
            measured_data, source_estimate_updated, dx, dt, c0, ...
            source_offset, grid_expansion, use_cpp_code);
        
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
        source_estimate_matrix(:, :, :, update_index) = source_estimate;
    end    
    
    % store step size
    step_size_history(update_index) = step_size;
    
    % increase step size
    step_size = step_size .* step_size_incr;
    disp(['Error reduced, step size increased to ' num2str(step_size)]);
    
    % compute error in data
    linf_err_conv(update_index) = 100 * max(abs(measured_data(:) - modelled_data(:))) ./ max(abs(measured_data(:)));
    l2_err_conv(update_index)   = 100 * sqrt( ...
        sum((measured_data(:) - modelled_data(:)).^2) ./ sum(measured_data(:).^2) ...
        );
    
    % update plot
    if plot_updates
        updatePlot(f1, measured_data, modelled_data, linf_err_conv, l2_err_conv, step_size_history, update_index, source_estimate);
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

function [error, gradient, modelled_data] = errorFunctional(measured_data, source_estimate, dx, dt, c0, source_offset, grid_expansion, use_cpp_code)
% author: Bradley Treeby
% date: 20th February 2018
% last update: 20th March 2018

% grid size parameters
Nx = size(measured_data, 1) + 2 * grid_expansion;
Ny = size(measured_data, 2) + 2 * grid_expansion;
Nt = size(measured_data, 3);

% set indices for grid expansion
x1 = 1 + grid_expansion;
x2 = Nx - grid_expansion;
y1 = 1 + grid_expansion;
y2 = Ny - grid_expansion;

% use source plane to compute hologram plane
adjoint_flag = false;
modelled_data = kWaveProjection(source_estimate, dx, dt, c0, source_offset, adjoint_flag, use_cpp_code);
modelled_data = modelled_data(1 + grid_expansion:end - grid_expansion, 1 + grid_expansion:end - grid_expansion, :);

% compute difference between complex fields
data_error = modelled_data - measured_data;

% compute error
error = sum(data_error(:).^ 2);

% assign input plane
input_plane = zeros(Nx, Ny, Nt);
input_plane(x1:x2, y1:y2, :) = data_error;

% compute gradient using adjoint
adjoint_flag = true;
gradient = kWaveProjection(input_plane, dx, dt, c0, source_offset, adjoint_flag, use_cpp_code);

% =========================================================================

function output_plane = kWaveProjection(input_plane, dx, dt, c0, source_offset, adjoint, use_cpp_code)
% author: Bradley Treeby
% date: 20th February 2018
% last update: 25th February 2018

% time reverse if computing adjoint
if adjoint
    input_plane = flip(input_plane, 3);
end

% create k-Wave grid
kgrid = kWaveGrid(size(input_plane, 1), dx, size(input_plane, 2), dx, 1 + source_offset, dx);

% assign medium
medium.sound_speed = c0;

% work out how long to run the projection d = vt
diag_dist = sqrt(kgrid.x_size.^2 + kgrid.y_size.^2 + kgrid.z_size.^2);
Nt_projection = ceil(diag_dist / (dt * c0));
Nt = Nt_projection + size(input_plane, 3);

% create time array
kgrid.setTime(Nt, dt);

% create source mask
source.p_mask = zeros(kgrid.Nx, kgrid.Ny, kgrid.Nz);
if adjoint
    source.p_mask(:, :, end) = 1;
else
    source.p_mask(:, :, 1) = 1;
end

% assign source input
source.p = reshape(input_plane, kgrid.Nx * kgrid.Ny, []);

% create sensor mask
sensor.mask = zeros(kgrid.Nx, kgrid.Ny, kgrid.Nz);
if adjoint
    sensor.mask(:, :, 1) = 1;
else
    sensor.mask(:, :, end) = 1;
end

% set sensible PML size
pml_size = getOptimalPMLSize(kgrid);

% run simulation
switch use_cpp_code
    case 0
        
        % MATLAB code
        output_plane = kspaceFirstOrder3D(kgrid, medium, source, sensor, ...
            'PMLInside', false, ...
            'PMLSize', pml_size, ...
            'DataCast', 'single', ...
            'PlotSim', false);        
        
    case 1
        
        % C++ CPU code
        output_plane = kspaceFirstOrder3DC(kgrid, medium, source, sensor, ...
            'PMLInside', false, ...
            'PMLSize', pml_size);
        
    case 2
        
        % C++ GPU code
        output_plane = kspaceFirstOrder3DG(kgrid, medium, source, sensor, ...
            'PMLInside', false, ...
            'PMLSize', pml_size);
        
end

% reshape sensor data
output_plane = reshape(output_plane, kgrid.Nx, kgrid.Ny, []);

% work out where to trim the data (distance between planes)
Nt_trim = round(kgrid.z_size / (dt * c0));

% trim the time domain data to remove the propagation time from plane to
% plane
output_plane = output_plane(:, :, Nt_trim + 1:Nt_trim + size(input_plane, 3));

% time reverse if computing adjoint
if adjoint
    output_plane = flip(output_plane, 3);
end

% =========================================================================

function updatePlot(figure_handle, measured_data, modelled_data, linf_err_conv, l2_err_conv, step_size_history, update_index, source_estimate)
% author: Bradley Treeby
% date: 26th February 2018
% last update: 9th April 2018

% plot data
figure(figure_handle);
subplot(3, 3, 1);
[~, scale, prefix] = scaleSI(max(abs(measured_data(:))));
imagesc(scale * max(measured_data, [], 3));
axis image;
cb = colorbar;
ylabel(cb, ['[' prefix 'Pa]']);
title('Measured Data MIP');

subplot(3, 3, 4);
[~, scale, prefix] = scaleSI(max(abs(modelled_data(:))));
imagesc(scale * max(modelled_data, [], 3));
axis image;
cb = colorbar;
ylabel(cb, ['[' prefix 'Pa]']);
title('Modelled Data MIP');

subplot(3, 3, 7);
[~, scale, prefix] = scaleSI(max(abs(measured_data(:) - modelled_data(:))));
imagesc(scale * max(abs(measured_data - modelled_data), [], 3));
axis image;
cb = colorbar;
ylabel(cb, ['[' prefix 'Pa]']);
title('Error MIP');

subplot(3, 3, 2);
semilogy(linf_err_conv(1:update_index));
xlabel('Iteration');
title('L\infty Error [%]');

subplot(3, 3, 5);
semilogy(l2_err_conv(1:update_index));
xlabel('Iteration');
title('L2 Error [%]');

subplot(3, 3, 8);
plot(step_size_history(1:update_index));
xlabel('Iteration');
ylabel('Step Size');
title('Step Size');

subplot(3, 3, 3);
[~, scale, prefix] = scaleSI(max(abs(source_estimate)));
imagesc(scale * max(source_estimate, [], 3));
axis image;
cb = colorbar;
ylabel(cb, ['[' prefix 'Pa]']);
title('Source Estimate MIP');

subplot(3, 3, 6);
trace_measd = squeeze(measured_data(floor(end/2), floor(end/2), :));
trace_model = squeeze(modelled_data(floor(end/2), floor(end/2), :));
[~, scale, prefix] = scaleSI(max(abs(trace_measd)));
plot(scale * trace_measd);
hold on;
plot(scale * trace_model, '--');
hold off;
ylabel(['Pressure [' prefix 'Pa]']);
xlabel('Time Index');
title('Data Central Trace');
legend('Measured', 'Modelled', 'Location', 'Best');

subplot(3, 3, 9);
trace_src = squeeze(source_estimate(floor(end/2), floor(end/2), :));
[~, scale, prefix] = scaleSI(max(abs(trace_src)));
plot(scale * trace_src);
ylabel(['Pressure [' prefix 'Pa]']);
xlabel('Time Index');
title('Source Central Trace');

drawnow;
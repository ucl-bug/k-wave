function test_pass = angularSpectrumCW_data_cast(plot_comparisons, ~)
% DESCRIPTION:
%     Unit test to compare a projections with the datacast input.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 8th January 2019
%     last update - 19th February 2019
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2019- Bradley Treeby

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
comparison_threshold_single = 1e-6;
comparison_threshold_double = 1e-15;

% set pass variable
test_pass = true;

% =========================================================================
% RUN TEST
% =========================================================================

% simulation settings
Nx = 128;
Ny = 128;
Nz = 32;
dx = 1e-4;
f0 = 2e6;
c0 = 1500;
a0 = 2;
y  = 1.5;

% generate input plane as a disc
input_plane = makeDisc(Nx, Ny, Nx/2, Ny/2, Nx/4);

% run test with and without absorption
for ind = 1:2
    
    % assign medium
    if ind == 1
        medium = c0;
    else
        clear medium;
        medium.sound_speed = c0;
        medium.alpha_power = y;
        medium.alpha_coeff = a0;
    end

    % run simulation with different data cast flags
    output_data_double = angularSpectrumCW(input_plane, dx, (0:(Nz - 1)) * dx, f0, medium);
    output_data_single = angularSpectrumCW(input_plane, dx, (0:(Nz - 1)) * dx, f0, medium, ...
        'DataCast', 'single', ...
        'DataRecast', true);

    % take amplitude
    output_data_double = abs(output_data_double);
    output_data_single = abs(output_data_single);

    % compute relative L_inf error
    err_single = max(abs(output_data_double(:) - output_data_single(:))) / max(abs(output_data_double(:)));


    % try GPU implementations
    use_gpu = false;
    try

        % double precision on GPU
        output_data_gpu_double = angularSpectrumCW(input_plane, dx, (0:(Nz - 1)) * dx, f0, medium, ...
            'DataCast', 'gpuArray-double', ...
            'DataRecast', true);

        % single precision on GPU
        output_data_gpu_single = angularSpectrumCW(input_plane, dx, (0:(Nz - 1)) * dx, f0, medium, ...
            'DataCast', 'gpuArray-single', ...
            'DataRecast', true);

        % take amplitude
        output_data_gpu_double = abs(output_data_gpu_double);
        output_data_gpu_single = abs(output_data_gpu_single);

        % compute relative L_inf error
        err_gpu_double = max(abs(output_data_double(:) - output_data_gpu_double(:))) / max(abs(output_data_double(:)));
        err_gpu_single = max(abs(output_data_double(:) - output_data_gpu_single(:))) / max(abs(output_data_double(:)));

        % set flag
        use_gpu = true;

    catch

        % GPU code couldn't execute, so ignore errors
        err_gpu_double = 0;
        err_gpu_single = 0;

    end

    % check for test pass
    if (err_single > comparison_threshold_single) || ...
            (err_gpu_single > comparison_threshold_single) || ...
            (err_gpu_double > comparison_threshold_double)
        test_pass = false;
    end

    % plot maximum pressure
    if plot_comparisons

        figure;

        subplot(2, 2, 1);
        imagesc(max(output_data_double, [], 3));
        axis image;
        colorbar;
        title('Reference');

        subplot(2, 2, 2);
        imagesc(max(abs(output_data_single - output_data_double), [], 3) ./ max(abs(output_data_double(:))));
        axis image;
        colorbar;
        title('Error - Single Precision');

        if use_gpu

            subplot(2, 2, 3);
            imagesc(max(abs(output_data_gpu_single - output_data_double), [], 3) ./ max(abs(output_data_double(:))));
            axis image;
            colorbar;
            title('Error - GPU Single Precision');

            subplot(2, 2, 4);
            imagesc(max(abs(output_data_gpu_double - output_data_double), [], 3) ./ max(abs(output_data_double(:))));
            axis image;
            colorbar;
            title('Error - GPU Double Precision');        

        end

    end
    
end
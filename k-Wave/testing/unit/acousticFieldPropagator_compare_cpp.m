function test_pass = acousticFieldPropagator_compare_cpp(plot_comparisons, ~)
% DESCRIPTION:
%     Unit test to compare acousticFieldPropagator with a k-Wave simulation
%     using kspaceFirstOrder2D.
%
% ABOUT:
%     author        - Bradley Treeby and Jakub Budisky
%     date          - 28th March 2017
%     last update   - 25th April 2019
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2017-2019 Bradley Treeby

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
amplitude_max_error = 1e-4;
phase_max_error = 0.05;

% set pass variable
test_pass = true;

% set default name for the binary
if isunix
    binary_name = 'acousticFieldPropagator-OMP';
else
    binary_name = 'acousticFieldPropagator-OMP.exe';
end

% check for the C++ binary
if ~exist([getkWavePath('binaries'), binary_name], 'file') || ismac
    disp('WARNING: No binary found for acousticFieldPropagator. Exiting test.');
    return
end

% =========================================================================
% RUN TEST
% =========================================================================

% input grid size
Nx = 64;
Ny = 48;
Nz = 26;
dx = 0.1e-3;    % [m]

% expanded grid size
Nx_ex = [216 217 217 217 217 216 216 216];
Ny_ex = [192 192 193 192 193 193 193 192];
Nz_ex = [180 180 180 181 181 180 181 181];

% define the properties of the propagation medium
c0 = 1510;       % [m/s] 

% set the source frequency
f0 = 2e6;        % [Hz]

% set the amplitude
width = 9;
amp_in = zeros(Nx, Ny, Nz);
amp_in(1, Ny/2 - width:Ny/2 + width, Nz/2 - width:Nz/2 + width) = 1;

% loop through the propagators
for prop_index = 4
    
    % set the phase
    phase_in = zeros(Nx, Ny, Nz);
    switch prop_index
        case 1
            phase_in(amp_in == 1) = 0;
        case 2
            phase_in(amp_in == 1) = pi/2;
        case 3
            phase_in(amp_in == 1) = pi/4;
        case 4
            phase_in = pi * rand(Nx, Ny, Nz);
            amp_in = rand(Nx, Ny, Nz);
    end
            
    % loop through expanded sizes
    for index = 1:length(Nx_ex)

        % set the expanded grid size
        sz_ex = [Nx_ex(index), Ny_ex(index), Nz_ex(index)];

        % calculate beam pattern
        [amp_out, phase_out] = acousticFieldPropagator(amp_in, phase_in, dx, f0, c0, 'ExpandedGridSize', sz_ex);

        % run C++ code
        [amp_out_cpp, phase_out_cpp] = acousticFieldPropagatorC(amp_in, phase_in, dx, f0, c0, 'ExpandedGridSize', sz_ex);

        % compute the am,plitude error
        L_inf_amp = max(abs(amp_out_cpp(:) - amp_out(:))) / max(abs(amp_out(:))); 

        % compute the phase error accounting for phase wrapping
        angle_diff = @(ref, cmp) (mod(ref - cmp + pi, 2*pi) - pi);
        L_inf_phase = max(abs(angle_diff(phase_out(:), phase_out_cpp(:)))) / (2*pi);

        % print errors
        disp(['Expanded grid size = ' num2str(sz_ex)]);
        disp(['Amplitude error    = ' num2str(L_inf_amp)]);
        disp(['Phase error        = ' num2str(L_inf_phase)]);

        if (L_inf_amp < amplitude_max_error) && (L_inf_phase < phase_max_error)
            disp('Test Passed');
        else
            disp('Test Failed');
            test_pass = false;
        end

        if plot_comparisons

            % extract a plane
            xy_amp_matlab = squeeze(amp_out(:, :, Nz/2));
            xy_amp_cpp = squeeze(amp_out_cpp(:, :, Nz/2));
            xy_phase_matlab = squeeze(phase_out(:, :, Nz/2));
            xy_phase_cpp = squeeze(phase_out_cpp(:, :, Nz/2));

            % plot the field
            figure;
            subplot(2, 3, 1);
            imagesc(xy_amp_matlab);
            axis image;
            colorbar;
            title('MATLAB');

            subplot(2, 3, 2);
            imagesc(xy_amp_cpp);
            axis image;
            colorbar;
            title('C++');

            subplot(2, 3, 3);
            imagesc(abs(xy_amp_matlab - xy_amp_cpp) ./ max(abs(amp_out(:))));
            axis image;
            colorbar;
            title('Error');

            % plot the phase
            subplot(2, 3, 4);
            imagesc(xy_phase_matlab);
            axis image;
            colorbar;
            title('MATLAB');

            subplot(2, 3, 5);
            imagesc(xy_phase_cpp);
            axis image;
            colorbar;
            title('C++');

            subplot(2, 3, 6);
            imagesc(abs(xy_phase_matlab - xy_phase_cpp) ./ (2*pi));
            axis image;
            colorbar;
            title('Error');

            drawnow;

        end

    end
    
end

% display status
if test_pass
    disp('--------------------------------------------------');
    disp('acousticFieldPropagator C++ comparison - TEST PASS');
    disp('--------------------------------------------------');
else
    disp('--------------------------------------------------');
    disp('acousticFieldPropagator C++ comparison - TEST FAIL');
    disp('--------------------------------------------------');
end
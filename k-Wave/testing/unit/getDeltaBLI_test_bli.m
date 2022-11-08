function test_pass = getDeltaBLI_test_bli(plot_comparisons, ~)
%GETDELTABLI_TEST_BLI Check off-grid delta function's BLI.
%
% DESCRIPTION:
%     getDeltaBLI_test_bli samples an off-grid delta function, interpolates
%     using zero-padding, and compares with the exact expression.
%
% INPUTS:
%     plot_comparisons - Boolean controlling whether any comparisons (e.g.,
%                        error norms) are plotted 
%
% OUTPUTS:
%     test_pass        - Boolean denoting whether the test passed
%
% ABOUT:
%     author           - Elliott Wise
%     date             - 31st May 2018
%     last update      - 4th July 2019
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2018-2019 Elliott Wise

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
COMPARISON_THRESH   = 1e-13;

% set the number of upsampled points
NN = 10002;

% include the imaginary component of the off-grid delta function
include_imag = true;

% run test for both even and odd off-grid deltas
for loop_index = 1:2
    
    % specify the computational grid
    Nx = 30 + loop_index;
    dx = 2/Nx;
    x_grid = linspace(-1, 1, Nx+1);
    x_grid = x_grid(1:end-1);
    
    % create a list of delta function positions
    positions = linspace(0, dx, 7);
        
    % loop through positions
    for position_index = 1:length(positions)
        
        % get current position
        position = positions(position_index);
    
        % generate samples
        f_grid = getDeltaBLI(Nx, dx, x_grid, position, include_imag);
        
        % upsample in Fourier-space
        f_grid_k = fft(f_grid);
        nyqst = ceil((Nx+1)/2);
        f_upsampled_k = [f_grid_k(1:nyqst), zeros(1, NN-Nx), f_grid_k(nyqst+1:Nx)];
        f_upsampled = ifft(f_upsampled_k) * NN / Nx;
        
        % evaluate on upsampled grid
        x_upsampled = linspace(-1, 1, NN+1);
        x_upsampled = x_upsampled(1:end-1);
        f_eval = getDeltaBLI(Nx, dx, x_upsampled, position, include_imag);
        
        % plot
        if plot_comparisons
            figure
            subplot(2, 1, 1)
            plot(x_grid, real(f_grid), '.')
            hold on
            plot(x_upsampled, real(f_upsampled), '-')
            plot(x_upsampled, real(f_eval), '--')
            subplot(2, 1, 2)
            plot(x_grid, imag(f_grid), '.')
            hold on
            plot(x_upsampled, imag(f_upsampled), '-')
            plot(x_upsampled, imag(f_eval), '--')
        end
        
        % compare exact BLI expression with Fourier-space upsampling
        if max(abs(f_upsampled - f_eval)) > COMPARISON_THRESH
            test_pass = false;
        end
        
    end
    
end
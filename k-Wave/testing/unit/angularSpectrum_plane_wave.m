function test_pass = angularSpectrum_plane_wave(plot_comparisons, plot_simulations)
% DESCRIPTION:
%     Unit test to propagate a plane wave, checking the input and output
%     signals are the same.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 6th March 2018
%     last update - 13th February 2019
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2018- Bradley Treeby

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

% set comparison threshold
COMPARISON_THRESH = 1e-14;

% set pass variable
test_pass = true;

% simulation settings
Nx = 64;
Ny = 64;
Nz = 12;
dx = 1e-4;
f0 = 2e6;
c0 = 1500;

% calculate the time step using an integer number of points per period
cfl = 0.3;
ppw = c0 / (f0 * dx);     % points per wavelength
ppp = ceil(ppw / cfl);    % points per period
T   = 1 / f0;             % period [s]
dt  = T / ppp;            % time step [s]

% set number of points
Nt = ppp;

% create input signals
t_array = (0:Nt - 1) * dt;
sig = sin(2 * pi * f0 * t_array);

% create input plane
source_plane_time = repmat(reshape(sig, 1, 1, Nt), Nx, Ny, 1);

% run simulation
[~, output_data_tv] = angularSpectrum(source_plane_time, dx, dt, Nz * dx, ...
    c0, 'Plot', plot_simulations, 'FFTLength', Nx);

% compute error
err = max(output_data_tv(:) - source_plane_time(:));

% check for test pass
if err > COMPARISON_THRESH
    test_pass = false;
end

% plot comparison
if plot_comparisons
    
    output_sig = squeeze(output_data_tv(Nx/2, Ny/2, :)).';
    
    figure;
    subplot(2, 1, 1);
    plot(sig);
    hold on;
    plot(output_sig, '--')
    legend('Input', 'Output');
    title('Signals');
    
    subplot(2, 1, 2);
    plot(abs(sig - output_sig));
    title('Error');
    
end
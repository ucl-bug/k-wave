function test_pass = angularSpectrumCW_plane_wave(plot_comparisons, ~)
% DESCRIPTION:
%     Unit test to propagate a plane wave, checking the input and output
%     signals are the same.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 6th March 2018
%     last update - 19th February 2019
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2018-2019 Bradley Treeby

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

% generate input plane
input_plane = ones(Nx, Ny);

% run simulation
pressure = abs(angularSpectrumCW(input_plane, dx, Nz * dx, f0, c0, 'FFTLength', Nx));

% compute error
err = max(pressure(:) - 1);

% check for test pass
if err > COMPARISON_THRESH
    test_pass = false;
end

% plot comparison
if plot_comparisons
    
    figure;
    imagesc(pressure - 1);
    colorbar
    title('Error in output plane');
    
end